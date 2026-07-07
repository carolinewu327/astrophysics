#!/usr/bin/env python
"""Compare BOSS observed kappa stacks with BigMDPL simulation stacks."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


logger = logging.getLogger(__name__)


RPERP_CONFIG = {
    "rperp5": {
        "center": 5.0,
        "boss_suffix": "5_frac100",
        "note": "BOSS 5 h^-1 Mpc, full random",
    },
    "rperp10": {
        "center": 10.0,
        "boss_suffix": "10_frac100",
        "note": "BOSS 10 h^-1 Mpc, full random",
    },
    "rperp20": {
        "center": 20.0,
        "boss_suffix": "20_rpar10_frac100",
        "note": "BOSS 20 h^-1 Mpc, rpar10 full random",
    },
}


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_map(path: Path) -> np.ndarray:
    arr = pd.read_csv(path, index_col=0).to_numpy(dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D CSV map at {path}, got shape {arr.shape}")
    if not np.isfinite(arr).all():
        raise ValueError(f"Non-finite values found in {path}")
    return arr


def axis_for_map(arr: np.ndarray, box_size_hmpc: float = 100.0) -> np.ndarray:
    n = arr.shape[0]
    if n % 2 == 1:
        return np.linspace(-0.5 * box_size_hmpc, 0.5 * box_size_hmpc, n)
    cell = box_size_hmpc / n
    return np.linspace(-0.5 * box_size_hmpc + 0.5 * cell, 0.5 * box_size_hmpc - 0.5 * cell, n)


def interpolate_single_template(
    single: np.ndarray,
    target_axis: np.ndarray,
    rperp_center: float,
) -> np.ndarray:
    try:
        from scipy.interpolate import RegularGridInterpolator
    except ImportError as exc:
        raise RuntimeError("scipy is required to interpolate the single-halo template") from exc

    source_axis = axis_for_map(single)
    target_x, target_y = np.meshgrid(target_axis, target_axis)
    interpolator = RegularGridInterpolator(
        (source_axis, source_axis),
        single,
        bounds_error=False,
        fill_value=0.0,
    )
    halo_offset = 0.5 * rperp_center
    left_points = np.column_stack([target_y.ravel(), (target_x + halo_offset).ravel()])
    right_points = np.column_stack([target_y.ravel(), (target_x - halo_offset).ravel()])
    template = interpolator(left_points) + interpolator(right_points)
    return template.reshape(target_x.shape)


def resample_map(
    arr: np.ndarray,
    source_axis: np.ndarray,
    target_axis: np.ndarray,
) -> np.ndarray:
    try:
        from scipy.interpolate import RegularGridInterpolator
    except ImportError as exc:
        raise RuntimeError("scipy is required to resample comparison maps") from exc

    target_x, target_y = np.meshgrid(target_axis, target_axis)
    interpolator = RegularGridInterpolator(
        (source_axis, source_axis),
        arr,
        bounds_error=False,
        fill_value=np.nan,
    )
    points = np.column_stack([target_y.ravel(), target_x.ravel()])
    out = interpolator(points).reshape(target_x.shape)
    if not np.isfinite(out).all():
        raise ValueError("Resampled map contains non-finite values")
    return out


def bridge_masks(axis: np.ndarray, rperp_center: float) -> tuple[np.ndarray, np.ndarray]:
    x_grid, y_grid = np.meshgrid(axis, axis)
    bridge = (np.abs(x_grid) <= 0.35 * rperp_center) & (np.abs(y_grid) <= 0.15 * rperp_center)
    side = (
        (np.abs(x_grid) <= 0.35 * rperp_center)
        & (np.abs(y_grid) >= 0.45 * rperp_center)
        & (np.abs(y_grid) <= 0.85 * rperp_center)
    )
    return bridge, side


def map_stats(
    arr: np.ndarray,
    axis: np.ndarray,
    rperp_center: float,
    prefix: str,
) -> dict[str, float]:
    bridge, side = bridge_masks(axis, rperp_center)
    center_idx = len(axis) // 2
    bridge_mean = float(arr[bridge].mean())
    side_mean = float(arr[side].mean())
    return {
        f"{prefix}_midpoint_kappa": float(arr[center_idx, center_idx]),
        f"{prefix}_bridge_mean_kappa": bridge_mean,
        f"{prefix}_sideband_mean_kappa": side_mean,
        f"{prefix}_bridge_excess_kappa": bridge_mean - side_mean,
    }


def count_pairs(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        return max(sum(1 for _ in f) - 1, 0)


def save_map_csv(path: Path, arr: np.ndarray, axis: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(arr, index=axis, columns=axis).to_csv(path)
    logger.info("Saved %s", path)


def common_norm(arrays: list[np.ndarray]):
    vmax = float(max(np.nanpercentile(np.abs(arr), 99.5) for arr in arrays))
    if vmax == 0.0:
        vmax = 1e-12
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def plot_map_comparison(rows: list[dict], output: Path) -> None:
    obs_norm = common_norm([row["obs_filament"] for row in rows])
    sim_norm = common_norm([row["sim_residual"] for row in rows])
    diff_norm = common_norm([row["obs_minus_sim_residual"] for row in rows])

    fig, axes = plt.subplots(
        len(rows),
        3,
        figsize=(12.5, 11.0),
        constrained_layout=True,
    )
    for row_i, row in enumerate(rows):
        rperp_center = row["rperp_center"]
        halo_offset = 0.5 * rperp_center
        panels = [
            ("BOSS filament", row["obs_filament"], obs_norm),
            ("Simulation residual", row["sim_residual"], sim_norm),
            ("BOSS - simulation", row["obs_minus_sim_residual"], diff_norm),
        ]
        for col_i, (title, arr, norm) in enumerate(panels):
            ax = axes[row_i, col_i]
            im = ax.imshow(
                arr,
                origin="lower",
                extent=(-50.0, 50.0, -50.0, 50.0),
                interpolation="nearest",
                cmap="RdBu_r",
                norm=norm,
                aspect="equal",
            )
            ax.axhline(0.0, color="black", lw=0.7, alpha=0.35)
            ax.axvline(0.0, color="black", lw=0.7, alpha=0.35)
            ax.scatter([-halo_offset, halo_offset], [0.0, 0.0], s=42, marker="x", c="black", linewidths=1.2)
            if row_i == 0:
                ax.set_title(title)
            ax.set_xlabel(r"$X$ [$h^{-1}$ Mpc]")
            ax.set_ylabel(f"{row['label']}\n" + r"$Y$ [$h^{-1}$ Mpc]" if col_i == 0 else r"$Y$ [$h^{-1}$ Mpc]")
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cbar.set_label(r"$\kappa$")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def plot_bridge_comparison(stats: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), constrained_layout=True)
    x = stats["rperp_center_hmpc"].to_numpy(dtype=float)

    axes[0].plot(x, stats["obs_filament_bridge_excess_kappa"], marker="o", lw=1.8, label="BOSS filament")
    axes[0].plot(x, stats["sim_residual_bridge_excess_kappa"], marker="o", lw=1.8, label="simulation residual")
    axes[0].axhline(0.0, color="black", lw=0.8, alpha=0.35)
    axes[0].set_xlabel(r"$r_\perp$ [$h^{-1}$ Mpc]")
    axes[0].set_ylabel(r"bridge excess $\kappa$")
    axes[0].set_title("Residual/Filament Bridge Excess")
    axes[0].legend(frameon=False)

    axes[1].plot(x, stats["obs_raw_bridge_excess_kappa"], marker="o", lw=1.8, label="BOSS corrected pair")
    axes[1].plot(x, stats["sim_raw_bridge_excess_kappa"], marker="o", lw=1.8, label="simulation pair")
    axes[1].axhline(0.0, color="black", lw=0.8, alpha=0.35)
    axes[1].set_xlabel(r"$r_\perp$ [$h^{-1}$ Mpc]")
    axes[1].set_ylabel(r"bridge excess $\kappa$")
    axes[1].set_title("Raw Pair Bridge Excess")
    axes[1].legend(frameon=False)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare BOSS observed stacks against simulation stacks.")
    parser.add_argument("--boss-results-dir", default="analysis/boss/results")
    parser.add_argument("--sim-results-dir", default="analysis/sim/results")
    parser.add_argument("--output-dir", default="analysis/sim/results/observed_comparison")
    parser.add_argument("--mass-label", default="mass13")
    parser.add_argument("--format", default="png", choices=["png", "pdf"])
    parser.add_argument("--rperp-labels", nargs="+", default=["rperp5", "rperp10", "rperp20"])
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()

    boss_dir = Path(args.boss_results_dir)
    sim_dir = Path(args.sim_results_dir)
    out_dir = Path(args.output_dir)

    sim_single = load_map(sim_dir / f"kappa_single_sim_{args.mass_label}.csv")

    rows = []
    stats_rows = []
    for label in args.rperp_labels:
        if label not in RPERP_CONFIG:
            raise ValueError(f"Unknown rperp label {label!r}; expected one of {sorted(RPERP_CONFIG)}")
        cfg = RPERP_CONFIG[label]
        rperp_center = float(cfg["center"])
        boss_suffix = cfg["boss_suffix"]

        obs_raw = load_map(boss_dir / f"kappa_corrected_pairs_{boss_suffix}_BOSS.csv")
        obs_filament = load_map(boss_dir / f"kappa_filament_{boss_suffix}_BOSS.csv")
        sim_raw = load_map(sim_dir / f"kappa_pairs_sim_{args.mass_label}_{label}.csv")

        sim_axis = axis_for_map(sim_raw)
        sim_template = interpolate_single_template(sim_single, sim_axis, rperp_center)
        sim_residual = sim_raw - sim_template
        obs_axis = axis_for_map(obs_filament)

        sim_residual_path = out_dir / f"kappa_pairs_sim_{args.mass_label}_{label}_fixed_residual.csv"
        save_map_csv(sim_residual_path, sim_residual, sim_axis)

        # Resample simulation 101x101 maps to the 100x100 BOSS grid only for
        # direct image differencing. Bridge statistics use native axes.
        sim_for_difference = resample_map(sim_residual, sim_axis, obs_axis)
        obs_minus_sim = obs_filament - sim_for_difference

        row_stats = {
            "rperp_label": label,
            "rperp_center_hmpc": rperp_center,
            "boss_product_note": cfg["note"],
            "boss_pair_count": None,
            "sim_pair_count": count_pairs(sim_dir / f"pairs_{args.mass_label}_{label}.csv"),
        }
        row_stats.update(map_stats(obs_raw, axis_for_map(obs_raw), rperp_center, "obs_raw"))
        row_stats.update(map_stats(obs_filament, obs_axis, rperp_center, "obs_filament"))
        row_stats.update(map_stats(sim_raw, sim_axis, rperp_center, "sim_raw"))
        row_stats.update(map_stats(sim_residual, sim_axis, rperp_center, "sim_residual"))
        stats_rows.append(row_stats)

        rows.append(
            {
                "label": label,
                "rperp_center": rperp_center,
                "obs_filament": obs_filament,
                "sim_residual": sim_for_difference,
                "obs_minus_sim_residual": obs_minus_sim,
            }
        )

    stats = pd.DataFrame(stats_rows)
    stats_path = out_dir / f"boss_vs_sim_bridge_stats_{args.mass_label}.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    stats.to_csv(stats_path, index=False)
    logger.info("Saved %s", stats_path)

    plot_bridge_comparison(
        stats,
        out_dir / f"boss_vs_sim_bridge_stats_{args.mass_label}.{args.format}",
    )
    plot_map_comparison(
        rows,
        out_dir / f"boss_vs_sim_map_comparison_{args.mass_label}.{args.format}",
    )


if __name__ == "__main__":
    main()
