#!/usr/bin/env python
"""Compare simulation-only kappa stacks across pair-separation bins."""

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


RPERP_CENTERS = {
    "rperp5": 5.0,
    "rperp10": 10.0,
    "rperp20": 20.0,
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


def make_single_template(
    single: np.ndarray,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    rperp_center: float,
) -> np.ndarray:
    try:
        from scipy.interpolate import RegularGridInterpolator
    except ImportError as exc:
        raise RuntimeError("scipy is required to interpolate the single-halo template") from exc

    single_axis = np.linspace(-49.5, 49.5, single.shape[0])
    interpolator = RegularGridInterpolator(
        (single_axis, single_axis),
        single,
        bounds_error=False,
        fill_value=0.0,
    )

    halo_offset = 0.5 * rperp_center
    left_points = np.column_stack([y_grid.ravel(), (x_grid + halo_offset).ravel()])
    right_points = np.column_stack([y_grid.ravel(), (x_grid - halo_offset).ravel()])
    template = interpolator(left_points) + interpolator(right_points)
    return template.reshape(x_grid.shape)


def bridge_masks(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_grid, y_grid = np.meshgrid(axis, axis)
    bridge = (np.abs(x_grid) <= 0.35) & (np.abs(y_grid) <= 0.15)
    side = (np.abs(x_grid) <= 0.35) & (np.abs(y_grid) >= 0.45) & (np.abs(y_grid) <= 0.85)
    axis_strip = np.abs(y_grid) <= 0.15
    return bridge, side, axis_strip


def summarize_bridge(
    label: str,
    rperp_center: float,
    pair_count: int,
    pair_norm: np.ndarray,
    residual_norm: np.ndarray,
    norm_axis: np.ndarray,
) -> dict[str, float | int | str]:
    bridge, side, _ = bridge_masks(norm_axis)
    center_idx = len(norm_axis) // 2

    raw_bridge = float(pair_norm[bridge].mean())
    raw_side = float(pair_norm[side].mean())
    residual_bridge = float(residual_norm[bridge].mean())
    residual_side = float(residual_norm[side].mean())

    return {
        "rperp_label": label,
        "rperp_center_hmpc": rperp_center,
        "pair_count": pair_count,
        "raw_midpoint_kappa": float(pair_norm[center_idx, center_idx]),
        "raw_bridge_mean_kappa": raw_bridge,
        "raw_sideband_mean_kappa": raw_side,
        "raw_bridge_excess_kappa": raw_bridge - raw_side,
        "residual_midpoint_kappa": float(residual_norm[center_idx, center_idx]),
        "residual_bridge_mean_kappa": residual_bridge,
        "residual_sideband_mean_kappa": residual_side,
        "residual_bridge_excess_kappa": residual_bridge - residual_side,
    }


def count_pairs(path: Path) -> int:
    with path.open("r", encoding="utf-8") as f:
        return max(sum(1 for _ in f) - 1, 0)


def common_norm(arrays: list[np.ndarray], diverging: bool = False):
    if diverging:
        vmax = float(max(np.nanpercentile(np.abs(arr), 99.5) for arr in arrays))
        return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax), "RdBu_r"

    vmin = float(min(np.nanmin(arr) for arr in arrays))
    vmax = float(max(np.nanpercentile(arr, 99.8) for arr in arrays))
    if vmin < 0.0:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax), "RdBu_r"
    return None, "magma"


def save_residual_csv(path: Path, residual: np.ndarray, axis: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(residual, index=axis, columns=axis).to_csv(path)
    logger.info("Saved %s", path)


def plot_comparison(
    rows: list[dict],
    stats: pd.DataFrame,
    output: Path,
    fmt: str,
) -> None:
    fixed_norm, fixed_cmap = common_norm([row["pair_fixed"] for row in rows])
    normed_norm, normed_cmap = common_norm([row["pair_norm"] for row in rows])
    residual_norm, residual_cmap = common_norm([row["residual_norm"] for row in rows], diverging=True)

    fig, axes = plt.subplots(
        len(rows),
        4,
        figsize=(16.0, 11.0),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.0, 1.12]},
    )

    for row_i, row in enumerate(rows):
        label = row["label"]
        rperp_center = row["rperp_center"]
        halo_offset = 0.5 * rperp_center
        stat = stats.loc[stats["rperp_label"] == label].iloc[0]

        panels = [
            (
                row["pair_fixed"],
                (-50.0, 50.0, -50.0, 50.0),
                fixed_norm,
                fixed_cmap,
                r"fixed pair",
                [(-halo_offset, 0.0), (halo_offset, 0.0)],
            ),
            (
                row["pair_norm"],
                (-2.5, 2.5, -2.5, 2.5),
                normed_norm,
                normed_cmap,
                r"normalized pair",
                [(-0.5, 0.0), (0.5, 0.0)],
            ),
            (
                row["residual_norm"],
                (-2.5, 2.5, -2.5, 2.5),
                residual_norm,
                residual_cmap,
                r"normalized residual",
                [(-0.5, 0.0), (0.5, 0.0)],
            ),
        ]

        for col_i, (arr, extent, norm, cmap, title, markers) in enumerate(panels):
            ax = axes[row_i, col_i]
            im = ax.imshow(
                arr,
                origin="lower",
                extent=extent,
                interpolation="nearest",
                cmap=cmap,
                norm=norm,
                aspect="equal",
            )
            ax.axhline(0.0, color="black", lw=0.6, alpha=0.35)
            ax.axvline(0.0, color="black", lw=0.6, alpha=0.35)
            for marker_x, marker_y in markers:
                ax.scatter([marker_x], [marker_y], s=38, marker="x", c="black", linewidths=1.1)
            if row_i == 0:
                ax.set_title(title)
            if col_i == 0:
                ax.set_ylabel(f"{label}\nY")
            else:
                ax.set_ylabel("Y")
            ax.set_xlabel("X")
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cbar.set_label(r"$\kappa$")

        ax = axes[row_i, 3]
        axis = row["norm_axis"]
        _, _, axis_strip = bridge_masks(axis)
        x_grid, _ = np.meshgrid(axis, axis)
        raw_profile = np.array([row["pair_norm"][axis_strip & np.isclose(x_grid, x)].mean() for x in axis])
        residual_profile = np.array([row["residual_norm"][axis_strip & np.isclose(x_grid, x)].mean() for x in axis])
        ax.plot(axis, raw_profile, color="#4c78a8", lw=1.7, label="raw")
        ax.plot(axis, residual_profile, color="#d62728", lw=1.7, label="residual")
        ax.axhline(0.0, color="black", lw=0.7, alpha=0.35)
        ax.axvline(-0.5, color="black", lw=0.7, ls=":", alpha=0.55)
        ax.axvline(0.5, color="black", lw=0.7, ls=":", alpha=0.55)
        ax.set_xlim(-1.1, 1.1)
        ax.set_xlabel(r"$X/r_\perp$")
        ax.set_ylabel(r"mean $\kappa$, $|Y/r_\perp|\leq0.15$")
        if row_i == 0:
            ax.set_title("axis-strip profile")
            ax.legend(frameon=False, fontsize=9, loc="upper right")
        ax.text(
            0.02,
            0.04,
            f"resid. bridge excess = {stat['residual_bridge_excess_kappa']:.2e}",
            transform=ax.transAxes,
            fontsize=9,
            va="bottom",
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def plot_residual_heatmap(row: dict, output: Path) -> None:
    arr = row["residual_norm"]
    vmax = float(np.nanpercentile(np.abs(arr), 99.5))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(6.2, 5.5), constrained_layout=True)
    im = ax.imshow(
        arr,
        origin="lower",
        extent=(-2.5, 2.5, -2.5, 2.5),
        interpolation="nearest",
        cmap="RdBu_r",
        norm=norm,
        aspect="equal",
    )
    ax.axhline(0.0, color="black", lw=0.7, alpha=0.35)
    ax.axvline(0.0, color="black", lw=0.7, alpha=0.35)
    ax.scatter([-0.5, 0.5], [0.0, 0.0], s=48, marker="x", c="black", linewidths=1.2)
    ax.set_title(f"Residual normalized pair stack, {row['label']}")
    ax.set_xlabel(r"$X/r_\perp$")
    ax.set_ylabel(r"$Y/r_\perp$")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\kappa_{\rm pair} - \kappa_{\rm two halos}$")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def plot_bridge_stats(stats: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)

    x = stats["rperp_center_hmpc"].to_numpy(dtype=float)
    axes[0].plot(x, stats["raw_bridge_excess_kappa"], marker="o", lw=1.8, label="raw")
    axes[0].plot(x, stats["residual_bridge_excess_kappa"], marker="o", lw=1.8, label="residual")
    axes[0].axhline(0.0, color="black", lw=0.8, alpha=0.35)
    axes[0].set_xlabel(r"$r_\perp$ [$h^{-1}$ Mpc]")
    axes[0].set_ylabel(r"bridge excess $\kappa$")
    axes[0].set_title("Bridge Excess")
    axes[0].legend(frameon=False)

    axes[1].plot(x, stats["raw_midpoint_kappa"], marker="o", lw=1.8, label="raw midpoint")
    axes[1].plot(x, stats["residual_midpoint_kappa"], marker="o", lw=1.8, label="residual midpoint")
    axes[1].axhline(0.0, color="black", lw=0.8, alpha=0.35)
    axes[1].set_xlabel(r"$r_\perp$ [$h^{-1}$ Mpc]")
    axes[1].set_ylabel(r"midpoint $\kappa$")
    axes[1].set_title("Midpoint Signal")
    axes[1].legend(frameon=False)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare simulation stacks across pair-separation bins.")
    parser.add_argument("--results-dir", default="analysis/sim/results")
    parser.add_argument("--mass-label", default="mass13")
    parser.add_argument("--format", default="png", choices=["png", "pdf"])
    parser.add_argument("--rperp-labels", nargs="+", default=["rperp5", "rperp10", "rperp20"])
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    results = Path(args.results_dir)

    single = load_map(results / f"kappa_single_sim_{args.mass_label}.csv")
    fixed_axis = np.linspace(-50.0, 50.0, 101)
    norm_axis = np.linspace(-2.5, 2.5, 101)
    norm_x, norm_y = np.meshgrid(norm_axis, norm_axis)

    rows = []
    stats_rows = []
    for label in args.rperp_labels:
        if label not in RPERP_CENTERS:
            raise ValueError(f"Unknown rperp label {label!r}; expected one of {sorted(RPERP_CENTERS)}")

        rperp_center = RPERP_CENTERS[label]
        pair_fixed = load_map(results / f"kappa_pairs_sim_{args.mass_label}_{label}.csv")
        pair_norm = load_map(results / f"kappa_pairs_sim_{args.mass_label}_{label}_normalized.csv")
        template_norm = make_single_template(
            single,
            x_grid=norm_x * rperp_center,
            y_grid=norm_y * rperp_center,
            rperp_center=rperp_center,
        )
        residual_norm = pair_norm - template_norm

        residual_path = results / f"kappa_pairs_sim_{args.mass_label}_{label}_normalized_residual.csv"
        save_residual_csv(residual_path, residual_norm, norm_axis)

        pair_count = count_pairs(results / f"pairs_{args.mass_label}_{label}.csv")
        rows.append(
            {
                "label": label,
                "rperp_center": rperp_center,
                "pair_fixed": pair_fixed,
                "pair_norm": pair_norm,
                "residual_norm": residual_norm,
                "norm_axis": norm_axis,
                "fixed_axis": fixed_axis,
            }
        )
        stats_rows.append(
            summarize_bridge(
                label=label,
                rperp_center=rperp_center,
                pair_count=pair_count,
                pair_norm=pair_norm,
                residual_norm=residual_norm,
                norm_axis=norm_axis,
            )
        )

    stats = pd.DataFrame(stats_rows)
    stats_path = results / f"sim_bridge_stats_{args.mass_label}.csv"
    stats.to_csv(stats_path, index=False)
    logger.info("Saved %s", stats_path)

    plot_comparison(
        rows=rows,
        stats=stats,
        output=results / f"sim_separation_comparison_{args.mass_label}.{args.format}",
        fmt=args.format,
    )
    plot_bridge_stats(
        stats=stats,
        output=results / f"sim_bridge_stats_{args.mass_label}.{args.format}",
    )
    for row in rows:
        plot_residual_heatmap(
            row=row,
            output=results / f"kappa_pairs_sim_{args.mass_label}_{row['label']}_normalized_residual_heatmap.{args.format}",
        )


if __name__ == "__main__":
    main()
