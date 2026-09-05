#!/usr/bin/env python
"""Check BOSS 5/10/20 products for cut consistency and simple error bars.
.. warning::

   This script scores with **separation-scaled** Y bands
   (``|Y| <= 0.15 r_perp``, off-centre ``0.45-0.85 r_perp``), not the fixed
   physical bands (``|Y| <= 1.5``, ``1.5-10.5 h^-1 Mpc``) that
   ``lib/geometry.py`` defines and that the BOSS estimator and the paper use.
   The two differ by 15-44 per cent and the ratio itself varies with
   separation, so numbers from here are **not** comparable with BOSS numbers
   and must not be quoted beside them.  Kept as a legacy sensitivity axis --
   scaled bands answer a different question.  Outputs carry a ``_scaledband``
   suffix and a ``band_convention`` column so they cannot be mistaken for the
   production statistic.  See "Scoring the mock with the fixed physical bands"
   in ``paper/README.md``.
"""

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

from sim_utils import radial_symmetrize_map


logger = logging.getLogger(__name__)


BOX_SIZE_HMPC = 100.0
RPERP_CONFIG = {
    "rperp5": {
        "center": 5.0,
        "pair_tag": "5",
        "suffix": "5_frac100",
        "uses_rpar10": False,
        "note": "5 h^-1 Mpc, full-random product",
    },
    "rperp10": {
        "center": 10.0,
        "pair_tag": "10",
        "suffix": "10_frac100",
        "uses_rpar10": False,
        "note": "10 h^-1 Mpc, full-random product",
    },
    "rperp20": {
        "center": 20.0,
        "pair_tag": "20_rpar10",
        "suffix": "20_rpar10_frac100",
        "uses_rpar10": True,
        "note": "20 h^-1 Mpc, rpar10 full-random product",
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


def axis_for_map(arr: np.ndarray) -> np.ndarray:
    n = arr.shape[0]
    if n % 2 == 1:
        return np.linspace(-0.5 * BOX_SIZE_HMPC, 0.5 * BOX_SIZE_HMPC, n)
    cell = BOX_SIZE_HMPC / n
    return np.linspace(-0.5 * BOX_SIZE_HMPC + 0.5 * cell, 0.5 * BOX_SIZE_HMPC - 0.5 * cell, n)


def reconcile_shapes(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if a.shape == b.shape:
        return a, b
    rows = min(a.shape[0], b.shape[0])
    cols = min(a.shape[1], b.shape[1])

    def trim(arr: np.ndarray) -> np.ndarray:
        row0 = (arr.shape[0] - rows) // 2
        col0 = (arr.shape[1] - cols) // 2
        return arr[row0 : row0 + rows, col0 : col0 + cols]

    return trim(a), trim(b)


def build_control_pair_map(single_map: np.ndarray, separation_hmpc: float) -> np.ndarray:
    grid_size = single_map.shape[1]
    cell_size = BOX_SIZE_HMPC / grid_size
    shift_pixels = separation_hmpc / 2.0 / cell_size

    if np.isclose(shift_pixels, np.round(shift_pixels)):
        shift_int = int(np.round(shift_pixels))
        right = np.roll(single_map, shift_int, axis=1)
        left = np.roll(single_map, -shift_int, axis=1)
    else:
        from scipy.ndimage import shift as ndimage_shift

        right = ndimage_shift(single_map, [0, shift_pixels], order=3, mode="wrap")
        left = ndimage_shift(single_map, [0, -shift_pixels], order=3, mode="wrap")
    return right + left


def bridge_masks(axis: np.ndarray, rperp_center: float) -> tuple[np.ndarray, np.ndarray]:
    """LEGACY separation-scaled bands -- see the module warning."""
    x_grid, y_grid = np.meshgrid(axis, axis)
    bridge = (np.abs(x_grid) <= 0.35 * rperp_center) & (np.abs(y_grid) <= 0.15 * rperp_center)
    side = (
        (np.abs(x_grid) <= 0.35 * rperp_center)
        & (np.abs(y_grid) >= 0.45 * rperp_center)
        & (np.abs(y_grid) <= 0.85 * rperp_center)
    )
    return bridge, side


def bridge_stats(arr: np.ndarray, rperp_center: float, prefix: str) -> dict[str, float]:
    axis = axis_for_map(arr)
    bridge, side = bridge_masks(axis, rperp_center)
    center = len(axis) // 2
    bridge_mean = float(arr[bridge].mean())
    side_mean = float(arr[side].mean())
    return {
        f"{prefix}_midpoint_kappa": float(arr[center, center]),
        f"{prefix}_bridge_mean_kappa": bridge_mean,
        f"{prefix}_sideband_mean_kappa": side_mean,
        f"{prefix}_bridge_excess_kappa": bridge_mean - side_mean,
    }


def ns_error(values: list[float]) -> tuple[float, float]:
    if len(values) < 2:
        return np.nan, np.nan
    arr = np.asarray(values, dtype=float)
    return float(arr.std(ddof=1) / np.sqrt(len(arr))), float(abs(arr[0] - arr[1]) / 2.0)


def region_maps(results_dir: Path, cfg: dict, region: str) -> tuple[np.ndarray, np.ndarray]:
    single_gal = load_map(results_dir / f"kappa_single_galaxy_BOSS_{region}.csv")
    single_random = load_map(results_dir / f"kappa_single_random_frac100_BOSS_{region}.csv")
    pair_gal = load_map(results_dir / f"kappa_pairs_galaxy_{cfg['pair_tag']}_BOSS_{region}.csv")
    pair_random = load_map(results_dir / f"kappa_pairs_random_{cfg['suffix']}_BOSS_{region}.csv")

    corrected_single = radial_symmetrize_map(single_gal - single_random)
    control = build_control_pair_map(corrected_single, float(cfg["center"]))
    control = 0.5 * (control + control[:, ::-1])
    corrected_pair = pair_gal - pair_random
    corrected_pair, control = reconcile_shapes(corrected_pair, control)
    filament = corrected_pair - control
    return corrected_pair, filament


def summarize_boss(results_dir: Path, labels: list[str]) -> pd.DataFrame:
    rows = []
    for label in labels:
        cfg = RPERP_CONFIG[label]
        suffix = cfg["suffix"]
        center = float(cfg["center"])

        corrected_pair = load_map(results_dir / f"kappa_corrected_pairs_{suffix}_BOSS.csv")
        filament = load_map(results_dir / f"kappa_filament_{suffix}_BOSS.csv")
        row = {
            "rperp_label": label,
            "rperp_center_hmpc": center,
            "boss_suffix": suffix,
            "boss_pair_tag": cfg["pair_tag"],
            "uses_rpar10": bool(cfg["uses_rpar10"]),
            "cut_consistency_note": cfg["note"],
            "consistent_with_5_10_full_random_cut": not bool(cfg["uses_rpar10"]),
        }
        row.update(bridge_stats(corrected_pair, center, "combined_corrected_pair"))
        row.update(bridge_stats(filament, center, "combined_filament"))

        region_values = {"corrected_pair": [], "filament": []}
        for region in ("North", "South"):
            try:
                region_pair, region_filament = region_maps(results_dir, cfg, region)
            except FileNotFoundError as exc:
                logger.warning("Missing region file for %s/%s: %s", label, region, exc)
                continue
            pair_stat = bridge_stats(region_pair, center, f"{region.lower()}_corrected_pair")
            filament_stat = bridge_stats(region_filament, center, f"{region.lower()}_filament")
            row.update(pair_stat)
            row.update(filament_stat)
            region_values["corrected_pair"].append(pair_stat[f"{region.lower()}_corrected_pair_bridge_excess_kappa"])
            region_values["filament"].append(filament_stat[f"{region.lower()}_filament_bridge_excess_kappa"])

        pair_sem, pair_halfdiff = ns_error(region_values["corrected_pair"])
        filament_sem, filament_halfdiff = ns_error(region_values["filament"])
        row["ns_region_count"] = len(region_values["filament"])
        row["corrected_pair_ns_sem_kappa"] = pair_sem
        row["corrected_pair_ns_halfdiff_kappa"] = pair_halfdiff
        row["filament_ns_sem_kappa"] = filament_sem
        row["filament_ns_halfdiff_kappa"] = filament_halfdiff
        rows.append(row)
    return pd.DataFrame(rows)


def plot_boss_summary(stats: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), constrained_layout=True)
    x = stats["rperp_center_hmpc"].to_numpy(dtype=float)
    labels = stats["rperp_label"].str.replace("rperp", "", regex=False).to_list()

    axes[0].errorbar(
        x,
        stats["combined_filament_bridge_excess_kappa"].to_numpy(dtype=float),
        yerr=stats["filament_ns_sem_kappa"].to_numpy(dtype=float),
        marker="o",
        lw=1.8,
        capsize=3,
        label="filament",
    )
    axes[1].errorbar(
        x,
        stats["combined_corrected_pair_bridge_excess_kappa"].to_numpy(dtype=float),
        yerr=stats["corrected_pair_ns_sem_kappa"].to_numpy(dtype=float),
        marker="o",
        lw=1.8,
        capsize=3,
        label="corrected pair",
    )

    for ax, title in zip(axes, ["BOSS Filament", "BOSS Corrected Pair"]):
        ax.axhline(0.0, color="black", lw=0.8, alpha=0.35)
        ax.set_xticks(x, labels)
        ax.set_xlabel(r"$r_\perp$ [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"bridge excess $\kappa$")
        ax.set_title(title)
        ax.legend(frameon=False)
    axes[0].text(
        0.02,
        0.02,
        "Error bars: North/South split only; 20 uses rpar10 full-random product.",
        transform=axes[0].transAxes,
        fontsize=8,
        va="bottom",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check BOSS product consistency and simple error bars.")
    parser.add_argument("--results-dir", type=Path, default=Path("analysis/boss/results"))
    parser.add_argument("--output", type=Path, default=Path("analysis/sim/results/observed_comparison/boss_consistency_stats_scaledband.csv"))
    parser.add_argument("--plot", type=Path, default=Path("analysis/sim/results/observed_comparison/boss_consistency_stats_scaledband.png"))
    parser.add_argument("--rperp-labels", nargs="+", default=["rperp5", "rperp10", "rperp20"])
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    stats = summarize_boss(args.results_dir, args.rperp_labels)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stats["band_convention"] = "scaled_0.15/0.45/0.85_x_rperp"
    stats.to_csv(args.output, index=False)
    logger.info("Saved %s", args.output)
    plot_boss_summary(stats, args.plot)


if __name__ == "__main__":
    main()
