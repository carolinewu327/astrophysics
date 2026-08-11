#!/usr/bin/env python
"""Summarize completed simulation sensitivity variants from a manifest."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sim_utils import make_two_halo_template, radial_symmetrize_map


logger = logging.getLogger(__name__)


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


def count_pairs(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        return max(sum(1 for _ in f) - 1, 0)


def make_single_template(
    single: np.ndarray,
    norm_axis: np.ndarray,
    rperp_center: float,
) -> np.ndarray:
    # radial-profile construction: exact for the symmetrized single map,
    # symmetric by construction, and covers radii out to the map corners
    x_grid, y_grid = np.meshgrid(norm_axis * rperp_center, norm_axis * rperp_center)
    return make_two_halo_template(single, x_grid, y_grid, rperp_center)


def bridge_masks(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_grid, y_grid = np.meshgrid(axis, axis)
    bridge = (np.abs(x_grid) <= 0.35) & (np.abs(y_grid) <= 0.15)
    side = (np.abs(x_grid) <= 0.35) & (np.abs(y_grid) >= 0.45) & (np.abs(y_grid) <= 0.85)
    return bridge, side


def map_stats(arr: np.ndarray, axis: np.ndarray, prefix: str) -> dict[str, float]:
    bridge, side = bridge_masks(axis)
    center = len(axis) // 2
    bridge_mean = float(arr[bridge].mean())
    side_mean = float(arr[side].mean())
    return {
        f"{prefix}_midpoint_kappa": float(arr[center, center]),
        f"{prefix}_bridge_mean_kappa": bridge_mean,
        f"{prefix}_sideband_mean_kappa": side_mean,
        f"{prefix}_bridge_excess_kappa": bridge_mean - side_mean,
    }


def summarize_variant(variant: dict) -> list[dict]:
    single_path = Path(variant["single_path"])
    if not single_path.exists():
        logger.warning("Skipping %s: missing %s", variant["name"], single_path)
        return []

    single = radial_symmetrize_map(load_map(single_path))
    rows = []
    for rperp_label, pair_cfg in variant["pairs"].items():
        normalized_path = Path(pair_cfg["normalized_path"])
        if not normalized_path.exists():
            logger.warning("Skipping %s/%s: missing %s", variant["name"], rperp_label, normalized_path)
            continue

        pair_norm = load_map(normalized_path)
        norm_axis = np.linspace(-2.5, 2.5, pair_norm.shape[0])
        rperp_center = float(pair_cfg["rperp_center_hmpc"])
        template = make_single_template(single, norm_axis, rperp_center)
        residual = pair_norm - template

        row = {
            "variant": variant["name"],
            "kind": variant["kind"],
            "rperp_label": rperp_label,
            "rperp_center_hmpc": rperp_center,
            "pair_count": count_pairs(Path(pair_cfg["pair_path"])),
        }
        for key, value in variant.get("parameters", {}).items():
            row[f"param_{key}"] = value
        row.update(map_stats(pair_norm, norm_axis, "raw"))
        row.update(map_stats(residual, norm_axis, "residual"))
        rows.append(row)
    return rows


def plot_summary(stats: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), constrained_layout=True)
    for kind, kind_df in stats.groupby("kind", sort=False):
        for variant, df in kind_df.groupby("variant", sort=False):
            x = df.sort_values("rperp_center_hmpc")["rperp_center_hmpc"].to_numpy(dtype=float)
            residual = df.sort_values("rperp_center_hmpc")["residual_bridge_excess_kappa"].to_numpy(dtype=float)
            raw = df.sort_values("rperp_center_hmpc")["raw_bridge_excess_kappa"].to_numpy(dtype=float)
            label = f"{kind}: {variant}" if kind != "baseline" else "baseline"
            axes[0].plot(x, residual, marker="o", lw=1.4, label=label)
            axes[1].plot(x, raw, marker="o", lw=1.4, label=label)

    for ax, title in zip(axes, ["Residual Bridge Excess", "Raw Pair Bridge Excess"]):
        ax.axhline(0.0, color="black", lw=0.8, alpha=0.35)
        ax.set_xlabel(r"$r_\perp$ [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"bridge excess $\kappa$")
        ax.set_title(title)
    axes[1].legend(frameon=False, fontsize=8, loc="best")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize completed simulation sensitivity variants.")
    parser.add_argument("--manifest", type=Path, default=Path("analysis/sim/results/sensitivity/manifest.json"))
    parser.add_argument("--output", type=Path, default=Path("analysis/sim/results/sensitivity/sim_sensitivity_summary.csv"))
    parser.add_argument("--plot", type=Path, default=Path("analysis/sim/results/sensitivity/sim_sensitivity_summary.png"))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))

    rows = []
    for variant in manifest["variants"]:
        rows.extend(summarize_variant(variant))
    if not rows:
        raise RuntimeError("No completed sensitivity products found to summarize")

    stats = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(args.output, index=False)
    logger.info("Saved %s", args.output)
    plot_summary(stats, args.plot)


if __name__ == "__main__":
    main()
