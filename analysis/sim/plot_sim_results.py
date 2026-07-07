#!/usr/bin/env python
"""Plot simulation kappa stack products from CSV outputs."""

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


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_map(path: Path) -> np.ndarray:
    arr = pd.read_csv(path, index_col=0).to_numpy()
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D CSV map at {path}, got shape {arr.shape}")
    if not np.isfinite(arr).all():
        raise ValueError(f"Non-finite values found in {path}")
    return arr


def norm_and_cmap(arr: np.ndarray):
    vmin = float(np.nanmin(arr))
    vmax = float(np.nanmax(arr))
    if vmin < 0.0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax), "RdBu_r"
    return None, "magma"


def plot_map(
    arr: np.ndarray,
    output: Path,
    title: str,
    extent: tuple[float, float, float, float],
    xlabel: str,
    ylabel: str,
    markers: list[tuple[float, float, str]] | None = None,
) -> None:
    norm, cmap = norm_and_cmap(arr)
    fig, ax = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    im = ax.imshow(
        arr,
        origin="lower",
        extent=extent,
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
    )
    ax.axhline(0.0, color="black", lw=0.7, alpha=0.35)
    ax.axvline(0.0, color="black", lw=0.7, alpha=0.35)
    if markers:
        for x, y, label in markers:
            ax.scatter([x], [y], s=55, marker="x", c="black", linewidths=1.4)
            ax.text(x, y, f" {label}", color="black", fontsize=9, va="bottom")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\kappa$")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def plot_overview(
    single: np.ndarray,
    pair: np.ndarray,
    pair_norm: np.ndarray,
    output: Path,
    fixed_halo_offset: float,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True)
    panels = [
        (
            axes[0],
            single,
            "Single halos",
            (-50.0, 50.0, -50.0, 50.0),
            [],
        ),
        (
            axes[1],
            pair,
            "Pairs, fixed separation",
            (-50.0, 50.0, -50.0, 50.0),
            [(-fixed_halo_offset, 0.0), (fixed_halo_offset, 0.0)],
        ),
        (
            axes[2],
            pair_norm,
            "Pairs, normalized separation",
            (-2.5, 2.5, -2.5, 2.5),
            [(-0.5, 0.0), (0.5, 0.0)],
        ),
    ]
    for ax, arr, title, extent, markers in panels:
        norm, cmap = norm_and_cmap(arr)
        im = ax.imshow(
            arr,
            origin="lower",
            extent=extent,
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            aspect="equal",
        )
        ax.axhline(0.0, color="black", lw=0.6, alpha=0.35)
        ax.axvline(0.0, color="black", lw=0.6, alpha=0.35)
        for x, y in markers:
            ax.scatter([x], [y], s=40, marker="x", c="black", linewidths=1.2)
        ax.set_title(title)
        ax.set_xlabel(r"$X$ [$h^{-1}$ Mpc]" if "normalized" not in title else r"$X/r_\perp$")
        ax.set_ylabel(r"$Y$ [$h^{-1}$ Mpc]" if "normalized" not in title else r"$Y/r_\perp$")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"$\kappa$")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    logger.info("Saved %s", output)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot simulation kappa stack heatmaps.")
    parser.add_argument("--results-dir", default="analysis/sim/results")
    parser.add_argument(
        "--mass-label",
        default="mass13",
        help="Mass label used in simulation result filenames.",
    )
    parser.add_argument(
        "--rperp-label",
        default="rperp5",
        choices=["rperp5", "rperp10", "rperp20"],
        help="Pair-separation label used in simulation result filenames.",
    )
    parser.add_argument(
        "--rperp-center",
        type=float,
        default=None,
        help="Central projected pair separation in h^-1 Mpc. Defaults from --rperp-label.",
    )
    parser.add_argument("--format", default="png", choices=["png", "pdf"])
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    results = Path(args.results_dir)
    fmt = args.format
    default_centers = {"rperp5": 5.0, "rperp10": 10.0, "rperp20": 20.0}
    rperp_center = args.rperp_center or default_centers[args.rperp_label]
    halo_offset = 0.5 * rperp_center

    single_path = results / f"kappa_single_sim_{args.mass_label}.csv"
    pair_path = results / f"kappa_pairs_sim_{args.mass_label}_{args.rperp_label}.csv"
    pair_norm_path = results / f"kappa_pairs_sim_{args.mass_label}_{args.rperp_label}_normalized.csv"

    single = load_map(single_path)
    pair = load_map(pair_path)
    pair_norm = load_map(pair_norm_path)

    plot_map(
        single,
        results / f"kappa_single_sim_{args.mass_label}_heatmap.{fmt}",
        f"Single-halo kappa stack, {args.mass_label}",
        (-50.0, 50.0, -50.0, 50.0),
        r"$X$ [$h^{-1}$ Mpc]",
        r"$Y$ [$h^{-1}$ Mpc]",
        markers=[(0.0, 0.0, "halo")],
    )
    plot_map(
        pair,
        results / f"kappa_pairs_sim_{args.mass_label}_{args.rperp_label}_heatmap.{fmt}",
        rf"Pair kappa stack, {args.rperp_label}",
        (-50.0, 50.0, -50.0, 50.0),
        r"$X$ [$h^{-1}$ Mpc]",
        r"$Y$ [$h^{-1}$ Mpc]",
        markers=[(-halo_offset, 0.0, "halo"), (halo_offset, 0.0, "halo")],
    )
    plot_map(
        pair_norm,
        results / f"kappa_pairs_sim_{args.mass_label}_{args.rperp_label}_normalized_heatmap.{fmt}",
        rf"Normalized pair kappa stack, {args.rperp_label}",
        (-2.5, 2.5, -2.5, 2.5),
        r"$X/r_\perp$",
        r"$Y/r_\perp$",
        markers=[(-0.5, 0.0, "halo"), (0.5, 0.0, "halo")],
    )
    plot_overview(
        single,
        pair,
        pair_norm,
        results / f"kappa_sim_{args.mass_label}_{args.rperp_label}_overview.{fmt}",
        fixed_halo_offset=halo_offset,
    )


if __name__ == "__main__":
    main()
