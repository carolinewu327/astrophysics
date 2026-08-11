#!/usr/bin/env python
"""Does loosening the Planck smoothing improve filament S/N?

8 arcmin was chosen to suppress Planck's lensing reconstruction noise.  Loosening
it keeps more small-scale signal but admits more of that noise, so the net effect
on S/N has to be measured rather than argued.

Layout: one panel per separation showing the filament profile at each available
smoothing, plus a summary panel of bridge excess and S/N.

Smoothings are discovered from the accumulators on disk, so adding a 4 arcmin run
requires no change here -- just stack with ``--fwhm-arcmin 4 --label fwhm4`` (and
``galaxy_<sep>_fwhm4`` for the pairs) and re-run.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/plot_smoothing_scan.py
"""

from __future__ import annotations

import argparse
import logging
import os

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from catalog import setup_logging
from jackknife import jackknife_error
from plot_separation_summary import BRIDGE_X_FRAC, band_profile, build_terms
from weighting_sensitivity import SEPARATIONS

logger = logging.getLogger(__name__)

# Categorical slots in fixed order; coarsest smoothing first.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e3e2df"

# (label, single-stack tag, pair-accumulator suffix), coarsest first.
CANDIDATES = [
    ("8 arcmin", "_scw", ""),
    ("6 arcmin", "_scw_fwhm6", "_fwhm6"),
    ("4 arcmin", "_scw_fwhm4", "_fwhm4"),
    ("2 arcmin", "_scw_fwhm2", "_fwhm2"),
]

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def available(results_dir: str, dataset: str, region_label: str) -> list[tuple[str, str, str]]:
    jk = os.path.join(results_dir, "jk")
    found = []
    for label, tag, suffix in CANDIDATES:
        gal = os.path.join(jk, f"acc_single_galaxy{tag}_{dataset}_{region_label}.npz")
        if not os.path.exists(gal):
            logger.info("Skipping %s (no single-galaxy accumulator)", label)
            continue
        # A smoothing is included as soon as its single stack exists; separations
        # whose pair accumulator is missing are simply absent from that panel,
        # so a partial scan still plots what it has.
        have = [k for k in SEPARATIONS if os.path.exists(os.path.join(
            jk, f"acc_pairs_galaxy_{SEPARATIONS[k]['galaxy_label']}{suffix}_"
                f"{dataset}_{region_label}.npz"))]
        if not have:
            logger.info("Skipping %s (no pair accumulators)", label)
            continue
        if len(have) < len(SEPARATIONS):
            logger.warning("%s has only %s of %d separations: %s", label, len(have),
                           len(SEPARATIONS), ", ".join(have))
        found.append((label, tag, suffix))
    return found


def bridge_stat(t: dict) -> tuple[float, float, np.ndarray, np.ndarray]:
    axis, sep = t["axis"], t["sep"]
    window = np.abs(axis) <= BRIDGE_X_FRAC * sep
    prof = band_profile(t["filament"], axis)
    loo = np.array([band_profile(m.reshape(t["filament"].shape), axis)
                    for m in t["loo"]["filament"]])
    _, prof_err = jackknife_error(loo)
    _, bridge_err = jackknife_error(loo[:, window].mean(axis=1)[:, None])
    return float(prof[window].mean()), float(bridge_err[0]), prof, prof_err


def main(argv=None):
    parser = argparse.ArgumentParser(description="Filament S/N versus Planck smoothing.")
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument("--single-random-tag", default="_scw_frac100")
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--output-dir", default="output/plots")
    args = parser.parse_args(argv)
    regions = [r.strip() for r in args.regions.split(",")]
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    smoothings = available(args.results_dir, args.dataset, "_".join(regions))
    if not smoothings:
        raise SystemExit("No smoothing runs found.")
    logger.info("Comparing: %s", ", ".join(s[0] for s in smoothings))

    keys = list(SEPARATIONS)
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.6), constrained_layout=True)
    panels = [axes[0, 0], axes[0, 1], axes[1, 0]]
    rows = []

    for ax, key in zip(panels, keys):
        sep = SEPARATIONS[key]["center"]
        for i, (label, tag, suffix) in enumerate(smoothings):
            acc = os.path.join(
                args.results_dir, "jk",
                f"acc_pairs_galaxy_{SEPARATIONS[key]['galaxy_label']}{suffix}_"
                f"{args.dataset}_{'_'.join(regions)}.npz")
            if not os.path.exists(acc):
                continue
            t = build_terms(args.results_dir, args.dataset, regions, key, tag,
                            args.single_random_tag, pair_suffix=suffix)
            value, err, prof, prof_err = bridge_stat(t)
            rows.append({"smoothing": label, "separation_hmpc": sep,
                         "filament_bridge_excess": value, "filament_bridge_err": err,
                         "snr": value / err if err else np.nan,
                         "n_pairs": t["n_pairs"]})
            color = PALETTE[i % len(PALETTE)]
            axis = t["axis"]
            ax.plot(axis, prof * 1e4, lw=2, color=color, label=label, zorder=3 + i)
            ax.fill_between(axis, (prof - prof_err) * 1e4, (prof + prof_err) * 1e4,
                            color=color, alpha=0.13, lw=0, zorder=1)
        ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6)
        for xx in (-sep / 2, sep / 2):
            ax.axvline(xx, color=INK_MUTED, lw=0.8, ls="--", alpha=0.6)
        ax.axvspan(-BRIDGE_X_FRAC * sep, BRIDGE_X_FRAC * sep, color=INK_MUTED,
                   alpha=0.07, lw=0, zorder=0)
        ax.set_xlim(-2.5 * sep, 2.5 * sep)
        ax.set_title(f"Filament profile, $r_\\perp$ = {sep:g} $h^{{-1}}$ Mpc",
                     loc="left", color=INK, pad=5)
        ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$\kappa$  [$10^{-4}$]")
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if ax is panels[0]:
            ax.legend(loc="upper right")

    out = pd.DataFrame(rows)
    path = os.path.join(args.results_dir, f"smoothing_scan_{args.dataset}.csv")
    out.to_csv(path, index=False)

    ax = axes[1, 1]
    width = 0.8 / len(smoothings)
    x = np.arange(len(keys))
    for i, (label, _, _) in enumerate(smoothings):
        sub = out[out.smoothing == label].set_index("separation_hmpc")
        vals = [sub.loc[SEPARATIONS[k]["center"], "snr"]
                if SEPARATIONS[k]["center"] in sub.index else np.nan for k in keys]
        ax.bar(x + (i - (len(smoothings) - 1) / 2) * width, vals, width * 0.9,
               color=PALETTE[i % len(PALETTE)], label=label, zorder=3)
    ax.axhline(0, color=INK_MUTED, lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{SEPARATIONS[k]['center']:g}" for k in keys])
    ax.set_xlabel(r"pair separation  [$h^{-1}$ Mpc]")
    ax.set_ylabel("filament S/N  ($\\sigma$)")
    ax.set_title("Bridge-excess significance by smoothing", loc="left", color=INK, pad=5)
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.suptitle(
        "Does loosening the Planck smoothing help? — BOSS filament, joint N+S\n"
        "Shaded bands are jackknife errors over 287 cells; grey window is the bridge "
        "region averaged for the S/N panel.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    base = os.path.join(args.output_dir, "filament_snr_vs_smoothing")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s.png / .pdf and %s", base, path)

    logger.info("\nFilament bridge excess (kappa x 1e4):")
    logger.info("  sep   " + "".join(f"{s[0]:>24}" for s in smoothings))
    for k in keys:
        sep = SEPARATIONS[k]["center"]
        cells = []
        for label, _, _ in smoothings:
            m = out[(out.smoothing == label) & (out.separation_hmpc == sep)]
            if m.empty:
                cells.append("not run")
                continue
            r = m.iloc[0]
            cells.append(f"{r.filament_bridge_excess * 1e4:8.2f} ± {r.filament_bridge_err * 1e4:5.2f}"
                         f" ({r.snr:5.2f}s)")
        logger.info("  %4.0f  %s", sep, "".join(f"{c:>24}" for c in cells))


if __name__ == "__main__":
    main()
