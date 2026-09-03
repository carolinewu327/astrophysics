#!/usr/bin/env python
"""Where does the compressed band fit's amplitude actually come from?

The kappa_band fit returns A = 1.30 / 0.84 / 1.51 at 3.9 / 2.1 / 3.6 sigma,
which reads as a detection until you ask which part of the profile produced it.

Writing u = B s / (s^T B s), the estimator is the linear form A = u . d and
u . s = 1 identically, so the per-bin products u_i s_i partition the fit's
response into fractions summing to one.  Plotting that partition answers the
question directly, and the answer is that the halo bin dominates: the fit is
largely measuring the two galaxies, not the bridge between them.

This is the third estimator to fail this way -- the 2-D pixel fit put 80-98% of
its weight on the perpendicular deficit, the 1-D profile fit put 68-81% on the
halo peaks -- so it is a property of whole-template fitting here, not of any one
compression.

Reads the per-bin CSVs written by fit_band_covariance.py.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts python \\
        analysis/boss/scripts/plot_band_response.py
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
from geometry import band_response_zones
from weighting_sensitivity import SEPARATIONS

logger = logging.getLogger(__name__)

INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e3e2df"
BLUE, ORANGE, NEUTRAL = "#2a78d6", "#eb6834", "#b8b7b3"

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def main(argv=None):
    p = argparse.ArgumentParser(description="Plot the band fit's response split.")
    p.add_argument("--dataset", default="BOSS")
    p.add_argument("--separations", default="5,10,20")
    p.add_argument("--x-bin", type=float, default=4.0)
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--output-dir", default="output/plots")
    args = p.parse_args(argv)
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    keys = [s.strip() for s in args.separations.split(",")]
    fig, axes = plt.subplots(1, len(keys), figsize=(3.4 * len(keys), 3.9),
                             constrained_layout=True, sharey=True)
    axes = np.atleast_1d(axes)

    peak_frac = 0.0
    for ax, key in zip(axes, keys):
        sep = SEPARATIONS[key]["center"]
        df = pd.read_csv(os.path.join(
            args.results_dir, f"band_fit_response_sep{key}_{args.dataset}.csv"))
        x = df.x_hmpc.to_numpy()
        frac = df.response_fraction.to_numpy() * 100.0
        bridge, halo, outer = band_response_zones(x, sep, args.x_bin)
        peak_frac = max(peak_frac, float(frac.max()))

        # Colour carries the zone, but every bar is also labelled by position on
        # the X axis and the two named zones are annotated directly, so identity
        # never rests on hue alone.
        colors = np.where(bridge, BLUE, np.where(halo, ORANGE, NEUTRAL))
        ax.bar(x, frac, width=args.x_bin * 0.85, color=list(colors),
               edgecolor="white", linewidth=0.6, zorder=3)
        ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.7, zorder=2)

        for mask, color, name in ((bridge, BLUE, "bridge"), (halo, ORANGE, "halo")):
            if not mask.any():
                continue
            total = frac[mask].sum()
            peak = x[mask][np.argmax(np.abs(frac[mask]))]
            ax.annotate(f"{name}\n{total:+.0f}%", (peak, frac[mask].max()),
                        xytext=(0, 6), textcoords="offset points", ha="center",
                        color=color, fontsize=8, fontweight="bold", zorder=5)

        if not halo.any():
            ax.annotate("bridge and halo\nfall in one bin —\nnot separable here",
                        (0.97, 0.55), xycoords="axes fraction", ha="right",
                        va="top", color=INK_2, fontsize=7.5, style="italic")

        ax.set_title(rf"separation {sep:g} $h^{{-1}}$ Mpc", loc="left",
                     color=INK, pad=5)
        ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
        ax.grid(True, axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("share of the fitted amplitude  [%]")
    # Headroom for the zone labels, which sit just above the tallest bar and
    # would otherwise ride up into the panel titles.
    axes[0].set_ylim(top=peak_frac * 1.32)

    fig.suptitle(
        "The band fit is driven by the halo peaks, not the bridge\n"
        "Each bar is one $X$ bin's share of $\\hat{A}$; the shares sum to 100%. "
        "Bars can exceed 100% or go negative where zones partly cancel.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    base = os.path.join(args.output_dir, "band_fit_response")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s.png / .pdf", base)


if __name__ == "__main__":
    main()
