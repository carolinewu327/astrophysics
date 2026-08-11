#!/usr/bin/env python
"""Filament signal versus line-of-sight cut, from rpar_scan.py.

Two panels, deliberately not one panel with two y-axes: the signal amplitude and
the predicted-S/N figure of merit are different quantities on different scales.

  (a) simulation bridge signal S(C) with jackknife bands
  (b) predicted BOSS figure of merit S(C) * sqrt(N(C)), normalised per r_perp to
      the current production cut

The headline is the difference in *shape* between the two r_perp values: at
r_perp = 5 the signal rises before it falls, at r_perp = 10 it only falls.

Usage
-----
    PYTHONPATH=lib python analysis/sim/plot_rpar_scan.py
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

logger = logging.getLogger(__name__)

# Repo house style, matching the BOSS summary sheets so the figures sit together.
# Blue/orange is the canonical CVD-safe two-hue pair; distinct markers and direct
# labels carry identity as well, so nothing depends on colour alone.
SERIES = {5.0: ("#2a78d6", "o", r"$r_\perp$ = 5"),
          10.0: ("#eb6834", "s", r"$r_\perp$ = 10")}
PRODUCTION_CUT = {5.0: 5.0, 10.0: 10.0}
INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e3e2df"

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def tidy(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main(argv=None):
    p = argparse.ArgumentParser(description="Plot the LOS-cut scan.")
    p.add_argument("--scan", default="analysis/sim/results/rpar_scan.csv")
    p.add_argument("--output-dir", default="output/plots")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.scan)
    fig, (ax_s, ax_f) = plt.subplots(1, 2, figsize=(9.6, 4.0), constrained_layout=True)

    for rperp, grp in df.groupby("rperp_center_hmpc"):
        color, marker, label = SERIES[rperp]
        cuts = grp.rpar_cut_hmpc.to_numpy()
        sig, err = grp.bridge_excess.to_numpy() * 1e4, grp.bridge_excess_err.to_numpy() * 1e4

        ax_s.plot(cuts, sig, lw=2, color=color, marker=marker, ms=5, label=label, zorder=3)
        ax_s.fill_between(cuts, sig - err, sig + err, color=color, alpha=0.15, lw=0, zorder=1)
        ax_f.plot(cuts, grp.fom.to_numpy(), lw=2, color=color, marker=marker, ms=5,
                  label=label, zorder=3)

        # Direct label at the right edge: identity without relying on colour.
        ax_s.annotate(label, (cuts[-1], sig[-1]), xytext=(4, 0),
                      textcoords="offset points", color=color, fontsize=8,
                      va="center", ha="left")

        for ax in (ax_s, ax_f):
            ax.axvline(PRODUCTION_CUT[rperp], color=color, lw=1.1, ls=":", alpha=0.8, zorder=0)

    ax_s.set_xlabel(r"line-of-sight cut  $r_{\parallel,\mathrm{rsd}} \leq C$  [$h^{-1}$ Mpc]")
    ax_s.set_ylabel(r"bridge signal,  $\kappa$  [$10^{-4}$]")
    ax_s.set_title("(a) Simulation signal versus LOS cut", loc="left", color=INK, pad=5)
    ax_s.set_xlim(0, 56)
    tidy(ax_s)

    ax_f.axhline(1.0, color=INK_MUTED, lw=0.8, alpha=0.6)
    ax_f.set_xlabel(r"line-of-sight cut  $r_{\parallel,\mathrm{rsd}} \leq C$  [$h^{-1}$ Mpc]")
    ax_f.set_ylabel(r"$S(C)\,\sqrt{N(C)}$  /  value at current cut")
    ax_f.set_title("(b) Predicted BOSS figure of merit", loc="left", color=INK, pad=5)
    ax_f.set_xlim(0, 56)
    # Centre-right is the one empty band in this panel: the blue curve sits above
    # 1.5 and the orange below 1.15 across the whole right-hand side.
    ax_f.legend(loc="center right")
    tidy(ax_f)

    fig.suptitle(
        "Does a wider line-of-sight cut help? — density-matched mock, 8 arcmin-equivalent smoothing\n"
        "Dotted verticals mark the cut currently used at each $r_\\perp$.  Panel (b) assumes added "
        "pairs are independent, which overstates the gain at wide cuts.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    base = os.path.join(args.output_dir, "rpar_scan")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s.png / .pdf", base)


if __name__ == "__main__":
    main()
