#!/usr/bin/env python
"""Diagnostic: how much of the profile_c / panel-(d) divergence is the method?

The archived ``profile_c_filament_*.pdf`` curves and panel (d) of
``boss_vs_sim_filament_sep*.png`` are both "the filament profile", but they
disagree badly -- correlation 0.65 on the same map, with sign flips.  Two things
could explain it: the two figures were made three months apart across the
half-pixel template fix, or the two collapsing methods simply measure different
quantities.

This isolates the second by applying *both* collapsing methods to the *same*
current maps:

  profile_c method   sum the 5 rows with |Y| <= 2.5, no off-band subtraction
  panel (d) method   mean of 3 central rows minus mean of 18 off-centre rows

If the two panels disagree, the method alone accounts for it and the vintage
difference is a second, separate effect.

Diagnostic only -- the panel-(d) statistic is the production one, because a
filament is an excess *along* the axis relative to just off it, and only the
difference cancels a map-wide offset.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts python \\
        analysis/boss/scripts/plot_collapse_comparison.py --separation 20
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
from geometry import BRIDGE_HALF_X_FRAC, band_profile
from jackknife import jackknife_error
from plot_obs_vs_sim_filament import sim_filament
from plot_separation_summary import build_terms
from weighting_sensitivity import SEPARATIONS

logger = logging.getLogger(__name__)

INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e3e2df"
BLUE, ORANGE = "#2a78d6", "#eb6834"

# plot_results.PROFILE_HALF_WIDTH, the band profile_c collapses over.
PROFILE_HALF_WIDTH = 2.5

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def profile_c_collapse(arr: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """plot_results.extract_profile: sum rows within +/- PROFILE_HALF_WIDTH of Y=0."""
    return arr[np.abs(axis) <= PROFILE_HALF_WIDTH, :].sum(axis=0)


def main(argv=None):
    p = argparse.ArgumentParser(description="Compare the two profile collapsing methods.")
    p.add_argument("--dataset", default="BOSS")
    p.add_argument("--regions", default="North,South")
    p.add_argument("--separation", default="20")
    p.add_argument("--single-tag", default="_scw")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin.csv")
    p.add_argument("--output-dir", default="output/plots")
    args = p.parse_args(argv)
    regions = [r.strip() for r in args.regions.split(",")]
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    key = args.separation
    sep = SEPARATIONS[key]["center"]
    t = build_terms(args.results_dir, args.dataset, regions, key,
                    args.single_tag, args.single_random_tag)
    obs, axis, loo = t["filament"], t["axis"], t["loo"]["filament"]
    shape = obs.shape
    sim, _ = sim_filament(args.sim_dir, key, args.sim_single)

    methods = {
        "profile_c": (profile_c_collapse,
                      rf"(a) profile_c method — sum of rows $|Y| \leq {PROFILE_HALF_WIDTH:g}$",
                      r"summed $\kappa$  [$10^{-4}$]"),
        "panel_d": (band_profile,
                    "(b) panel (d) method — central band − off-centre band",
                    r"$\kappa$  [$10^{-4}$]"),
    }

    zoom = 2.5 * sep
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), constrained_layout=True,
                             sharex=True)
    out = {"X_hmpc": axis}

    for ax, (name, (fn, title, ylabel)) in zip(axes, methods.items()):
        o = fn(obs, axis) * 1e4
        s = fn(sim, axis) * 1e4
        _, err = jackknife_error(np.array([fn(m.reshape(shape), axis) for m in loo]))
        err = err * 1e4

        ax.plot(axis, o, lw=2, color=BLUE, label="BOSS", zorder=3)
        ax.fill_between(axis, o - err, o + err, color=BLUE, alpha=0.15, lw=0, zorder=1)
        ax.plot(axis, s, lw=2, color=ORANGE, label="Simulation", zorder=4)

        ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6)
        for xx in (-sep / 2, sep / 2):
            ax.axvline(xx, color=INK_MUTED, lw=0.8, ls="--", alpha=0.7)
        ax.axvspan(-BRIDGE_HALF_X_FRAC * sep, BRIDGE_HALF_X_FRAC * sep,
                   color=INK_MUTED, alpha=0.07, lw=0, zorder=0)
        ax.set_xlim(-zoom, zoom)
        ax.set_title(title, loc="left", color=INK, pad=5)
        ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(ylabel)
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.legend(loc="upper right")

        out[f"boss_{name}"] = fn(obs, axis)
        out[f"boss_{name}_err"] = err / 1e4
        out[f"sim_{name}"] = fn(sim, axis)

    inner = np.abs(axis) <= 30
    r = np.corrcoef(out["boss_profile_c"][inner], out["boss_panel_d"][inner])[0, 1]

    fig.suptitle(
        f"Same filament map, two ways of collapsing it — separation {sep:g} "
        f"$h^{{-1}}$ Mpc\n"
        f"Identical input maps and vintage, so every difference here is the method. "
        f"BOSS curves correlate at r = {r:.2f} over $|X| \\leq 30$.  Dashed lines "
        f"mark the galaxies; shading is the bridge window.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    base = os.path.join(args.output_dir, f"collapse_comparison_sep{key}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight")
    plt.close(fig)

    csv = os.path.join(args.results_dir, f"collapse_comparison_sep{key}_{args.dataset}.csv")
    pd.DataFrame(out).to_csv(csv, index=False)

    logger.info("BOSS curves correlate at r = %.4f over |X| <= 30", r)
    logger.info("   X    profile_c(sum)   panel_d(diff)")
    for xv in (0, 5, 10, 15, 20, 25):
        i = int(np.argmin(np.abs(axis - xv)))
        logger.info("  %4.0f   %12.2f   %12.2f", axis[i],
                    out["boss_profile_c"][i] * 1e4, out["boss_panel_d"][i] * 1e4)
    logger.info("Saved %s.png / .pdf and %s", base, csv)


if __name__ == "__main__":
    main()
