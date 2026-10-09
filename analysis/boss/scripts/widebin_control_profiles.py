#!/usr/bin/env python
"""What the control looks like, next to the pair stack it is subtracted from.

The band-profile figures only ever show the residual, pair stack minus control.
This plots the two terms themselves along the pair axis, in the central band
(|Y| <= 1.5 h^-1 Mpc, the strip through both galaxies), for

    archived BOSS   4-6, 9-11, 18-22; nominal single-separation control
    wide BOSS       4-6, 6-15, 15-25; control averaged over the pairs' actual
                    separations (the nominal control dashed, for reference)
    wide mock       same bins; the mock pair stack and its own averaged control
                    (the nominal control dashed, as for BOSS)

In the mock the averaged control matches the central band to ~1e-4 and the
single-R control misses the halo peak by up to ~1e-3, so the dashed curve is a
reference for what averaging changes, not an alternative model.

Writes widebin_rpar5/comparison/control_profiles_archived_vs_wide.{png,pdf}.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python \\
        analysis/boss/scripts/widebin_control_profiles.py
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

from catalog import setup_logging
from combine_filament_jackknife import (BRIDGE_CENTERS, load_grouping_pairs, load_random_subbins,
                                        loo_and_total, pooled_random_map)
from combine_jackknife import load_accumulator
from geometry import (BRIDGE_HALF_X_FRAC, band_profiles, reflect_symmetrize_map,
                      separation_averaged_template, separation_mix, symmetrize_map, two_halo_template)
from jackknife import jackknife_error
from plot_separation_summary import build_terms
from weighting_sensitivity import load_map

logger = logging.getLogger(__name__)

ROWS = [("5", "g4_6", "4-6", "4-6"), ("10", "g6_15", "9-11", "6-15"), ("20", "g15_25", "18-22", "15-25")]
BLUE, ORANGE, INK, MUTED, GREEN = "#2a78d6", "#eb6834", "#0b0b0b", "#8a8985", "#3a9b5c"
AX = np.linspace(-50, 50, 101)
X, Y = np.meshgrid(AX, AX)


def central(arr):
    return band_profiles(arr, AX)[0]


def central_err(loo_flat):
    return jackknife_error(np.array([central(m.reshape(101, 101)) for m in loo_flat]))[1]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--wide-dir", default="analysis/boss/results/widebin_rpar5")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--random-label", default="random_w3_25_rpar5_scw_frac100")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    setup_logging()
    out = os.path.join(args.wide_dir, "comparison", "control_profiles_archived_vs_wide")
    if os.path.exists(out + ".png") and not args.overwrite:
        raise FileExistsError(f"{out}.png exists (use --overwrite)")

    gal = load_accumulator(os.path.join(args.results_dir, "jk", "acc_single_galaxy_scw_BOSS_North_South.npz"))
    rnd = load_accumulator(os.path.join(args.results_dir, "jk",
                                        "acc_single_random_scw_frac100_BOSS_North_South.npz"))
    single = symmetrize_map((loo_and_total(gal)[0] - loo_and_total(rnd)[0]).reshape(100, 100))
    sim_single = load_map(os.path.join(args.sim_dir, "kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv"))

    fig, axes = plt.subplots(3, 3, figsize=(12.5, 9.5), sharey="row", constrained_layout=True)
    m = (AX >= 0) & (AX <= 25)
    for r, (key, name, old_bin, new_bin) in enumerate(ROWS):
        # archived BOSS
        t = build_terms(args.results_dir, "BOSS", ["North", "South"], key, "_scw", "_scw_frac100")
        sep = t["sep"]
        cols = [(axes[r, 0], f"archived BOSS  r$_\\perp$ = {old_bin}", t["physical"], central_err(t["loo"]["physical"]),
                 [(t["control"], "control (at R)", "-")],
                 [(sep / 2, sep / 2)])]
        # wide BOSS
        lo, hi = (float(v) for v in name[1:].split("_"))
        pair = load_grouping_pairs(os.path.join(args.wide_dir, "jk",
                                                f"acc_pairs_galaxy_{name}_rpar5_scw_BOSS_North_South.npz"), lo, hi)
        rand = pooled_random_map(load_random_subbins(args.wide_dir, args.random_label, "BOSS",
                                                     ["North", "South"], lo, hi))
        pt, pl = loo_and_total(pair)
        phys = reflect_symmetrize_map(pt.reshape(101, 101)) - rand
        phys_err = central_err(np.array([(reflect_symmetrize_map(q.reshape(101, 101)) - rand).ravel() for q in pl]))
        s, w, c = separation_mix(pair["fine_w"], pair["fine_wr"], pair["fine_cov"])
        avg = reflect_symmetrize_map(separation_averaged_template(single, X, Y, s, w, coverage=c))
        cols.append((axes[r, 1], f"wide BOSS  r$_\\perp$ = {new_bin}", phys, phys_err,
                     [(avg, "control (averaged over pairs)", "-"),
                      (two_halo_template(single, X, Y, BRIDGE_CENTERS[name]), "control at R only", "--")],
                     [(lo / 2, hi / 2)]))
        # wide mock
        base = os.path.join(args.sim_dir, "hod_pairs", "widebin_rpar5", f"stack_{name}_rpar5_avgctrl")
        blocks = np.load(base + "_blocks.npz")
        sums, counts = blocks["sums"], blocks["counts"]
        mock = load_map(base + "_map.csv")
        loo_mock = [reflect_symmetrize_map(((sums.sum(0) - sums[b]) / (counts.sum() - counts[b])).astype(np.float32))
                    .astype(np.float64).ravel() for b in range(len(counts))]
        ms, mw, _ = separation_mix(blocks["fine_w"], blocks["fine_wr"])
        cols.append((axes[r, 2], f"wide mock  r$_\\perp$ = {new_bin}", mock, central_err(np.array(loo_mock)),
                     [(separation_averaged_template(sim_single, X, Y, ms, mw), "control (averaged over pairs)", "-"),
                      (two_halo_template(sim_single, X, Y, BRIDGE_CENTERS[name]), "control at R only", "--")],
                     [(lo / 2, hi / 2)]))

        for ax, title, data, err, controls, gal_pos in cols:
            ax.axvspan(0, BRIDGE_HALF_X_FRAC * sep, color="#e3e2df", zorder=0)
            for a, b in gal_pos:
                if a == b:
                    ax.axvline(a, color=GREEN, ls=":", lw=1)
                else:
                    ax.axvspan(a, b, color=BLUE, alpha=0.06, zorder=0)
            ax.errorbar(AX[m], central(data)[m] * 1e4, yerr=err[m] * 1e4, fmt="o", ms=3.2, color=BLUE,
                        elinewidth=1, capsize=1.5, label="pair stack", zorder=3)
            for ctrl, label, ls in controls:
                ax.plot(AX[m], central(ctrl)[m] * 1e4, color=ORANGE if ls != ":" else INK, ls=ls, lw=1.8,
                        alpha=0.9 if ls == "-" else 0.7, label=label)
            ax.axhline(0, color=MUTED, lw=0.7)
            ax.set_title(title, loc="left", fontsize=9.5)
            ax.set_xlabel(r"$X$ along pair axis [$h^{-1}$Mpc]")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if r == 0:
                ax.legend(frameon=False, fontsize=7.5, loc="upper right")
        axes[r, 0].set_ylabel(r"central band $\kappa$ [$10^{-4}$]")
        logger.info("row %s done", key)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight")
    logger.info("Saved -> %s.png", out)


if __name__ == "__main__":
    main()
