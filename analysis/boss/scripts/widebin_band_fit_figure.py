#!/usr/bin/env python
"""Candidate figure: the 22-element band vectors and their amplitude fits, wide bins.

For each wide r_perp bin (4-6, 6-15, 15-25) draws the two halves of the
compressed data vector of the text -- the central band (|Y| <= 1.5) and the
off-centre band (1.5 <= |Y| <= 10.5) in 4 h^-1 Mpc bins of X >= 0 -- with
jackknife errors (square root of the covariance diagonal), the mock residual
profile at unit amplitude, and the best fit A * mock from the full-covariance
fit over all 22 elements.  The fit is fit_band_covariance.run_fit, the same
call widebin_vs_archived.py makes, so A, sigma_A and chi^2 are the numbers in
comparison/bridge_study_archived_vs_wide.csv; the script checks that they agree
and stops if they do not.

Not yet in the manuscript.  Writes to widebin_rpar5/comparison/paper_figures/
(PNG and PDF, drawn at print size); paper/sts/Figure/ is not touched.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python \\
        analysis/boss/scripts/widebin_band_fit_figure.py
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

from catalog import setup_logging
from combine_filament_jackknife import loo_and_total
from combine_jackknife import load_accumulator
from fit_band_covariance import run_fit
from geometry import BRIDGE_HALF_X_FRAC
from weighting_sensitivity import load_map

import widebin_paper_figures as wpf
from deconvolve_pair_profile import BLUE, GREEN, INK, INK_MUTED, ORANGE, tidy

logger = logging.getLogger(__name__)

X_MAX, X_BIN = 40.0, 4.0   # the text's compression: 11 bins per band, 22 elements
STEM = "fig_band_fit_22"


def fit_bin(s):
    fit, (x, d, m, cov, _) = run_fit(s["res"], s["res_loo"], s["sim_res"], s["axis"], X_MAX, X_BIN, s["sep"])
    return fit, x, d, m, np.sqrt(np.diag(cov))


def check(fits, table):
    """The plotted fit must be the one the text quotes."""
    ref = pd.read_csv(table)
    ref = ref[ref.version == "wide"].set_index("bin")
    for b, fit in fits.items():
        for ours, theirs in (("A", "A_full"), ("sigma_A", "A_full_err"), ("chi2_A0", "chi2_profile"),
                             ("chi2_min", "chi2_fit")):
            if not np.isclose(fit[ours], ref.loc[b, theirs], rtol=1e-6):
                raise RuntimeError(f"{b}: {ours} = {fit[ours]:.6g}, table has {ref.loc[b, theirs]:.6g}")
    logger.info("Fits match %s", table)


def figure(S, fits, out):
    fig, axes = plt.subplots(2, 3, figsize=(wpf.PRINT_WIDTH, 3.6), sharex="col", sharey="col",
                             constrained_layout=True)
    for c, (key, *_rest) in enumerate(wpf.ROWS):
        s = S[key]
        fit, x, d, m, err = fits[s["bin"]]
        n = len(x)
        lo, hi = (float(t) for t in s["bin"].split("-"))
        for r, (sl, band) in enumerate(((slice(0, n), "central"), (slice(n, 2 * n), "off-center"))):
            ax = axes[r, c]
            ax.axhline(0, color=INK, lw=0.9)
            ax.axvspan(0, BRIDGE_HALF_X_FRAC * s["sep"], color=INK_MUTED, alpha=wpf.SHADE, lw=0, zorder=0)
            ax.axvspan(lo / 2, hi / 2, color=GREEN, alpha=wpf.SHADE, lw=0, zorder=0)
            ax.plot(x, m[sl] * 1e4, lw=1.0, ls="--", color=ORANGE, label="mock", zorder=2)
            ax.plot(x, fit["A"] * m[sl] * 1e4, lw=1.6, color=ORANGE, label=r"$\hat A\,\times$ mock", zorder=2)
            ax.errorbar(x, d[sl] * 1e4, yerr=err[sl] * 1e4, fmt="o", ms=2.6, lw=0, elinewidth=0.8, capsize=0,
                        color=BLUE, ecolor=BLUE, label="BOSS", zorder=3)
            if r == 0:
                wpf.title(ax, rf"$r_\perp$ = {s['bin'].replace('-', '–')} $h^{{-1}}$Mpc" "\n"
                          rf"$\hat A = {fit['A']:.2f} \pm {fit['sigma_A']:.2f}$,  "
                          rf"$\chi^2 = {fit['chi2_min']:.1f}$ / {fit['dof']}")
            else:
                ax.set_xlabel(r"$X$  [$h^{-1}$ Mpc]")
            if c == 0:
                ax.set_ylabel(rf"$\langle F \rangle_{{\rm {'cen' if r == 0 else 'off'}}}$  [$10^{{-4}}$]")
            ax.set_xlim(0, X_MAX + X_BIN / 2)
            tidy(ax)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside upper center", ncol=3, frameon=False)
    out.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(out / f"{STEM}.{e}", bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("Saved %s", out / STEM)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--wide-dir", default="analysis/boss/results/widebin_rpar5")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv")
    p.add_argument("--random-label", default="random_w3_25_rpar5_scw_frac100")
    p.add_argument("--single-tag", default="_scw")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    args = p.parse_args(argv)
    setup_logging()
    wpf.PAPER = True               # bare labels, print-size fonts
    plt.rcParams.update(wpf.PAPER_RC)

    wide_dir, sim_dir = Path(args.wide_dir), Path(args.sim_dir)
    out = wide_dir / "comparison" / "paper_figures"
    sim_single = load_map(sim_dir / args.sim_single)
    gal = load_accumulator(os.path.join(args.results_dir, "jk",
                                        f"acc_single_galaxy{args.single_tag}_BOSS_North_South.npz"))
    rnd = load_accumulator(os.path.join(args.results_dir, "jk",
                                        f"acc_single_random{args.single_random_tag}_BOSS_North_South.npz"))
    (gt, gl), (rt, rl) = loo_and_total(gal), loo_and_total(rnd)
    singles = (gt - rt, gl - rl, gal["grid_size"])

    S = {key: wpf.wide_sample(wide_dir, sim_dir, sim_single, singles, name, new_bin, args.random_label)
         for key, name, _old, new_bin in wpf.ROWS}
    fits = {S[key]["bin"]: fit_bin(S[key]) for key, *_rest in wpf.ROWS}
    check({b: f[0] for b, f in fits.items()},
          wide_dir / "comparison" / "bridge_study_archived_vs_wide.csv")
    logger.info("X bin centres: %s", np.round(fits[S["5"]["bin"]][1], 2))
    figure(S, fits, out)


if __name__ == "__main__":
    main()
