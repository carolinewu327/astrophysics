#!/usr/bin/env python
"""The manuscript's figures for the wide bins, alone or beside the archived bins.

This is the source of the figures in paper/sts/Figure/ (run with --paper), and
of the archived-vs-wide comparison figures in widebin_rpar5/comparison/.

Redraws Figures 2-7 of paper/styau (the submitted draft) twice, once per
binning, in one image per figure so the two can be read side by side:

    archived   r_perp 4-6 (r_par<=5), 9-11 and 18-22 (r_par<=10); nominal
               single-separation control -- reproduces the manuscript panels
    wide       r_perp 4-6, 6-15, 15-25 (all r_par<=5); every stack weighted by
               1/Sigma_crit^2; control averaged over the pairs' actual
               separations, in the data and in the mock alike

  fig2  mock pair stack beside its control, R = 10
  fig3  mock residual maps, R = 5, 10, 20
  fig4  mock central band, off-centre band and their difference
  fig5  BOSS single stack and the three pair stacks
  fig6  observed beside mock residual maps
  fig7  B(X), BOSS against the mock, with the residual strip

Figure 1 (footprint, n(z)) does not depend on the binning.  Figure 8 (smoothing
scan) would need the wide stacks re-run at 6 and 4 arcmin and is not redrawn.
With --wide-only, figures 2-6 show the wide bins alone (files *_wide.*) and the
reader compares them with the manuscript; figure 7 keeps both binnings, since
B(X) is the measurement itself.

With --paper, all six figures show the wide bins alone, panel titles are bare
labels, and the files keep the names the manuscript's \\includegraphics lines
use.  They are written as PDF only into paper/sts/Figure/ (--output-dir to
change).  paper/styau/ is the submitted version and is never written.

Also writes the residual at the galaxy positions (the "2e-3 vs 0.9e-3" sentence
of section 5.1) for both binnings.

Outputs go to widebin_rpar5/comparison/paper_figures/, or paper/sts/Figure/
under --paper; nothing in the archive or paper/styau/ is written.

Usage
-----
    export PYTHONPATH=lib:analysis/boss/scripts:analysis/sim
    S=analysis/boss/scripts/widebin_paper_figures.py
    python $S --paper       # figures 2-7 for the manuscript -> paper/sts/Figure/
    python $S               # archived beside wide           -> comparison/paper_figures/
    python $S --wide-only   # wide alone, annotated titles   -> comparison/paper_figures/*_wide.*
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
from matplotlib.colors import LogNorm

from catalog import setup_logging
from combine_filament_jackknife import (BRIDGE_CENTERS, grouping_filament, load_grouping_pairs,
                                        load_random_subbins, loo_and_total, pooled_random_map)
from combine_jackknife import load_accumulator
from geometry import (BRIDGE_HALF_X_FRAC, band_profile, band_profiles, reflect_symmetrize_map,
                      separation_averaged_template, separation_mix, two_halo_template)
from jackknife import jackknife_error
from plot_observed_figures import SIM_STACK, boss_single, sim_band_error, sim_diff
from plot_separation_summary import build_terms
from weighting_sensitivity import axis_for_map, load_map

from deconvolve_pair_profile import BLUE, GREEN, GRID, INK, INK_MUTED, ORANGE, tidy

logger = logging.getLogger(__name__)

ROWS = [("5", "g4_6", "4-6", "4-6"), ("10", "g6_15", "9-11", "6-15"), ("20", "g15_25", "18-22", "15-25")]
VERSIONS = ("archived", "wide")
SHOW = VERSIONS      # binnings drawn in figures 2-6; ("wide",) under --wide-only
SUFFIX = ""          # "_wide" under --wide-only
FIG7_SHOW = VERSIONS  # binnings drawn in figure 7; ("wide",) under --paper
PAPER = False        # --paper: manuscript file names and bare panel labels
STEMS = {2: "fig2_mock_pair_vs_control", 3: "fig3_mock_residual_maps", 4: "fig4_mock_band_profiles",
         5: "fig5_boss_single_and_pair_maps", 6: "fig6_obs_vs_mock_residual_maps",
         7: "fig7_band_profiles_obs_vs_mock"}
# Names of the manuscript's own figures (figure manifests in paper/README.md and paper/sts/README.md).
PAPER_STEMS = {2: "sim_pair_vs_control_maps_sep10", 3: "sim_quadrupole_diff_maps",
               4: "sim_band_profiles_all_seps", 5: "boss_single_and_pair_maps",
               6: "obs_vs_sim_quadrupole_maps", 7: "band_profiles_obs_vs_sim"}
ZOOM_FIG5 = 26.0


def zoom(sep):
    return float(min(max(2.2 * sep, 16.0), 32.0))


def crop(arr, axis, half):
    m = np.abs(axis) <= half
    return arr[np.ix_(m, m)], [-half, half, -half, half]


def mark(ax, sep, color=INK):
    ax.plot([-0.5 * sep, 0.5 * sep], [0, 0], ls="none", marker="+", ms=10, mew=1.6, color=color)


def label_axes(ax):
    ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
    ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")


def title(ax, text):
    ax.set_title(text, loc="left", color=INK, fontsize=10.5, pad=5)


# --------------------------------------------------------------------------
# One sample = everything the figures need for one bin of one binning.
def archived_sample(results_dir, sim_dir, sim_single, key, bin_text, single_tag, single_random_tag):
    t = build_terms(results_dir, "BOSS", ["North", "South"], key, single_tag, single_random_tag)
    sim_res, sim_axis = sim_diff(key)
    sim_pair = load_map(sim_dir / "hod_pairs" / f"{SIM_STACK[key]}.csv")
    gx, gy = np.meshgrid(sim_axis, sim_axis)
    sim_ctl = two_halo_template(sim_single, gx, gy, float(key))
    loo = np.asarray([band_profile(m.reshape(t["filament"].shape), t["axis"]) for m in t["loo"]["filament"]])
    return {"bin": bin_text, "sep": t["sep"], "n_pairs": t["n_pairs"], "axis": t["axis"],
            "pair": t["physical"], "res": t["filament"], "B_err": jackknife_error(loo)[1],
            "sim_axis": sim_axis, "sim_pair": sim_pair, "sim_ctl": sim_ctl, "sim_res": sim_res,
            "sim_B_err": sim_band_error(key, sim_axis, sim_ctl)}


def wide_sample(wide_dir, sim_dir, sim_single, singles, name, bin_text, random_label):
    st, sl, g = singles
    lo, hi = (float(v) for v in name[1:].split("_"))
    center = BRIDGE_CENTERS[name]
    pair = load_grouping_pairs(wide_dir / "jk" / f"acc_pairs_galaxy_{name}_rpar5_scw_BOSS_North_South.npz",
                               lo, hi)
    rand = pooled_random_map(load_random_subbins(str(wide_dir), random_label, "BOSS",
                                                 ["North", "South"], lo, hi))
    pt, pl = loo_and_total(pair)
    fine = (pair["fine_w"], pair["fine_wr"], pair["fine_cov"])

    def fil(p, s, mix):
        return grouping_filament(p, s, rand, mix, center, "separation_averaged", g, symmetrize=True)[0]

    res = fil(pt, st, separation_mix(*fine))
    axis = axis_for_map(res)
    loo = np.asarray([band_profile(fil(pl[k], sl[k], separation_mix(*fine, exclude_region=k)), axis)
                      for k in range(sl.shape[0])])
    phys = reflect_symmetrize_map(pt.reshape(res.shape)) - rand

    base = sim_dir / "hod_pairs" / "widebin_rpar5" / f"stack_{name}_rpar5_avgctrl"
    blocks = np.load(f"{base}_blocks.npz")
    sim_pair = load_map(f"{base}_map.csv")
    sim_axis = axis_for_map(sim_pair)
    gx, gy = np.meshgrid(sim_axis, sim_axis)
    seps, wts, _ = separation_mix(blocks["fine_w"], blocks["fine_wr"])
    sim_ctl = separation_averaged_template(sim_single, gx, gy, seps, wts)
    sums, counts = blocks["sums"], blocks["counts"]
    sim_loo = [band_profile(reflect_symmetrize_map(((sums.sum(0) - sums[b]) / (counts.sum() - counts[b]))
                                                   .astype(np.float32)) - sim_ctl, sim_axis)
               for b in range(len(counts)) if counts.sum() - counts[b] > 0]
    return {"bin": bin_text, "sep": center, "n_pairs": int(pair["n_pairs"].sum()), "axis": axis,
            "pair": phys, "res": res, "B_err": jackknife_error(loo)[1],
            "sim_axis": sim_axis, "sim_pair": sim_pair, "sim_ctl": sim_ctl, "sim_res": sim_pair - sim_ctl,
            "sim_B_err": jackknife_error(np.asarray(sim_loo))[1]}


def stem(n):
    if PAPER:
        return PAPER_STEMS[n]
    return STEMS[n] + ("" if n == 7 else SUFFIX)


def row_label(version, s):
    if PAPER:
        return rf"$r_\perp$ = {s['bin'].replace('-', '–')} $h^{{-1}}$Mpc"
    return f"{version}  r$_\\perp$ = {s['bin']}"


# --------------------------------------------------------------------------
def fig2(S, out):
    """Mock pair stack beside its control at R = 10, one row per binning."""
    fig, axes = plt.subplots(len(SHOW), 2, figsize=(9.6, 4.3 * len(SHOW)), constrained_layout=True,
                             squeeze=False)
    for r, v in enumerate(SHOW):
        s = S[v]["10"]
        half = zoom(s["sep"])
        a, ext = crop(s["sim_pair"], s["sim_axis"], half)
        b, _ = crop(s["sim_ctl"], s["sim_axis"], half)
        norm = LogNorm(vmin=max(min(a.min(), b.min()), 1e-5) * 1e4, vmax=max(a.max(), b.max()) * 1e4)
        ctl_name = "control at R" if v == "archived" else "control averaged over pairs"
        for c, (arr, lab) in enumerate(((a, "true pair stack"), (b, ctl_name))):
            text = f"({'ab'[c]}) {lab}" if PAPER else f"mock, {row_label(v, s)}: {lab}"
            ax = axes[r, c]
            im = ax.imshow(arr * 1e4, origin="lower", extent=ext, cmap="magma_r", norm=norm)
            mark(ax, s["sep"], GREEN)
            label_axes(ax)
            title(ax, text)
        fig.colorbar(im, ax=axes[r], label=r"$\kappa$  [$10^{-4}$]", shrink=0.9)
    save(fig, out / stem(2))


def fig3(S, out):
    """Mock residual maps, one row per binning, one colour scale per column."""
    fig, axes = plt.subplots(len(SHOW), 3, figsize=(14.7, 4.4 * len(SHOW)), constrained_layout=True,
                             squeeze=False)
    for c, (key, *_rest) in enumerate(ROWS):
        half = zoom(S["archived"][key]["sep"])
        crops = {v: crop(S[v][key]["sim_res"], S[v][key]["sim_axis"], half) for v in SHOW}
        lim = max(float(np.nanpercentile(np.abs(crops[v][0]), 99.0)) for v in SHOW) * 1e4
        for r, v in enumerate(SHOW):
            ax = axes[r, c]
            im = ax.imshow(crops[v][0] * 1e4, origin="lower", extent=crops[v][1], cmap="RdBu_r",
                           vmin=-lim, vmax=lim)
            mark(ax, S[v][key]["sep"])
            label_axes(ax)
            title(ax, row_label(v, S[v][key]) if PAPER else f"mock, {row_label(v, S[v][key])}")
        fig.colorbar(im, ax=axes[:, c], label=r"$\kappa$  [$10^{-4}$]", shrink=0.6 if len(SHOW) > 1 else 0.88)
    save(fig, out / stem(3))


def fig4(S, out):
    """Mock central band, off-centre band and difference."""
    fig, axes = plt.subplots(len(SHOW), 3, figsize=(14.7, 3.8 * len(SHOW)), constrained_layout=True,
                             sharey="col", squeeze=False)
    for c, (key, *_rest) in enumerate(ROWS):
        for r, v in enumerate(SHOW):
            s, ax = S[v][key], axes[r, c]
            cen, off = band_profiles(s["sim_res"], s["sim_axis"])
            half = zoom(s["sep"])
            m = np.abs(s["sim_axis"]) <= half
            x = s["sim_axis"][m]
            ax.plot(x, cen[m] * 1e4, lw=2.0, color=BLUE, label=r"central band  $|Y| \leq 1.5$")
            ax.plot(x, off[m] * 1e4, lw=2.0, color=ORANGE, label=r"off-centre band  $1.5 \leq |Y| \leq 10.5$")
            ax.plot(x, (cen - off)[m] * 1e4, lw=2.3, color=INK, label="difference")
            ax.axhline(0, color=INK_MUTED, lw=0.9)
            ax.axvspan(-BRIDGE_HALF_X_FRAC * s["sep"], BRIDGE_HALF_X_FRAC * s["sep"],
                       color=INK_MUTED, alpha=0.09, lw=0)
            for xx in (-0.5 * s["sep"], 0.5 * s["sep"]):
                ax.axvline(xx, color=GREEN, lw=1.1, ls=":")
            ax.set_xlim(-half, half)
            ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
            ax.set_ylabel(r"$\kappa$  [$10^{-4}$]")
            title(ax, row_label(v, s) if PAPER else f"mock, {row_label(v, s)}")
            if r == 0 and c == 0:
                ax.legend(loc="lower left", fontsize=8, frameon=False)
            tidy(ax)
    save(fig, out / stem(4))


def fig5(single, S, out):
    """Observed single stack and the pair stacks; one linear scale throughout."""
    panels = [(r, c + 1, (f"({'bcd'[c]}) " if PAPER else "") + f"pairs, {row_label(v, S[v][key])}",
               S[v][key]["pair"], S[v][key]["axis"],
               S[v][key]["sep"]) for r, v in enumerate(SHOW) for c, (key, *_rest) in enumerate(ROWS)]
    s_axis = axis_for_map(single)
    vmax = max([float(crop(a, ax_, ZOOM_FIG5)[0].max()) for *_, a, ax_, _ in panels]
               + [float(crop(single, s_axis, ZOOM_FIG5)[0].max())]) * 1e4
    fig, axes = plt.subplots(len(SHOW), 4, figsize=(17.2, 4.3 * len(SHOW)), constrained_layout=True,
                             squeeze=False)
    a, ext = crop(single, s_axis, ZOOM_FIG5)
    im = axes[0, 0].imshow(a * 1e4, origin="lower", extent=ext, cmap="magma_r", vmin=0, vmax=vmax,
                           interpolation="nearest")
    axes[0, 0].plot([0], [0], ls="none", marker="+", ms=10, mew=1.6, color=GREEN)
    label_axes(axes[0, 0])
    title(axes[0, 0], ("(a) " if PAPER else "") + "single galaxies" + (" (both binnings)" if len(SHOW) > 1 else ""))
    for r in range(1, len(SHOW)):
        axes[r, 0].axis("off")
    for r, c, lab, arr, axis, sep in panels:
        a, ext = crop(arr, axis, ZOOM_FIG5)
        axes[r, c].imshow(a * 1e4, origin="lower", extent=ext, cmap="magma_r", vmin=0, vmax=vmax,
                          interpolation="nearest")
        mark(axes[r, c], sep, GREEN)
        label_axes(axes[r, c])
        title(axes[r, c], lab)
    fig.colorbar(im, ax=axes, label=r"$\kappa$  [$10^{-4}$]", shrink=0.6 if len(SHOW) > 1 else 0.9)
    save(fig, out / stem(5))


def fig6(S, out):
    """Observed and mock residual maps; each row (separation) shares one scale."""
    fig, axes = plt.subplots(3, 2 * len(SHOW), figsize=(4.75 * 2 * len(SHOW), 13.2), constrained_layout=True)
    for r, (key, *_rest) in enumerate(ROWS):
        half = zoom(S["archived"][key]["sep"])
        lim = max(float(np.nanpercentile(np.abs(crop(S[v][key]["res"], S[v][key]["axis"], half)[0]), 99.0))
                  for v in SHOW) * 1e4
        for c, (v, which) in enumerate([(v, w) for v in SHOW for w in ("observed", "mock")]):
            s = S[v][key]
            arr, axis = (s["res"], s["axis"]) if which == "observed" else (s["sim_res"], s["sim_axis"])
            a, ext = crop(arr, axis, half)
            im = axes[r, c].imshow(a * 1e4, origin="lower", extent=ext, cmap="RdBu_r", vmin=-lim, vmax=lim,
                                   interpolation="nearest")
            mark(axes[r, c], s["sep"])
            label_axes(axes[r, c])
            title(axes[r, c], f"{row_label(v, s)} — {which}" if PAPER else f"{row_label(v, s)}, {which}")
        fig.colorbar(im, ax=axes[r], label=r"$\kappa$  [$10^{-4}$]", shrink=0.88)
    save(fig, out / stem(6))


def fig7(S, out):
    """B(X), BOSS against the mock, with (BOSS - mock)/sigma strips."""
    n = len(FIG7_SHOW)
    fig, axes = plt.subplots(2 * n, 3, figsize=(13.8, 5.4 * n), height_ratios=[3, 1] * n,
                             constrained_layout=True, squeeze=False)
    for c, (key, *_rest) in enumerate(ROWS):
        for r, v in enumerate(FIG7_SHOW):
            s = S[v][key]
            axis, sep = s["axis"], s["sep"]
            prof, err = band_profile(s["res"], axis), s["B_err"]
            sim_axis = s["sim_axis"]
            sim_prof, sim_err = band_profile(s["sim_res"], sim_axis), s["sim_B_err"]
            half = zoom(sep)
            m, ms = np.abs(axis) <= half, np.abs(sim_axis) <= half

            top, bot = axes[2 * r, c], axes[2 * r + 1, c]
            top.sharex(bot)
            top.axhline(0, color=INK, lw=0.9)
            top.axvspan(-BRIDGE_HALF_X_FRAC * sep, BRIDGE_HALF_X_FRAC * sep, color=INK_MUTED, alpha=0.09,
                        lw=0, zorder=0)
            if v == "archived":
                for xx in (-0.5 * sep, 0.5 * sep):
                    top.axvline(xx, color=GREEN, lw=1.1, ls=":")
            else:
                lo, hi = (float(t) for t in s["bin"].split("-"))
                for sgn in (-1, 1):
                    top.axvspan(sgn * lo / 2, sgn * hi / 2, color=GREEN, alpha=0.08, lw=0, zorder=0)
            top.fill_between(sim_axis[ms], (sim_prof - sim_err)[ms] * 1e4, (sim_prof + sim_err)[ms] * 1e4,
                             color=ORANGE, alpha=0.22, lw=0, zorder=1)
            top.plot(sim_axis[ms], sim_prof[ms] * 1e4, lw=2.0, color=ORANGE, label="mock", zorder=2)
            top.errorbar(axis[m], prof[m] * 1e4, yerr=err[m] * 1e4, fmt="o", ms=3.0, lw=0, elinewidth=1.0,
                         capsize=0, color=BLUE, ecolor=BLUE, label="BOSS", zorder=3)
            top.set_ylabel(r"$B(X)$  [$10^{-4}$]")
            title(top, f"{row_label(v, s)}  ({s['n_pairs']:,} pairs)")
            if r == 0 and c == 0:
                top.legend(loc="upper right", fontsize=8.4, frameon=False)
            tidy(top)
            plt.setp(top.get_xticklabels(), visible=False)

            resid = (prof[m] - np.interp(axis[m], sim_axis, sim_prof)) / err[m]
            bot.axhline(0, color=INK, lw=0.9)
            for lev in (-2.0, 2.0):
                bot.axhline(lev, color=INK_MUTED, lw=0.7, ls="--")
            bot.plot(axis[m], resid, lw=0, marker="o", ms=2.6, color=BLUE)
            bot.set_ylim(-4.0, 4.0)
            bot.set_xlim(-half, half)
            bot.set_ylabel(r"$\Delta/\sigma$")
            bot.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
            bot.grid(True, axis="y", color=GRID, linewidth=0.6)
            tidy(bot)
        # one y-range per column so the two binnings compare directly
        lims = [axes[2 * r, c].get_ylim() for r in range(n)]
        for r in range(n):
            axes[2 * r, c].set_ylim(min(l[0] for l in lims), max(l[1] for l in lims))
    save(fig, out / stem(7))


def galaxy_residuals(S):
    """Residual at the galaxy positions (X = +/-R/2, Y = 0) and the map maximum, data and mock."""
    rows = []
    for v in VERSIONS:
        for key, *_rest in ROWS:
            s = S[v][key]
            i0 = int(np.argmin(np.abs(s["axis"])))
            ig = int(np.argmin(np.abs(s["axis"] - s["sep"] / 2)))
            j0 = int(np.argmin(np.abs(s["sim_axis"])))
            jg = int(np.argmin(np.abs(s["sim_axis"] - s["sep"] / 2)))
            rows.append({"version": v, "bin": s["bin"], "obs_at_galaxy": s["res"][i0, ig],
                         "mock_at_galaxy": s["sim_res"][j0, jg],
                         "obs_max_within_R": crop(s["res"], s["axis"], s["sep"])[0].max(),
                         "mock_max_within_R": crop(s["sim_res"], s["sim_axis"], s["sep"])[0].max()})
    return pd.DataFrame(rows)


def save(fig, base):
    base.parent.mkdir(parents=True, exist_ok=True)
    for e in (("pdf",) if PAPER else ("png", "pdf")):
        fig.savefig(f"{base}.{e}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Saved %s", base)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--wide-dir", default="analysis/boss/results/widebin_rpar5")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv")
    p.add_argument("--random-label", default="random_w3_25_rpar5_scw_frac100")
    p.add_argument("--single-tag", default="_scw")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    p.add_argument("--wide-only", action="store_true",
                   help="Figures 2-6 for the wide bins alone (*_wide files); figure 7 keeps both.")
    p.add_argument("--paper", action="store_true",
                   help="Manuscript versions: wide bins only in all six figures, bare panel labels, "
                        "the manuscript's file names, PDF only, written to paper/sts/Figure/.")
    p.add_argument("--output-dir", default=None,
                   help="Default: <wide-dir>/comparison/paper_figures, or paper/sts/Figure "
                        "under --paper.")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    setup_logging()

    global SHOW, SUFFIX, FIG7_SHOW, PAPER
    if args.wide_only:
        SHOW, SUFFIX = ("wide",), "_wide"
    if args.paper:
        PAPER, SHOW, FIG7_SHOW = True, ("wide",), ("wide",)
    wide_dir, sim_dir = Path(args.wide_dir), Path(args.sim_dir)
    default_out = Path("paper/sts/Figure") if PAPER else wide_dir / "comparison" / "paper_figures"
    out = Path(args.output_dir) if args.output_dir else default_out
    targets = [out / f"{stem(n)}.pdf" for n in ((2, 3, 4, 5, 6, 7) if PAPER else (2, 3, 4, 5, 6))]
    if any(t.exists() for t in targets) and not args.overwrite:
        raise FileExistsError(f"Figures exist in {out} (use --overwrite)")
    out.mkdir(parents=True, exist_ok=True)

    sim_single = load_map(sim_dir / args.sim_single)
    gal = load_accumulator(os.path.join(args.results_dir, "jk",
                                        f"acc_single_galaxy{args.single_tag}_BOSS_North_South.npz"))
    rnd = load_accumulator(os.path.join(args.results_dir, "jk",
                                        f"acc_single_random{args.single_random_tag}_BOSS_North_South.npz"))
    (gt, gl), (rt, rl) = loo_and_total(gal), loo_and_total(rnd)
    singles = (gt - rt, gl - rl, gal["grid_size"])

    S = {"archived": {}, "wide": {}}
    for key, name, old_bin, new_bin in ROWS:
        S["archived"][key] = archived_sample(args.results_dir, sim_dir, sim_single, key, old_bin,
                                             args.single_tag, args.single_random_tag)
        S["wide"][key] = wide_sample(wide_dir, sim_dir, sim_single, singles, name, new_bin, args.random_label)
        logger.info("R=%s: archived %d pairs, wide %d pairs", key, S["archived"][key]["n_pairs"],
                    S["wide"][key]["n_pairs"])

    single = boss_single(args.results_dir, "BOSS", ["North", "South"], args.single_tag, args.single_random_tag)
    fig2(S, out)
    fig3(S, out)
    fig4(S, out)
    fig5(single, S, out)
    fig6(S, out)
    fig7(S, out)
    if not PAPER:
        galaxy_residuals(S).to_csv(out / "residual_at_galaxies.csv", index=False)
    logger.info("Saved -> %s", out)


if __name__ == "__main__":
    main()
