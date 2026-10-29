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

Also writes the residual at the mean galaxy position, with its jackknife error,
for both binnings (residual_at_galaxies.csv; --residuals-only writes just that).

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
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter

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
# Under --paper each figure is drawn at the width it prints, so font sizes and line widths below are
# printed sizes: AASTeX text width 7.1 in, Figures 2-5 and 7 included at 0.9\linewidth, Figure 6 at
# 0.7\linewidth.  (The STS report's 6.5 in text width prints them 8 per cent smaller.)
PRINT_WIDTH = 0.9 * 7.1
PAPER_RC = {"font.size": 8.5, "axes.labelsize": 9, "axes.titlesize": 9,
            "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8}
# Figure 3's residual-map colours (RdBu_r): the central band is the positive, red strip and the
# off-center band its negative, blue surroundings, so Figure 4 draws the two bands in these colours.
RED, DEEP_BLUE = "#b2182b", "#2166ac"
SHADE = 0.2  # opacity of the bridge-window (grey) and galaxy-range (green) shading


def zoom(sep):
    return float(min(max(2.2 * sep, 16.0), 32.0))


def crop(arr, axis, half):
    m = np.abs(axis) <= half
    return arr[np.ix_(m, m)], [-half, half, -half, half]


def mark(ax, sep, color=INK):
    xs = [0.0] if sep == 0 else [-0.5 * sep, 0.5 * sep]
    ax.plot(xs, [0] * len(xs), ls="none", marker="+", ms=7 if PAPER else 10, mew=1.3 if PAPER else 1.6,
            color=color)


def label_axes(ax, x=True, y=True):
    """Under --paper, side-by-side panels share the Y label and stacked ones the X label."""
    if x or not PAPER:
        ax.set_xlabel(r"$X$  [$h^{-1}$ Mpc]")
    if y or not PAPER:
        ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")


def log_ticks(cbar):
    """Label a log colour bar at 1, 2, 5 x 10^n so values can be read off it."""
    cbar.ax.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    cbar.ax.yaxis.set_minor_formatter(NullFormatter())


def at_galaxy(arr, axis, x):
    """Value on the pair axis (Y = 0) at X = x, linearly interpolated between grid points."""
    return float(np.interp(x, axis, arr[int(np.argmin(np.abs(axis)))]))


def title(ax, text):
    # --paper centres titles: on a narrow map a left-aligned title runs into the next colour bar
    ax.set_title(text, loc="center" if PAPER else "left", color=INK, pad=4,
                 **({} if PAPER else {"fontsize": 10.5}))


# --------------------------------------------------------------------------
# One sample = everything the figures need for one bin of one binning.
def archived_sample(results_dir, sim_dir, sim_single, key, bin_text, single_tag, single_random_tag):
    t = build_terms(results_dir, "BOSS", ["North", "South"], key, single_tag, single_random_tag)
    sim_res, sim_axis = sim_diff(key)
    sim_pair = load_map(sim_dir / "hod_pairs" / f"{SIM_STACK[key]}.csv")
    gx, gy = np.meshgrid(sim_axis, sim_axis)
    sim_ctl = two_halo_template(sim_single, gx, gy, float(key))
    loo = np.asarray([band_profile(m.reshape(t["filament"].shape), t["axis"]) for m in t["loo"]["filament"]])
    gal_x = t["sep"] / 2     # narrow symmetric bin: the mean galaxy position is the bin centre
    gal_loo = np.asarray([[at_galaxy(m.reshape(t["filament"].shape), t["axis"], x) for x in (gal_x, 0.0)]
                          for m in t["loo"]["filament"]])
    gal_err, mid_err = jackknife_error(gal_loo)[1]
    return {"bin": bin_text, "sep": t["sep"], "n_pairs": t["n_pairs"], "axis": t["axis"],
            "gal_x": gal_x, "gal_err": float(gal_err), "mid_err": float(mid_err),
            "sim_gal_err": np.nan, "sim_mid_err": np.nan,
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
    loo_maps = [fil(pl[k], sl[k], separation_mix(*fine, exclude_region=k)) for k in range(sl.shape[0])]
    loo = np.asarray([band_profile(m, axis) for m in loo_maps])
    gal_x = float(pair["fine_wr"].sum() / pair["fine_w"].sum()) / 2   # half the mean pair separation
    gal_err, mid_err = jackknife_error(np.asarray([[at_galaxy(m, axis, x) for x in (gal_x, 0.0)]
                                                   for m in loo_maps]))[1]
    phys = reflect_symmetrize_map(pt.reshape(res.shape)) - rand

    base = sim_dir / "hod_pairs" / "widebin_rpar5" / f"stack_{name}_rpar5_avgctrl"
    blocks = np.load(f"{base}_blocks.npz")
    sim_pair = load_map(f"{base}_map.csv")
    sim_axis = axis_for_map(sim_pair)
    gx, gy = np.meshgrid(sim_axis, sim_axis)
    seps, wts, _ = separation_mix(blocks["fine_w"], blocks["fine_wr"])
    sim_ctl = separation_averaged_template(sim_single, gx, gy, seps, wts)
    sums, counts = blocks["sums"], blocks["counts"]
    sim_loo_maps = [reflect_symmetrize_map(((sums.sum(0) - sums[b]) / (counts.sum() - counts[b]))
                                           .astype(np.float32)) - sim_ctl
                    for b in range(len(counts)) if counts.sum() - counts[b] > 0]
    sim_loo = [band_profile(m, sim_axis) for m in sim_loo_maps]
    sim_gal_err, sim_mid_err = jackknife_error(np.asarray(
        [[at_galaxy(m, sim_axis, x) for x in (gal_x, 0.0)] for m in sim_loo_maps]))[1]
    return {"bin": bin_text, "sep": center, "n_pairs": int(pair["n_pairs"].sum()), "axis": axis,
            "gal_x": gal_x, "gal_err": float(gal_err), "mid_err": float(mid_err),
            "sim_gal_err": float(sim_gal_err), "sim_mid_err": float(sim_mid_err),
            "pair": phys, "res": res, "res_loo": np.asarray(loo_maps), "B_err": jackknife_error(loo)[1],
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
    size = (PRINT_WIDTH, 2.95) if PAPER else (9.6, 4.3 * len(SHOW))
    fig, axes = plt.subplots(len(SHOW), 2, figsize=size, constrained_layout=True, squeeze=False)
    for r, v in enumerate(SHOW):
        s = S[v]["10"]
        half = zoom(s["sep"])
        a, ext = crop(s["sim_pair"], s["sim_axis"], half)
        b, _ = crop(s["sim_ctl"], s["sim_axis"], half)
        norm = LogNorm(vmin=max(min(a.min(), b.min()), 1e-5) * 1e4, vmax=max(a.max(), b.max()) * 1e4)
        ctl_name = "control at R" if v == "archived" else "control pair stack"
        for c, (arr, lab) in enumerate(((a, "true pair stack"), (b, ctl_name))):
            text = f"({'ab'[c]}) {lab}" if PAPER else f"mock, {row_label(v, s)}: {lab}"
            ax = axes[r, c]
            im = ax.imshow(arr * 1e4, origin="lower", extent=ext, cmap="magma_r", norm=norm)
            mark(ax, s["sep"], GREEN)
            label_axes(ax, y=c == 0)
            title(ax, text)
        log_ticks(fig.colorbar(im, ax=axes[r], label=r"$\kappa$  [$10^{-4}$]", shrink=0.9))
    save(fig, out / stem(2))


def fig3(S, out):
    """Mock residual maps, one row per binning, one colour scale per column."""
    size = (PRINT_WIDTH, 1.85) if PAPER else (14.7, 4.4 * len(SHOW))
    fig, axes = plt.subplots(len(SHOW), 3, figsize=size, constrained_layout=True, squeeze=False)
    for c, (key, *_rest) in enumerate(ROWS):
        half = zoom(S["archived"][key]["sep"])
        crops = {v: crop(S[v][key]["sim_res"], S[v][key]["sim_axis"], half) for v in SHOW}
        lim = max(float(np.nanpercentile(np.abs(crops[v][0]), 99.0)) for v in SHOW) * 1e4
        for r, v in enumerate(SHOW):
            ax = axes[r, c]
            im = ax.imshow(crops[v][0] * 1e4, origin="lower", extent=crops[v][1], cmap="RdBu_r",
                           vmin=-lim, vmax=lim)
            mark(ax, S[v][key]["sep"])
            label_axes(ax, y=c == 0)
            title(ax, row_label(v, S[v][key]) if PAPER else f"mock, {row_label(v, S[v][key])}")
        # under --paper the three bars share one unit, so only the last is labelled
        fig.colorbar(im, ax=axes[:, c], label=r"$\kappa$  [$10^{-4}$]" if c == 2 or not PAPER else None,
                     shrink=0.6 if len(SHOW) > 1 else (0.8 if PAPER else 0.88), pad=0.06 if PAPER else 0.05)
    save(fig, out / stem(3))


def fig4(S, out):
    """Mock central band, off-centre band and difference."""
    size = (PRINT_WIDTH, 2.4) if PAPER else (14.7, 3.8 * len(SHOW))
    fig, axes = plt.subplots(len(SHOW), 3, figsize=size, constrained_layout=True, sharey="col", squeeze=False)
    lw = 1.6 if PAPER else 2.0
    for c, (key, *_rest) in enumerate(ROWS):
        for r, v in enumerate(SHOW):
            s, ax = S[v][key], axes[r, c]
            cen, off = band_profiles(s["sim_res"], s["sim_axis"])
            half = zoom(s["sep"])
            m = np.abs(s["sim_axis"]) <= half
            x = s["sim_axis"][m]
            ax.plot(x, cen[m] * 1e4, lw=lw, color=RED, label=r"central band  $|Y| \leq 1.5$")
            ax.plot(x, off[m] * 1e4, lw=lw, color=DEEP_BLUE, label=r"off-center band  $1.5 \leq |Y| \leq 10.5$")
            ax.plot(x, (cen - off)[m] * 1e4, lw=lw + 0.3, color=INK, label="difference")
            ax.axhline(0, color=INK_MUTED, lw=0.9)
            ax.axvspan(-BRIDGE_HALF_X_FRAC * s["sep"], BRIDGE_HALF_X_FRAC * s["sep"],
                       color=INK_MUTED, alpha=SHADE, lw=0)
            for xx in (-0.5 * s["sep"], 0.5 * s["sep"]):
                ax.axvline(xx, color=GREEN, lw=1.8, ls=":")
            ax.set_xlim(-half, half)
            ax.set_xlabel(r"$X$  [$h^{-1}$ Mpc]")
            if c == 0 or not PAPER:
                ax.set_ylabel(r"$\kappa$  [$10^{-4}$]")
            title(ax, row_label(v, s) if PAPER else f"mock, {row_label(v, s)}")
            if r == 0 and c == 0 and not PAPER:
                ax.legend(loc="lower left", frameon=False)
            tidy(ax)
    if PAPER:  # above the panels, where it cannot cover a curve
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside upper center", ncol=3, frameon=False)
    save(fig, out / stem(4))


def fig5(single, S, out):
    """Observed single stack and the pair stacks; one linear scale throughout."""
    panels = [(r, c + 1, (f"({'bcd'[c]}) pairs\n" if PAPER else "pairs, ") + row_label(v, S[v][key]),
               S[v][key]["pair"], S[v][key]["axis"],
               S[v][key]["sep"]) for r, v in enumerate(SHOW) for c, (key, *_rest) in enumerate(ROWS)]
    s_axis = axis_for_map(single)
    vmax = max([float(crop(a, ax_, ZOOM_FIG5)[0].max()) for *_, a, ax_, _ in panels]
               + [float(crop(single, s_axis, ZOOM_FIG5)[0].max())]) * 1e4
    size = (PRINT_WIDTH, 1.9) if PAPER else (17.2, 4.3 * len(SHOW))
    fig, axes = plt.subplots(len(SHOW), 4, figsize=size, constrained_layout=True, squeeze=False)
    a, ext = crop(single, s_axis, ZOOM_FIG5)
    im = axes[0, 0].imshow(a * 1e4, origin="lower", extent=ext, cmap="magma_r", vmin=0, vmax=vmax,
                           interpolation="nearest")
    mark(axes[0, 0], 0, GREEN)
    label_axes(axes[0, 0])
    title(axes[0, 0], ("(a) single\ngalaxies" if PAPER else "single galaxies")
          + (" (both binnings)" if len(SHOW) > 1 else ""))
    for r in range(1, len(SHOW)):
        axes[r, 0].axis("off")
    for r, c, lab, arr, axis, sep in panels:
        a, ext = crop(arr, axis, ZOOM_FIG5)
        axes[r, c].imshow(a * 1e4, origin="lower", extent=ext, cmap="magma_r", vmin=0, vmax=vmax,
                          interpolation="nearest")
        mark(axes[r, c], sep, GREEN)
        label_axes(axes[r, c], y=False)
        title(axes[r, c], lab)
    fig.colorbar(im, ax=axes, label=r"$\kappa$  [$10^{-4}$]",
                 shrink=0.6 if len(SHOW) > 1 else (0.72 if PAPER else 0.9), pad=0.04 if PAPER else 0.05)
    save(fig, out / stem(5))


def fig6(S, out):
    """Observed and mock residual maps; each row (separation) shares one scale."""
    size = (0.7 * 7.1, 6.9) if PAPER else (4.75 * 2 * len(SHOW), 13.2)
    fig, axes = plt.subplots(3, 2 * len(SHOW), figsize=size, constrained_layout=True)
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
            label_axes(axes[r, c], x=r == len(ROWS) - 1, y=c == 0)
            title(axes[r, c], f"{which}, {row_label(v, s)}" if PAPER else f"{row_label(v, s)}, {which}")
        fig.colorbar(im, ax=axes[r], label=r"$\kappa$  [$10^{-4}$]", shrink=0.88)
    save(fig, out / stem(6))


def fig7(S, out):
    """B(X), BOSS against the mock, with (BOSS - mock)/sigma strips."""
    n = len(FIG7_SHOW)
    fig, axes = plt.subplots(2 * n, 3, figsize=(PRINT_WIDTH, 3.0) if PAPER else (13.8, 5.4 * n),
                             height_ratios=[3, 1] * n,
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
            top.axvspan(-BRIDGE_HALF_X_FRAC * sep, BRIDGE_HALF_X_FRAC * sep, color=INK_MUTED, alpha=SHADE,
                        lw=0, zorder=0)
            if v == "archived":
                for xx in (-0.5 * sep, 0.5 * sep):
                    top.axvline(xx, color=GREEN, lw=1.8, ls=":")
            else:
                lo, hi = (float(t) for t in s["bin"].split("-"))
                for sgn in (-1, 1):
                    top.axvspan(sgn * lo / 2, sgn * hi / 2, color=GREEN, alpha=SHADE, lw=0, zorder=0)
            top.fill_between(sim_axis[ms], (sim_prof - sim_err)[ms] * 1e4, (sim_prof + sim_err)[ms] * 1e4,
                             color=ORANGE, alpha=0.22, lw=0, zorder=1)
            top.plot(sim_axis[ms], sim_prof[ms] * 1e4, lw=1.6 if PAPER else 2.0, color=ORANGE, label="mock",
                     zorder=2)
            top.errorbar(axis[m], prof[m] * 1e4, yerr=err[m] * 1e4, fmt="o", ms=2.2 if PAPER else 3.0, lw=0,
                         elinewidth=0.8 if PAPER else 1.0, capsize=0, color=BLUE, ecolor=BLUE, label="BOSS",
                         zorder=3)
            if c == 0 or not PAPER:
                top.set_ylabel(r"$B(X)$  [$10^{-4}$]")
            title(top, f"{row_label(v, s)}\n{s['n_pairs']:,} pairs" if PAPER
                  else f"{row_label(v, s)}  ({s['n_pairs']:,} pairs)")
            if r == 0 and c == 0 and not PAPER:
                top.legend(loc="upper right", frameon=False)
            tidy(top)
            plt.setp(top.get_xticklabels(), visible=False)

            resid = (prof[m] - np.interp(axis[m], sim_axis, sim_prof)) / err[m]
            bot.axhline(0, color=INK, lw=0.9)
            for lev in (-2.0, 2.0):
                bot.axhline(lev, color=INK_MUTED, lw=0.7, ls="--")
            bot.plot(axis[m], resid, lw=0, marker="o", ms=2.0 if PAPER else 2.6, color=BLUE)
            bot.set_ylim(-4.0, 4.0)
            bot.set_yticks([-2, 0, 2])
            bot.set_xlim(-half, half)
            if c == 0 or not PAPER:
                bot.set_ylabel(r"$\Delta/\sigma$")
            bot.set_xlabel(r"$X$  [$h^{-1}$ Mpc]")
            bot.grid(True, axis="y", color=GRID, linewidth=0.6)
            tidy(bot)
        # one y-range per column so the two binnings compare directly
        lims = [axes[2 * r, c].get_ylim() for r in range(n)]
        for r in range(n):
            axes[2 * r, c].set_ylim(min(l[0] for l in lims), max(l[1] for l in lims))
    if PAPER:  # above the panels, where it cannot cover a data point
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside upper center", ncol=2, frameon=False)
    save(fig, out / stem(7))


def galaxy_residuals(S):
    """Residual on the pair axis at the mean galaxy position, data and mock.

    The maps are reflection-symmetrized, so the value is the mean of the two
    galaxies' pixels; an unsymmetrized map (the kappa_filament_* CSVs) differs
    between +X and -X by about one per-pixel error.  The position is half the
    mean pair separation, interpolated between grid points, and the error is the
    jackknife error of that interpolated value.  In a wide bin the galaxies span
    a range of X, so the peak along the axis is reported as well.  The same is
    given at the pair midpoint (X = 0), with the pair-stack convergence there for
    scale.  Mock errors are the 25-block jackknife with the control held fixed;
    they are not computed for the archived bins.
    """
    rows = []
    for v in VERSIONS:
        for key, *_rest in ROWS:
            s = S[v][key]
            obs = at_galaxy(s["res"], s["axis"], s["gal_x"])
            m = (s["axis"] >= 0) & (s["axis"] <= s["sep"])
            row = s["res"][int(np.argmin(np.abs(s["axis"])))]
            rows.append({"version": v, "bin": s["bin"], "x_galaxy_hmpc": s["gal_x"],
                         "obs_at_galaxy": obs, "obs_at_galaxy_err": s["gal_err"],
                         "obs_at_galaxy_sigma": obs / s["gal_err"],
                         "mock_at_galaxy": at_galaxy(s["sim_res"], s["sim_axis"], s["gal_x"]),
                         "mock_at_galaxy_err": s["sim_gal_err"],
                         "obs_at_midpoint": at_galaxy(s["res"], s["axis"], 0.0),
                         "obs_at_midpoint_err": s["mid_err"],
                         "mock_at_midpoint": at_galaxy(s["sim_res"], s["sim_axis"], 0.0),
                         "mock_at_midpoint_err": s["sim_mid_err"],
                         "obs_pair_kappa_at_midpoint": at_galaxy(s["pair"], s["axis"], 0.0),
                         "mock_pair_kappa_at_midpoint": at_galaxy(s["sim_pair"], s["sim_axis"], 0.0),
                         "obs_axis_peak": float(row[m].max()),
                         "obs_axis_peak_x_hmpc": float(s["axis"][m][int(np.argmax(row[m]))]),
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
    p.add_argument("--residuals-only", action="store_true",
                   help="Write residual_at_galaxies.csv and draw no figure.")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    setup_logging()

    global SHOW, SUFFIX, FIG7_SHOW, PAPER
    if args.wide_only:
        SHOW, SUFFIX = ("wide",), "_wide"
    if args.paper:
        PAPER, SHOW, FIG7_SHOW = True, ("wide",), ("wide",)
        plt.rcParams.update(PAPER_RC)
    wide_dir, sim_dir = Path(args.wide_dir), Path(args.sim_dir)
    default_out = Path("paper/sts/Figure") if PAPER else wide_dir / "comparison" / "paper_figures"
    out = Path(args.output_dir) if args.output_dir else default_out
    targets = [out / f"{stem(n)}.pdf" for n in ((2, 3, 4, 5, 6, 7) if PAPER else (2, 3, 4, 5, 6))]
    if args.residuals_only:
        targets = [out / "residual_at_galaxies.csv"]
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

    if args.residuals_only:
        galaxy_residuals(S).to_csv(out / "residual_at_galaxies.csv", index=False)
        logger.info("Saved %s", out / "residual_at_galaxies.csv")
        return
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
