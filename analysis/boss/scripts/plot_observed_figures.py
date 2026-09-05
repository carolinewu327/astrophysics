#!/usr/bin/env python
"""The observed-side figures for section 4.2 of the paper.

  fig1  the stacked symmetrized single, then the 5, 10, 20 pair kappa maps
  fig2  the (a - b) difference maps, observed beside theoretical, one row per
        separation

Both are BOSS CMASS x Planck, North and South stacked jointly, galaxy stack
minus random stack throughout.  The control (b) is the plain superposed single
-- two copies of the measured single-galaxy profile at +/- a/2, no fitted
parameters -- built natively on the pair grid via geometry.two_halo_template.
That is the unfitted control on purpose: these figures exist to show what is
left over *before* any fitting, which is the motivation for the core+tail model
rather than a result of it.

The r_perp = 5 sample carries r_par <= 5 while 10 and 20 carry r_par <= 10.
The matching r_par <= 5 random-pair maps are the ones that exist; using the
r_par <= 10 randoms there would subtract a differently-cut background.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python \\
        analysis/boss/scripts/plot_observed_figures.py
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, "analysis/sim")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from catalog import setup_logging
from combine_jackknife import load_accumulator, total_map
from geometry import (BRIDGE_HALF_X_FRAC, band_profile, reflect_symmetrize_map,
                      symmetrize_map, two_halo_template)
from jackknife import jackknife_error
from plot_separation_summary import assert_comparable, build_terms
from weighting_sensitivity import SEPARATIONS

from deconvolve_pair_profile import (BLUE, GREEN, GRID, INK, INK_MUTED,
                                     ORANGE, axis_for, load_map, tidy)

logger = logging.getLogger(__name__)

SIM_DIR = Path("analysis/sim/results")
SIM_SINGLE = "kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv"
KEYS = ("5", "10", "20")
ZOOM_FIG1 = 26.0

# Simulation counterpart for each BOSS separation, chosen to match that
# sample's line-of-sight cut.  This deliberately does NOT reuse
# deconvolve_pair_profile.PAIR_STACK: that maps "5" to the r_par <= 10 stack,
# which is right for the simulation-only work but wrong here, because the BOSS
# "5" sample is r_par <= 5.  Comparing the two would put a differently-cut
# theory panel next to the observation and call it the same pipeline.  The cut
# is read back out of each stack's sidecar and asserted below.
SIM_STACK = {"5": "stack_rperp5_matched",
             "10": "stack_rperp10_matched",
             "20": "stack_rperp20_matched"}

# Per-block accumulators for the mock's box jackknife.  Named separately
# because jackknife_pair_stack.py writes them without the "_matched" suffix
# the stack CSV carries, so they cannot be derived from SIM_STACK by string
# surgery.  sim_band_error() checks each one's pair count against the stack's
# sidecar, so a file from a different cut cannot be paired with a stack.
SIM_BLOCKS = {"5": "stack_rperp5_blocks.npz",
              "10": "stack_rperp10_blocks.npz",
              "20": "stack_rperp20_blocks.npz"}

# Set by --paper.  Module-level rather than threaded through every fig*()
# because it affects only presentation, never a number.
PAPER = False

# Set by --allow-rperp-mismatch.  Off by default: a mock stack built
# from a different transverse bin than BOSS is a different sample, and
# silently plotting the two together is what put a 19-21 mock panel
# beside BOSS's 18-22 in a committed figure.
ALLOW_RPERP_MISMATCH = False


def caption(fig, text):
    """Explanatory paragraph above a figure.  Suppressed under ``--paper``.

    Diagnostic runs keep it so a figure read on its own still explains itself.
    The paper carries the same wording in the LaTeX caption instead, because a
    figure with the text baked into the image cannot be re-worded during review
    and duplicates the caption underneath it -- see ``paper/README.md``.
    """
    if PAPER:
        return
    fig.suptitle(text, fontsize=10.5, color=INK, ha="left", x=0.01)


def axis_of(arr):
    return axis_for(arr.shape[0])


def crop(arr, axis, half):
    m = np.abs(axis) <= half
    return arr[np.ix_(m, m)], [-half, half, -half, half]


def boss_single(results_dir, dataset, regions, tag, random_tag):
    """Corrected, radially symmetrized single stack -- the control's parent."""
    jk = os.path.join(results_dir, "jk")
    label = "_".join(regions)
    gal = load_accumulator(os.path.join(
        jk, f"acc_single_galaxy{tag}_{dataset}_{label}.npz"))
    rnd = load_accumulator(os.path.join(
        jk, f"acc_single_random{random_tag}_{dataset}_{label}.npz"))
    g = gal["grid_size"]
    return symmetrize_map((total_map(gal["sum_wk"], gal["sum_w"])
                           - total_map(rnd["sum_wk"], rnd["sum_w"])).reshape(g, g))


def boss_rpar_max(key):
    """The BOSS sample's line-of-sight cut, from its pair-catalog stem."""
    return float(SEPARATIONS[key]["pair_catalog"].split("_")[0])


def sim_diff(key):
    """Simulation (a - b): pair stack minus superposed singles.

    Refuses to proceed if the simulation stack was cut at a different r_par
    from the BOSS sample it is about to be plotted against.
    """
    stem = SIM_STACK[key]
    with open(SIM_DIR / "hod_pairs" / f"{stem}.json", encoding="utf-8") as fh:
        meta = json.load(fh)
    want, got = boss_rpar_max(key), float(meta["rpar_max_hmpc"])
    if not np.isclose(want, got):
        raise ValueError(
            f"r_perp = {key}: BOSS uses r_par <= {want:g} but {stem} was cut at "
            f"r_par <= {got:g}. Pick the matching simulation stack.")
    assert_comparable(key, meta, stem,
                      allow_rperp_mismatch=ALLOW_RPERP_MISMATCH)
    logger.info("sep %-3s: sim %s, r_par <= %g, r_perp %g-%g, %d pairs", key,
                stem, got, meta["rperp_min_hmpc"], meta["rperp_max_hmpc"],
                int(meta["n_pairs"]))
    single = load_map(SIM_DIR / SIM_SINGLE)
    pair = load_map(SIM_DIR / "hod_pairs" / f"{stem}.csv")
    axis = axis_of(pair)
    gx, gy = np.meshgrid(axis, axis)
    return pair - two_halo_template(single, gx, gy, float(key)), axis


# --------------------------------------------------------------------------
def fig1(single, terms, out_dir):
    """Single stack and the three pair stacks, shared linear scale."""
    panels = [("(a) single galaxies", single, axis_of(single), None)]
    for key in KEYS:
        t = terms[key]
        panels.append((rf"({'bcd'[KEYS.index(key)]}) pairs, "
                       rf"$r_\perp$ = {key} $h^{{-1}}$Mpc",
                       t["physical"], t["axis"], t["sep"]))

    vmax = max(float(np.nanmax(crop(a, ax, ZOOM_FIG1)[0]))
               for _, a, ax, _ in panels) * 1e4

    fig, axes = plt.subplots(1, 4, figsize=(17.2, 4.4), constrained_layout=True)
    for ax, (lab, arr, axis, sep) in zip(axes, panels):
        a, ext = crop(arr, axis, ZOOM_FIG1)
        im = ax.imshow(a * 1e4, origin="lower", extent=ext, cmap="magma_r",
                       vmin=0.0, vmax=vmax, interpolation="nearest")
        if sep is None:
            ax.plot([0], [0], ls="none", marker="+", ms=11, mew=1.8, color=GREEN)
        else:
            ax.plot([-0.5 * sep, 0.5 * sep], [0, 0], ls="none", marker="+",
                    ms=11, mew=1.8, color=GREEN)
        ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")
        ax.set_title(lab, loc="left", color=INK, fontsize=11, pad=5)
    fig.colorbar(im, ax=axes, label=r"$\kappa$  [$10^{-4}$]", shrink=0.9)
    caption(
        fig,
        "Observed stacked convergence — BOSS CMASS galaxies on the Planck "
        "lensing map, North and South stacked jointly\n"
        "Every panel is the galaxy stack minus the matching random stack.  The "
        "single stack is radially symmetrized; the pair stacks are "
        "reflection-symmetrized about both axes.\n"
        "Green crosses mark the galaxies.  All four panels share one linear "
        "colour scale, so the single halo and the pairs can be compared "
        "directly.",
    )
    save(fig, out_dir / "boss_single_and_pair_maps")


def fig2(terms, out_dir):
    """(a - b) observed beside theoretical, one row per separation."""
    fig, axes = plt.subplots(len(KEYS), 2, figsize=(9.8, 4.3 * len(KEYS)),
                             constrained_layout=True, squeeze=False)
    for row, key in enumerate(KEYS):
        t = terms[key]
        sep = t["sep"]
        half = float(min(max(2.2 * sep, 16.0), 32.0))
        obs, ext = crop(t["filament"], t["axis"], half)
        d_sim, sim_axis = sim_diff(key)
        sim, _ = crop(d_sim, sim_axis, half)

        # One scale per row, set by the observed panel: the comparison is
        # obs against theory at fixed separation, and BOSS is the noisier of
        # the two, so letting it set the range keeps the comparison honest
        # rather than flattering the simulation.
        v = float(np.nanpercentile(np.abs(obs), 99.0)) * 1e4
        for col, (arr, lab) in enumerate((
                (obs, "observed"), (sim, "theoretical"))):
            ax = axes[row][col]
            im = ax.imshow(arr * 1e4, origin="lower", extent=ext, cmap="RdBu_r",
                           vmin=-v, vmax=v, interpolation="nearest")
            ax.plot([-0.5 * sep, 0.5 * sep], [0, 0], ls="none", marker="+",
                    ms=11, mew=1.8, color=INK)
            ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
            ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")
            ax.set_title(rf"$r_\perp$ = {key} $h^{{-1}}$Mpc — {lab}",
                         loc="left", color=INK, fontsize=11, pad=5)
        fig.colorbar(im, ax=axes[row], label=r"$\kappa$  [$10^{-4}$]",
                     shrink=0.88)
        logger.info("sep %-3s: obs 99th pct %.3f, sim 99th pct %.3f  [1e-4]",
                    key, v, float(np.nanpercentile(np.abs(sim), 99.0)) * 1e4)

    caption(
        fig,
        r"Difference between the pair stack and the single-based control, "
        r"$(a-b)$ — observation beside theory""\n"
        "Left: BOSS CMASS on Planck.  Right: the matched simulation mock, run "
        "through the identical pipeline at the same line-of-sight cut "
        r"($r_\parallel \leq$ 5, 10, 10 $h^{-1}$Mpc by row).""\n"
        "The control is the unfitted superposed single in both.  "
        "Each row shares one colour scale, set from the observed panel, so the "
        "two are directly comparable — the simulation is far less noisy, and "
        "looks correspondingly flat.",
    )
    save(fig, out_dir / "obs_vs_sim_quadrupole_maps")


def zoom(sep):
    """Half-width of the plotted X range -- the same rule fig2 uses."""
    return float(min(max(2.2 * sep, 16.0), 32.0))


def sim_band_error(key, axis, control):
    """Box-jackknife error on the mock band profile, from the per-block stacks.

    The control is passed in and held fixed across realizations rather than
    refitted, matching what jackknife_pair_stack.py does: it is built from the
    *single* stack, which these blocks do not resample, so refitting it here
    would inject a variance the data never had.

    Each leave-one-out map is rebuilt as (sum_total - sum_b) / (n_total - n_b)
    and reflection-symmetrized, because that is what the saved stack is -- the
    blocks reconstruct it to 3e-8 of peak only after that step.
    """
    path = SIM_DIR / "hod_pairs" / SIM_BLOCKS[key]
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Re-run jackknife_pair_stack.py for this "
            "separation with --blocks-output; see paper/README.md step 1.")
    d = np.load(path)
    sums, counts = d["sums"], d["counts"]

    # The blocks must belong to the same run as the stack they will be
    # differenced against -- stack_rperp5_rpar10_blocks.npz sits in the same
    # directory and is a different line-of-sight cut of the same separation.
    with open(SIM_DIR / "hod_pairs" / f"{SIM_STACK[key]}.json",
              encoding="utf-8") as fh:
        meta = json.load(fh)
    if int(d["n_pairs"]) != int(meta["n_pairs"]):
        raise ValueError(
            f"r_perp = {key}: {path.name} holds {int(d['n_pairs'])} pairs but "
            f"{SIM_STACK[key]}.json records {int(meta['n_pairs'])}. These are "
            "different runs; regenerate the blocks for this stack.")
    tot_s, tot_c = sums.sum(axis=0), counts.sum()
    loo = []
    for b in range(len(counts)):
        n = tot_c - counts[b]
        if n <= 0:
            continue
        m = reflect_symmetrize_map(((tot_s - sums[b]) / n).astype(np.float32))
        loo.append(band_profile(m - control, axis))
    logger.info("sep %-3s: mock band error from %d of %d blocks", key,
                len(loo), len(counts))
    return jackknife_error(np.asarray(loo))[1]


def fig3(terms, out_dir):
    """B(X) observed against the mock, with a residual strip beneath."""
    fig, axes = plt.subplots(2, len(KEYS), figsize=(4.6 * len(KEYS), 5.4),
                             sharex="col", height_ratios=[3, 1],
                             constrained_layout=True, squeeze=False)
    for col, key in enumerate(KEYS):
        t = terms[key]
        sep, axis = t["sep"], t["axis"]
        boss_fil = t["filament"]

        # Never build this error from per-pixel errors: the 8 arcmin beam is
        # ~3.3 h^-1 Mpc against a 1 h^-1 Mpc pixel, so neighbouring pixels are
        # correlated and the independent-bin approximation understates the
        # band error by 2.2-2.6x (see the paper's Section 3.6).
        boss_prof = band_profile(boss_fil, axis)
        loo = np.asarray([band_profile(m.reshape(boss_fil.shape), axis)
                          for m in t["loo"]["filament"]])
        boss_err = jackknife_error(loo)[1]

        sim_fil, sim_axis = sim_diff(key)
        single = load_map(SIM_DIR / SIM_SINGLE)
        gx, gy = np.meshgrid(sim_axis, sim_axis)
        control = two_halo_template(single, gx, gy, float(key))
        sim_prof = band_profile(sim_fil, sim_axis)
        sim_err = sim_band_error(key, sim_axis, control)

        half = zoom(sep)
        m = np.abs(axis) <= half
        ms = np.abs(sim_axis) <= half

        top = axes[0][col]
        top.axhline(0, color=INK, lw=0.9)
        top.axvspan(-BRIDGE_HALF_X_FRAC * sep, BRIDGE_HALF_X_FRAC * sep,
                    color=INK_MUTED, alpha=0.09, lw=0, zorder=0)
        for xx in (-0.5 * sep, 0.5 * sep):
            top.axvline(xx, color=GREEN, lw=1.1, ls=":")
        top.fill_between(sim_axis[ms], (sim_prof[ms] - sim_err[ms]) * 1e4,
                         (sim_prof[ms] + sim_err[ms]) * 1e4,
                         color=ORANGE, alpha=0.22, lw=0, zorder=1)
        top.plot(sim_axis[ms], sim_prof[ms] * 1e4, lw=2.0, color=ORANGE,
                 label="mock", zorder=2)
        top.errorbar(axis[m], boss_prof[m] * 1e4, yerr=boss_err[m] * 1e4,
                     fmt="o", ms=3.0, lw=0, elinewidth=1.0, capsize=0,
                     color=BLUE, ecolor=BLUE, label="BOSS", zorder=3)
        top.set_ylabel(r"$B(X)$  [$10^{-4}$]")
        top.set_title(rf"$r_\perp$ = {key} $h^{{-1}}$Mpc  ({t['n_pairs']:,} pairs)",
                      loc="left", color=INK, fontsize=11, pad=5)
        if col == 0:
            top.legend(loc="upper right", fontsize=8.6, frameon=False)
        tidy(top)

        # The mock is on the pair grid and BOSS may not be, so interpolate the
        # mock onto the observed axis rather than assuming they line up.
        bot = axes[1][col]
        resid = (boss_prof[m] - np.interp(axis[m], sim_axis, sim_prof)) / boss_err[m]
        bot.axhline(0, color=INK, lw=0.9)
        for lev in (-2.0, 2.0):
            bot.axhline(lev, color=INK_MUTED, lw=0.7, ls="--")
        bot.plot(axis[m], resid, lw=0, marker="o", ms=2.6, color=BLUE)
        bot.set_ylim(-4.0, 4.0)
        bot.set_xlim(-half, half)
        bot.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
        bot.set_ylabel(r"$\Delta/\sigma$")
        bot.grid(True, axis="y", color=GRID, linewidth=0.6)
        bot.set_axisbelow(True)
        tidy(bot)

        n_out = int(np.sum(np.abs(resid) > 2.0))
        logger.info("sep %-3s: %d of %d residual bins beyond 2 sigma", key,
                    n_out, resid.size)

    caption(
        fig,
        r"Band profile $B(X)$ — BOSS against the mock""\n"
        "Points are the observed central-minus-off-centre band with delete-one "
        "jackknife errors over 287 sky patches; the shaded band is the mock's "
        "25-block box jackknife.\n"
        "Grey marks the bridge window, dotted green the galaxies.  Lower "
        r"strip: (BOSS $-$ mock) in units of the observed error.""\n"
        "Mock pairs are selected in the same transverse bins as BOSS at all "
        "three separations; the mock band is statistical only.",
    )
    save(fig, out_dir / "band_profiles_obs_vs_sim")


def save(fig, base):
    base.parent.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(f"{base}.{e}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Saved %s.png / .pdf", base)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Section 4.2 observed figures.")
    ap.add_argument("--results-dir", default="analysis/boss/results")
    ap.add_argument("--dataset", default="BOSS")
    ap.add_argument("--regions", default="North,South")
    ap.add_argument("--single-tag", default="_scw")
    ap.add_argument("--single-random-tag", default="_scw_frac100")
    ap.add_argument("--output-dir", default="output/plots")
    ap.add_argument("--only", default="", help="subset of 1,2,3")
    ap.add_argument("--allow-rperp-mismatch", action="store_true",
                    help="Plot a mock stack whose transverse bin differs from "
                         "the BOSS sample's. Provisional use only, while the "
                         "corrected pair catalogue is being generated; the "
                         "caption must carry the same caveat.")
    ap.add_argument("--paper", action="store_true",
                    help="Drop the explanatory paragraph above each figure. "
                         "Paper figures carry no text inside the image -- the "
                         "wording lives in the LaTeX caption instead. Panel "
                         "titles are labels and are kept either way.")
    args = ap.parse_args(argv)
    setup_logging()

    global PAPER, ALLOW_RPERP_MISMATCH
    PAPER = args.paper
    ALLOW_RPERP_MISMATCH = args.allow_rperp_mismatch
    regions = [r.strip() for r in args.regions.split(",")]
    out = Path(args.output_dir)

    single = boss_single(args.results_dir, args.dataset, regions,
                         args.single_tag, args.single_random_tag)
    terms = {}
    for key in KEYS:
        terms[key] = build_terms(args.results_dir, args.dataset, regions, key,
                                 args.single_tag, args.single_random_tag)
        logger.info("%s: %d pairs, a = %.1f", key, terms[key]["n_pairs"],
                    terms[key]["sep"])

    want = {s.strip() for s in args.only.split(",") if s.strip()} or {"1", "2", "3"}
    if "1" in want:
        fig1(single, terms, out)
    if "2" in want:
        fig2(terms, out)
    if "3" in want:
        fig3(terms, out)


if __name__ == "__main__":
    main()
