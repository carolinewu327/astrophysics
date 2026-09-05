#!/usr/bin/env python
"""The four theoretical-expectation figures for section 3.2 of the paper.

Zheng's list, one function each:

  fig1  true pair kappa map (a) beside the single-based control (b), at
        r_perp = 10 as the worked example
  fig2  the difference (a - b) at 5, 10, 20, to show the quadrupole
  fig3  central band, off-centre band, and their difference, taken from the
        (a - b) map, at 5, 10, 20
  fig4  the central band profile of the superposed singles, decomposed into
        core and tail

Everything is simulation.  The single stack is read on the *pair* grid
(101 px): reading it on the 100-px grid samples the band 2-4% low, because its
rows sit at |Y| = 0.5, 1.5 against the pair grid's 0, 1 and kappa falls off the
axis.  That bias is small on a map but shifts the fitted core amplitude by
0.02-0.04, so the matched-grid stack is the input here.

The control in fig1/fig2 is the plain superposed single -- two copies of the
measured single profile at +/- a/2, no fitted parameters.  That is what "single
-based control" means, and it is deliberately the *unfitted* control: fig2's
job is to show what is left over before any fitting, which is the motivation
for the core+tail model.

Usage
-----
    PYTHONPATH=lib:analysis/sim python analysis/sim/plot_theory_figures.py
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
from matplotlib.colors import LogNorm

from geometry import band_profiles, two_halo_template
from deconvolve_pair_profile import (GREEN, INK, INK_MUTED, PAIR_STACK,
                                     axis_for, load_map, tidy)
from stretched_single_control import HOD, SIM
from band_core_tail_control import (HALO_HALF, R_CORE_DEFAULT, band_spline,
                                    decompose_band, single_band_profile,
                                    superpose)

logger = logging.getLogger(__name__)

SINGLE = "kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv"
SEPS = ("5", "10", "20")

# Which mock stack stands for each separation *in the paper*.  This is NOT
# deconvolve_pair_profile.PAIR_STACK: that maps "5" to the r_par <= 10 stack,
# which is the right sample for the deconvolution work but not for the paper,
# where the BOSS 5 h^-1 Mpc sample is cut at r_par <= 5.  Section 4 states the
# LCDM expectation for the samples Section 5 measures, so the two must be the
# same stacks -- otherwise the figure and the number quoted beside it disagree
# (measured: 7.64e-4 plotted against 6.32e-4 quoted, 21 per cent apart).
PAPER_STACK = {"5": "stack_rperp5_matched",
               "10": "stack_rperp10_matched",
               "20": "stack_rperp20_matched"}
BLUE, ORANGE, PURPLE = "#2a78d6", "#eb6834", "#7d5ba6"

# Set by --paper.  Module-level rather than threaded through every fig*()
# because it affects only presentation, never a number.
PAPER = False


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


def zoom_for(sep):
    """Half-width of the plotted region.

    Capped well inside the 50 h^-1 Mpc box half-width: a shifted control
    template reaches the map edge before the map does, so the outermost few
    h^-1 Mpc are edge effects rather than signal and only distract.
    """
    return float(min(max(2.2 * sep, 16.0), 32.0))


def load_pair(key):
    return load_map(HOD / f"{PAPER_STACK[key]}.csv")


def control_map(single, axis, sep):
    """Plain superposed singles on the pair grid -- no fitted parameters."""
    gx, gy = np.meshgrid(axis, axis)
    return two_halo_template(single, gx, gy, sep)


def crop(arr, axis, half):
    m = np.abs(axis) <= half
    return arr[np.ix_(m, m)], [-half, half, -half, half]


# --------------------------------------------------------------------------
def fig1(single, axis, out_dir, key="10"):
    """(a) true pair map and (b) single-based control, one separation."""
    sep = float(key)
    pair = load_pair(key)
    ctl = control_map(single, axis, sep)
    half = zoom_for(sep)
    a, ext = crop(pair, axis, half)
    b, _ = crop(ctl, axis, half)

    lo = max(min(np.nanmin(a), np.nanmin(b)), 1e-5)
    hi = max(np.nanmax(a), np.nanmax(b))
    norm = LogNorm(vmin=lo * 1e4, vmax=hi * 1e4)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.5), constrained_layout=True)
    for ax, arr, lab in ((axes[0], a, "(a) true pair stack"),
                         (axes[1], b, "(b) single-based control")):
        im = ax.imshow(arr * 1e4, origin="lower", extent=ext, cmap="magma_r",
                       norm=norm)
        ax.plot([-0.5 * sep, 0.5 * sep], [0, 0], ls="none", marker="+",
                ms=11, mew=1.8, color=GREEN)
        ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")
        ax.set_title(lab, loc="left", color=INK, fontsize=11, pad=5)
    fig.colorbar(im, ax=axes, label=r"$\kappa$  [$10^{-4}$]", shrink=0.9)
    caption(
        fig,
        rf"Simulated pair stack and its single-based control, $r_\perp$ = {key} "
        r"$h^{-1}$Mpc""\n"
        "The control is two copies of the measured single-galaxy profile placed "
        r"at $\pm a/2$, with no fitted parameters.""\n"
        "Green crosses mark the galaxies.  Shared logarithmic colour scale, so "
        "the two panels are directly comparable — and nearly indistinguishable "
        "by eye, which is the point: the filament is a small residual of two "
        "large maps.",
    )
    save(fig, out_dir / f"sim_pair_vs_control_maps_sep{key}")


def fig2(single, axis, out_dir):
    """(a - b) at all three separations: the quadrupole."""
    fig, axes = plt.subplots(1, len(SEPS), figsize=(4.9 * len(SEPS), 4.4),
                             constrained_layout=True)
    for ax, key in zip(np.atleast_1d(axes), SEPS):
        sep = float(key)
        d = load_pair(key) - control_map(single, axis, sep)
        arr, ext = crop(d, axis, zoom_for(sep))
        v = float(np.nanpercentile(np.abs(arr), 99.0)) * 1e4
        im = ax.imshow(arr * 1e4, origin="lower", extent=ext, cmap="RdBu_r",
                       vmin=-v, vmax=v)
        ax.plot([-0.5 * sep, 0.5 * sep], [0, 0], ls="none", marker="+",
                ms=11, mew=1.8, color=INK)
        ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")
        ax.set_title(rf"$r_\perp$ = {key} $h^{{-1}}$Mpc", loc="left",
                     color=INK, fontsize=11, pad=5)
        fig.colorbar(im, ax=ax, label=r"$\kappa$  [$10^{-4}$]", shrink=0.88)
    caption(
        fig,
        "Difference between the pair stack and the single-based control, "
        r"$(a-b)$ — simulation""\n"
        "Positive along the pair axis, negative above and below it: the "
        "quadrupole left behind when two isolated halo profiles are subtracted "
        "from a real pair.\n"
        "Crosses mark the galaxies.  Colour scale is symmetric about zero and "
        "clipped at the 99th percentile of each panel.",
    )
    save(fig, out_dir / "sim_quadrupole_diff_maps")


def fig3(single, axis, out_dir):
    """Central band, off-centre band, and difference, from the (a - b) map."""
    fig, axes = plt.subplots(1, len(SEPS), figsize=(4.9 * len(SEPS), 4.0),
                             constrained_layout=True, sharey=False)
    for ax, key in zip(np.atleast_1d(axes), SEPS):
        sep = float(key)
        d = load_pair(key) - control_map(single, axis, sep)
        cen, off = band_profiles(d, axis)
        half = zoom_for(sep)
        m = np.abs(axis) <= half
        ax.plot(axis[m], cen[m] * 1e4, lw=2.1, color=BLUE,
                label=r"central band  $|Y| \leq 1.5$")
        ax.plot(axis[m], off[m] * 1e4, lw=2.1, color=ORANGE,
                label=r"off-centre band  $1.5 \leq |Y| \leq 10.5$")
        ax.plot(axis[m], (cen - off)[m] * 1e4, lw=2.4, color=INK,
                label="difference")
        ax.axhline(0, color=INK_MUTED, lw=0.9)
        for xx in (-0.5 * sep, 0.5 * sep):
            ax.axvline(xx, color=GREEN, lw=1.2, ls=":")
        ax.axvspan(-0.5 * sep, 0.5 * sep, color=INK_MUTED, alpha=0.055, lw=0)
        ax.set_xlim(-half, half)
        ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$\kappa$  [$10^{-4}$]")
        ax.set_title(rf"$r_\perp$ = {key} $h^{{-1}}$Mpc", loc="left",
                     color=INK, fontsize=11, pad=5)
        if key == SEPS[0]:
            ax.legend(loc="lower left", fontsize=8.2, frameon=False)
        tidy(ax)
    caption(
        fig,
        r"Band compression of the $(a-b)$ difference map — simulation""\n"
        "The central band runs along the pair axis where a filament would lie; "
        "the off-centre band brackets it above and below.\n"
        "Their difference cancels anything that varies slowly with $Y$, "
        "including the large-scale offset the two bands share.",
    )
    save(fig, out_dir / "sim_band_profiles_all_seps")


def fig4(single, axis, out_dir, r_c=R_CORE_DEFAULT):
    """Band-space core/tail split, and what it looks like superposed."""
    x_s, b_s = single_band_profile(single)
    bfun, x_max = band_spline(x_s, b_s)
    core, tail, sigma, _ = decompose_band(bfun, r_c)
    logger.info("band decomposition at r_c = %.1f: sigma = %.3f h^-1 Mpc",
                r_c, sigma)

    fig, axes = plt.subplots(1, 1 + len(SEPS),
                             figsize=(4.4 * (1 + len(SEPS)), 4.1),
                             constrained_layout=True)

    q = np.linspace(0.0, 30.0, 601)
    ax = axes[0]
    ax.plot(q, bfun(q) * 1e4, lw=2.4, color=INK, label="single, central band")
    ax.plot(q, core(q) * 1e4, lw=2.0, color=PURPLE, ls=(0, (5, 2)),
            label=r"core  ($r \leq r_c$, then Gaussian)")
    ax.plot(q, tail(q) * 1e4, lw=2.0, color=ORANGE, label="tail  (remainder)")
    ax.axvline(r_c, color=GREEN, lw=1.2, ls=":")
    ax.annotate(rf"$r_c$ = {r_c:g}", xy=(r_c, 0.55),
                xycoords=("data", "axes fraction"), xytext=(5, 0),
                textcoords="offset points", fontsize=8.6, color=INK_MUTED)
    ax.axhline(0, color=INK_MUTED, lw=0.9)
    ax.set_xlim(0, 30)
    ax.set_xlabel(r"$|X|$  [$h^{-1}$ Mpc]")
    ax.set_ylabel(r"central band $\kappa$  [$10^{-4}$]")
    ax.set_title("one galaxy, decomposed", loc="left", color=INK,
                 fontsize=11, pad=5)
    ax.legend(loc="upper right", fontsize=8.2, frameon=False)
    tidy(ax)

    for ax, key in zip(axes[1:], SEPS):
        sep = float(key)
        c0, t0 = float(core(0.5 * sep)), float(tail(0.5 * sep))
        share = 100.0 * c0 / (c0 + t0) if (c0 + t0) else np.nan
        logger.info("a = %4.1f: at the bridge centre the core supplies "
                    "%5.1f%% of the control", sep, share)
        seps, wts = np.array([sep]), np.array([1.0])
        tot = superpose(bfun, axis, seps, wts)
        c2 = superpose(core, axis, seps, wts)
        t2 = superpose(tail, axis, seps, wts)
        half = zoom_for(sep)
        m = (np.abs(axis) <= half) & np.isfinite(tot)
        ax.plot(axis[m], tot[m] * 1e4, lw=2.4, color=INK, label="superposed")
        ax.plot(axis[m], c2[m] * 1e4, lw=2.0, color=PURPLE, ls=(0, (5, 2)),
                label="core part")
        ax.plot(axis[m], t2[m] * 1e4, lw=2.0, color=ORANGE, label="tail part")
        ax.axhline(0, color=INK_MUTED, lw=0.9)
        for xx in (-0.5 * sep, 0.5 * sep):
            ax.axvline(xx, color=GREEN, lw=1.2, ls=":")
        for sgn in (-1, 1):
            ax.axvspan(sgn * (0.5 * sep - HALO_HALF),
                       sgn * (0.5 * sep + HALO_HALF),
                       color=GREEN, alpha=0.10, lw=0)
        ax.set_xlim(-half, half)
        ax.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"central band $\kappa$  [$10^{-4}$]")
        ax.set_title(rf"superposed, $r_\perp$ = {key} $h^{{-1}}$Mpc"
                     "\n"rf"core supplies {share:.0f}% at $X=0$",
                     loc="left", color=INK, fontsize=10.5, pad=5)
        if key == SEPS[0]:
            ax.legend(loc="upper right", fontsize=8.2, frameon=False)
        tidy(ax)

    caption(
        fig,
        "Central band profile of the single-based control, and its core/tail "
        f"split (matched core, $r_c$ = {r_c:g}, giving $\\sigma$ = {sigma:.2f} "
        r"$h^{-1}$Mpc)""\n"
        "Left: one galaxy's central band profile, taken from the symmetrized "
        "single map, split into a compact core and the remaining tail.\n"
        "Right: each piece placed at both galaxy positions.  The core is "
        "compact, so at 20 it has died away before mid-point and the bridge is "
        r"the tail's; at 5 the galaxies are closer together than the core is""\n"
        "wide, so the core supplies the whole of it.  Since the core amplitude "
        "is pinned by the halo peaks, nothing is left free to describe the "
        r"bridge at $r_\perp$ = 5 — which is why that separation fails.",
    )
    save(fig, out_dir / "band_core_tail_decomposition")


def save(fig, base):
    base.parent.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(f"{base}.{e}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Saved %s.png / .pdf", base)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Section 3.2 figures.")
    ap.add_argument("--single", default=SINGLE)
    ap.add_argument("--example-sep", default="10")
    ap.add_argument("--r-core", type=float, default=R_CORE_DEFAULT)
    ap.add_argument("--output-dir", default="output/plots")
    ap.add_argument("--only", default="",
                    help="comma-separated subset of 1,2,3,4")
    ap.add_argument("--paper", action="store_true",
                    help="Drop the explanatory paragraph above each figure. "
                         "Paper figures carry no text inside the image -- the "
                         "wording lives in the LaTeX caption instead. Panel "
                         "titles are labels and are kept either way.")
    args = ap.parse_args(argv)

    global PAPER
    PAPER = args.paper
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    single = load_map(SIM / args.single)
    axis = axis_for(single.shape[0])
    out = Path(args.output_dir)
    logger.info("single %s: grid %d", args.single, single.shape[0])
    if single.shape[0] != 101:
        logger.warning("single is not on the 101-px pair grid -- the control "
                       "will carry the band-sampling bias")

    want = {s.strip() for s in args.only.split(",") if s.strip()} or \
        {"1", "2", "3", "4"}
    if "1" in want:
        fig1(single, axis, out, args.example_sep)
    if "2" in want:
        fig2(single, axis, out)
    if "3" in want:
        fig3(single, axis, out)
    if "4" in want:
        fig4(single, axis, out, args.r_core)


if __name__ == "__main__":
    main()
