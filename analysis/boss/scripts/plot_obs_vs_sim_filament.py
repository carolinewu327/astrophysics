#!/usr/bin/env python
"""Observation vs simulation filament comparison, one sheet per separation.

  (a) BOSS filament map
  (b) simulation filament map -- matched HOD mock, matched smoothing, matched LOS cut
  (c) BOSS minus simulation
  (d) filament profile, BOSS vs simulation (central band minus far band, along the
      pair axis -- the same statistic as panel (d) of the BOSS summary sheets)

Matching, which is the whole point of the comparison:
  * mock          density-matched HOD (hodnmatch), not a halo-mass cut
  * smoothing     the simulation kappa map is pre-smoothed to 3.31 h^-1 Mpc FWHM,
                  the equivalent of Planck's 8 arcmin at z = 0.55
  * LOS cut       the simulation is re-stacked with the *redshift-space* cut that
                  BOSS applies (r_par,rsd <= 5 at sep 5, <= 10 at sep 10 and 20),
                  since a redshift survey only ever sees the RSD-displaced
                  separation.  The archived deep50/truth20 stacks used wider or
                  real-space cuts and are not matched.
  * control       both sides subtract two copies of their own single stack, built
                  through geometry.two_halo_template on the pair grid

The simulation curve carries no error band: its box is 2500 h^-1 Mpc and the
bridge excess is detected at ~20 sigma, so its uncertainty is small next to the
BOSS errors.  Its jackknife error is printed in the annotation.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/plot_obs_vs_sim_filament.py
"""

from __future__ import annotations

import argparse
import json
import logging
import os

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from catalog import setup_logging
from geometry import two_halo_template
from jackknife import jackknife_error
from plot_separation_summary import (
    BRIDGE_X_FRAC,
    band_profile,
    build_terms,
    draw_map,
)
from weighting_sensitivity import SEPARATIONS, axis_for_map, load_map

logger = logging.getLogger(__name__)

BLUE, ORANGE, INK, INK_2, INK_MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8985"
GRID = "#e3e2df"

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def sim_filament(sim_dir: str, sep_key: str, single_name: str) -> tuple[np.ndarray, dict]:
    """Simulation filament map: matched pair stack minus its own two-halo control."""
    sep = SEPARATIONS[sep_key]["center"]
    stack_path = os.path.join(sim_dir, "hod_pairs", f"stack_rperp{sep_key}_matched.csv")
    meta_path = os.path.join(sim_dir, "hod_pairs", f"stack_rperp{sep_key}_matched.json")
    if not os.path.exists(stack_path):
        raise FileNotFoundError(
            f"{stack_path} not found. Re-stack the simulation at the matched LOS cut "
            "with jackknife_pair_stack.py before running this script.")

    pair_stack = load_map(stack_path)
    # Simulation singles are symmetrized when produced.  Loading them verbatim
    # lets two_halo_template reject archived half-pixel-centered products rather
    # than silently re-symmetrizing already-corrupted radial bins.
    single = load_map(os.path.join(sim_dir, single_name))
    axis = axis_for_map(pair_stack)
    gx, gy = np.meshgrid(axis, axis)
    control = two_halo_template(single, gx, gy, sep)

    with open(meta_path, encoding="utf-8") as fh:
        meta = json.load(fh)
    return pair_stack - control, meta


def make_sheet(t: dict, sim_fil: np.ndarray, sim_meta: dict, out_base: str) -> None:
    sep, axis = t["sep"], t["axis"]
    zoom = 2.5 * sep
    boss_fil = t["filament"]

    if sim_fil.shape != boss_fil.shape:
        raise ValueError(
            f"Grid mismatch: BOSS {boss_fil.shape} vs simulation {sim_fil.shape}. "
            "Both should be the 101-pixel pair grid over a 100 h^-1 Mpc box.")
    diff = boss_fil - sim_fil

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.0), constrained_layout=True)
    shared = max(np.nanpercentile(np.abs(boss_fil), 99.5),
                 np.nanpercentile(np.abs(sim_fil), 99.5), 1e-12)
    norm = TwoSlopeNorm(vmin=-shared, vcenter=0.0, vmax=shared)

    for ax, arr, title in (
        (axes[0, 0], boss_fil, "(a) BOSS filament map"),
        (axes[0, 1], sim_fil, "(b) Simulation filament map (matched mock, cut, smoothing)"),
        (axes[1, 0], diff, "(c) BOSS − simulation"),
    ):
        im = draw_map(ax, arr, axis, sep, title, norm, zoom,
                      show_bands=(ax is axes[1, 0]))
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(r"$\kappa$")
        cb.ax.tick_params(labelsize=7)

    ax = axes[1, 1]
    boss_prof = band_profile(boss_fil, axis) * 1e4
    loo = np.array([band_profile(m.reshape(boss_fil.shape), axis)
                    for m in t["loo"]["filament"]])
    _, boss_err = jackknife_error(loo)
    boss_err = boss_err * 1e4
    sim_prof = band_profile(sim_fil, axis) * 1e4

    ax.plot(axis, boss_prof, lw=2, color=BLUE, label="BOSS", zorder=3)
    ax.fill_between(axis, boss_prof - boss_err, boss_prof + boss_err, color=BLUE,
                    alpha=0.15, lw=0, zorder=1)
    ax.plot(axis, sim_prof, lw=2, color=ORANGE, label="Simulation", zorder=4)

    ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6)
    for xx in (-sep / 2, sep / 2):
        ax.axvline(xx, color=INK_MUTED, lw=0.8, ls="--", alpha=0.7)
    ax.axvspan(-BRIDGE_X_FRAC * sep, BRIDGE_X_FRAC * sep, color=INK_MUTED, alpha=0.07,
               lw=0, zorder=0)
    ax.set_xlim(-zoom, zoom)
    ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
    ax.set_ylabel(r"central band − off-centre band,  $\kappa$  [$10^{-4}$]")
    ax.set_title("(d) Filament profile, observation vs simulation", loc="left",
                 color=INK, pad=5)
    ax.legend(loc="upper right")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    window = np.abs(axis) <= BRIDGE_X_FRAC * sep
    boss_bridge = float(band_profile(boss_fil, axis)[window].mean())
    sim_bridge = float(band_profile(sim_fil, axis)[window].mean())
    # Jackknife the *windowed average* directly. Taking the RMS of the per-x
    # errors instead would ignore that averaging over x reduces the error, and
    # would not match the bridge excess error in the filament table.
    loo_bridge = np.array([band_profile(m.reshape(boss_fil.shape), axis)[window].mean()
                           for m in t["loo"]["filament"]])
    _, boss_bridge_err = jackknife_error(loo_bridge[:, None])
    boss_bridge_err = float(boss_bridge_err[0])
    # The archived simulation statistic still uses the older separation-scaled Y
    # bands (analysis/sim/ has not been migrated to the fixed physical bands), so
    # it is NOT the same measurement as sim_bridge above and must not be read as
    # a 2D-vs-1D discrepancy. It is quoted only for its jackknife error, which
    # sets the scale of the simulation's own uncertainty.
    sim_stat = sim_meta["stats"]["residual_bridge_excess_kappa"]
    ax.text(0.02, 0.03,
            f"bridge excess, $|X|\\leq${BRIDGE_X_FRAC * sep:.1f}:\n"
            f"  BOSS  {boss_bridge * 1e4:+.2f} ± {boss_bridge_err * 1e4:.2f}\n"
            f"  sim   {sim_bridge * 1e4:+.2f}  "
            f"(2D box, sep-scaled bands: {sim_stat['value'] * 1e4:+.2f} "
            f"± {sim_stat['jackknife_error'] * 1e4:.2f})",
            transform=ax.transAxes, fontsize=7, color=INK_MUTED, va="bottom",
            family="monospace")

    fig.suptitle(
        f"BOSS vs simulation filament — separation {sep:g} $h^{{-1}}$ Mpc\n"
        f"BOSS: {t['n_pairs']:,} pairs, jackknife over 287 cells.   "
        f"Simulation: {sim_meta['n_pairs']:,} pairs, "
        f"$r_{{\\rm par,rsd}}\\leq${sim_meta['rpar_max_hmpc']:g} $h^{{-1}}$Mpc, "
        f"density-matched HOD, 8 arcmin-equivalent smoothing.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s.png / .pdf  (BOSS bridge %+.2f, sim %+.2f)",
                out_base, boss_bridge * 1e4, sim_bridge * 1e4)
    return {"separation_hmpc": sep, "boss_bridge_excess": boss_bridge,
            "boss_bridge_err": boss_bridge_err, "sim_bridge_excess": sim_bridge,
            "sim_bridge_err": sim_stat["jackknife_error"],
            "boss_minus_sim": boss_bridge - sim_bridge,
            "boss_over_sim": boss_bridge / sim_bridge if sim_bridge else np.nan,
            "boss_pairs": t["n_pairs"], "sim_pairs": sim_meta["n_pairs"]}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Observation vs simulation filament maps and profiles.")
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument("--separations", default="5,10,20")
    parser.add_argument("--single-tag", default="_scw")
    parser.add_argument("--single-random-tag", default="_scw_frac100")
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--sim-dir", default="analysis/sim/results")
    parser.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin.csv")
    parser.add_argument("--output-dir", default="output/plots")
    args = parser.parse_args(argv)
    args.regions = [r.strip() for r in args.regions.split(",")]
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    rows = []
    for sep_key in (s.strip() for s in args.separations.split(",")):
        logger.info("=== separation %s h^-1 Mpc ===", sep_key)
        t = build_terms(args.results_dir, args.dataset, args.regions, sep_key,
                        args.single_tag, args.single_random_tag)
        sim_fil, sim_meta = sim_filament(args.sim_dir, sep_key, args.sim_single)
        rows.append(make_sheet(
            t, sim_fil, sim_meta,
            os.path.join(args.output_dir, f"boss_vs_sim_filament_sep{sep_key}")))

    out = pd.DataFrame(rows)
    path = os.path.join(args.results_dir, f"boss_vs_sim_filament_{args.dataset}.csv")
    out.to_csv(path, index=False)
    logger.info("Saved -> %s", path)
    logger.info("\nBridge excess, kappa x 1e4:")
    logger.info("  sep     BOSS      err      sim      err   BOSS-sim   BOSS/sim")
    for _, r in out.iterrows():
        logger.info("  %4.0f  %7.2f  %7.2f  %7.2f  %7.2f  %9.2f  %9.2f",
                    r.separation_hmpc, r.boss_bridge_excess * 1e4, r.boss_bridge_err * 1e4,
                    r.sim_bridge_excess * 1e4, r.sim_bridge_err * 1e4,
                    r.boss_minus_sim * 1e4, r.boss_over_sim)


if __name__ == "__main__":
    main()
