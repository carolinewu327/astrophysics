#!/usr/bin/env python
"""Two figures for the joint-jackknife results, plus the log-bin cross-check table.

Figure 1 -- corrected kappa profile, BOSS vs the density-matched mock, with the
ratio underneath.  This is the "lensing is low" result: BOSS sits below the mock
at every radius above a few h^-1 Mpc.  Both binnings are drawn so the log-bin
cross-check is visible against the 1 h^-1 Mpc bins it is meant to check.

Figure 2 -- where the 20 h^-1 Mpc dip came from.  North, South, the old 50/50
average and the joint stack on one axis.  The dip exists only in the average.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/plot_jackknife_results.py
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
from combine_jackknife import radial_profiles

logger = logging.getLogger(__name__)

# Reference categorical palette, slots in fixed order (light mode).
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK_2, INK_MUTED = "#0b0b0b", "#52514e", "#8a8985"
GRID = "#e3e2df"

plt.rcParams.update({
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "font.size": 9,
    "axes.labelsize": 9.5,
    "axes.titlesize": 10,
    "axes.edgecolor": INK_MUTED,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.8,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def style_axes(ax) -> None:
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def figure_profile(R: str, mock_path: str, out_base: str) -> None:
    lin = pd.read_csv(
        os.path.join(R, "profile_corrected_scw_rand_scw_frac100_BOSS_North_South_joint_linear.csv"))
    log = pd.read_csv(
        os.path.join(R, "profile_corrected_scw_rand_scw_frac100_BOSS_North_South_joint_log10.csv"))
    mock_map = pd.read_csv(mock_path, index_col=0).to_numpy(dtype=float)

    def mock_on(edges):
        prof, _ = radial_profiles(mock_map.reshape(1, -1), edges, mock_map.shape[0])
        return prof[0]

    lin_edges = np.append(lin.r_lo_hmpc.to_numpy(), lin.r_hi_hmpc.to_numpy()[-1])
    log_edges = np.append(log.r_lo_hmpc.to_numpy(), log.r_hi_hmpc.to_numpy()[-1])
    mock_lin, mock_log = mock_on(lin_edges), mock_on(log_edges)

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.2, 5.4), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.08})

    sel = lin.r_center_hmpc >= 2
    ax.errorbar(lin.r_center_hmpc[sel], lin.kappa[sel] * 1e4,
                yerr=lin.kappa_err[sel] * 1e4, fmt="o", ms=3, lw=0,
                elinewidth=0.8, color=BLUE, alpha=0.30, zorder=2)
    ax.errorbar(log.r_center_hmpc, log.kappa * 1e4, yerr=log.kappa_err * 1e4,
                fmt="o", ms=8, lw=0, elinewidth=2, capsize=0, color=BLUE,
                markeredgecolor="white", markeredgewidth=1.2, zorder=4,
                label="BOSS CMASS (log bins)")
    ax.plot(lin.r_center_hmpc[sel], mock_lin[sel] * 1e4, lw=2, color=ORANGE,
            zorder=3, label="Density-matched mock")

    ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.5, zorder=1)
    ax.set_xscale("log")
    ax.set_ylabel(r"corrected $\kappa$  [$10^{-4}$]")
    ax.set_title("BOSS CMASS convergence sits below the mock above a few $h^{-1}$Mpc",
                 loc="left", color=INK, pad=8)
    ax.legend(loc="upper right")
    ax.text(0.985, 0.62, "faint points: 1 $h^{-1}$Mpc bins",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK_MUTED)
    style_axes(ax)

    ratio = log.kappa.to_numpy() / mock_log
    ratio_err = log.kappa_err.to_numpy() / mock_log
    axr.axhline(1.0, color=ORANGE, lw=2, zorder=2)
    axr.errorbar(log.r_center_hmpc, ratio, yerr=ratio_err, fmt="o", ms=8, lw=0,
                 elinewidth=2, color=BLUE, markeredgecolor="white",
                 markeredgewidth=1.2, zorder=3)
    axr.set_xscale("log")
    axr.set_ylim(0, 1.6)
    axr.set_xlabel(r"$r$  [$h^{-1}$ Mpc]")
    axr.set_ylabel("BOSS / mock")
    axr.set_xlim(1.8, 56)
    axr.set_xticks([2, 3, 5, 10, 20, 30, 50])
    axr.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    style_axes(axr)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s.png / .pdf", out_base)


def figure_dip(R: str, out_base: str) -> None:
    d = pd.read_csv(os.path.join(R, "profile_region_split_scw_BOSS_North_South.csv"))
    win = (d.r_center_hmpc >= 10) & (d.r_center_hmpc <= 35)
    d = d[win]
    r = d.r_center_hmpc.to_numpy()

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    series = [
        ("North", d.North.to_numpy() * 1e4, BLUE, "-"),
        ("South", d.South.to_numpy() * 1e4, ORANGE, "-"),
        ("Old 50/50 average", d.unweighted_avg.to_numpy() * 1e4, YELLOW, "--"),
        ("Joint stack", d.joint.to_numpy() * 1e4, AQUA, "-"),
    ]
    # Error bands on the two combination rules that matter: the South (whose
    # noise drives the artifact) and the joint stack (the corrected answer).
    ax.fill_between(r, (d.South - d.South_err) * 1e4, (d.South + d.South_err) * 1e4,
                    color=ORANGE, alpha=0.12, lw=0, zorder=1)
    ax.fill_between(r, (d.joint - d.joint_err) * 1e4, (d.joint + d.joint_err) * 1e4,
                    color=AQUA, alpha=0.16, lw=0, zorder=1)
    for name, y, color, ls in series:
        ax.plot(r, y, lw=2, color=color, linestyle=ls, zorder=3, label=name)
    # No direct labels here: the four curves converge on the right edge, so they
    # would collide. The legend carries identity, and the dashed 50/50 line adds
    # a non-color channel.

    ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6, zorder=2)
    ax.axvspan(19, 24, color=INK_MUTED, alpha=0.07, lw=0, zorder=0)

    i = int(np.argmin(np.where((r >= 18) & (r <= 25), d.South.to_numpy(), np.inf)))
    ax.annotate(
        f"South: {d.South.to_numpy()[i] * 1e4:.1f} ± {d.South_err.to_numpy()[i] * 1e4:.1f}"
        f"\n({d.South.to_numpy()[i] / d.South_err.to_numpy()[i]:.1f}$\\sigma$ — ordinary noise)",
        (r[i], d.South.to_numpy()[i] * 1e4), xytext=(24, -35),
        textcoords="offset points", fontsize=8, color=INK_2,
        arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8))

    ax.set_xlim(10, 35)
    ax.set_xlabel(r"$r$  [$h^{-1}$ Mpc]")
    ax.set_ylabel(r"corrected $\kappa$  [$10^{-4}$]")
    ax.set_title("The 20 $h^{-1}$Mpc dip appears only in the 50/50 average",
                 loc="left", color=INK, pad=8)
    ax.legend(loc="upper right", ncol=2)
    style_axes(ax)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s.png / .pdf", out_base)


def log_bin_table(R: str, mock_path: str) -> pd.DataFrame:
    mock_map = pd.read_csv(mock_path, index_col=0).to_numpy(dtype=float)
    frames = []
    for label in ("log6", "log10"):
        p = os.path.join(
            R, f"profile_corrected_scw_rand_scw_frac100_BOSS_North_South_joint_{label}.csv")
        d = pd.read_csv(p)
        edges = np.append(d.r_lo_hmpc.to_numpy(), d.r_hi_hmpc.to_numpy()[-1])
        mock, _ = radial_profiles(mock_map.reshape(1, -1), edges, mock_map.shape[0])
        d = d.assign(binning=label, mock_kappa=mock[0],
                     boss_over_mock=d.kappa / mock[0])
        frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(os.path.join(R, "log_bin_cross_check_BOSS.csv"), index=False)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description="Figures for the joint jackknife results.")
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--mock",
                        default="analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin.csv")
    parser.add_argument("--output-dir", default="output/plots")
    args = parser.parse_args(argv)
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    figure_profile(args.results_dir, args.mock,
                   os.path.join(args.output_dir, "jk_profile_boss_vs_mock"))
    figure_dip(args.results_dir, os.path.join(args.output_dir, "jk_dip_region_split"))

    table = log_bin_table(args.results_dir, args.mock)
    for label, grp in table.groupby("binning", sort=False):
        logger.info("\n%s bins (2-50 h^-1 Mpc), corrected kappa x 1e4:", label)
        logger.info("   r_lo   r_hi   r_cen     BOSS      err     mock   BOSS/mock")
        for _, x in grp.iterrows():
            logger.info("  %5.2f  %5.2f  %6.2f  %7.3f  %7.3f  %7.3f  %9.2f",
                        x.r_lo_hmpc, x.r_hi_hmpc, x.r_center_hmpc, x.kappa * 1e4,
                        x.kappa_err * 1e4, x.mock_kappa * 1e4, x.boss_over_mock)


if __name__ == "__main__":
    main()
