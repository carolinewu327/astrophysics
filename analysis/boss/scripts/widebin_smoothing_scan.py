#!/usr/bin/env python
"""Filament S/N versus Planck smoothing, for the wide bins (manuscript Figure 8).

The wide-bin counterpart of plot_smoothing_scan.py, and the source of
paper/sts/Figure/filament_snr_vs_smoothing.pdf.  One panel per bin with B(X) at
each smoothing, plus a panel of bridge-excess S/N.

Per smoothing, the galaxy single stack and the galaxy pair stacks are re-made
on the re-smoothed map:

    single   jk/acc_single_galaxy_scw{,_fwhm6,_fwhm4}_BOSS_North_South.npz
    pairs    widebin_rpar5/jk/acc_pairs_galaxy_g*_rpar5_scw{,_fwhm6,_fwhm4}_*.npz
             (stack_pairs_jk.py --fwhm-arcmin F --label galaxy_g*_rpar5_scw_fwhmF)

The random single stack and the random-pair sub-bin run are the 8 arcmin ones at
every smoothing, as in the submitted figure: the random-pair map moves the
bridge by ~0.3e-4 at 8 arcmin, far below the errors here.  The control is the
separation-averaged one throughout.

Writes the figure (PDF only) and widebin_rpar5/smoothing_scan_wide_BOSS.csv.
The 8 arcmin row must reproduce filament_jackknife_wide_avgctrl_*.csv.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python \\
        analysis/boss/scripts/widebin_smoothing_scan.py
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
from combine_filament_jackknife import (BRIDGE_CENTERS, grouping_filament, load_grouping_pairs,
                                        load_random_subbins, loo_and_total, pooled_random_map)
from combine_jackknife import load_accumulator
from geometry import BRIDGE_HALF_X_FRAC, band_profile, separation_mix
from jackknife import jackknife_error
from plot_smoothing_scan import GRID, INK, INK_MUTED, PALETTE
from weighting_sensitivity import axis_for_map

logger = logging.getLogger(__name__)

BINS = [("g4_6", "4–6"), ("g6_15", "6–15"), ("g15_25", "15–25")]
# (label, tag on the single and pair accumulators), coarsest first.
SMOOTHINGS = [("8 arcmin", ""), ("6 arcmin", "_fwhm6"), ("4 arcmin", "_fwhm4")]
# Legend text per smoothing; the CSV keeps the plain labels above.
LEGEND = {"8 arcmin": r"FWHM = $8'$", "6 arcmin": r"FWHM = $6'$", "4 arcmin": r"FWHM = $4'$"}
# Drawn at the width it prints (\linewidth, AASTeX text width 7.1 in), so the font sizes set in
# plot_smoothing_scan.py are the printed sizes.
FIGSIZE = (7.1, 5.7)


def measure(wide_dir, singles, rand, name, tag):
    """B(X), its error, and the bridge excess with its error for one bin at one smoothing."""
    st, sl, g = singles
    lo, hi = (float(v) for v in name[1:].split("_"))
    center = BRIDGE_CENTERS[name]
    pair = load_grouping_pairs(wide_dir / "jk" / f"acc_pairs_galaxy_{name}_rpar5_scw{tag}_BOSS_North_South.npz",
                               lo, hi)
    pt, pl = loo_and_total(pair)
    fine = (pair["fine_w"], pair["fine_wr"], pair["fine_cov"])

    def fil(p, s, mix):
        return grouping_filament(p, s, rand, mix, center, "separation_averaged", g, symmetrize=True)[0]

    res = fil(pt, st, separation_mix(*fine))
    axis = axis_for_map(res)
    loo = np.asarray([band_profile(fil(pl[k], sl[k], separation_mix(*fine, exclude_region=k)), axis)
                      for k in range(sl.shape[0])])
    prof = band_profile(res, axis)
    window = np.abs(axis) <= BRIDGE_HALF_X_FRAC * center
    bridge_err = jackknife_error(loo[:, window].mean(axis=1)[:, None])[1][0]
    return {"axis": axis, "prof": prof, "prof_err": jackknife_error(loo)[1], "center": center,
            "bridge": float(prof[window].mean()), "bridge_err": float(bridge_err),
            "n_pairs": int(pair["n_pairs"].sum())}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--wide-dir", default="analysis/boss/results/widebin_rpar5")
    p.add_argument("--random-label", default="random_w3_25_rpar5_scw_frac100")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    p.add_argument("--output-dir", default="paper/sts/Figure")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    setup_logging()

    wide_dir, out_dir = Path(args.wide_dir), Path(args.output_dir)
    fig_path = out_dir / "filament_snr_vs_smoothing.pdf"
    csv_path = wide_dir / "smoothing_scan_wide_BOSS.csv"
    existing = [str(f) for f in (fig_path, csv_path) if f.exists()]
    if existing and not args.overwrite:
        raise FileExistsError("Outputs exist (use --overwrite): " + ", ".join(existing))
    out_dir.mkdir(parents=True, exist_ok=True)

    jk = os.path.join(args.results_dir, "jk")
    rt, rl = loo_and_total(load_accumulator(
        os.path.join(jk, f"acc_single_random{args.single_random_tag}_BOSS_North_South.npz")))

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE, constrained_layout=True)
    panels = [axes[0, 0], axes[0, 1], axes[1, 0]]
    rows = []
    for ax, (name, bin_text) in zip(panels, BINS):
        lo, hi = (float(v) for v in name[1:].split("_"))
        rand = pooled_random_map(load_random_subbins(str(wide_dir), args.random_label, "BOSS",
                                                     ["North", "South"], lo, hi))
        for i, (label, tag) in enumerate(SMOOTHINGS):
            gal_path = os.path.join(jk, f"acc_single_galaxy_scw{tag}_BOSS_North_South.npz")
            pair_path = wide_dir / "jk" / f"acc_pairs_galaxy_{name}_rpar5_scw{tag}_BOSS_North_South.npz"
            if not (os.path.exists(gal_path) and pair_path.exists()):
                logger.warning("%s at %s: accumulators missing, skipped", bin_text, label)
                continue
            gal = load_accumulator(gal_path)
            gt, gl = loo_and_total(gal)
            m = measure(wide_dir, (gt - rt, gl - rl, gal["grid_size"]), rand, name, tag)
            rows.append({"smoothing": label, "bin": name, "bridge_center_hmpc": m["center"],
                         "filament_bridge_excess": m["bridge"], "filament_bridge_err": m["bridge_err"],
                         "snr": m["bridge"] / m["bridge_err"], "n_pairs": m["n_pairs"]})
            logger.info("%s, %s: %.2f +/- %.2f e-4 (%.2f sigma)", bin_text, label, m["bridge"] * 1e4,
                        m["bridge_err"] * 1e4, m["bridge"] / m["bridge_err"])
            # 4 arcmin at the bottom and 8 arcmin on top, so the noisiest curve does not hide the others;
            # the legend keeps the 8, 6, 4 order
            color, rank = PALETTE[i], len(SMOOTHINGS) - i
            ax.plot(m["axis"], m["prof"] * 1e4, lw=1.6, color=color, label=LEGEND[label], zorder=3 + rank)
            ax.fill_between(m["axis"], (m["prof"] - m["prof_err"]) * 1e4, (m["prof"] + m["prof_err"]) * 1e4,
                            color=color, alpha=0.13, lw=0, zorder=1 + 0.1 * rank)
        center = BRIDGE_CENTERS[name]
        ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6)
        for xx in (-center / 2, center / 2):
            ax.axvline(xx, color=INK_MUTED, lw=0.8, ls="--", alpha=0.6)
        ax.axvspan(-BRIDGE_HALF_X_FRAC * center, BRIDGE_HALF_X_FRAC * center, color=INK_MUTED,
                   alpha=0.2, lw=0, zorder=0)
        ax.set_xlim(-2.5 * center, 2.5 * center)
        ax.set_title(rf"Filament profile, $r_\perp$ = {bin_text} $h^{{-1}}$ Mpc", loc="left", color=INK, pad=5)
        ax.set_xlabel(r"$X$  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$B(X)$  [$10^{-4}$]")
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if ax is panels[0]:
            ax.legend(loc="upper right")

    out = pd.DataFrame(rows)
    out.to_csv(csv_path, index=False)

    ax = axes[1, 1]
    width = 0.8 / len(SMOOTHINGS)
    x = np.arange(len(BINS))
    for i, (label, _) in enumerate(SMOOTHINGS):
        sub = out[out.smoothing == label].set_index("bin")
        vals = [sub.loc[name, "snr"] if name in sub.index else np.nan for name, _ in BINS]
        ax.bar(x + (i - (len(SMOOTHINGS) - 1) / 2) * width, vals, width * 0.9, color=PALETTE[i],
               label=LEGEND[label], zorder=3)
    ax.axhline(0, color=INK_MUTED, lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([t for _, t in BINS])
    ax.set_xlabel(r"$r_\perp$ bin  [$h^{-1}$ Mpc]")
    ax.set_ylabel("filament S/N  ($\\sigma$)")
    ax.set_title("Bridge-excess significance by smoothing", loc="left", color=INK, pad=5)
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s and %s", fig_path, csv_path)


if __name__ == "__main__":
    main()
