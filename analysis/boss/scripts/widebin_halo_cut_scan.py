#!/usr/bin/env python
"""Halo-cut scan: does the bridge excess survive when the window moves away from the galaxies?

The bridge excess is the mean of B(X) over |X| <= 0.35 R.  In a wide bin the
closest pairs put their galaxies at the edge of that window, so a small
mismatch between the pair stack and the control at the haloes leaks into the
"bridge".  This shrinks the window one grid step at a time, from the default
down to the pair centre alone, and recomputes the statistic with its full
jackknife error each time.  A filament keeps its amplitude as the window
narrows; halo leakage fades.

Done for both binnings through one code path, with the mock scored in the same
windows:

    paper   r_perp 4-6 (r_par<=5), 9-11 and 18-22 (r_par<=10); nominal control
    wide    r_perp 4-6, 6-15, 15-25 and the 6-10 / 10-15 split (r_par<=5);
            1/Sigma_crit^2 weights; separation-averaged control

No stacking: it reads the jackknife accumulators and mock stacks the other
wide-bin scripts read.  The first row of each bin is the default window and
must reproduce the bridge excess in filament_jackknife_*.csv.

Writes widebin_rpar5/validation/halo_cut/halo_cut_scan.csv.  It supersedes
bridge_window_scan.csv in the same folder, which was made by hand for the wide
bins only; rows the two share agree exactly.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python \\
        analysis/boss/scripts/widebin_halo_cut_scan.py
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import numpy as np
import pandas as pd

from catalog import setup_logging
from combine_filament_jackknife import (BRIDGE_CENTERS, grouping_filament, load_grouping_pairs,
                                        load_random_subbins, loo_and_total, pooled_random_map)
from combine_jackknife import load_accumulator
from geometry import BRIDGE_HALF_X_FRAC, band_profile, separation_averaged_template, separation_mix
from jackknife import jackknife_error
from plot_observed_figures import sim_diff
from plot_separation_summary import build_terms
from weighting_sensitivity import axis_for_map, load_map

logger = logging.getLogger(__name__)

PAPER_BINS = [("5", 4.0, 6.0), ("10", 9.0, 11.0), ("20", 18.0, 22.0)]
WIDE_BINS = ["g4_6", "g6_15", "g15_25", "g6_10", "g10_15"]


def windows(center):
    """Default half-width, then every whole grid step below it down to the pair centre."""
    default = BRIDGE_HALF_X_FRAC * center
    return [default] + [float(w) for w in range(int(np.floor(default + 1e-9)), -1, -1)
                        if w < default - 1e-9]


def scan(version, bin_text, lo, center, mean_rperp, n_pairs, axis, prof, loo_prof, mock_axis, mock_prof):
    rows = []
    for i, w in enumerate(windows(center)):
        m = np.abs(axis) <= w + 1e-9
        value = float(prof[m].mean())
        err = float(jackknife_error(loo_prof[:, m].mean(axis=1)[:, None])[1][0])
        mock = float(mock_prof[np.abs(mock_axis) <= w + 1e-9].mean())
        rows.append({"version": version, "bin": bin_text, "bridge_center": center, "n_pairs": n_pairs,
                     "window_half": round(w, 3), "default_window": i == 0, "n_x_points": int(m.sum()),
                     "gap_to_nearest_galaxy": round(lo / 2 - w, 3),
                     "gap_to_mean_galaxy": round(mean_rperp / 2 - w, 3),
                     "bridge": value, "bridge_err": err, "bridge_snr": value / err,
                     "mock_bridge": mock, "obs_minus_mock_sigma": (value - mock) / err})
    return rows


def paper_bin(results_dir, key, lo, hi, single_tag, single_random_tag):
    t = build_terms(results_dir, "BOSS", ["North", "South"], key, single_tag, single_random_tag)
    axis = t["axis"]
    loo = np.asarray([band_profile(m.reshape(t["filament"].shape), axis) for m in t["loo"]["filament"]])
    sim_res, sim_axis = sim_diff(key)
    # Narrow symmetric bins: the mean separation is taken as the bin centre.
    return scan("paper", f"{lo:g}-{hi:g}", lo, t["sep"], t["sep"], t["n_pairs"], axis,
                band_profile(t["filament"], axis), loo, sim_axis, band_profile(sim_res, sim_axis))


def wide_bin(wide_dir, sim_dir, sim_single, singles, name, random_label):
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

    base = sim_dir / "hod_pairs" / "widebin_rpar5" / f"stack_{name}_rpar5_avgctrl"
    blocks = np.load(f"{base}_blocks.npz")
    sim_pair = load_map(f"{base}_map.csv")
    sim_axis = axis_for_map(sim_pair)
    gx, gy = np.meshgrid(sim_axis, sim_axis)
    seps, wts, _ = separation_mix(blocks["fine_w"], blocks["fine_wr"])
    sim_res = sim_pair - separation_averaged_template(sim_single, gx, gy, seps, wts)
    mean_rperp = float(pair["fine_wr"].sum() / pair["fine_w"].sum())
    return scan("wide", f"{lo:g}-{hi:g}", lo, center, mean_rperp, int(pair["n_pairs"].sum()), axis,
                band_profile(res, axis), loo, sim_axis, band_profile(sim_res, sim_axis))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--wide-dir", default="analysis/boss/results/widebin_rpar5")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv")
    p.add_argument("--random-label", default="random_w3_25_rpar5_scw_frac100")
    p.add_argument("--single-tag", default="_scw")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    setup_logging()

    wide_dir, sim_dir = Path(args.wide_dir), Path(args.sim_dir)
    out = wide_dir / "validation" / "halo_cut" / "halo_cut_scan.csv"
    if out.exists() and not args.overwrite:
        raise FileExistsError(f"{out} exists (use --overwrite)")
    out.parent.mkdir(parents=True, exist_ok=True)

    jk = os.path.join(args.results_dir, "jk")
    gal = load_accumulator(os.path.join(jk, f"acc_single_galaxy{args.single_tag}_BOSS_North_South.npz"))
    rnd = load_accumulator(os.path.join(jk, f"acc_single_random{args.single_random_tag}_BOSS_North_South.npz"))
    (gt, gl), (rt, rl) = loo_and_total(gal), loo_and_total(rnd)
    singles = (gt - rt, gl - rl, gal["grid_size"])
    sim_single = load_map(sim_dir / args.sim_single)

    rows = []
    for key, lo, hi in PAPER_BINS:
        rows += paper_bin(args.results_dir, key, lo, hi, args.single_tag, args.single_random_tag)
    for name in WIDE_BINS:
        rows += wide_bin(wide_dir, sim_dir, sim_single, singles, name, args.random_label)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)

    logger.info("\nBridge excess (kappa x 1e4) as the window shrinks; gap = nearest galaxy minus window edge")
    for (version, bin_text), sub in df.groupby(["version", "bin"], sort=False):
        logger.info("%s %s (%d pairs)", version, bin_text, sub.n_pairs.iloc[0])
        for r in sub.itertuples():
            logger.info("   |X| <= %-5g gap %5.2f   BOSS %6.2f +/- %5.2f (%+.2f sigma)   mock %5.2f",
                        r.window_half, r.gap_to_nearest_galaxy, r.bridge * 1e4, r.bridge_err * 1e4,
                        r.bridge_snr, r.mock_bridge * 1e4)
    logger.info("Saved %s", out)


if __name__ == "__main__":
    main()
