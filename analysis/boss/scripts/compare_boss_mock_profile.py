#!/usr/bin/env python
"""BOSS-vs-mock deficit significance, using the joint jackknife covariance.

Replaces the retracted 3.8/5.6/2.9 sigma deficit numbers.  Those combined an
over-weighted South with an error built as sqrt(sigma_N^2 + sigma_S^2) -- the
error on the sum rather than the average -- and a pixel-to-profile conversion
that assumed radial bins decorrelate.  This script uses the joint North+South
stack and the measured bin-to-bin covariance instead.

The mock is treated as noiseless.  BigMDPL is a 2500 h^-1 Mpc box, so its
single-halo profile is determined far better than the BOSS measurement; the
deficit error is therefore taken to be the BOSS error alone.  That is an
assumption, not a measurement -- it makes the quoted significances mild upper
bounds on what a full error budget would give.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/compare_boss_mock_profile.py \\
        --mock analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin.csv
"""

from __future__ import annotations

import argparse
import logging
import os

import numpy as np
import pandas as pd

from catalog import setup_logging
from combine_jackknife import DEFAULT_BANDS, band_average, grid_radius_hmpc, radial_profiles

logger = logging.getLogger(__name__)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Deficit of BOSS corrected kappa relative to a density-matched mock.")
    parser.add_argument("--mock", required=True, help="Simulation single-stack CSV.")
    parser.add_argument("--boss-profile", default=None)
    parser.add_argument("--boss-cov", default=None)
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--stem",
                        default="_scw_rand_scw_frac100_BOSS_North_South_joint_linear")
    parser.add_argument("--mock-label", default="hodnmatch_8arcmin")
    parser.add_argument("--validate-against", default=None,
                        help="Archived sim radial-profile CSV; checks that this "
                             "script's binning reproduces it.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    setup_logging()
    R = args.results_dir

    prof_path = args.boss_profile or os.path.join(R, f"profile_corrected{args.stem}.csv")
    cov_path = args.boss_cov or os.path.join(R, f"profile_cov_corrected{args.stem}.csv")
    boss = pd.read_csv(prof_path)
    cov = pd.read_csv(cov_path, index_col=0).to_numpy(dtype=float)
    edges = np.append(boss.r_lo_hmpc.to_numpy(), boss.r_hi_hmpc.to_numpy()[-1])
    counts = boss.n_pixels.to_numpy()
    logger.info("BOSS profile: %s (%d bins)", prof_path, len(boss))

    mock_map = pd.read_csv(args.mock, index_col=0).to_numpy(dtype=float)
    if mock_map.shape[0] != mock_map.shape[1]:
        raise ValueError(f"Mock map is not square: {mock_map.shape}")
    logger.info("Mock: %s %s", args.mock, mock_map.shape)

    mock_prof, mock_counts = radial_profiles(
        mock_map.reshape(1, -1), edges, mock_map.shape[0])
    mock_prof = mock_prof[0]

    if not np.array_equal(mock_counts, counts):
        raise ValueError(
            "Mock and BOSS maps give different pixel counts per radial bin; the "
            "two stacks are not on the same grid and cannot be differenced bin-wise.")

    if args.validate_against:
        ref = pd.read_csv(args.validate_against)
        merged = pd.DataFrame({"r": boss.r_center_hmpc, "mine": mock_prof}).merge(
            ref.rename(columns={"radius_hmpc": "r", "kappa": "archived"}), on="r", how="inner")
        if merged.empty:
            logger.warning("No overlapping radii with %s; skipping validation.",
                           args.validate_against)
        else:
            dev = np.abs(merged.mine - merged.archived).max()
            logger.info(
                "Binning validation vs %s: %d shared radii, max |diff| = %.3e "
                "(%.2f%% of peak)", args.validate_against, len(merged), dev,
                100 * dev / np.abs(merged.archived).max())

    deficit = boss.kappa.to_numpy() - mock_prof

    rows = []
    for lo, hi in DEFAULT_BANDS:
        boss_val, boss_err, n_bins = band_average(
            boss.kappa.to_numpy(), cov, edges, counts, lo, hi)
        mock_val, _, _ = band_average(mock_prof, cov, edges, counts, lo, hi)
        diff = boss_val - mock_val
        rows.append({
            "band_lo_hmpc": lo, "band_hi_hmpc": hi, "n_bins": n_bins,
            "boss_kappa": boss_val, "mock_kappa": mock_val,
            "deficit": diff, "deficit_err": boss_err,
            "deficit_significance": diff / boss_err if boss_err > 0 else np.nan,
            "boss_over_mock": boss_val / mock_val if mock_val != 0 else np.nan,
        })

    out = pd.DataFrame(rows)
    path = os.path.join(R, f"boss_vs_mock_deficit_{args.mock_label}.csv")
    out.to_csv(path, index=False)

    profile_out = pd.DataFrame({
        "r_center_hmpc": boss.r_center_hmpc, "boss_kappa": boss.kappa,
        "boss_err": boss.kappa_err, "mock_kappa": mock_prof, "deficit": deficit,
    })
    profile_out.to_csv(
        os.path.join(R, f"boss_vs_mock_profile_{args.mock_label}.csv"), index=False)

    logger.info("Saved -> %s", path)
    logger.info("\nBOSS vs %s (kappa x 1e4; negative deficit = BOSS below mock):",
                args.mock_label)
    logger.info("  band        BOSS    mock   deficit    err   sigma   BOSS/mock")
    for _, r in out.iterrows():
        logger.info("  %2.0f-%2.0f   %7.2f %7.2f %9.2f %6.2f %7.1f %11.2f",
                    r.band_lo_hmpc, r.band_hi_hmpc, r.boss_kappa * 1e4,
                    r.mock_kappa * 1e4, r.deficit * 1e4, r.deficit_err * 1e4,
                    r.deficit_significance, r.boss_over_mock)


if __name__ == "__main__":
    main()
