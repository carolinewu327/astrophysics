#!/usr/bin/env python
"""How much do the pair/filament products move under correct N/S weighting?

``plot_results.average_region_maps`` combines the North and South maps with an
unweighted mean, giving the South -- 27% of the galaxies -- half the weight.
On the single-galaxy profile that error manufactured an apparent dip at
r ~ 20 h^-1 Mpc out of a -1.4 sigma noise excursion in the South
(``region_split_check.py``).  The pair, control and filament maps are built
through the same averaging call, so they inherit the same defect.

This script rebuilds the derived maps under both combination rules from the
archived per-region CSVs and reports how far the bridge statistics move.  It
answers one question: is rebuilding the pair pipeline with per-region
accumulators necessary, or is the effect negligible?

The count-weighted rule used here is an approximation to the exact joint
estimator.  Exactly, each *pixel* should be weighted by its ``sum(w)`` in each
region; the archived CSVs store only the mean map, so a single scalar weight
per region stands in.  ``sum(w)`` varies smoothly across the map, so the two
agree closely -- but a genuinely joint pair stack is the only way to be exact,
and to get jackknife errors on the filament signal.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/weighting_sensitivity.py
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os

import numpy as np
import pandas as pd

from catalog import setup_logging
from geometry import BRIDGE_HALF_X_FRAC
from geometry import bridge_excess as _bridge_excess

logger = logging.getLogger(__name__)

BOX_SIZE_HMPC = 100.0

# Bridge/sideband geometry now lives in lib/geometry.py so that the BOSS and
# simulation pipelines cannot drift apart; see the band-geometry comment there.
BRIDGE_HALF_X = BRIDGE_HALF_X_FRAC

# Pair catalogs are named {r_par}_{r_perp_min}_{r_perp_max}: the separation is
# set by r_perp, while r_par is the line-of-sight cut. The 20 h^-1 Mpc analysis
# uses the r_par=10 variant, so its catalog is 10.0_18.0_22.0, NOT 20.0_18.0_22.0.
# Verified against the rpar/rperp fields in the random-pair .meta.json sidecars.
SEPARATIONS = {
    "5": {"center": 5.0, "galaxy_label": "5", "random_label": "5_frac100",
          "pair_catalog": "5.0_4.0_6.0"},
    "10": {"center": 10.0, "galaxy_label": "10", "random_label": "10_frac100",
           "pair_catalog": "10.0_9.0_11.0"},
    "20": {"center": 20.0, "galaxy_label": "20_rpar10", "random_label": "20_rpar10_frac100",
           "pair_catalog": "10.0_18.0_22.0"},
    # Zheng asked for one LOS cut across all separations, so the 5 h^-1 Mpc
    # sample moves to r_par <= 10 like the other two.  Added as a *separate* key
    # rather than by editing "5": every published sep-5 number was measured at
    # r_par <= 5, and silently repointing the old key would restate them at a cut
    # they were never computed with.
    #
    # The galaxy accumulator (acc_pairs_galaxy_5_rpar10_*, 344,306 pairs) and the
    # simulation stack (stack_rperp5_rpar10_matched) both exist.  The frac100
    # random-pair maps do NOT yet -- that is the ~15 h job.  This key therefore
    # raises FileNotFoundError until they land, which is the intended behaviour:
    # falling back to the r_par <= 5 randoms would subtract a differently-cut
    # background from the galaxy pairs and quietly bias the result.
    "5_rpar10": {"center": 5.0, "galaxy_label": "5_rpar10",
                 "random_label": "5_rpar10_frac100", "pair_catalog": "10.0_4.0_6.0"},
}


def load_map(path: str) -> np.ndarray:
    return pd.read_csv(path, index_col=0).to_numpy(dtype=float)


def axis_for_map(arr: np.ndarray) -> np.ndarray:
    n = arr.shape[0]
    if n % 2 == 1:
        return np.linspace(-0.5 * BOX_SIZE_HMPC, 0.5 * BOX_SIZE_HMPC, n)
    cell = BOX_SIZE_HMPC / n
    return np.linspace(-0.5 * BOX_SIZE_HMPC + 0.5 * cell,
                       0.5 * BOX_SIZE_HMPC - 0.5 * cell, n)


def reconcile_shapes(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if a.shape == b.shape:
        return a, b
    rows = min(a.shape[0], b.shape[0])
    cols = min(a.shape[1], b.shape[1])

    def trim(arr):
        r0 = (arr.shape[0] - rows) // 2
        c0 = (arr.shape[1] - cols) // 2
        return arr[r0:r0 + rows, c0:c0 + cols]

    return trim(a), trim(b)


def build_control_pair_map(single_map: np.ndarray, separation_hmpc: float,
                           target_axis: np.ndarray | None = None) -> np.ndarray:
    """Superposed-singles control, built on *target_axis* via the radial profile.

    Rolling/shifting the single map's pixel grid and trimming to a common shape
    carried a half-pixel error (the pair grid is 101 px on integer coordinates,
    the single grid 100 px on half-integers). Building the control on the target
    grid removes it -- see geometry.two_halo_template.
    """
    from geometry import two_halo_template
    if target_axis is None:
        target_axis = axis_for_map(single_map)
    gx, gy = np.meshgrid(target_axis, target_axis)
    return two_halo_template(single_map, gx, gy, separation_hmpc)


def bridge_excess(arr: np.ndarray, rperp_center: float) -> float:
    """Bridge statistic on a map whose axis follows from its own grid size."""
    return _bridge_excess(arr, axis_for_map(arr), rperp_center)


def count_gz_rows(path: str) -> int:
    """Rows in a gzipped CSV, excluding the header."""
    total = 0
    with gzip.open(path, "rb") as handle:
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            total += block.count(b"\n")
    return max(total - 1, 0)


def combine(maps: dict[str, np.ndarray], weights: dict[str, float] | None) -> np.ndarray:
    names = list(maps)
    if weights is None:
        return np.mean([maps[n] for n in names], axis=0)
    total = sum(weights[n] for n in names)
    return np.sum([maps[n] * (weights[n] / total) for n in names], axis=0)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Measure how far pair/filament products move under correct "
                    "North/South weighting.")
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--paircatalog-dir", default="data/paircatalogs/BOSS")
    parser.add_argument("--separations", default="5,10,20")
    args = parser.parse_args(argv)
    args.regions = [r.strip() for r in args.regions.split(",")]
    args.separations = [s.strip() for s in args.separations.split(",")]
    return args


def main(argv=None):
    args = parse_args(argv)
    setup_logging()
    R = args.results_dir

    # Galaxy counts after the CMASS redshift cut, as stacked.
    gal_counts = {"North": 579089.0, "South": 213205.0}
    logger.info("Single-galaxy weights: unweighted 0.500/0.500 vs count-weighted "
                "%.3f/%.3f", gal_counts["North"] / sum(gal_counts.values()),
                gal_counts["South"] / sum(gal_counts.values()))

    single_gal = {r: load_map(os.path.join(
        R, f"kappa_single_galaxy_{args.dataset}_{r}_nojk.csv")) for r in args.regions}
    single_rand = {r: load_map(os.path.join(
        R, f"kappa_single_random_frac100_{args.dataset}_{r}.csv")) for r in args.regions}

    rows = []
    for sep in args.separations:
        cfg = SEPARATIONS[sep]
        center = cfg["center"]

        pair_gal = {r: load_map(os.path.join(
            R, f"kappa_pairs_galaxy_{cfg['galaxy_label']}_{args.dataset}_{r}.csv"))
            for r in args.regions}
        pair_rand = {r: load_map(os.path.join(
            R, f"kappa_pairs_random_{cfg['random_label']}_{args.dataset}_{r}.csv"))
            for r in args.regions}

        gal_pair_counts = {}
        for r in args.regions:
            path = os.path.join(
                args.paircatalog_dir,
                f"galaxy_pairs_{args.dataset}_{r}_{cfg['pair_catalog']}hmpc.csv.gz")
            gal_pair_counts[r] = float(count_gz_rows(path))
        rand_pair_counts = {}
        for r in args.regions:
            meta = os.path.join(
                R, f"kappa_pairs_random_{cfg['random_label']}_{args.dataset}_{r}.meta.json")
            with open(meta, encoding="utf-8") as handle:
                rand_pair_counts[r] = float(json.load(handle)["total_pairs"])

        logger.info(
            "sep=%s: galaxy pairs N/S = %d/%d (South share %.1f%%), "
            "random pairs South share %.1f%%",
            sep, gal_pair_counts["North"], gal_pair_counts["South"],
            100 * gal_pair_counts["South"] / sum(gal_pair_counts.values()),
            100 * rand_pair_counts["South"] / sum(rand_pair_counts.values()),
        )

        for rule, w_single, w_gpair, w_rpair in (
            ("unweighted (current)", None, None, None),
            ("count-weighted", gal_counts, gal_pair_counts, rand_pair_counts),
        ):
            m1 = combine(single_gal, w_single)
            m3 = combine(single_rand, w_single)
            m2 = combine(pair_gal, w_gpair)
            m4 = combine(pair_rand, w_rpair)

            m1t, m3t = reconcile_shapes(m1, m3)
            corrected_single = m1t - m3t
            m2t, m4t = reconcile_shapes(m2, m4)
            corrected_pairs = m2t - m4t

            control = build_control_pair_map(
                corrected_single, center, target_axis=axis_for_map(corrected_pairs))
            cp, ct = reconcile_shapes(corrected_pairs, control)
            filament = cp - ct

            rows.append({
                "separation_hmpc": center,
                "rule": rule,
                "corrected_pairs_bridge_excess": bridge_excess(corrected_pairs, center),
                "control_bridge_excess": bridge_excess(control, center),
                "filament_bridge_excess": bridge_excess(filament, center),
            })

    out = pd.DataFrame(rows)
    pivot = out.pivot(index="separation_hmpc", columns="rule")
    path = os.path.join(R, f"weighting_sensitivity_{args.dataset}.csv")
    out.to_csv(path, index=False)
    logger.info("Saved -> %s", path)

    logger.info("\nBridge excess in kappa x 1e4, by combination rule:")
    for stat in ("corrected_pairs", "control", "filament"):
        col = f"{stat}_bridge_excess"
        a = pivot[col]["unweighted (current)"] * 1e4
        b = pivot[col]["count-weighted"] * 1e4
        logger.info("  %s:", stat)
        for sep in a.index:
            shift = (b[sep] - a[sep]) / abs(a[sep]) * 100 if a[sep] != 0 else np.nan
            logger.info("    sep=%4.0f  current %8.3f -> weighted %8.3f  (%+.1f%%)",
                        sep, a[sep], b[sep], shift)


if __name__ == "__main__":
    main()
