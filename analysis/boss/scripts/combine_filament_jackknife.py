#!/usr/bin/env python
"""Filament map and bridge excess with jackknife errors.

Chains the joint accumulators into the derived-map sequence that
``plot_results.py`` builds, but does it once per leave-one-out region so the
filament signal gets an uncertainty:

    corrected pairs  = galaxy pairs - random pairs
    corrected single = single galaxies - single randoms
    control          = corrected single superposed at +/- sep/2
    filament         = corrected pairs - control

Every term that carries galaxy shot noise -- the galaxy pairs and both single
stacks -- is deleted region by region *together*, using the shared
tessellation.  That coherence is the point: the control is built from the
single stack, so a jackknife that dropped a region from the pair term but not
from the control would mis-estimate the very cancellation the filament
measurement depends on.

The random-pair map is held fixed rather than jackknifed, because per-region
random-pair accumulators would require re-running the 23-36 hour full-random
pair jobs.  That approximation is quantified rather than assumed: the measured
bridge excess of the frac100 random-pair maps is reported, and since randoms
carry no filament signal its departure from zero is a direct estimate of the
noise this term contributes.  Empirically it is ~0.5e-4, well below the
galaxy-side error.

The frac10 random-pair maps are *not* usable here.  Subsampling the random
catalog to 10% leaves 1% of the pairs -- pair counts go as the square of the
fraction -- so their noise is ~10x the frac100 maps.  Re-measured under the
fixed Y bands (2026-08-06), swapping frac10 for frac100 shifts the per-region
random bridge excess by -8.8 and +24.2e-4 at sep 5 and by -4.8 and -6.3e-4 at
sep 10.  Jointly that is a +5.2e-4 shift in the sep-10 filament against a 5.1e-4
error: a full 1 sigma of pure bookkeeping.  Sep 20 survives the swap (0.2e-4)
only because it has 16x more random pairs than sep 5, so do not generalize from
it -- any new separation needs its own frac100 run.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/combine_filament_jackknife.py \\
        --separations 5,10,20
"""

from __future__ import annotations

import argparse
import json
import logging
import os

import numpy as np
import pandas as pd

from catalog import setup_logging
from combine_jackknife import leave_one_out_maps, load_accumulator, total_map
from geometry import reflect_symmetrize_map, symmetrize_map, two_halo_template
from jackknife import jackknife_error
from weighting_sensitivity import (
    SEPARATIONS,
    axis_for_map,
    bridge_excess,
    build_control_pair_map,
    combine,
    load_map,
    reconcile_shapes,
)

logger = logging.getLogger(__name__)

GALAXY_COUNTS = {"North": 579089.0, "South": 213205.0}


def loo_and_total(acc: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return (total_map, leave_one_out_maps) as flattened arrays."""
    return (
        total_map(acc["sum_wk"], acc["sum_w"]),
        leave_one_out_maps(acc["sum_wk"], acc["sum_w"]),
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Filament bridge excess with jackknife errors.")
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument("--separations", default="5,10,20")
    parser.add_argument("--single-tag", default="_scw")
    parser.add_argument("--single-random-tag", default=None,
                        help="Defaults to --single-tag; use e.g. _scw_frac100.")
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--no-symmetrize-single", action="store_true",
                        help="Build the control from the unsymmetrized single map. "
                             "Off by default: the archived pipeline symmetrizes, and "
                             "not doing so leaks azimuthal noise into the filament.")
    parser.add_argument("--symmetrize", action="store_true",
                        help="Reflection-symmetrize each pair realization before "
                             "deriving statistics (matches the archived maps).")
    args = parser.parse_args(argv)
    args.regions = [r.strip() for r in args.regions.split(",")]
    args.separations = [s.strip() for s in args.separations.split(",")]
    if args.single_random_tag is None:
        args.single_random_tag = args.single_tag
    return args


def main(argv=None):
    args = parse_args(argv)
    setup_logging()
    R = args.results_dir
    jk_dir = os.path.join(R, "jk")
    region_label = "_".join(args.regions)

    # ---- single-galaxy side: corrected single, per leave-one-out region ----
    gal = load_accumulator(os.path.join(
        jk_dir, f"acc_single_galaxy{args.single_tag}_{args.dataset}_{region_label}.npz"))
    rnd = load_accumulator(os.path.join(
        jk_dir, f"acc_single_random{args.single_random_tag}_{args.dataset}_{region_label}.npz"))
    if gal["jk_digest"] != rnd["jk_digest"]:
        raise ValueError("Single galaxy/random accumulators use different tessellations.")

    gal_total, gal_loo = loo_and_total(gal)
    rnd_total, rnd_loo = loo_and_total(rnd)
    single_total = gal_total - rnd_total
    single_loo = gal_loo - rnd_loo
    n_regions = single_loo.shape[0]
    g = gal["grid_size"]

    rows = []
    for sep in args.separations:
        cfg = SEPARATIONS[sep]
        center = cfg["center"]

        pair_acc_path = os.path.join(
            jk_dir, f"acc_pairs_galaxy_{cfg['galaxy_label']}_{args.dataset}_{region_label}.npz")
        if not os.path.exists(pair_acc_path):
            logger.warning("Missing %s; skipping sep=%s", pair_acc_path, sep)
            continue
        pair = load_accumulator(pair_acc_path)
        if pair["jk_digest"] != gal["jk_digest"]:
            raise ValueError(
                f"Pair accumulator for sep={sep} uses tessellation {pair['jk_digest']} "
                f"but the single stacks use {gal['jk_digest']}; the filament jackknife "
                "would delete different sky from the two terms.")
        gp = pair["grid_size"]
        pair_total, pair_loo = loo_and_total(pair)

        # Random pairs: fixed, count-weighted across regions (no accumulators).
        rand_pair_maps, rand_pair_counts = {}, {}
        for reg in args.regions:
            rand_pair_maps[reg] = load_map(os.path.join(
                R, f"kappa_pairs_random_{cfg['random_label']}_{args.dataset}_{reg}.csv"))
            with open(os.path.join(
                R, f"kappa_pairs_random_{cfg['random_label']}_{args.dataset}_{reg}.meta.json"),
                encoding="utf-8") as handle:
                rand_pair_counts[reg] = float(json.load(handle)["total_pairs"])
        rand_pair = combine(rand_pair_maps, rand_pair_counts)
        rand_pair_floor = bridge_excess(rand_pair, center)

        def filament_from(single_flat: np.ndarray, pair_flat: np.ndarray) -> np.ndarray:
            pair_map = pair_flat.reshape(gp, gp)
            if args.symmetrize:
                pair_map = reflect_symmetrize_map(pair_map)
            cp, rp = reconcile_shapes(pair_map, rand_pair)
            corrected_pairs = cp - rp
            # Radially symmetrize the single map before superposing it: the
            # single-galaxy signal is radially symmetric by construction, so this
            # suppresses azimuthal noise without biasing the control. Skipping it
            # leaks that noise straight into the filament residual.
            single_map = single_flat.reshape(g, g)
            if not args.no_symmetrize_single:
                single_map = symmetrize_map(single_map)
            # Build the control natively on the pair grid via the radial profile,
            # rather than shifting a pixel array and trimming to a common shape.
            # The latter left a half-pixel misalignment (101-pixel pair grid on
            # integers vs 100-pixel single grid on half-integers) that produced
            # false residuals at the halo peaks. See geometry.two_halo_template.
            pair_axis = axis_for_map(corrected_pairs)
            gx, gy = np.meshgrid(pair_axis, pair_axis)
            # The symmetry check exists to catch archived mis-centered stacks.
            # Under --no-symmetrize-single the caller is handing in a raw map on
            # purpose, so the check would reject the very thing being tested.
            control = two_halo_template(single_map, gx, gy, center,
                                        validate=not args.no_symmetrize_single)
            return corrected_pairs - control, corrected_pairs, control

        fil_total, cpairs_total, control_total = filament_from(single_total, pair_total)
        be_total = bridge_excess(fil_total, center)

        be_loo = np.empty(n_regions)
        fil_loo = np.empty((n_regions, fil_total.size))
        for k in range(n_regions):
            fil_k, _, _ = filament_from(single_loo[k], pair_loo[k])
            fil_loo[k] = fil_k.ravel()
            be_loo[k] = bridge_excess(fil_k, center)

        _, be_err = jackknife_error(be_loo)
        _, fil_err = jackknife_error(fil_loo)

        rows.append({
            "separation_hmpc": center,
            "n_pairs": int(pair["n_pairs"].sum()),
            "corrected_pairs_bridge_excess": bridge_excess(cpairs_total, center),
            "control_bridge_excess": bridge_excess(control_total, center),
            "filament_bridge_excess": be_total,
            "filament_bridge_err": be_err,
            "significance": be_total / be_err if be_err > 0 else np.nan,
            "random_pair_noise_floor": rand_pair_floor,
        })

        pd.DataFrame(fil_total).to_csv(os.path.join(
            R, f"kappa_filament_{cfg['galaxy_label']}_{args.dataset}_{region_label}_joint.csv"))
        pd.DataFrame(fil_err.reshape(fil_total.shape)).to_csv(os.path.join(
            R, f"error_filament_{cfg['galaxy_label']}_{args.dataset}_{region_label}_joint.csv"))
        logger.info("sep=%s: filament bridge excess %.3fe-4 +/- %.3fe-4 (%.1f sigma)",
                    sep, be_total * 1e4, be_err * 1e4,
                    be_total / be_err if be_err > 0 else np.nan)

    if not rows:
        logger.warning("No separations produced results.")
        return

    out = pd.DataFrame(rows)
    path = os.path.join(R, f"filament_jackknife_{args.dataset}_{region_label}.csv")
    out.to_csv(path, index=False)
    logger.info("Saved -> %s", path)

    logger.info("\nJoint filament measurement (kappa x 1e4):")
    logger.info("  sep   corrected_pairs   control   filament   error   sigma   rand-floor")
    for _, r in out.iterrows():
        logger.info("  %4.0f  %14.2f  %8.2f  %9.2f  %6.2f  %5.1f  %9.2f",
                    r.separation_hmpc, r.corrected_pairs_bridge_excess * 1e4,
                    r.control_bridge_excess * 1e4, r.filament_bridge_excess * 1e4,
                    r.filament_bridge_err * 1e4, r.significance,
                    r.random_pair_noise_floor * 1e4)


if __name__ == "__main__":
    main()
