#!/usr/bin/env python
"""Diagnose the ~20 h^-1 Mpc dip: is it physical, or an artifact of averaging N and S?

The old pipeline formed the combined single-galaxy map as an unweighted mean of
the North and South maps (``plot_results.average_region_maps``), which gives the
South -- 27% of the galaxies -- the same weight as the North.  This script
re-derives the North-only, South-only, unweighted-average, and joint profiles
from the *same* per-region accumulators, so the only thing that varies between
them is the combination rule.

Because the NGC and SGC footprints are disjoint on the sky, every jackknife
cell belongs to exactly one hemisphere.  A hemisphere's profile and its
jackknife error therefore follow from the joint accumulators by restricting to
that hemisphere's regions -- no re-stacking, and the errors stay internally
consistent with the joint measurement.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/region_split_check.py \\
        --dataset BOSS --tag _scw --random-tag _scw_frac10
"""

from __future__ import annotations

import argparse
import logging
import os

import numpy as np
import pandas as pd

from catalog import load_catalog_lightweight, resolve_catalog_path, setup_logging
from combine_jackknife import (
    grid_radius_hmpc,
    leave_one_out_maps,
    load_accumulator,
    radial_profiles,
    total_map,
)
from jackknife import JackknifeRegions, jackknife_error

logger = logging.getLogger(__name__)


def hemisphere_regions(
    regions: JackknifeRegions, dataset: str, data_dir: str, survey_regions: list[str],
    fraction: float, seed: int,
) -> dict[str, np.ndarray]:
    """Map each survey region name to the jackknife region indices it occupies."""
    z_min, z_max = (0.4, 0.7) if dataset == "BOSS" else (0.0, 10000.0)
    occupied = {}
    for name in survey_regions:
        path = resolve_catalog_path(data_dir, dataset, name, "random")
        data = load_catalog_lightweight(
            path, columns=("RA", "DEC", "Z"), fraction=fraction,
            z_min=z_min, z_max=z_max, seed=seed,
        )
        labels = regions.assign(data["RA"], data["DEC"])
        occupied[name] = np.unique(labels)
        logger.info("%s occupies %d jackknife regions", name, len(occupied[name]))

    names = list(survey_regions)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            shared = np.intersect1d(occupied[a], occupied[b])
            if shared.size:
                raise ValueError(
                    f"{shared.size} jackknife region(s) straddle {a} and {b}: "
                    f"{shared.tolist()}. A hemisphere split is not clean at this "
                    "nside; lower --jk-nside or split by survey region instead."
                )
    return occupied


def subset_profile(
    gal_acc: dict, rand_acc: dict, region_idx: np.ndarray, edges: np.ndarray,
    grid_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Corrected profile and jackknife error using only ``region_idx`` regions."""
    gal_wk, gal_w = gal_acc["sum_wk"][region_idx], gal_acc["sum_w"][region_idx]
    rand_wk, rand_w = rand_acc["sum_wk"][region_idx], rand_acc["sum_w"][region_idx]

    corrected = total_map(gal_wk, gal_w) - total_map(rand_wk, rand_w)
    loo = leave_one_out_maps(gal_wk, gal_w) - leave_one_out_maps(rand_wk, rand_w)

    profiles, _ = radial_profiles(
        np.vstack([corrected[None, :], loo]), edges, grid_size
    )
    _, sigma = jackknife_error(profiles[1:])
    return profiles[0], sigma


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare N-only, S-only, unweighted-average and joint profiles."
    )
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument("--tag", default="_scw")
    parser.add_argument("--random-tag", default=None)
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--jk-fraction", type=float, default=0.10)
    parser.add_argument("--jk-seed", type=int, default=20260727)
    parser.add_argument("--r-max", type=float, default=50.0)
    parser.add_argument("--bin-width", type=float, default=1.0)
    args = parser.parse_args(argv)
    args.regions = [r.strip() for r in args.regions.split(",") if r.strip()]
    if args.random_tag is None:
        args.random_tag = args.tag
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()

    region_label = "_".join(args.regions)
    jk_dir = os.path.join(args.results_dir, "jk")
    gal = load_accumulator(
        os.path.join(jk_dir, f"acc_single_galaxy{args.tag}_{args.dataset}_{region_label}.npz")
    )
    rand = load_accumulator(
        os.path.join(jk_dir, f"acc_single_random{args.random_tag}_{args.dataset}_{region_label}.npz")
    )
    if gal["jk_digest"] != rand["jk_digest"]:
        raise ValueError("Galaxy and random accumulators use different tessellations.")

    regions = JackknifeRegions.load(
        os.path.join(jk_dir, f"regions_{args.dataset}_{region_label}_nside{gal['jk_nside']}.npz")
    )
    if regions.digest != gal["jk_digest"]:
        raise ValueError(
            f"Cached regions digest {regions.digest} does not match the accumulators "
            f"({gal['jk_digest']}); the tessellation was rebuilt after stacking."
        )

    grid_size = gal["grid_size"]
    edges = np.arange(0.0, args.r_max + args.bin_width, args.bin_width)
    occupied = hemisphere_regions(
        regions, args.dataset, args.data_dir, args.regions, args.jk_fraction, args.jk_seed
    )

    columns = {}
    for name, idx in occupied.items():
        p, s = subset_profile(gal, rand, idx, edges, grid_size)
        columns[name] = p
        columns[f"{name}_err"] = s

    all_idx = np.arange(gal["sum_wk"].shape[0])
    joint, joint_err = subset_profile(gal, rand, all_idx, edges, grid_size)
    columns["joint"] = joint
    columns["joint_err"] = joint_err

    # The rule the old pipeline used: mean of the per-region maps, weight 1/2 each.
    columns["unweighted_avg"] = np.mean([columns[n] for n in occupied], axis=0)

    n_gal = {n: int(gal["n_objects"][idx].sum()) for n, idx in occupied.items()}
    total_gal = sum(n_gal.values())
    logger.info("Galaxies per hemisphere: %s (South is %.1f%% of the sample)",
                n_gal, 100.0 * min(n_gal.values()) / total_gal)
    columns["count_weighted_avg"] = np.sum(
        [columns[n] * (n_gal[n] / total_gal) for n in occupied], axis=0
    )

    out = pd.DataFrame({"r_center_hmpc": 0.5 * (edges[:-1] + edges[1:]), **columns})
    path = os.path.join(
        args.results_dir,
        f"profile_region_split{args.tag}_{args.dataset}_{region_label}.csv",
    )
    out.to_csv(path, index=False)
    logger.info("Saved region-split comparison -> %s", path)

    window = (out.r_center_hmpc >= 15) & (out.r_center_hmpc <= 30)
    view = out.loc[window, ["r_center_hmpc", *occupied, "unweighted_avg",
                            "count_weighted_avg", "joint"]]
    logger.info("Corrected kappa x 1e4, r = 15-30 h^-1 Mpc:\n%s",
                (view.set_index("r_center_hmpc") * 1e4).round(2).to_string())

    for name in occupied:
        p, s = out[name].to_numpy(), out[f"{name}_err"].to_numpy()
        seg = window.to_numpy()
        i = int(np.argmin(np.where(seg, p, np.inf)))
        logger.info(
            "%s: minimum in 15-30 is %.2fe-4 at r=%.1f, which is %.1f sigma from zero "
            "(sigma = %.2fe-4)",
            name, p[i] * 1e4, out.r_center_hmpc[i], p[i] / s[i], s[i] * 1e4,
        )
    logger.info(
        "Joint at that radius: %.2fe-4 +/- %.2fe-4",
        out.joint[i] * 1e4, out.joint_err[i] * 1e4,
    )


if __name__ == "__main__":
    main()
