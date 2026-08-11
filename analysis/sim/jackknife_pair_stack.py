#!/usr/bin/env python
"""Pair stack with spatial-jackknife errors on the residual bridge statistics.

Splits the pairs into K = blocks_per_side^2 spatial blocks by pair-center
position, stacks each block once (total cost = one full stack), and combines
leave-one-out means to estimate jackknife uncertainties of the bridge
statistics of (pair stack - matched two-halo template).

The template is built from the supplied single stack (same map/smoothing as
the pair stack: matched-template rule) and held fixed across jackknife
regions -- its own sampling noise is negligible when the single stack uses
the full galaxy sample.

Supports nested LOS cuts via --rpar-max (filters r_parallel_rsd), so one
deep pair catalog serves every candidate cut. Maps are expected pre-smoothed;
stacking is always fixed-separation (normalized stacking is a Phase 4
sensitivity axis).
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from sim_utils import (
    make_two_halo_template,
    open_kappa_memmap,
    radial_symmetrize_map,
    reflect_symmetrize_map,
    save_map_csv,
    setup_logging,
)
from stack_pairs_sim import stack_pairs
from summarize_sim_sensitivity import load_map, map_stats


logger = logging.getLogger(__name__)


def jackknife_stats(theta: np.ndarray) -> tuple[float, float]:
    """Jackknife mean and error from leave-one-out estimates."""
    k = len(theta)
    mean = float(np.mean(theta))
    var = (k - 1) / k * float(np.sum((theta - mean) ** 2))
    return mean, float(np.sqrt(var))


def run(args: argparse.Namespace) -> dict:
    pairs = pd.read_csv(args.pairs)
    n_total = len(pairs)
    rpar_col = "r_parallel_real" if args.rpar_space == "real" else "r_parallel_rsd"
    if args.rpar_min is not None or args.rpar_max is not None:
        if rpar_col not in pairs.columns:
            raise ValueError(
                f"Pair catalog lacks {rpar_col!r}; it predates the real/rsd split. "
                "Regenerate it with find_pairs_sim.py.")
        values = pairs[rpar_col].to_numpy()
        keep = np.ones(len(pairs), dtype=bool)
        if args.rpar_min is not None:
            keep &= values >= args.rpar_min
        if args.rpar_max is not None:
            keep &= values <= args.rpar_max
        pairs = pairs[keep]
    if args.max_pairs is not None:
        pairs = pairs.iloc[: args.max_pairs]
    logger.info("Pairs: %d of %d after cuts (%s in [%s, %s])", len(pairs), n_total,
                rpar_col, args.rpar_min, args.rpar_max)

    kappa_map, info = open_kappa_memmap(args.kappa_map)
    box = info.box_size_hmpc
    cell = box / args.blocks_per_side
    bx = np.clip((pairs["pair_center_x"].to_numpy() // cell).astype(int), 0, args.blocks_per_side - 1)
    by = np.clip((pairs["pair_center_y"].to_numpy() // cell).astype(int), 0, args.blocks_per_side - 1)
    block = bx * args.blocks_per_side + by

    grid = args.grid_size
    sums = []
    counts = []
    for b in range(args.blocks_per_side ** 2):
        sub = pairs[block == b]
        if len(sub) == 0:
            continue
        mean_map = stack_pairs(
            pairs=sub,
            kappa_map=kappa_map,
            pixel_size_hmpc=info.pixel_size_hmpc,
            map_box_size_hmpc=box,
            stack_box_size_hmpc=args.box_size,
            grid_size=grid,
            normalize_separation=False,
            normalized_half_size=2.5,
        ).astype(np.float64)
        sums.append(mean_map * len(sub))
        counts.append(len(sub))
    sums = np.array(sums)
    counts = np.array(counts, dtype=np.float64)
    k = len(counts)
    logger.info("Stacked %d non-empty blocks (%d pairs)", k, int(counts.sum()))

    single = radial_symmetrize_map(load_map(Path(args.single)))
    axis = np.linspace(-0.5 * args.box_size, 0.5 * args.box_size, grid)
    x_grid, y_grid = np.meshgrid(axis, axis)
    template = make_two_halo_template(single, x_grid, y_grid, args.rperp_center)
    norm_axis = axis / args.rperp_center

    full_map = sums.sum(axis=0) / counts.sum()
    full_stats = {}
    full_stats.update(map_stats(full_map, norm_axis, "raw"))
    full_stats.update(map_stats(full_map - template, norm_axis, "residual"))

    loo_keys = None
    loo_values = []
    total_sum = sums.sum(axis=0)
    total_count = counts.sum()
    for b in range(k):
        loo_map = (total_sum - sums[b]) / (total_count - counts[b])
        stats = {}
        stats.update(map_stats(loo_map, norm_axis, "raw"))
        stats.update(map_stats(loo_map - template, norm_axis, "residual"))
        if loo_keys is None:
            loo_keys = list(stats)
        loo_values.append([stats[key] for key in loo_keys])
    loo_values = np.array(loo_values)

    result = {
        "pairs": str(args.pairs),
        "kappa_map": str(args.kappa_map),
        "single": str(args.single),
        "rperp_center_hmpc": args.rperp_center,
        "rpar_max_hmpc": args.rpar_max,
        "rpar_min_hmpc": args.rpar_min,
        "rpar_space": args.rpar_space,
        "n_pairs": int(counts.sum()),
        "n_blocks": k,
        "blocks_per_side": args.blocks_per_side,
        "stats": {},
    }
    for j, key in enumerate(loo_keys):
        _, err = jackknife_stats(loo_values[:, j])
        result["stats"][key] = {"value": full_stats[key], "jackknife_error": err}

    if args.stack_output:
        save_map_csv(args.stack_output, reflect_symmetrize_map(full_map.astype(np.float32)))
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pair stack with spatial-jackknife bridge errors.")
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--kappa-map", required=True, help="Pre-smoothed map; no smoothing applied here.")
    parser.add_argument("--single", required=True,
                        help="Matched single stack CSV (same map/smoothing as --kappa-map).")
    parser.add_argument("--rperp-center", type=float, required=True)
    parser.add_argument("--rpar-max", type=float, default=None,
                        help="Optional upper LOS cut, applied in --rpar-space.")
    parser.add_argument("--rpar-min", type=float, default=None,
                        help="Optional lower LOS cut. Combined with --rpar-max and "
                             "--rpar-space real this selects NON-PHYSICAL pairs: close "
                             "in projection but far apart along the line of sight, so "
                             "they are not physically associated. Their residual bridge "
                             "excess should be zero if the two-halo template is an "
                             "unbiased control.")
    parser.add_argument("--rpar-space", choices=["redshift", "real"], default="redshift",
                        help="Which LOS separation the cuts apply to (default: redshift, "
                             "matching the production/BOSS-like selection).")
    parser.add_argument("--blocks-per-side", type=int, default=5)
    parser.add_argument("--grid-size", type=int, default=101)
    parser.add_argument("--box-size", type=float, default=100.0)
    parser.add_argument("--max-pairs", type=int, default=None, help="Pilot cap.")
    parser.add_argument("--output", required=True, help="Output JSON with stats +/- jackknife errors.")
    parser.add_argument("--stack-output", default=None, help="Optional full stacked-map CSV.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    result = run(args)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
    r = result["stats"]["residual_bridge_excess_kappa"]
    logger.info(
        "residual bridge excess = %.4e +/- %.4e (S/N %.2f) -> %s",
        r["value"], r["jackknife_error"],
        r["value"] / r["jackknife_error"] if r["jackknife_error"] > 0 else float("nan"),
        out,
    )


if __name__ == "__main__":
    main()
