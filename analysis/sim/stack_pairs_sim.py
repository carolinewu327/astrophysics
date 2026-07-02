#!/usr/bin/env python
"""Stack a simulated kappa map around periodic BigMDPL halo pairs."""

from __future__ import annotations

import argparse
import logging
import os
import time

import numpy as np
import pandas as pd
from tqdm import tqdm

from sim_utils import (
    OBS_BOX_SIZE_HMPC,
    OBS_PAIR_GRID_SIZE,
    apply_periodic_gaussian_smoothing,
    open_kappa_memmap,
    pair_stack_offsets,
    periodic_bilinear_sample,
    reflect_symmetrize_map,
    save_map_csv,
    setup_logging,
)


logger = logging.getLogger(__name__)


def stack_pairs(
    pairs: pd.DataFrame,
    kappa_map: np.ndarray,
    pixel_size_hmpc: float,
    map_box_size_hmpc: float,
    stack_box_size_hmpc: float,
    grid_size: int,
    normalize_separation: bool,
    normalized_half_size: float,
) -> np.ndarray:
    if normalize_separation:
        axis = np.linspace(-normalized_half_size, normalized_half_size, grid_size)
        x_grid, y_grid = np.meshgrid(axis, axis)
    else:
        _, x_grid, y_grid = pair_stack_offsets(stack_box_size_hmpc, grid_size)

    stack_sum = np.zeros((grid_size, grid_size), dtype=np.float64)
    n_used = 0

    center_x = pairs["pair_center_x"].to_numpy(dtype=np.float64, copy=False)
    center_y = pairs["pair_center_y"].to_numpy(dtype=np.float64, copy=False)
    cos_t = pairs["cos_theta"].to_numpy(dtype=np.float64, copy=False)
    sin_t = pairs["sin_theta"].to_numpy(dtype=np.float64, copy=False)
    rperp = pairs["r_perp"].to_numpy(dtype=np.float64, copy=False)

    for i in tqdm(range(len(pairs)), desc="Stacking halo pairs"):
        if normalize_separation:
            local_x = x_grid * rperp[i]
            local_y = y_grid * rperp[i]
        else:
            local_x = x_grid
            local_y = y_grid

        sim_x = center_x[i] + cos_t[i] * local_x - sin_t[i] * local_y
        sim_y = center_y[i] + sin_t[i] * local_x + cos_t[i] * local_y
        vals = periodic_bilinear_sample(
            kappa_map,
            sim_x,
            sim_y,
            pixel_size_hmpc=pixel_size_hmpc,
            box_size_hmpc=map_box_size_hmpc,
        )
        if not np.all(np.isfinite(vals)):
            continue
        stack_sum += vals
        n_used += 1

    if n_used == 0:
        logger.warning("No valid pairs were stacked")
        return np.zeros((grid_size, grid_size), dtype=np.float32)

    mean_map = (stack_sum / n_used).astype(np.float32)
    return reflect_symmetrize_map(mean_map).astype(np.float32)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stack simulated kappa around halo pairs.")
    parser.add_argument("--pairs", default="analysis/sim/results/pairs_mass13_rperp5.csv")
    parser.add_argument("--kappa-map", default="analysis/sim/results/kappa_map_l0p5.float32")
    parser.add_argument("--output", default="analysis/sim/results/kappa_pairs_sim_mass13_rperp5.csv")
    parser.add_argument("--box-size", type=float, default=OBS_BOX_SIZE_HMPC)
    parser.add_argument("--grid-size", type=int, default=OBS_PAIR_GRID_SIZE)
    parser.add_argument("--normalize-separation", action="store_true")
    parser.add_argument(
        "--normalized-half-size",
        type=float,
        default=2.5,
        help="Half-width in X/r_perp and Y/r_perp when --normalize-separation is set.",
    )
    parser.add_argument("--smooth", choices=["none", "2arcmin", "4arcmin", "8arcmin"], default="none")
    parser.add_argument("--max-pairs", type=int, default=None, help="Optional pair cap for smoke tests.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    t0 = time.time()

    if os.path.exists(args.output) and not args.overwrite:
        logger.info("Output exists: %s (use --overwrite)", args.output)
        return

    pairs = pd.read_csv(args.pairs)
    if args.max_pairs is not None:
        pairs = pairs.iloc[: args.max_pairs].copy()
    logger.info("Loaded %d pairs", len(pairs))

    kappa_map, info = open_kappa_memmap(args.kappa_map)
    map_for_stack = apply_periodic_gaussian_smoothing(kappa_map, info.pixel_size_hmpc, args.smooth)

    stacked = stack_pairs(
        pairs=pairs,
        kappa_map=map_for_stack,
        pixel_size_hmpc=info.pixel_size_hmpc,
        map_box_size_hmpc=info.box_size_hmpc,
        stack_box_size_hmpc=args.box_size,
        grid_size=args.grid_size,
        normalize_separation=args.normalize_separation,
        normalized_half_size=args.normalized_half_size,
    )

    save_map_csv(args.output, stacked)
    logger.info("Saved pair stack -> %s in %.1f s", args.output, time.time() - t0)


if __name__ == "__main__":
    main()

