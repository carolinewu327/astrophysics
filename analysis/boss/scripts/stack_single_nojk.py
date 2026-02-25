#!/usr/bin/env python
"""
stack_single_nojk.py -- Stack single galaxies or random points against the
Planck kappa (convergence) map, WITHOUT jackknife resampling.

Uses the full-sample stack for galaxies (consistent with pair stacking),
avoiding the jackknife valid_mask that can drop galaxies outside good
HEALPix regions.

Usage examples:
    python stack_single_nojk.py --dataset BOSS --region North --catalog-type galaxy
    python stack_single_nojk.py --dataset BOSS --region North --catalog-type random --fraction 0.10
"""

import argparse
import logging
import os
import time

import numpy as np
import pandas as pd
import healpy as hp
from tqdm import tqdm

from catalog import (
    load_catalog, load_catalog_lightweight, load_kappa_map,
    preprocess_catalog_galactic,
    resolve_catalog_path, resolve_planck_paths,
)
from constants import BOX_SIZE_HMPC, GRID_SIZE, NSIDE
from geometry import symmetrize_map


# ===========================================================================
# Derived constants
# ===========================================================================
CELL_SIZE_HMPC = BOX_SIZE_HMPC / GRID_SIZE
HALF_BOX_HMPC = BOX_SIZE_HMPC / 2

# Pre-compute the grid of offsets (in h^-1 Mpc) -- same for every galaxy
OFFSETS = np.linspace(
    -HALF_BOX_HMPC + CELL_SIZE_HMPC / 2,
    HALF_BOX_HMPC - CELL_SIZE_HMPC / 2,
    GRID_SIZE,
)
OFF_X, OFF_Y = np.meshgrid(OFFSETS, OFFSETS)
OFF_X = OFF_X.ravel()
OFF_Y = OFF_Y.ravel()

logger = logging.getLogger(__name__)


# ===========================================================================
# Core stacking function
# ===========================================================================
def stack_kappa(data, weights, kmap, mask, label="galaxies"):
    """Stack the kappa map at galaxy positions on a projected grid.

    For each galaxy the code projects a 100x100 grid of physical offsets
    (in h^-1 Mpc) onto the sky, converts to HEALPix pixels, reads the
    convergence value, and accumulates weighted averages.

    Parameters
    ----------
    data : structured ndarray
        Galaxy catalog (must contain RA, DEC, Z columns).
    weights : ndarray
        Per-galaxy weights (same length as *data*).
    kmap : ndarray
        HEALPix convergence map.
    mask : ndarray
        HEALPix lensing mask.
    label : str
        Label used in log messages and progress bar.

    Returns
    -------
    kappa_mean : 2-D ndarray  (GRID_SIZE x GRID_SIZE)
        Weighted mean convergence map.
    kappa_sigma : 2-D ndarray
        Per-pixel standard error of the mean.
    kappa_sn : 2-D ndarray
        Signal-to-noise map (kappa_mean / kappa_sigma).
    """
    sz = GRID_SIZE ** 2
    sum_wk = np.zeros(sz)
    sum_wk2 = np.zeros(sz)
    sum_w = np.zeros(sz)

    # Convert catalog to Galactic coordinates + comoving distances
    l_arr, b_arr, D_arr, data_valid, weights_valid = preprocess_catalog_galactic(data, weights)

    n_gal = len(data_valid)
    logger.info("Stacking %d %s", n_gal, label)

    for i in tqdm(range(n_gal), desc=f"Stacking {label}"):
        l_gal = l_arr[i]
        b_gal = b_arr[i]
        D_gal = D_arr[i]

        # Project physical offsets onto angular offsets (degrees)
        cosb = np.cos(np.radians(b_gal))
        dl = (OFF_X / D_gal) * (180.0 / np.pi) / np.clip(cosb, 1e-6, None)
        db = (OFF_Y / D_gal) * (180.0 / np.pi)

        l_grid = l_gal + dl
        b_grid = b_gal + db

        # Convert to HEALPix theta/phi (colatitude, longitude in radians)
        theta = np.radians(90.0 - b_grid)
        theta = np.clip(theta, 0.0, np.pi)
        phi = np.radians(l_grid)

        pix = hp.ang2pix(NSIDE, theta, phi)

        # Weight = galaxy weight * mask value at each grid pixel
        w = weights_valid[i] * mask[pix]
        if np.sum(w) == 0:
            continue

        kappa_vals = kmap[pix]
        sum_wk += w * kappa_vals
        sum_wk2 += w * kappa_vals ** 2
        sum_w += w

        # Progress logging every 1000 galaxies
        if (i + 1) % 1000 == 0:
            logger.info("  stacked %d / %d %s", i + 1, n_gal, label)

    # Compute weighted mean, variance, and S/N
    valid = sum_w > 0
    kappa_mean = np.zeros_like(sum_w)
    kappa_var = np.zeros_like(sum_w)
    kappa_mean[valid] = sum_wk[valid] / sum_w[valid]
    kappa_var[valid] = sum_wk2[valid] / sum_w[valid] - kappa_mean[valid] ** 2

    kappa_sigma = np.zeros_like(kappa_var)
    kappa_sigma[valid] = np.sqrt(kappa_var[valid]) / np.sqrt(sum_w[valid])

    kappa_sn = np.zeros_like(kappa_mean)
    pos = valid & (kappa_sigma > 0)
    kappa_sn[pos] = kappa_mean[pos] / kappa_sigma[pos]

    return (
        kappa_mean.reshape(GRID_SIZE, GRID_SIZE),
        kappa_sigma.reshape(GRID_SIZE, GRID_SIZE),
        kappa_sn.reshape(GRID_SIZE, GRID_SIZE),
    )


# ===========================================================================
# Output helpers
# ===========================================================================
def output_kappa_path(output_dir, catalog_type, dataset, region):
    return os.path.join(output_dir, f"kappa_single_{catalog_type}_{dataset}_{region}_nojk.csv")


# ===========================================================================
# Main entry point
# ===========================================================================
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Stack single galaxies (or randoms) against the Planck kappa map (no jackknife)."
    )
    parser.add_argument(
        "--dataset", default="BOSS", choices=["BOSS", "eBOSS"],
        help="Galaxy survey dataset (default: BOSS)."
    )
    parser.add_argument(
        "--region", required=True,
        help="Survey region (e.g. North, South for BOSS; NGC, SGC for eBOSS)."
    )
    parser.add_argument(
        "--catalog-type", required=True, choices=["galaxy", "random"],
        help="Type of catalog to stack: 'galaxy' or 'random'."
    )
    parser.add_argument(
        "--fraction", type=float, default=1.0,
        help="Fraction of catalog to use (default: 1.0; use 0.10 for randoms)."
    )
    parser.add_argument(
        "--output-dir", default="analysis/boss/results",
        help="Directory for output CSV files (default: analysis/boss/results)."
    )
    parser.add_argument(
        "--data-dir", default="data",
        help="Root data directory (default: data)."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite existing output files."
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # ------------------------------------------------------------------
    # Configure logging
    # ------------------------------------------------------------------
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger.info("=== stack_single_nojk.py ===")
    logger.info("Dataset   : %s", args.dataset)
    logger.info("Region    : %s", args.region)
    logger.info("Catalog   : %s", args.catalog_type)
    logger.info("Fraction  : %.4f", args.fraction)
    logger.info("Output dir: %s", args.output_dir)
    logger.info("Data dir  : %s", args.data_dir)
    logger.info("Overwrite : %s", args.overwrite)

    # ------------------------------------------------------------------
    # Check output paths early
    # ------------------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)

    kappa_out = output_kappa_path(args.output_dir, args.catalog_type,
                                  args.dataset, args.region)

    if not args.overwrite:
        if os.path.exists(kappa_out):
            logger.info("Output already exists: %s  (use --overwrite to rerun)", kappa_out)
            return

    # ------------------------------------------------------------------
    # Resolve file paths
    # ------------------------------------------------------------------
    catalog_path = resolve_catalog_path(
        args.data_dir, args.dataset, args.region, args.catalog_type
    )
    alm_path, mask_path = resolve_planck_paths(args.data_dir)

    logger.info("Catalog file : %s", catalog_path)
    logger.info("Planck ALM   : %s", alm_path)
    logger.info("Planck mask  : %s", mask_path)

    # ------------------------------------------------------------------
    # Dataset-specific settings
    # ------------------------------------------------------------------
    if args.dataset == "BOSS":
        z_min, z_max = 0.4, 0.7
        weight_scheme = "CMASS" if args.catalog_type == "galaxy" else False
    else:  # eBOSS
        z_min, z_max = 0.0, 10000.0
        weight_scheme = True if args.catalog_type == "galaxy" else False

    # For random catalogs, apply the --fraction argument
    random_fraction = args.fraction if args.catalog_type == "random" else None
    # For galaxy catalogs with fraction < 1, also sub-sample
    if args.catalog_type == "galaxy" and args.fraction < 1.0:
        random_fraction = args.fraction

    # ------------------------------------------------------------------
    # Load Planck kappa map + mask
    # ------------------------------------------------------------------
    t0 = time.time()
    kmap, planck_mask, _ = load_kappa_map(alm_file=alm_path, mask_file=mask_path)
    logger.info("Planck map loaded in %.1f s", time.time() - t0)

    # ------------------------------------------------------------------
    # Load catalog
    # ------------------------------------------------------------------
    t0 = time.time()
    logger.info("Loading catalog: %s", catalog_path)
    if args.catalog_type == "random":
        data = load_catalog_lightweight(
            catalog_path, columns=("RA", "DEC", "Z"),
            fraction=random_fraction,
            z_min=z_min, z_max=z_max,
        )
        weights = np.ones(len(data))
    else:
        data, weights = load_catalog(
            catalog_path,
            weights=weight_scheme,
            random_fraction=random_fraction,
            z_min=z_min,
            z_max=z_max,
        )
    logger.info("Loaded %d objects in %.1f s", len(data), time.time() - t0)

    # ------------------------------------------------------------------
    # Stack (full sample, no jackknife)
    # ------------------------------------------------------------------
    t_start = time.time()

    kappa_map, _, _ = stack_kappa(
        data, weights, kmap, planck_mask,
        label=args.catalog_type,
    )

    elapsed = time.time() - t_start
    logger.info("Stacking completed in %.1f s (%.1f min)", elapsed, elapsed / 60)

    # ------------------------------------------------------------------
    # Symmetrize
    # ------------------------------------------------------------------
    kappa_map = symmetrize_map(kappa_map)
    logger.info("Applied radial symmetrization to kappa map")

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    pd.DataFrame(kappa_map).to_csv(kappa_out, index=True)
    logger.info("Saved kappa map -> %s", kappa_out)

    logger.info("Done.")


if __name__ == "__main__":
    main()
