#!/usr/bin/env python
"""Stack galaxy pairs (or random pairs) from a pair catalog against the Planck kappa map.

This is a standalone script version of the legacy notebook
``legacy/stack_pairs_from_pair_catalog.ipynb``.  It reads a pair catalog CSV,
loads the Planck CMB lensing convergence map, orients a grid along each pair
axis, samples kappa at each grid position, and writes the reflection-
symmetrized mean stacked kappa map to a CSV file.

Usage examples
--------------
    python stack_pairs.py --pair-catalog data/paircatalogs/BOSS/galaxy_pairs_BOSS_South_20.0_18.0_22.0hmpc.csv --label galaxy_20 --region South
    python stack_pairs.py --pair-catalog data/paircatalogs/BOSS/galaxy_pairs_BOSS_North_20.0_18.0_22.0hmpc.csv --label galaxy_20 --region North
"""

from __future__ import annotations

import argparse
import logging
import os
import time

import healpy as hp
import numpy as np
import pandas as pd
from tqdm import tqdm

from catalog import load_kappa_map, resolve_planck_paths, setup_logging
from constants import BOX_SIZE_HMPC
from geometry import reflect_symmetrize_map

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Grid constants (pair stacking uses 101x101, NOT 100 like single stacking)
# ---------------------------------------------------------------------------
GRID_RES: int = 101
HALF_SIZE: float = BOX_SIZE_HMPC / 2.0


# ---------------------------------------------------------------------------
# Core stacking routine
# ---------------------------------------------------------------------------
def stack_pairs(
    pairs: pd.DataFrame,
    kmap: np.ndarray,
    mask: np.ndarray,
    nside: int,
) -> np.ndarray:
    """Stack all pairs from *pairs* onto a 2-D kappa grid.

    Parameters
    ----------
    pairs : DataFrame
        Must contain columns l1, b1, l2, b2 (degrees), w1, w2 (weights),
        and Dmid (midpoint comoving distance in h^-1 Mpc).
    kmap : array
        HEALPix convergence map.
    mask : array
        HEALPix mask (0 = masked).
    nside : int
        HEALPix nside of *kmap* and *mask*.

    Returns
    -------
    kappa_mean : ndarray, shape (GRID_RES, GRID_RES)
        Reflection-symmetrized mean stacked kappa map.
    """

    # --- Validate required columns ---
    required = {'l1', 'b1', 'l2', 'b2', 'w1', 'w2', 'Dmid'}
    missing = required - set(pairs.columns)
    if missing:
        raise KeyError(f"Pair catalog missing columns: {missing}")

    # --- Pre-extract columns as numpy arrays (avoids iterrows overhead) ---
    l1_arr = np.deg2rad(pairs['l1'].values)
    b1_arr = np.deg2rad(pairs['b1'].values)
    l2_arr = np.deg2rad(pairs['l2'].values)
    b2_arr = np.deg2rad(pairs['b2'].values)
    Dmid_arr = pairs['Dmid'].values
    w1_arr = pairs['w1'].values
    w2_arr = pairs['w2'].values

    X_vals = np.linspace(-HALF_SIZE, HALF_SIZE, GRID_RES)
    Y_vals = np.linspace(-HALF_SIZE, HALF_SIZE, GRID_RES)
    X_grid, Y_grid = np.meshgrid(X_vals, Y_vals)

    kappa_stack_sum = np.zeros((GRID_RES, GRID_RES))
    weight_stack_sum = np.zeros((GRID_RES, GRID_RES))

    n_pairs = len(pairs)
    n_skipped = 0

    for i in tqdm(range(n_pairs), desc="Stacking pairs"):
        l1 = l1_arr[i]
        b1 = b1_arr[i]
        l2 = l2_arr[i]
        b2 = b2_arr[i]

        # --- Enforce consistent ordering of longitudes ---
        dl_raw = (l2 - l1 + np.pi) % (2 * np.pi) - np.pi

        if dl_raw < 0:
            l1, l2 = l2, l1
            b1, b2 = b2, b1
            dl_raw = -dl_raw

        # --- Midpoint (wrap-safe) ---
        lc = l1 + 0.5 * dl_raw
        bc = 0.5 * (b1 + b2)
        Dc = Dmid_arr[i]

        # --- Rotation angle theta ---
        dl = dl_raw * np.cos(bc)
        db = b2 - b1
        norm = np.hypot(dl, db)
        cos_theta = dl / norm
        sin_theta = db / norm

        # Inverse transform (X, Y) -> (l, b) offsets
        dl_cosbc = (cos_theta * (X_grid / Dc)) - (sin_theta * (Y_grid / Dc))
        db_grid = (sin_theta * (X_grid / Dc)) + (cos_theta * (Y_grid / Dc))
        l_grid = lc + dl_cosbc / np.cos(bc)
        b_grid = bc + db_grid

        # Convert to degrees
        l_grid_deg = np.rad2deg(l_grid)
        b_grid_deg = np.rad2deg(b_grid)

        # Convert to HEALPix pixel coordinates
        l_grid_deg_wrapped = np.mod(l_grid_deg, 360)
        theta = np.radians(90.0 - b_grid_deg)
        phi = np.radians(l_grid_deg_wrapped)

        if np.min(theta) < 0:
            logger.warning(
                "Pair %d: theta < 0 encountered (min=%.6f). Skipping.", i, np.min(theta)
            )
            n_skipped += 1
            continue
        if np.max(theta) > np.pi:
            logger.warning(
                "Pair %d: theta > pi encountered (max=%.6f). Skipping.", i, np.max(theta)
            )
            n_skipped += 1
            continue

        pix = hp.ang2pix(nside, theta.ravel(), phi.ravel())

        # Sample kappa and mask
        kappa_vals = kmap[pix].reshape(GRID_RES, GRID_RES)
        mask_vals = mask[pix].reshape(GRID_RES, GRID_RES)

        valid_mask = (mask_vals != 0) & np.isfinite(kappa_vals)
        kappa_vals[~valid_mask] = 0.0

        weight = w1_arr[i] * w2_arr[i]

        kappa_stack_sum += weight * kappa_vals
        weight_stack_sum += weight * valid_mask.astype(float)

    # --- Mean map (avoid division by zero) ---
    nonzero = weight_stack_sum > 0
    kappa_stack_mean = np.zeros_like(kappa_stack_sum)
    kappa_stack_mean[nonzero] = kappa_stack_sum[nonzero] / weight_stack_sum[nonzero]

    logger.info(
        "Stacking complete. Pairs processed: %d, skipped: %d",
        len(pairs) - n_skipped,
        n_skipped,
    )

    # --- Apply reflection symmetry ---
    kappa_sym = reflect_symmetrize_map(kappa_stack_mean)

    return kappa_sym


# ---------------------------------------------------------------------------
# Auto-generate output filename
# ---------------------------------------------------------------------------
def make_output_path(output_dir: str, label: str, dataset: str, region: str) -> str:
    """Return the default output CSV path."""
    fname = f"kappa_pairs_{label}_{dataset}_{region}.csv"
    return os.path.join(output_dir, fname)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stack galaxy pairs from a pair catalog against the Planck kappa map.",
    )
    parser.add_argument(
        "--pair-catalog",
        required=True,
        help="Path to pair catalog CSV (must contain l1, b1, z1, l2, b2, z2 columns).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Explicit output path for the kappa CSV. Overrides auto-naming.",
    )
    parser.add_argument(
        "--output-dir",
        default="analysis/boss/results",
        help="Output directory (default: analysis/boss/results).",
    )
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Root data directory (default: data).",
    )
    parser.add_argument(
        "--label",
        default="galaxy",
        help='Label for auto-generated output filename, e.g. "galaxy_20" or "random_20".',
    )
    parser.add_argument(
        "--dataset",
        default="BOSS",
        choices=["BOSS", "eBOSS"],
        help="Dataset name used for auto-naming (default: BOSS).",
    )
    parser.add_argument(
        "--region",
        default="South",
        help="Region identifier for output naming, e.g. North/South (default: South).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file.",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()

    t0 = time.time()

    # ------------------------------------------------------------------
    # Resolve output path
    # ------------------------------------------------------------------
    if args.output:
        output_path = args.output
    else:
        output_path = make_output_path(args.output_dir, args.label, args.dataset, args.region)

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    if os.path.exists(output_path) and not args.overwrite:
        logger.info("Output already exists: %s  (use --overwrite to replace). Skipping.", output_path)
        return

    # ------------------------------------------------------------------
    # Load pair catalog
    # ------------------------------------------------------------------
    logger.info("Loading pair catalog: %s", args.pair_catalog)
    pairs = pd.read_csv(args.pair_catalog)
    logger.info("Loaded %d pairs.", len(pairs))

    # ------------------------------------------------------------------
    # Load Planck kappa map
    # ------------------------------------------------------------------
    alm_file, mask_file = resolve_planck_paths(args.data_dir)
    kmap, mask, nside = load_kappa_map(alm_file=alm_file, mask_file=mask_file)

    # ------------------------------------------------------------------
    # Stack
    # ------------------------------------------------------------------
    kappa_sym = stack_pairs(pairs, kmap, mask, nside)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    pd.DataFrame(kappa_sym).to_csv(output_path, index=True)
    logger.info("Saved stacked kappa map to %s", output_path)

    elapsed = time.time() - t0
    logger.info("Total elapsed time: %.1f s (%.1f min)", elapsed, elapsed / 60.0)


if __name__ == "__main__":
    main()
