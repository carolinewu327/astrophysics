#!/usr/bin/env python
"""
find_pairs.py - Standalone script to find galaxy (or random point) pairs
based on parallel and perpendicular separation criteria.

Ported from legacy/find_galaxy_pairs_optimized.ipynb.

Usage examples:
    # BOSS galaxy pairs, default separations
    python find_pairs.py --dataset BOSS --region North --catalog-type galaxy

    # BOSS random pairs with 10% subsample
    python find_pairs.py --dataset BOSS --region North --catalog-type random --fraction 0.10

    # Custom separation range
    python find_pairs.py --dataset BOSS --region South --catalog-type galaxy \
        --rpar 25 --rperp-min 10 --rperp-max 15
"""

import argparse
import logging
import os
import time

import numpy as np
import pandas as pd
from astropy.cosmology import Planck18 as cosmo
from functools import partial
from multiprocessing import Pool, cpu_count

from catalog import (
    load_catalog, load_catalog_lightweight,
    preprocess_catalog_galactic, resolve_catalog_path,
)
from geometry import angular_separation

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Module-level function (required for multiprocessing pickle)
# ---------------------------------------------------------------------------
def process_chunk_optimized(chunk_start_end, l, b, D, z, weights_valid, ids,
                            r_par_max, r_perp_min, r_perp_max,
                            sorted_indices, sorted_D):
    """
    Process a chunk of galaxies using the two-pass vectorized approach:
      Pass 1 (rough): approximate angular filter with generous margins
      Pass 2 (precise): exact arccos-based angular separation

    Parameters
    ----------
    chunk_start_end : tuple of (int, int)
        Start and end indices into the sorted_indices array for this chunk.
    l, b : ndarray
        Galactic longitude and latitude (degrees) for all galaxies.
    D : ndarray
        Comoving distance in h^-1 Mpc for all galaxies.
    z : ndarray
        Redshift for all galaxies.
    weights_valid : ndarray
        Weight for each galaxy.
    ids : ndarray
        Integer ID for each galaxy.
    r_par_max : float
        Maximum parallel (line-of-sight) distance in Mpc/h.
    r_perp_min : float
        Minimum perpendicular (transverse) distance in Mpc/h.
    r_perp_max : float
        Maximum perpendicular (transverse) distance in Mpc/h.
    sorted_indices : ndarray
        Indices that sort galaxies by comoving distance.
    sorted_D : ndarray
        Comoving distances sorted in ascending order.

    Returns
    -------
    pairs : list of dict
        Each dict has keys: l1, b1, z1, w1, Dc1, ID1, l2, b2, z2, w2, Dc2, ID2, Dmid
    """
    chunk_start, chunk_end = chunk_start_end

    # Collect candidate pairs (first pass: rough angular filter)
    i_list = []
    j_list = []

    start_time = time.time()
    last_print_time = start_time

    for idx in range(chunk_start, chunk_end):
        i = sorted_indices[idx]
        Dc1 = D[i]

        # Progress logging every 1000 galaxies within the chunk
        galaxies_done = idx - chunk_start
        if galaxies_done > 0 and galaxies_done % 1000 == 0:
            current_time = time.time()
            interval = current_time - last_print_time
            rate = 1000 / interval if interval > 0 else 0
            print(
                f"  Chunk [{chunk_start}-{chunk_end}]: "
                f"{galaxies_done}/{chunk_end - chunk_start} galaxies "
                f"({rate:.1f} gal/s, {len(i_list)} candidate pairs so far)"
            )
            last_print_time = current_time

        # Binary search: galaxies within r_par_max in comoving distance
        left_bound = np.searchsorted(sorted_D, Dc1 - r_par_max, side='left')
        right_bound = np.searchsorted(sorted_D, Dc1 + r_par_max, side='right')
        search_start = max(idx + 1, left_bound)  # avoid self-pairs and duplicates
        search_end = right_bound

        if search_end <= search_start:
            continue

        j_candidates = sorted_indices[search_start:search_end]
        Dc2_candidates = D[j_candidates]

        # --- First pass: rough angular filter ---
        l1_val, b1_val = l[i], b[i]
        l2_candidates = l[j_candidates]
        b2_candidates = b[j_candidates]

        dl = np.abs(l2_candidates - l1_val)
        db = np.abs(b2_candidates - b1_val)
        # Handle longitude wrapping at 360 degrees
        dl = np.minimum(dl, 360.0 - dl)

        D_avg = (Dc1 + Dc2_candidates) / 2.0
        min_angle = r_perp_min / D_avg
        max_angle = r_perp_max / D_avg

        # cosb correction: without this, rough_angle overestimates true angle
        # at high latitudes, potentially causing valid pairs to be missed
        cosb1 = np.cos(np.radians(b1_val))
        rough_angle = np.sqrt((np.radians(dl * cosb1)) ** 2 + (np.radians(db)) ** 2)

        # Generous margins: 0.7x on min, 2x on max
        angle_mask = (rough_angle >= min_angle * 0.7) & (rough_angle <= max_angle * 2.0)

        if not np.any(angle_mask):
            continue

        j_filtered = j_candidates[angle_mask]
        i_list.extend([i] * len(j_filtered))
        j_list.extend(j_filtered.tolist())

    if len(i_list) == 0:
        total_time = time.time() - start_time
        print(
            f"  Chunk [{chunk_start}-{chunk_end}] COMPLETE: "
            f"{chunk_end - chunk_start} galaxies in {total_time:.1f}s, found 0 pairs"
        )
        return []

    # --- Second pass: precise angular separation ---
    i_arr = np.array(i_list)
    j_arr = np.array(j_list)

    theta = angular_separation(l[i_arr], b[i_arr], l[j_arr], b[j_arr])

    Dc1_arr = D[i_arr]
    Dc2_arr = D[j_arr]
    D_avg_arr = (Dc1_arr + Dc2_arr) / 2.0
    r_perp_arr = D_avg_arr * theta

    final_mask = (r_perp_arr >= r_perp_min) & (r_perp_arr <= r_perp_max)

    i_final = i_arr[final_mask]
    j_final = j_arr[final_mask]

    # Build output pairs
    pairs = []
    for k in range(len(i_final)):
        ii = i_final[k]
        jj = j_final[k]
        Dmid = cosmo.comoving_distance((z[ii] + z[jj]) / 2.0).value * cosmo.h
        pairs.append({
            'l1': l[ii], 'b1': b[ii], 'z1': z[ii],
            'w1': weights_valid[ii], 'Dc1': D[ii], 'ID1': ids[ii],
            'l2': l[jj], 'b2': b[jj], 'z2': z[jj],
            'w2': weights_valid[jj], 'Dc2': D[jj], 'ID2': ids[jj],
            'Dmid': Dmid,
        })

    total_time = time.time() - start_time
    print(
        f"  Chunk [{chunk_start}-{chunk_end}] COMPLETE: "
        f"{chunk_end - chunk_start} galaxies in {total_time:.1f}s, found {len(pairs)} pairs"
    )
    return pairs


# ---------------------------------------------------------------------------
# Main pair-finding routine
# ---------------------------------------------------------------------------
def build_pair_catalog(data, weights, catalog_type, dataset, region,
                       r_par_max, r_perp_min, r_perp_max,
                       n_processes, chunk_size, output_dir):
    """
    Build a pair catalog from galaxy or random-point data.

    Parameters
    ----------
    data : FITS record array
        Loaded catalog data (must have RA, DEC, Z columns).
    weights : ndarray
        Per-object weights.
    catalog_type : str
        "galaxy" or "random".
    dataset : str
        "BOSS" or "eBOSS".
    region : str
        Region identifier (e.g. "North", "South").
    r_par_max : float
        Max parallel distance in Mpc/h.
    r_perp_min : float
        Min perpendicular distance in Mpc/h.
    r_perp_max : float
        Max perpendicular distance in Mpc/h.
    n_processes : int
        Number of worker processes for multiprocessing.
    chunk_size : int
        Galaxies per processing chunk.
    output_dir : str
        Directory to save pair catalog and checkpoints.

    Returns
    -------
    pairs_df : pd.DataFrame
        DataFrame with columns l1, b1, z1, w1, Dc1, ID1,
        l2, b2, z2, w2, Dc2, ID2, Dmid.
    """
    # Preprocess: convert RA/Dec to Galactic, compute comoving distances
    l, b, D, data_filtered, weights_valid = preprocess_catalog_galactic(data, weights)
    z = np.array(data_filtered['Z'])

    # Handle IDs: BOSS CMASS has 'ID' field; randoms may not
    try:
        ids = np.array(data_filtered['ID'])
    except (KeyError, ValueError):
        logger.info("No 'ID' column found in catalog; generating sequential IDs.")
        ids = np.arange(len(data_filtered))

    n_galaxies = len(data_filtered)
    logger.info(f"Catalog loaded: {n_galaxies} objects after preprocessing")
    logger.info(
        f"Search parameters: r_par_max={r_par_max}, "
        f"r_perp_min={r_perp_min}, r_perp_max={r_perp_max}"
    )

    # Sort by comoving distance for efficient binary search
    logger.info("Sorting by comoving distance...")
    sorted_indices = np.argsort(D)
    sorted_D = D[sorted_indices]

    avg_D = np.median(D)
    logger.info(f"Median comoving distance: {avg_D:.1f} h^-1 Mpc")

    # Split into chunks
    n_chunks = int(np.ceil(n_galaxies / chunk_size))
    chunk_ranges = [
        (i * chunk_size, min((i + 1) * chunk_size, n_galaxies))
        for i in range(n_chunks)
    ]
    logger.info(f"Processing {n_chunks} chunks of size {chunk_size} using {n_processes} process(es)")

    # Output file path
    output_file = os.path.join(
        output_dir,
        f"{catalog_type}_pairs_{dataset}_{region}"
        f"_{r_par_max}_{r_perp_min}_{r_perp_max}hmpc.csv"
    )

    # Create partial function with fixed parameters
    process_func = partial(
        process_chunk_optimized,
        l=l, b=b, D=D, z=z,
        weights_valid=weights_valid, ids=ids,
        r_par_max=r_par_max, r_perp_min=r_perp_min, r_perp_max=r_perp_max,
        sorted_indices=sorted_indices, sorted_D=sorted_D,
    )

    all_pairs = []
    global_start = time.time()

    if n_processes > 1:
        logger.info(
            f"Starting multiprocessing pool with {n_processes} workers..."
        )
        chunk_start_time = time.time()
        with Pool(processes=n_processes) as pool:
            for i, result in enumerate(pool.imap(process_func, chunk_ranges)):
                chunk_time = time.time() - chunk_start_time
                all_pairs.extend(result)

                chunk_s, chunk_e = chunk_ranges[i]
                n_in_chunk = len(result)
                total_so_far = len(all_pairs)
                gal_per_sec = (chunk_e - chunk_s) / chunk_time if chunk_time > 0 else 0
                logger.info(
                    f"Chunk {i + 1}/{n_chunks} [{chunk_s}-{chunk_e}]: "
                    f"{n_in_chunk:,} pairs in {chunk_time:.1f}s "
                    f"({gal_per_sec:.1f} gal/s) | Total: {total_so_far:,} pairs"
                )

                # Checkpoint save every 10 chunks
                if (i + 1) % 10 == 0:
                    _save_checkpoint(all_pairs, output_file, i + 1)

                chunk_start_time = time.time()
    else:
        logger.info("Running in single-process mode...")
        for i, chunk in enumerate(chunk_ranges):
            chunk_start_time = time.time()
            result = process_func(chunk)
            chunk_time = time.time() - chunk_start_time
            all_pairs.extend(result)

            chunk_s, chunk_e = chunk
            n_in_chunk = len(result)
            total_so_far = len(all_pairs)
            gal_per_sec = (chunk_e - chunk_s) / chunk_time if chunk_time > 0 else 0
            logger.info(
                f"Chunk {i + 1}/{n_chunks} [{chunk_s}-{chunk_e}]: "
                f"{n_in_chunk:,} pairs in {chunk_time:.1f}s "
                f"({gal_per_sec:.1f} gal/s) | Total: {total_so_far:,} pairs"
            )

            # Checkpoint save every 10 chunks
            if (i + 1) % 10 == 0:
                _save_checkpoint(all_pairs, output_file, i + 1)

    elapsed = time.time() - global_start
    logger.info(f"Pair finding complete: {len(all_pairs):,} total pairs in {elapsed:.1f}s")

    # Save final output
    if len(all_pairs) == 0:
        logger.warning("No pairs found. Writing empty CSV with header only.")
        pairs_df = pd.DataFrame(
            columns=['l1', 'b1', 'z1', 'w1', 'Dc1', 'ID1',
                     'l2', 'b2', 'z2', 'w2', 'Dc2', 'ID2', 'Dmid']
        )
    else:
        pairs_df = pd.DataFrame(all_pairs)

    pairs_df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(pairs_df):,} pairs to {output_file}")

    return pairs_df


def _save_checkpoint(all_pairs, output_file, chunk_num):
    """Save an intermediate checkpoint CSV."""
    if len(all_pairs) == 0:
        return
    temp_df = pd.DataFrame(all_pairs)
    temp_file = output_file.replace('.csv', f'_checkpoint_{chunk_num}.csv')
    temp_df.to_csv(temp_file, index=False)
    logger.info(f"Checkpoint: saved {len(temp_df):,} pairs to {temp_file}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Find galaxy (or random) pairs by parallel and perpendicular separation."
    )
    parser.add_argument(
        "--dataset", type=str, default="BOSS", choices=["BOSS", "eBOSS"],
        help="Survey dataset (default: BOSS)"
    )
    parser.add_argument(
        "--region", type=str, default="North",
        help="Catalog region, e.g. North/South for BOSS, NGC/SGC for eBOSS (default: North)"
    )
    parser.add_argument(
        "--catalog-type", type=str, default="galaxy", choices=["galaxy", "random"],
        help="Type of catalog: 'galaxy' or 'random' (default: galaxy)"
    )
    parser.add_argument(
        "--fraction", type=float, default=1.0,
        help="Fraction of catalog to use; use 0.10 for randoms (default: 1.0)"
    )
    parser.add_argument(
        "--rpar", type=float, default=20.0,
        help="Max parallel (line-of-sight) distance in Mpc/h (default: 20)"
    )
    parser.add_argument(
        "--rperp-min", type=float, default=18.0,
        help="Min perpendicular (transverse) distance in Mpc/h (default: 18)"
    )
    parser.add_argument(
        "--rperp-max", type=float, default=22.0,
        help="Max perpendicular (transverse) distance in Mpc/h (default: 22)"
    )
    parser.add_argument(
        "--n-processes", type=int, default=None,
        help="Number of parallel processes (default: auto = 75%% of cores, capped at 12)"
    )
    parser.add_argument(
        "--chunk-size", type=int, default=10000,
        help="Galaxies per processing chunk (default: 10000)"
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help="Output directory for pair catalogs (default: data/paircatalogs/{dataset})"
    )
    parser.add_argument(
        "--data-dir", type=str, default="data",
        help="Root data directory (default: data)"
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite existing output file if it exists"
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

    logger.info("=" * 60)
    logger.info("find_pairs.py - Galaxy/Random Pair Finder")
    logger.info("=" * 60)
    logger.info(f"Dataset:      {args.dataset}")
    logger.info(f"Region:       {args.region}")
    logger.info(f"Catalog type: {args.catalog_type}")
    logger.info(f"Fraction:     {args.fraction}")
    logger.info(f"r_par_max:    {args.rpar} Mpc/h")
    logger.info(f"r_perp_min:   {args.rperp_min} Mpc/h")
    logger.info(f"r_perp_max:   {args.rperp_max} Mpc/h")

    # ------------------------------------------------------------------
    # Resolve output directory (default: data/paircatalogs/{dataset})
    # ------------------------------------------------------------------
    if args.output_dir is None:
        args.output_dir = os.path.join("data", "paircatalogs", args.dataset)

    # ------------------------------------------------------------------
    # Determine number of processes
    # ------------------------------------------------------------------
    if args.n_processes is not None:
        n_processes = args.n_processes
    else:
        n_processes = max(1, int(cpu_count() * 0.75))
        n_processes = min(n_processes, 12)
    logger.info(f"Processes:    {n_processes}")
    logger.info(f"Chunk size:   {args.chunk_size}")

    # ------------------------------------------------------------------
    # Check output file
    # ------------------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(
        args.output_dir,
        f"{args.catalog_type}_pairs_{args.dataset}_{args.region}"
        f"_{args.rpar}_{args.rperp_min}_{args.rperp_max}hmpc.csv"
    )
    if os.path.exists(output_file) and not args.overwrite:
        logger.info(f"Output file already exists: {output_file}")
        logger.info("Use --overwrite to replace it. Exiting.")
        return

    # ------------------------------------------------------------------
    # Resolve data path and load catalog
    # ------------------------------------------------------------------
    fits_path = resolve_catalog_path(
        args.data_dir, args.dataset, args.region, args.catalog_type
    )
    logger.info(f"Loading catalog from: {fits_path}")

    # Determine weighting and redshift cuts
    if args.dataset == "BOSS":
        z_min, z_max = 0.4, 0.7
        if args.catalog_type == "galaxy":
            weight_scheme = "CMASS"
        else:
            weight_scheme = False
    else:  # eBOSS
        z_min, z_max = 0, 10000
        if args.catalog_type == "galaxy":
            weight_scheme = True
        else:
            weight_scheme = False

    # Apply fraction (subsample) -- use random_fraction for values < 1
    random_fraction = args.fraction if args.fraction < 1.0 else None

    if args.catalog_type == "random":
        # Lightweight loader: only read RA/DEC/Z, subsample before loading
        data = load_catalog_lightweight(
            fits_path, columns=("RA", "DEC", "Z"),
            fraction=random_fraction,
            z_min=z_min, z_max=z_max,
        )
        weights = np.ones(len(data))
    else:
        data, weights = load_catalog(
            fits_path, weights=weight_scheme,
            random_fraction=random_fraction,
            z_min=z_min, z_max=z_max,
        )
    logger.info(f"Loaded {len(data):,} objects (z in [{z_min}, {z_max}], fraction={args.fraction})")

    # ------------------------------------------------------------------
    # Run pair finding
    # ------------------------------------------------------------------
    t0 = time.time()
    pairs_df = build_pair_catalog(
        data=data,
        weights=weights,
        catalog_type=args.catalog_type,
        dataset=args.dataset,
        region=args.region,
        r_par_max=args.rpar,
        r_perp_min=args.rperp_min,
        r_perp_max=args.rperp_max,
        n_processes=n_processes,
        chunk_size=args.chunk_size,
        output_dir=args.output_dir,
    )
    elapsed = time.time() - t0

    logger.info("=" * 60)
    logger.info(f"DONE. {len(pairs_df):,} pairs found in {elapsed:.1f}s")
    logger.info(f"Output: {output_file}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
