#!/usr/bin/env python
"""
stack_single.py -- Standalone script to stack single galaxies or random points
against the Planck kappa (convergence) map.

Mirrors the stacking logic from the legacy notebook finalMODkmap.ipynb, but
structured as a command-line tool with argparse, logging, progress bars,
and jackknife error estimation.

New in this version: Optional Gaussian kernel density estimation (KDE) for
extracting smooth radial profiles that respect the Planck map resolution limit.

Usage examples:
    # Standard stacking (2D maps)
    python stack_single.py --dataset BOSS --region North --catalog-type galaxy

    # With kernel smoothing for smooth radial profiles
    python stack_single.py --dataset BOSS --region North --catalog-type galaxy \
        --extract-radial-kde --kernel-width 3.5

    # Random catalog
    python stack_single.py --dataset BOSS --region North --catalog-type random --fraction 0.10
"""

import argparse
import logging
import os
import time

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
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
# Kernel smoothing for radial profiles
# ===========================================================================
def extract_radial_profile_kde(kappa_map, kappa_err=None,
                               kernel_width_hmpc=3.0,
                               n_radial_bins=500,
                               r_max=None):
    """Extract smooth radial profile using Gaussian kernel density estimation.

    This function addresses the fundamental resolution limit of the Planck kappa
    map (~8 arcmin FWHM ≈ 3-4 Mpc at z~0.5) by applying kernel smoothing to the
    discrete grid rather than naively increasing grid resolution.

    Parameters
    ----------
    kappa_map : 2-D ndarray (GRID_SIZE x GRID_SIZE)
        Stacked convergence map from stack_kappa().
    kappa_err : 2-D ndarray, optional
        Error map (same shape as kappa_map). If provided, errors are propagated
        through the kernel smoothing.
    kernel_width_hmpc : float
        Gaussian kernel width (sigma) in h^-1 Mpc. Should be comparable to or
        larger than the Planck resolution (~3-4 Mpc). Default: 3.0.
    n_radial_bins : int
        Number of radial evaluation points for the smooth profile. Default: 500.
    r_max : float, optional
        Maximum radius in h^-1 Mpc. If None, uses half the box size.

    Returns
    -------
    r_profile : ndarray (n_radial_bins,)
        Radial distances in h^-1 Mpc.
    kappa_profile : ndarray (n_radial_bins,)
        Smoothed radial convergence profile.
    kappa_profile_err : ndarray (n_radial_bins,) or None
        Propagated errors if kappa_err was provided.

    Notes
    -----
    The kernel width is in h^-1 Mpc, matching the physical scale of the box.
    For reference:
    - Planck FWHM = 8 arcmin
    - At z=0.5: 8 arcmin ≈ 3.5 h^-1 Mpc
    - Recommended kernel_width_hmpc >= 3.0

    Examples
    --------
    >>> kappa_map, kappa_err, _ = stack_kappa(data, weights, kmap, mask)
    >>> r, profile, err = extract_radial_profile_kde(kappa_map, kappa_err,
    ...                                               kernel_width_hmpc=4.0)
    """
    if r_max is None:
        r_max = HALF_BOX_HMPC

    # Create coordinate grids (in h^-1 Mpc from center)
    x = np.linspace(-HALF_BOX_HMPC, HALF_BOX_HMPC, GRID_SIZE)
    y = np.linspace(-HALF_BOX_HMPC, HALF_BOX_HMPC, GRID_SIZE)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)

    # Flatten arrays
    r_flat = R.ravel()
    kappa_flat = kappa_map.ravel()

    # Remove NaN/inf values
    valid = np.isfinite(kappa_flat) & np.isfinite(r_flat)
    r_flat = r_flat[valid]
    kappa_flat = kappa_flat[valid]

    if kappa_err is not None:
        err_flat = kappa_err.ravel()[valid]

    # Sort by radius for efficient binning
    sort_idx = np.argsort(r_flat)
    r_sorted = r_flat[sort_idx]
    kappa_sorted = kappa_flat[sort_idx]
    if kappa_err is not None:
        err_sorted = err_flat[sort_idx]

    # Define evaluation radii
    r_profile = np.linspace(0, r_max, n_radial_bins)

    # Gaussian kernel smoothing
    # For each evaluation radius, compute weighted average using Gaussian kernel
    kappa_profile = np.zeros(n_radial_bins)
    kappa_profile_err = np.zeros(n_radial_bins) if kappa_err is not None else None

    sigma = kernel_width_hmpc  # Kernel width in h^-1 Mpc

    for i, r_eval in enumerate(r_profile):
        # Gaussian weights: exp(-0.5 * ((r - r_eval) / sigma)^2)
        dr = r_sorted - r_eval
        weights = np.exp(-0.5 * (dr / sigma)**2)

        # Normalize weights
        weight_sum = np.sum(weights)
        if weight_sum > 0:
            kappa_profile[i] = np.sum(weights * kappa_sorted) / weight_sum

            if kappa_err is not None:
                # Error propagation: σ_weighted = sqrt(sum(w^2 * σ^2)) / sum(w)
                kappa_profile_err[i] = np.sqrt(np.sum((weights * err_sorted)**2)) / weight_sum
        else:
            kappa_profile[i] = np.nan
            if kappa_err is not None:
                kappa_profile_err[i] = np.nan

    return r_profile, kappa_profile, kappa_profile_err


def extract_radial_profile_binned(kappa_map, kappa_err=None, n_bins=50):
    """Extract traditional binned radial profile (no kernel smoothing).

    This is the standard approach: azimuthally average kappa in radial annuli.
    Provided for comparison with the kernel-smoothed version.

    Parameters
    ----------
    kappa_map : 2-D ndarray
        Stacked convergence map.
    kappa_err : 2-D ndarray, optional
        Error map.
    n_bins : int
        Number of radial bins. Default: 50.

    Returns
    -------
    r_bins : ndarray
        Bin center radii in h^-1 Mpc.
    kappa_profile : ndarray
        Mean kappa in each radial bin.
    kappa_profile_err : ndarray or None
        Errors in each bin.
    """
    # Create coordinate grids
    x = np.linspace(-HALF_BOX_HMPC, HALF_BOX_HMPC, GRID_SIZE)
    y = np.linspace(-HALF_BOX_HMPC, HALF_BOX_HMPC, GRID_SIZE)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)

    # Define radial bins
    r_edges = np.linspace(0, HALF_BOX_HMPC, n_bins + 1)
    r_bins = 0.5 * (r_edges[:-1] + r_edges[1:])

    kappa_profile = np.zeros(n_bins)
    kappa_profile_err = np.zeros(n_bins) if kappa_err is not None else None

    for i in range(n_bins):
        mask = (R >= r_edges[i]) & (R < r_edges[i+1])
        if np.sum(mask) > 0:
            kappa_profile[i] = np.nanmean(kappa_map[mask])

            if kappa_err is not None:
                # Standard error of the mean in this bin
                kappa_profile_err[i] = np.sqrt(np.nanmean(kappa_err[mask]**2))
        else:
            kappa_profile[i] = np.nan
            if kappa_err is not None:
                kappa_profile_err[i] = np.nan

    return r_bins, kappa_profile, kappa_profile_err


# ===========================================================================
# Jackknife error estimation
# ===========================================================================
def assign_jackknife_regions_healpix(catalog, rand_cat, nside_jk=10, threshold=0.9):
    """Assign HEALPix-based jackknife regions to galaxies.

    Parameters
    ----------
    catalog : structured ndarray
        Galaxy catalog with RA, DEC columns.
    rand_cat : structured ndarray
        Random catalog used to define survey footprint.
    nside_jk : int
        HEALPix nside for the jackknife tessellation.
    threshold : float
        Minimum fraction of mean random density to consider a pixel as
        inside the survey footprint.

    Returns
    -------
    gal_pix : ndarray
        Jackknife region index for each *retained* galaxy.
    mask : bool ndarray
        Boolean mask into *catalog* indicating which galaxies fall inside
        valid jackknife regions.
    """
    rand_coord = SkyCoord(ra=rand_cat['RA'] * u.deg, dec=rand_cat['DEC'] * u.deg)
    rand_pix = hp.ang2pix(nside_jk,
                           np.radians(90.0 - rand_coord.dec.value),
                           np.radians(rand_coord.ra.value))

    unique, counts = np.unique(rand_pix, return_counts=True)
    mean_density = np.mean(counts)
    good_pix = unique[counts > threshold * mean_density]

    gal_coord = SkyCoord(ra=catalog['RA'] * u.deg, dec=catalog['DEC'] * u.deg)
    gal_pix = hp.ang2pix(nside_jk,
                          np.radians(90.0 - gal_coord.dec.value),
                          np.radians(gal_coord.ra.value))

    valid = np.isin(gal_pix, good_pix)
    return gal_pix[valid], valid


def jackknife_stack_healpix(data, weights, rand_cat, kmap, mask,
                            nside_jk=10, threshold=0.9):
    """Jackknife resampling of the single-galaxy stack.

    Drops one HEALPix jackknife region at a time and re-stacks, then
    computes the jackknife error estimate.

    Parameters
    ----------
    data : structured ndarray
        Galaxy catalog.
    weights : ndarray
        Per-galaxy weights.
    rand_cat : structured ndarray
        Random catalog (for defining jackknife regions).
    kmap : ndarray
        HEALPix convergence map.
    mask : ndarray
        HEALPix lensing mask.
    nside_jk : int
        HEALPix nside for jackknife regions.
    threshold : float
        Footprint threshold (see assign_jackknife_regions_healpix).

    Returns
    -------
    kappa_mean : 2-D ndarray
        Mean of jackknife realisations.
    kappa_err : 2-D ndarray
        Jackknife error estimate per pixel.
    """
    logger.info("Assigning jackknife regions (nside=%d) ...", nside_jk)
    jk_regions, valid_mask = assign_jackknife_regions_healpix(
        data, rand_cat, nside_jk, threshold
    )
    data = data[valid_mask]
    weights = weights[valid_mask]
    unique_regions = np.unique(jk_regions)
    n_regions = len(unique_regions)
    logger.info("Found %d jackknife regions", n_regions)

    kappa_jk = []
    for reg in tqdm(unique_regions, desc="Jackknife (HEALPix)"):
        mask_jk = jk_regions != reg
        data_jk = data[mask_jk]
        weights_jk = weights[mask_jk]
        kappa, _, _ = stack_kappa(data_jk, weights_jk, kmap, mask,
                                  label=f"JK-{reg}")
        kappa_jk.append(kappa)

    kappa_jk = np.array(kappa_jk)
    kappa_mean = np.mean(kappa_jk, axis=0)
    kappa_err = np.sqrt(
        (n_regions - 1) / n_regions * np.sum((kappa_jk - kappa_mean) ** 2, axis=0)
    )

    return kappa_mean, kappa_err


# ===========================================================================
# Output helpers
# ===========================================================================
def output_kappa_path(output_dir, catalog_type, dataset, region):
    return os.path.join(output_dir, f"kappa_single_{catalog_type}_{dataset}_{region}.csv")


def output_error_path(output_dir, dataset, region):
    return os.path.join(output_dir, f"error_single_{dataset}_{region}.csv")


# ===========================================================================
# Main entry point
# ===========================================================================
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Stack single galaxies (or randoms) against the Planck kappa map."
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
    parser.add_argument(
        "--jackknife-nside", type=int, default=10,
        help="HEALPix nside for jackknife region tessellation (default: 10)."
    )
    parser.add_argument(
        "--extract-radial-kde", action="store_true",
        help="Extract smooth radial profile using Gaussian kernel density estimation."
    )
    parser.add_argument(
        "--kernel-width", type=float, default=3.0,
        help="Gaussian kernel width (sigma) in h^-1 Mpc for KDE smoothing (default: 3.0). "
             "Should be >= Planck resolution (~3-4 Mpc at z~0.5)."
    )
    parser.add_argument(
        "--n-radial-bins", type=int, default=500,
        help="Number of radial evaluation points for smooth KDE profile (default: 500)."
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

    logger.info("=== stack_single.py ===")
    logger.info("Dataset   : %s", args.dataset)
    logger.info("Region    : %s", args.region)
    logger.info("Catalog   : %s", args.catalog_type)
    logger.info("Fraction  : %.4f", args.fraction)
    logger.info("Output dir: %s", args.output_dir)
    logger.info("Data dir  : %s", args.data_dir)
    logger.info("Overwrite : %s", args.overwrite)
    logger.info("JK nside  : %d", args.jackknife_nside)
    logger.info("Radial KDE: %s", args.extract_radial_kde)
    if args.extract_radial_kde:
        logger.info("  Kernel width   : %.2f h^-1 Mpc", args.kernel_width)
        logger.info("  Radial bins    : %d", args.n_radial_bins)

    # ------------------------------------------------------------------
    # Check output paths early
    # ------------------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)

    kappa_out = output_kappa_path(args.output_dir, args.catalog_type,
                                  args.dataset, args.region)
    error_out = output_error_path(args.output_dir, args.dataset, args.region)

    if not args.overwrite:
        # For galaxy stacks we also produce the error map; for randoms we
        # only produce the kappa map.
        if os.path.exists(kappa_out):
            logger.info("Output already exists: %s  (use --overwrite to rerun)", kappa_out)
            return
        if args.catalog_type == "galaxy" and os.path.exists(error_out):
            logger.info("Error map already exists: %s  (use --overwrite to rerun)", error_out)
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
        # Lightweight loader for large random catalogs (memmap + column selection)
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
    # Stack
    # ------------------------------------------------------------------
    t_start = time.time()

    if args.catalog_type == "galaxy":
        # For galaxy catalogs: run jackknife to get the error map.
        # We need a random catalog to define jackknife regions.
        rand_path = resolve_catalog_path(
            args.data_dir, args.dataset, args.region, "random"
        )
        logger.info("Loading random catalog for jackknife regions: %s", rand_path)
        rand_data = load_catalog_lightweight(
            rand_path, columns=("RA", "DEC", "Z"),
            fraction=0.10, z_min=z_min, z_max=z_max,
        )
        logger.info("Loaded %d randoms for jackknife assignment", len(rand_data))

        kappa_map, kappa_err = jackknife_stack_healpix(
            data, weights, rand_data, kmap, planck_mask,
            nside_jk=args.jackknife_nside,
        )
    else:
        # For random catalogs: simple stack (no jackknife needed)
        kappa_map, _, _ = stack_kappa(
            data, weights, kmap, planck_mask, label="random"
        )
        kappa_err = None

    elapsed = time.time() - t_start
    logger.info("Stacking completed in %.1f s (%.1f min)", elapsed, elapsed / 60)

    # ------------------------------------------------------------------
    # Symmetrize
    # ------------------------------------------------------------------
    kappa_map = symmetrize_map(kappa_map)
    logger.info("Applied radial symmetrization to kappa map")

    # ------------------------------------------------------------------
    # Optional: Extract smooth radial profile via kernel smoothing
    # ------------------------------------------------------------------
    if args.extract_radial_kde:
        logger.info("Extracting smooth radial profile with KDE (kernel_width=%.2f Mpc/h) ...",
                    args.kernel_width)

        r_profile, kappa_profile, kappa_profile_err = extract_radial_profile_kde(
            kappa_map,
            kappa_err=kappa_err if args.catalog_type == "galaxy" else None,
            kernel_width_hmpc=args.kernel_width,
            n_radial_bins=args.n_radial_bins,
        )

        # Also extract traditional binned profile for comparison
        r_binned, kappa_binned, kappa_binned_err = extract_radial_profile_binned(
            kappa_map,
            kappa_err=kappa_err if args.catalog_type == "galaxy" else None,
            n_bins=50,
        )

        # Save radial profiles
        radial_kde_path = os.path.join(
            args.output_dir,
            f"radial_kde_{args.catalog_type}_{args.dataset}_{args.region}.csv"
        )
        radial_binned_path = os.path.join(
            args.output_dir,
            f"radial_binned_{args.catalog_type}_{args.dataset}_{args.region}.csv"
        )

        # Save KDE profile
        df_kde = pd.DataFrame({
            'r_hmpc': r_profile,
            'kappa': kappa_profile,
        })
        if kappa_profile_err is not None:
            df_kde['kappa_err'] = kappa_profile_err
        df_kde.to_csv(radial_kde_path, index=False)
        logger.info("Saved KDE radial profile -> %s", radial_kde_path)

        # Save binned profile
        df_binned = pd.DataFrame({
            'r_hmpc': r_binned,
            'kappa': kappa_binned,
        })
        if kappa_binned_err is not None:
            df_binned['kappa_err'] = kappa_binned_err
        df_binned.to_csv(radial_binned_path, index=False)
        logger.info("Saved binned radial profile -> %s", radial_binned_path)

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    pd.DataFrame(kappa_map).to_csv(kappa_out, index=True)
    logger.info("Saved kappa map -> %s", kappa_out)

    if kappa_err is not None:
        pd.DataFrame(kappa_err).to_csv(error_out, index=True)
        logger.info("Saved error map -> %s", error_out)

    logger.info("Done.")


if __name__ == "__main__":
    main()
