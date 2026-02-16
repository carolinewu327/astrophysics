#!/usr/bin/env python
"""
plot_results.py -- Standalone script that loads stacked kappa map CSVs,
computes corrected and derived maps, and generates all plots from the
adviser's feedback.

This script works entirely from CSV files -- no FITS data is needed.

Usage examples:
    python plot_results.py --separation 20
    python plot_results.py --separation 20 --dataset BOSS --regions North,South
    python plot_results.py --separation 10 --dataset eBOSS --regions NGC,SGC --format png
    python plot_results.py --separation 20 --results-dir analysis/boss/results --output-dir output/plots
"""

import argparse
import logging
import os

import matplotlib
matplotlib.use("Agg")  # non-interactive backend; safe for scripts
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from constants import BOX_SIZE_HMPC, GRID_SIZE

logger = logging.getLogger(__name__)


# ===========================================================================
# Constants
# ===========================================================================
HALF_BOX = BOX_SIZE_HMPC / 2  # 50 Mpc/h
EXTENT = [-HALF_BOX, HALF_BOX, -HALF_BOX, HALF_BOX]
PROFILE_HALF_WIDTH = 5.0  # Mpc/h -- band for profile extraction


# ===========================================================================
# CSV loading helpers
# ===========================================================================
def load_map_csv(path):
    """Load a 2-D kappa map from a CSV file written by pandas.

    Parameters
    ----------
    path : str
        Path to the CSV file (expects ``index_col=0`` format).

    Returns
    -------
    ndarray
        2-D numpy array.
    """
    return pd.read_csv(path, index_col=0).values


def try_load(path, label):
    """Attempt to load a CSV map, logging the outcome.

    Returns the map array on success or ``None`` on failure.
    """
    if not os.path.isfile(path):
        logger.warning("MISSING  %s: %s", label, path)
        return None
    arr = load_map_csv(path)
    logger.info("Loaded   %s: %s  (shape %s)", label, path, arr.shape)
    return arr


# ===========================================================================
# Map arithmetic helpers
# ===========================================================================
def average_region_maps(maps):
    """Element-wise average of a list of 2-D arrays, skipping ``None`` entries.

    All arrays must have the same shape.
    """
    valid = [m for m in maps if m is not None]
    if len(valid) == 0:
        return None
    return np.mean(valid, axis=0)


def reconcile_shapes(map_a, map_b):
    """Return copies of *map_a* and *map_b* trimmed to the same shape.

    Single-galaxy stacking typically produces 100x100 grids while
    pair stacking produces 101x101.  This helper trims the larger map
    to the smaller size, cutting symmetrically from both sides when
    possible (i.e. trim the last row/column of the larger map).
    """
    if map_a is None or map_b is None:
        return map_a, map_b
    if map_a.shape == map_b.shape:
        return map_a, map_b

    min_rows = min(map_a.shape[0], map_b.shape[0])
    min_cols = min(map_a.shape[1], map_b.shape[1])

    def trim(arr, nr, nc):
        dr = arr.shape[0] - nr
        dc = arr.shape[1] - nc
        r0 = dr // 2
        c0 = dc // 2
        return arr[r0:r0 + nr, c0:c0 + nc]

    a_out = trim(map_a, min_rows, min_cols)
    b_out = trim(map_b, min_rows, min_cols)
    logger.info(
        "Reconciled shapes: %s & %s -> %s",
        map_a.shape, map_b.shape, a_out.shape,
    )
    return a_out, b_out


def build_control_pair_map(single_map, separation_hmpc):
    """Construct a control pair map by shifting the single-galaxy map.

    The control pair map places two copies of the single-galaxy halo at
    +/-(separation/2) along the X axis and averages them, emulating what
    a pair stack would look like if there were *no* filament connecting
    the galaxies.

    For sub-pixel shifts we use ``scipy.ndimage.shift``; for integer
    shifts we use ``numpy.roll`` (faster).

    Parameters
    ----------
    single_map : ndarray
        Corrected single-galaxy kappa map.
    separation_hmpc : float
        Pair separation in h^-1 Mpc.

    Returns
    -------
    ndarray
        Control pair map (same shape as *single_map*).
    """
    grid_size = single_map.shape[1]
    cell_size = BOX_SIZE_HMPC / grid_size
    shift_pixels = separation_hmpc / 2.0 / cell_size

    if np.isclose(shift_pixels, np.round(shift_pixels)):
        # Integer shift -- use numpy.roll (wraps around, which is fine
        # because the map edges are ~zero for well-centred stacks).
        shift_int = int(np.round(shift_pixels))
        shifted_right = np.roll(single_map, +shift_int, axis=1)
        shifted_left = np.roll(single_map, -shift_int, axis=1)
        logger.info(
            "Control pair map: integer shift of %d pixels (%.2f Mpc/h)",
            shift_int, shift_int * cell_size,
        )
    else:
        # Sub-pixel shift -- use scipy for interpolation
        from scipy.ndimage import shift as ndimage_shift
        shifted_right = ndimage_shift(single_map, [0, +shift_pixels],
                                      order=3, mode='wrap')
        shifted_left = ndimage_shift(single_map, [0, -shift_pixels],
                                     order=3, mode='wrap')
        logger.info(
            "Control pair map: sub-pixel shift of %.2f pixels (%.2f Mpc/h)",
            shift_pixels, shift_pixels * cell_size,
        )

    return 0.5 * (shifted_right + shifted_left)


# ===========================================================================
# Profile extraction
# ===========================================================================
def extract_profile(kappa_map, box_size=BOX_SIZE_HMPC, half_width=PROFILE_HALF_WIDTH):
    """Sum pixel rows within +/- *half_width* Mpc/h of the Y=0 axis.

    Parameters
    ----------
    kappa_map : ndarray  (ny, nx)
        2-D kappa map.
    box_size : float
        Physical box size in h^-1 Mpc.
    half_width : float
        Half-width of the summation band along Y in h^-1 Mpc.

    Returns
    -------
    x_coords : ndarray
        1-D array of X coordinates (Mpc/h) at pixel centres.
    profile : ndarray
        Summed kappa values within the band as a function of X.
    """
    ny, nx = kappa_map.shape
    cell_size = box_size / ny
    half_box = box_size / 2.0

    # Y coordinate of each row's centre
    y_centres = np.linspace(-half_box + cell_size / 2,
                            half_box - cell_size / 2, ny)

    # Select rows within the band
    row_mask = np.abs(y_centres) <= half_width
    if not np.any(row_mask):
        logger.warning("No rows within +/-%.1f Mpc/h -- returning zeros", half_width)
        x_coords = np.linspace(-half_box + cell_size / 2,
                               half_box - cell_size / 2, nx)
        return x_coords, np.zeros(nx)

    band = kappa_map[row_mask, :]
    profile = np.sum(band, axis=0)

    x_coords = np.linspace(-half_box + cell_size / 2,
                            half_box - cell_size / 2, nx)
    return x_coords, profile


# ===========================================================================
# Plotting routines
# ===========================================================================
def plot_single_map(kappa_map, title, cmap, output_path, fmt,
                    center_on_zero=False, extent=None):
    """Save a single 2-D kappa map as an image file.

    Parameters
    ----------
    kappa_map : ndarray
        2-D map to plot.
    title : str
        Figure title.
    cmap : str
        Matplotlib colormap name.
    output_path : str
        Output file path (without extension -- *fmt* is appended).
    fmt : str
        File format ('pdf' or 'png').
    center_on_zero : bool
        If True, use ``TwoSlopeNorm`` to centre the colorbar on zero.
    extent : list or None
        ``[xmin, xmax, ymin, ymax]`` for ``imshow``.
    """
    if extent is None:
        extent = EXTENT
    fig, ax = plt.subplots(figsize=(7, 6))

    if center_on_zero:
        vabs = max(abs(np.nanmin(kappa_map)), abs(np.nanmax(kappa_map)))
        if vabs == 0:
            vabs = 1e-10
        norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
        im = ax.imshow(kappa_map, origin='lower', extent=extent,
                       cmap=cmap, norm=norm)
    else:
        im = ax.imshow(kappa_map, origin='lower', extent=extent,
                       cmap=cmap)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r'$\kappa$')
    ax.set_xlabel(r'X ($h^{-1}\,$Mpc)')
    ax.set_ylabel(r'Y ($h^{-1}\,$Mpc)')
    ax.set_title(title)
    ax.axhline(0, color='grey', ls='--', lw=0.5, alpha=0.4)
    ax.axvline(0, color='grey', ls='--', lw=0.5, alpha=0.4)

    plt.tight_layout()
    full_path = f"{output_path}.{fmt}"
    fig.savefig(full_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info("Saved plot: %s", full_path)


def plot_profiles(profiles, title, output_path, fmt):
    """Plot one or more 1-D profiles on a single figure.

    Parameters
    ----------
    profiles : list of (x, y, label) tuples
        Each tuple gives X coordinates, Y values, and a legend label.
    title : str
        Figure title.
    output_path : str
        Output path (no extension).
    fmt : str
        File format.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    for x, y, label in profiles:
        ax.plot(x, y, label=label, lw=1.4)
    ax.set_xlabel(r'X ($h^{-1}\,$Mpc)')
    ax.set_ylabel(r'Summed $\kappa$')
    ax.set_title(title)
    ax.legend()
    ax.axhline(0, color='grey', ls='--', lw=0.5, alpha=0.4)
    ax.axvline(0, color='grey', ls='--', lw=0.5, alpha=0.4)
    plt.tight_layout()
    full_path = f"{output_path}.{fmt}"
    fig.savefig(full_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info("Saved plot: %s", full_path)


# ===========================================================================
# Derived-map CSV saving
# ===========================================================================
def save_derived_map(arr, name, results_dir):
    """Write a derived map to CSV (same format as the stacking scripts)."""
    path = os.path.join(results_dir, f"{name}.csv")
    pd.DataFrame(arr).to_csv(path, index=True)
    logger.info("Saved derived map: %s", path)


# ===========================================================================
# Argparse
# ===========================================================================
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Load stacked kappa map CSVs, compute corrected/derived maps, "
            "and generate analysis plots."
        ),
    )
    parser.add_argument(
        "--dataset", default="BOSS", choices=["BOSS", "eBOSS"],
        help="Galaxy survey dataset (default: BOSS).",
    )
    parser.add_argument(
        "--separation", type=float, required=True,
        help="Pair separation in Mpc/h (e.g. 20).",
    )
    parser.add_argument(
        "--results-dir", default="analysis/boss/results",
        help="Directory containing kappa CSV files (default: analysis/boss/results).",
    )
    parser.add_argument(
        "--output-dir", default="output/plots",
        help="Directory for output plot files (default: output/plots).",
    )
    parser.add_argument(
        "--format", default="pdf", choices=["pdf", "png"], dest="fmt",
        help="Plot file format (default: pdf).",
    )
    parser.add_argument(
        "--regions", default="North,South",
        help="Comma-separated list of regions to combine (default: North,South).",
    )
    return parser.parse_args(argv)


# ===========================================================================
# Main
# ===========================================================================
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

    sep = args.separation
    sep_label = f"{sep:g}"  # e.g. "20" or "10.5"
    regions = [r.strip() for r in args.regions.split(",")]

    logger.info("=" * 60)
    logger.info("plot_results.py -- Filament Analysis Plots")
    logger.info("=" * 60)
    logger.info("Dataset     : %s", args.dataset)
    logger.info("Separation  : %s Mpc/h", sep_label)
    logger.info("Regions     : %s", regions)
    logger.info("Results dir : %s", args.results_dir)
    logger.info("Output dir  : %s", args.output_dir)
    logger.info("Plot format : %s", args.fmt)

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.results_dir, exist_ok=True)

    # ==================================================================
    # 1. Load input CSV maps
    # ==================================================================
    logger.info("-" * 40)
    logger.info("Loading input maps ...")

    # --- Map (1): single galaxy maps per region, then average ---
    single_galaxy_per_region = []
    for reg in regions:
        path = os.path.join(
            args.results_dir,
            f"kappa_single_galaxy_{args.dataset}_{reg}.csv",
        )
        m = try_load(path, f"Map(1) single galaxy ({reg})")
        if m is not None:
            single_galaxy_per_region.append(m)
    map1 = average_region_maps(single_galaxy_per_region)

    # --- Map (2): galaxy pair maps per region, then average ---
    galaxy_pairs_per_region = []
    for reg in regions:
        path = os.path.join(
            args.results_dir,
            f"kappa_pairs_galaxy_{sep_label}_{args.dataset}_{reg}.csv",
        )
        m = try_load(path, f"Map(2) galaxy pairs ({reg})")
        if m is not None:
            galaxy_pairs_per_region.append(m)
    map2 = average_region_maps(galaxy_pairs_per_region)

    # --- Map (3): single random maps per region, then average ---
    single_random_per_region = []
    for reg in regions:
        path = os.path.join(
            args.results_dir,
            f"kappa_single_random_{args.dataset}_{reg}.csv",
        )
        m = try_load(path, f"Map(3) single random ({reg})")
        if m is not None:
            single_random_per_region.append(m)
    map3 = average_region_maps(single_random_per_region)

    # --- Map (4): random pair maps per region, then average ---
    random_pairs_per_region = []
    for reg in regions:
        path = os.path.join(
            args.results_dir,
            f"kappa_pairs_random_{sep_label}_{args.dataset}_{reg}.csv",
        )
        m = try_load(path, f"Map(4) random pairs ({reg})")
        if m is not None:
            random_pairs_per_region.append(m)
    map4 = average_region_maps(random_pairs_per_region)

    # ==================================================================
    # 2. Compute derived maps
    # ==================================================================
    logger.info("-" * 40)
    logger.info("Computing derived maps ...")

    # --- Map (5): corrected single = (1) - (3) ---
    map5 = None
    if map1 is not None and map3 is not None:
        map1_r, map3_r = reconcile_shapes(map1, map3)
        map5 = map1_r - map3_r
        logger.info(
            "Map(5) corrected single: range [%.6e, %.6e]",
            map5.min(), map5.max(),
        )
        save_derived_map(map5, f"kappa_corrected_single_{args.dataset}", args.results_dir)
    elif map1 is not None:
        logger.warning("Map(3) missing -- using Map(1) directly as corrected single")
        map5 = map1.copy()
    else:
        logger.warning("Cannot compute Map(5): Map(1) is missing")

    # --- Map (6): corrected pairs = (2) - (4) ---
    map6 = None
    if map2 is not None and map4 is not None:
        map2_r, map4_r = reconcile_shapes(map2, map4)
        map6 = map2_r - map4_r
        logger.info(
            "Map(6) corrected pairs: range [%.6e, %.6e]",
            map6.min(), map6.max(),
        )
        save_derived_map(
            map6,
            f"kappa_corrected_pairs_{sep_label}_{args.dataset}",
            args.results_dir,
        )
    elif map2 is not None:
        logger.warning("Map(4) missing -- using Map(2) directly as corrected pairs")
        map6 = map2.copy()
    else:
        logger.warning("Cannot compute Map(6): Map(2) is missing")

    # --- Map (7): control pair map from (5) ---
    map7 = None
    if map5 is not None:
        map7 = build_control_pair_map(map5, sep)
        logger.info(
            "Map(7) control pair: range [%.6e, %.6e]",
            map7.min(), map7.max(),
        )
        save_derived_map(
            map7,
            f"kappa_control_pair_{sep_label}_{args.dataset}",
            args.results_dir,
        )
    else:
        logger.warning("Cannot compute Map(7): Map(5) is missing")

    # --- Map (8): filament = (6) - (7) ---
    map8 = None
    if map6 is not None and map7 is not None:
        map6_r, map7_r = reconcile_shapes(map6, map7)
        map8 = map6_r - map7_r
        logger.info(
            "Map(8) filament: range [%.6e, %.6e]",
            map8.min(), map8.max(),
        )
        save_derived_map(
            map8,
            f"kappa_filament_{sep_label}_{args.dataset}",
            args.results_dir,
        )
    else:
        logger.warning("Cannot compute Map(8): Map(6) or Map(7) is missing")

    # ==================================================================
    # 3. Generate 8 individual map plots
    # ==================================================================
    logger.info("-" * 40)
    logger.info("Generating map plots ...")

    map_specs = [
        (map1, f"(1) Single Galaxies ({args.dataset})",
         'viridis', False, "map_1_single_galaxy"),
        (map2, f"(2) Galaxy Pairs ({sep_label} Mpc/h)",
         'viridis', False, "map_2_galaxy_pairs"),
        (map3, f"(3) Single Randoms ({args.dataset})",
         'viridis', False, "map_3_single_random"),
        (map4, f"(4) Random Pairs ({sep_label} Mpc/h)",
         'viridis', False, "map_4_random_pairs"),
        (map5, f"(5) Corrected Single = (1)-(3)",
         'coolwarm', True, "map_5_corrected_single"),
        (map6, f"(6) Corrected Pairs = (2)-(4)",
         'coolwarm', True, "map_6_corrected_pairs"),
        (map7, f"(7) Control Pair Map ({sep_label} Mpc/h)",
         'coolwarm', True, "map_7_control_pair"),
        (map8, f"(8) Filament = (6)-(7) [{sep_label} Mpc/h]",
         'coolwarm', True, "map_8_filament"),
    ]

    for m, title, cmap, centered, fname in map_specs:
        if m is None:
            logger.warning("Skipping plot %s -- map not available", fname)
            continue
        # Compute the correct extent based on actual map shape
        ny, nx = m.shape
        cell_x = BOX_SIZE_HMPC / nx
        cell_y = BOX_SIZE_HMPC / ny
        ext = [
            -HALF_BOX + cell_x / 2 - cell_x / 2,  # simplifies to -HALF_BOX
            HALF_BOX,
            -HALF_BOX,
            HALF_BOX,
        ]
        plot_single_map(
            m, title, cmap,
            os.path.join(args.output_dir, fname),
            args.fmt,
            center_on_zero=centered,
            extent=ext,
        )

    # ==================================================================
    # 4. Generate 3 profile plots
    # ==================================================================
    logger.info("-" * 40)
    logger.info("Generating profile plots ...")

    # Helper: extract profile or return None
    def safe_profile(kappa_map, label):
        if kappa_map is None:
            logger.warning("Profile skipped for %s -- map not available", label)
            return None
        return extract_profile(kappa_map)

    # --- Profile (a): maps (1), (3), (5) ---
    p1 = safe_profile(map1, "Map(1)")
    p3 = safe_profile(map3, "Map(3)")
    p5 = safe_profile(map5, "Map(5)")
    profs_a = []
    if p1 is not None:
        profs_a.append((*p1, "(1) Single Galaxies"))
    if p3 is not None:
        profs_a.append((*p3, "(3) Single Randoms"))
    if p5 is not None:
        profs_a.append((*p5, "(5) Corrected Single"))
    if profs_a:
        plot_profiles(
            profs_a,
            f"Single-Galaxy Profiles ({args.dataset})",
            os.path.join(args.output_dir, "profile_a_single"),
            args.fmt,
        )
    else:
        logger.warning("Skipping profile plot (a) -- no data available")

    # --- Profile (b): maps (2), (4), (6) ---
    p2 = safe_profile(map2, "Map(2)")
    p4 = safe_profile(map4, "Map(4)")
    p6 = safe_profile(map6, "Map(6)")
    profs_b = []
    if p2 is not None:
        profs_b.append((*p2, f"(2) Galaxy Pairs ({sep_label} Mpc/h)"))
    if p4 is not None:
        profs_b.append((*p4, f"(4) Random Pairs ({sep_label} Mpc/h)"))
    if p6 is not None:
        profs_b.append((*p6, f"(6) Corrected Pairs"))
    if profs_b:
        plot_profiles(
            profs_b,
            f"Pair-Stacked Profiles ({sep_label} Mpc/h, {args.dataset})",
            os.path.join(args.output_dir, "profile_b_pairs"),
            args.fmt,
        )
    else:
        logger.warning("Skipping profile plot (b) -- no data available")

    # --- Profile (c): maps (6), (7), (8) ---
    p7 = safe_profile(map7, "Map(7)")
    p8 = safe_profile(map8, "Map(8)")
    profs_c = []
    if p6 is not None:
        profs_c.append((*p6, "(6) Corrected Pairs"))
    if p7 is not None:
        profs_c.append((*p7, "(7) Control Pair"))
    if p8 is not None:
        profs_c.append((*p8, "(8) Filament"))
    if profs_c:
        plot_profiles(
            profs_c,
            f"Filament Extraction ({sep_label} Mpc/h, {args.dataset})",
            os.path.join(args.output_dir, "profile_c_filament"),
            args.fmt,
        )
    else:
        logger.warning("Skipping profile plot (c) -- no data available")

    # ==================================================================
    # Summary
    # ==================================================================
    logger.info("=" * 60)
    logger.info("All done.")
    available = sum(1 for m in [map1, map2, map3, map4, map5, map6, map7, map8]
                    if m is not None)
    logger.info("Maps available: %d / 8", available)
    logger.info("Plots saved to: %s", os.path.abspath(args.output_dir))
    logger.info("Derived CSVs in: %s", os.path.abspath(args.results_dir))
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
