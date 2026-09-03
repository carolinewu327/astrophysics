#!/usr/bin/env python
"""Combine galaxy and random jackknife accumulators into corrected maps and profiles.

Consumes the per-region accumulators written by ``stack_single_jk.py`` and
forms, for each of the K jackknife regions, the leave-one-out estimate of the
*corrected* map -- galaxies minus randoms with the same patch of sky deleted
from both.  Two things follow that were not available before:

* **Errors on the subtracted map.**  The old pipeline jackknifed galaxies only
  and carried that error onto ``galaxy - random``, ignoring the random stack's
  own sampling noise and its correlation with the galaxy stack over the same
  footprint.
* **Profile-level errors with a covariance.**  Radial bins are re-derived
  inside every leave-one-out realization, so the covariance between radii is
  measured rather than assumed.  Collapsing per-pixel errors into a profile
  by assuming pixels decorrelate is not defensible here: the 8 arcmin Planck
  beam is ~3.3 h^-1 Mpc at z = 0.55, more than three times the 1 h^-1 Mpc
  pixel, so neighbouring pixels and neighbouring radial bins are strongly
  correlated.

Radii are measured from the true stack center (physical r = 0, which falls
*between* pixels 49 and 50 on the 100-pixel grid), not from pixel index 50 as
``lib/geometry.radial_profile`` does.  That legacy half-pixel offset is 0.5
h^-1 Mpc -- half a bin width at the default binning.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/combine_jackknife.py \\
        --dataset BOSS --regions North,South

    # log bins, for comparison against published wide-bin measurements
    PYTHONPATH=lib python analysis/boss/scripts/combine_jackknife.py \\
        --dataset BOSS --regions North,South \\
        --binning log --r-min 0.5 --r-max 50 --n-bins 12
"""

from __future__ import annotations

import argparse
import json
import logging
import os

import numpy as np
import pandas as pd

from catalog import setup_logging
from constants import BOX_SIZE_HMPC, GRID_SIZE
from geometry import stack_axis, symmetrize_map
from jackknife import jackknife_covariance, jackknife_error

logger = logging.getLogger(__name__)

CELL_SIZE_HMPC = BOX_SIZE_HMPC / GRID_SIZE
HALF_BOX_HMPC = BOX_SIZE_HMPC / 2.0

# Bands used for the BOSS-vs-simulation significance quotes.
DEFAULT_BANDS = ((10.0, 20.0), (15.0, 25.0), (20.0, 40.0))


def grid_radius_hmpc(grid_size: int = GRID_SIZE) -> np.ndarray:
    """Physical radius of every grid pixel, measured from the stack center.

    Working in physical offsets rather than pixel indices sidesteps the
    question of which index "the center" is.  The offsets come from
    geometry.stack_axis, so an even grid puts the object between the two
    central pixels (+/-0.5, 1.5, ...) and an odd grid puts it on the middle
    pixel -- the same rule the pair stacks and the simulation use.
    """
    offsets, _ = stack_axis(BOX_SIZE_HMPC, grid_size)
    off_x, off_y = np.meshgrid(offsets, offsets)
    return np.sqrt(off_x**2 + off_y**2).ravel()


def load_accumulator(path: str) -> dict[str, object]:
    with np.load(path, allow_pickle=False) as data:
        # Single stacks count objects, pair stacks count pairs. Expose both
        # names so downstream code can treat the two interchangeably.
        count_key = "n_objects" if "n_objects" in data.files else "n_pairs"
        counts = data[count_key]
        acc = {
            "sum_wk": data["sum_wk"],
            "sum_w": data["sum_w"],
            "n_used": data["n_used"],
            "n_objects": counts,
            "n_pairs": counts,
            "seed_pix": data["seed_pix"],
            "jk_digest": str(data["jk_digest"]),
            "jk_nside": int(data["jk_nside"]),
            "grid_size": int(data["grid_size"]),
            "args": json.loads(str(data["args_json"])),
        }
    logger.info(
        "Loaded %s: %d regions, %d objects, digest=%s",
        path, acc["sum_wk"].shape[0], int(acc["n_objects"].sum()), acc["jk_digest"],
    )
    return acc


def leave_one_out_maps(sum_wk: np.ndarray, sum_w: np.ndarray) -> np.ndarray:
    """Weighted-mean map with each region deleted in turn -- shape (K, n_pix).

    Deleting a region is a subtraction from the totals, so the whole set costs
    one pass over an already-accumulated array rather than K re-stacks.
    """
    total_wk = sum_wk.sum(axis=0)
    total_w = sum_w.sum(axis=0)
    loo_wk = total_wk[None, :] - sum_wk
    loo_w = total_w[None, :] - sum_w
    out = np.zeros_like(loo_wk)
    good = loo_w > 0
    out[good] = loo_wk[good] / loo_w[good]
    return out


def total_map(sum_wk: np.ndarray, sum_w: np.ndarray) -> np.ndarray:
    total_wk = sum_wk.sum(axis=0)
    total_w = sum_w.sum(axis=0)
    out = np.zeros_like(total_wk)
    good = total_w > 0
    out[good] = total_wk[good] / total_w[good]
    return out


def bin_edges_from_args(args: argparse.Namespace) -> np.ndarray:
    if args.bin_edges:
        edges = np.array([float(x) for x in args.bin_edges.split(",")], dtype=float)
        if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0):
            raise ValueError("--bin-edges must be an increasing comma-separated list.")
        return edges
    if args.binning == "log":
        if args.r_min <= 0:
            raise ValueError("--r-min must be positive for log binning.")
        return np.geomspace(args.r_min, args.r_max, args.n_bins + 1)
    return np.arange(0.0, args.r_max + args.bin_width, args.bin_width)


def radial_profiles(
    maps: np.ndarray, edges: np.ndarray, grid_size: int = GRID_SIZE
) -> tuple[np.ndarray, np.ndarray]:
    """Bin a stack of flattened maps into radial bins.

    ``maps`` is (n_maps, n_pix).  Returns ``(profiles, n_pix_per_bin)`` where
    profiles is (n_maps, n_bins), the unweighted mean over pixels in each
    annulus.
    """
    radius = grid_radius_hmpc(grid_size)
    # np.digitize returns 0 for r < edges[0] and len(edges) for r >= edges[-1];
    # shifting by one puts valid bins at 0..n_bins-1.
    which = np.digitize(radius, edges) - 1
    n_bins = len(edges) - 1
    inside = (which >= 0) & (which < n_bins)

    counts = np.bincount(which[inside], minlength=n_bins).astype(float)
    if np.any(counts == 0):
        empty = np.flatnonzero(counts == 0)
        raise ValueError(
            f"Radial bin(s) {empty.tolist()} contain no pixels; widen the bins "
            f"(edges {edges[empty]}-{edges[empty + 1]} h^-1 Mpc are narrower "
            "than the 1 h^-1 Mpc grid)."
        )

    profiles = np.empty((maps.shape[0], n_bins), dtype=float)
    for i, flat in enumerate(maps):
        profiles[i] = np.bincount(which[inside], weights=flat[inside], minlength=n_bins) / counts
    return profiles, counts


def band_average(
    profile: np.ndarray,
    cov: np.ndarray,
    edges: np.ndarray,
    counts: np.ndarray,
    lo: float,
    hi: float,
) -> tuple[float, float, int]:
    """Pixel-count-weighted average of a profile over [lo, hi), with its error.

    The error uses the full covariance -- ``sqrt(w^T C w)`` -- because the
    bins inside a band are correlated by the beam, and treating them as
    independent is exactly the approximation this script exists to remove.
    """
    centers = 0.5 * (edges[:-1] + edges[1:])
    sel = (centers >= lo) & (centers < hi)
    if not sel.any():
        raise ValueError(f"No radial bins fall in band [{lo}, {hi}).")
    w = np.zeros_like(centers)
    w[sel] = counts[sel] / counts[sel].sum()
    value = float(w @ profile)
    sigma = float(np.sqrt(max(w @ cov @ w, 0.0)))
    return value, sigma, int(sel.sum())


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine galaxy/random jackknife accumulators into corrected "
        "maps, per-pixel errors, and profile-level errors with a covariance.",
    )
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument(
        "--tag", default="_scw",
        help="Filename tag of the accumulators to combine (default: _scw).",
    )
    parser.add_argument(
        "--random-tag", default=None,
        help="Tag of the random accumulator when it differs from --tag, e.g. "
             "'_scw_frac10' for a 10%% random stack. The tag is carried into "
             "the output names so a 10%% and a 100%% random subtraction never "
             "overwrite each other.",
    )
    parser.add_argument("--galaxy-acc", default=None, help="Explicit galaxy accumulator .npz.")
    parser.add_argument("--random-acc", default=None, help="Explicit random accumulator .npz.")
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument(
        "--binning", default="linear", choices=["linear", "log"],
        help="Radial binning scheme (default: linear).",
    )
    parser.add_argument("--bin-width", type=float, default=1.0,
                        help="Linear bin width in h^-1 Mpc (default: 1.0).")
    parser.add_argument("--r-min", type=float, default=0.5,
                        help="Inner radius for log binning (default: 0.5).")
    parser.add_argument("--r-max", type=float, default=50.0,
                        help="Outer radius in h^-1 Mpc (default: 50.0).")
    parser.add_argument("--n-bins", type=int, default=12,
                        help="Number of log bins (default: 12).")
    parser.add_argument("--bin-edges", default=None,
                        help="Explicit comma-separated bin edges, overriding --binning. "
                             "Use this to match a published binning.")
    parser.add_argument("--profile-label", default=None,
                        help="Suffix for profile outputs (default: the binning scheme).")
    parser.add_argument("--overwrite", action="store_true")

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
    gal_path = args.galaxy_acc or os.path.join(
        jk_dir, f"acc_single_galaxy{args.tag}_{args.dataset}_{region_label}.npz"
    )
    rand_path = args.random_acc or os.path.join(
        jk_dir, f"acc_single_random{args.random_tag}_{args.dataset}_{region_label}.npz"
    )
    # Outputs inherit the random tag when it differs, so a 10% and a 100%
    # random subtraction produce distinct files instead of clobbering.
    out_tag = args.tag if args.random_tag == args.tag else f"{args.tag}_rand{args.random_tag}"

    gal = load_accumulator(gal_path)
    rand = load_accumulator(rand_path)

    if gal["jk_digest"] != rand["jk_digest"]:
        raise ValueError(
            f"Tessellation mismatch: galaxy digest {gal['jk_digest']} vs random "
            f"digest {rand['jk_digest']}. The two stacks used different jackknife "
            "regions and cannot be combined leave-one-out."
        )
    if gal["grid_size"] != rand["grid_size"]:
        raise ValueError("Galaxy and random accumulators use different grid sizes.")
    grid_size = gal["grid_size"]
    n_regions = gal["sum_wk"].shape[0]
    logger.info("Combining %d jackknife regions on a %dx%d grid.",
                n_regions, grid_size, grid_size)

    # ------------------------------------------------------------------
    # Maps: total and leave-one-out
    # ------------------------------------------------------------------
    gal_total = total_map(gal["sum_wk"], gal["sum_w"])
    rand_total = total_map(rand["sum_wk"], rand["sum_w"])
    corrected_total = gal_total - rand_total

    gal_loo = leave_one_out_maps(gal["sum_wk"], gal["sum_w"])
    rand_loo = leave_one_out_maps(rand["sum_wk"], rand["sum_w"])
    corrected_loo = gal_loo - rand_loo

    _, gal_err = jackknife_error(gal_loo)
    _, corrected_err = jackknife_error(corrected_loo)

    logger.info(
        "Per-pixel jackknife error, median over the inner 20 h^-1 Mpc: "
        "galaxies-only %.3e, corrected %.3e (ratio %.2f)",
        np.median(gal_err[grid_radius_hmpc(grid_size) < 20.0]),
        np.median(corrected_err[grid_radius_hmpc(grid_size) < 20.0]),
        np.median(corrected_err[grid_radius_hmpc(grid_size) < 20.0])
        / max(np.median(gal_err[grid_radius_hmpc(grid_size) < 20.0]), 1e-30),
    )

    square = (grid_size, grid_size)
    outputs = {
        f"kappa_single_galaxy{out_tag}_{args.dataset}_{region_label}_joint.csv":
            symmetrize_map(gal_total.reshape(square)),
        f"kappa_single_random{out_tag}_{args.dataset}_{region_label}_joint.csv":
            symmetrize_map(rand_total.reshape(square)),
        f"kappa_corrected_single{out_tag}_{args.dataset}_{region_label}_joint.csv":
            symmetrize_map(corrected_total.reshape(square)),
        f"kappa_corrected_single{out_tag}_{args.dataset}_{region_label}_joint_raw.csv":
            corrected_total.reshape(square),
        f"error_single{out_tag}_{args.dataset}_{region_label}_joint.csv":
            corrected_err.reshape(square),
    }
    for name, arr in outputs.items():
        path = os.path.join(args.results_dir, name)
        if os.path.exists(path) and not args.overwrite:
            logger.warning("Exists, not overwriting: %s", path)
            continue
        pd.DataFrame(arr).to_csv(path, index=True)
        logger.info("Saved %s", path)

    # ------------------------------------------------------------------
    # Profile-level jackknife
    # ------------------------------------------------------------------
    edges = bin_edges_from_args(args)
    centers = 0.5 * (edges[:-1] + edges[1:])

    all_maps = np.vstack([corrected_total[None, :], corrected_loo])
    profiles, counts = radial_profiles(all_maps, edges, grid_size)
    profile_total, profile_loo = profiles[0], profiles[1:]

    profile_jk_mean, profile_sigma = jackknife_error(profile_loo)
    cov = jackknife_covariance(profile_loo)
    corr = cov / np.outer(profile_sigma, profile_sigma)

    label = args.profile_label or args.binning
    stem = f"{out_tag}_{args.dataset}_{region_label}_joint_{label}"

    profile_df = pd.DataFrame(
        {
            "r_lo_hmpc": edges[:-1],
            "r_hi_hmpc": edges[1:],
            "r_center_hmpc": centers,
            "n_pixels": counts,
            "kappa": profile_total,
            "kappa_jk_mean": profile_jk_mean,
            "kappa_err": profile_sigma,
        }
    )
    profile_path = os.path.join(args.results_dir, f"profile_corrected{stem}.csv")
    profile_df.to_csv(profile_path, index=False)
    logger.info("Saved radial profile -> %s", profile_path)

    cov_path = os.path.join(args.results_dir, f"profile_cov_corrected{stem}.csv")
    pd.DataFrame(cov, index=centers, columns=centers).to_csv(cov_path)
    logger.info("Saved profile covariance -> %s", cov_path)

    # ------------------------------------------------------------------
    # Band averages, with and without the independent-bin approximation
    # ------------------------------------------------------------------
    band_rows = []
    for lo, hi in DEFAULT_BANDS:
        if hi > edges[-1]:
            logger.warning("Band %.0f-%.0f exceeds r_max=%.0f; skipping.", lo, hi, edges[-1])
            continue
        value, sigma_full, n_used_bins = band_average(
            profile_total, cov, edges, counts, lo, hi
        )
        w = np.zeros_like(centers)
        sel = (centers >= lo) & (centers < hi)
        w[sel] = counts[sel] / counts[sel].sum()
        sigma_diag = float(np.sqrt(np.sum((w * profile_sigma) ** 2)))
        band_rows.append(
            {
                "band_lo_hmpc": lo,
                "band_hi_hmpc": hi,
                "n_bins": n_used_bins,
                "kappa": value,
                "kappa_err": sigma_full,
                "kappa_err_diag_only": sigma_diag,
                "significance": value / sigma_full if sigma_full > 0 else np.nan,
                "significance_diag_only": value / sigma_diag if sigma_diag > 0 else np.nan,
            }
        )
    if band_rows:
        bands = pd.DataFrame(band_rows)
        band_path = os.path.join(args.results_dir, f"profile_bands_corrected{stem}.csv")
        bands.to_csv(band_path, index=False)
        logger.info("Saved band averages -> %s", band_path)
        logger.info(
            "Band significances (full covariance vs diagonal-only):\n%s",
            bands[["band_lo_hmpc", "band_hi_hmpc", "kappa",
                   "kappa_err", "kappa_err_diag_only",
                   "significance", "significance_diag_only"]].to_string(index=False),
        )

    off_diagonal = corr[np.triu_indices_from(corr, k=1)]
    logger.info(
        "Bin-to-bin correlation: median |rho| = %.2f, max = %.2f. Treating radial "
        "bins as independent is only safe if these are near zero.",
        float(np.median(np.abs(off_diagonal))), float(np.max(off_diagonal)),
    )


if __name__ == "__main__":
    main()
