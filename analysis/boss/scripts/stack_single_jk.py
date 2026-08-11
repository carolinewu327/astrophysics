#!/usr/bin/env python
"""Joint-footprint single-object stack with per-region jackknife accumulators.

Replaces the per-region ``stack_single.py`` workflow for anything that needs
error bars.  Two things change:

1. **North and South are stacked as one region.**  There is no per-region
   profile to average afterwards, so the weighting of the two footprints is
   handled automatically by the stack itself (each object enters with its own
   weight) rather than by a 50/50 average of two maps.
2. **The jackknife costs one pass, not K passes.**  Instead of re-stacking the
   whole catalog once per deleted region, this script accumulates ``sum(w*k)``
   and ``sum(w)`` *per region* in a single pass.  Every leave-one-out estimate
   is then a subtraction of one region's accumulator from the total.  This is
   the pattern already used for simulation pairs in
   ``analysis/sim/jackknife_pair_stack.py``.

Galaxies and randoms are stacked by separate invocations that share one
tessellation (``lib/jackknife.py``), so ``combine_jackknife.py`` can form the
leave-one-out *corrected* map -- deleting the same patch of sky from both --
rather than propagating a galaxies-only error onto a subtracted map.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/stack_single_jk.py \\
        --dataset BOSS --regions North,South --catalog-type galaxy

    PYTHONPATH=lib python analysis/boss/scripts/stack_single_jk.py \\
        --dataset BOSS --regions North,South --catalog-type random --fraction 1.0 \\
        --n-processes 28 --chunk-size 20000
"""

from __future__ import annotations

import argparse
import json
import logging
import multiprocessing as mp
import os
import time

import healpy as hp
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18 as cosmo

from catalog import (
    load_catalog,
    load_catalog_lightweight,
    load_kappa_map,
    resolve_catalog_path,
    resolve_planck_paths,
    setup_logging,
    sigma_crit_weights,
)
from constants import BOX_SIZE_HMPC, FWHM_ARCMIN, GRID_SIZE, NSIDE
from geometry import fast_icrs_to_galactic, symmetrize_map
from jackknife import JackknifeRegions, build_jackknife_regions

logger = logging.getLogger(__name__)

CELL_SIZE_HMPC = BOX_SIZE_HMPC / GRID_SIZE
HALF_BOX_HMPC = BOX_SIZE_HMPC / 2.0
DEG_PER_RAD = 180.0 / np.pi

# Physical offsets of the stacking grid, identical to stack_single.py so the
# two pipelines produce directly comparable maps.
OFFSETS = np.linspace(
    -HALF_BOX_HMPC + CELL_SIZE_HMPC / 2,
    HALF_BOX_HMPC - CELL_SIZE_HMPC / 2,
    GRID_SIZE,
)
OFF_X, OFF_Y = (arr.ravel() for arr in np.meshgrid(OFFSETS, OFFSETS))

# Objects per vectorized ang2pix call.  Each object costs GRID_SIZE**2 sky
# lookups, so a batch of 64 holds ~640k pixels -- large enough to amortize
# Python overhead, small enough that 28 concurrent workers stay under ~1 GB
# of transient batch memory between them.
STACK_BATCH = 64

# Worker globals, populated in the parent before Pool() so fork() shares the
# kappa map copy-on-write instead of pickling 400 MB per worker.
KMAP: np.ndarray | None = None
MASK: np.ndarray | None = None
L_ARR: np.ndarray | None = None
B_ARR: np.ndarray | None = None
D_ARR: np.ndarray | None = None
W_ARR: np.ndarray | None = None
REGION_ARR: np.ndarray | None = None


def _set_runtime_globals(
    *,
    kmap: np.ndarray,
    mask: np.ndarray,
    l_arr: np.ndarray,
    b_arr: np.ndarray,
    d_arr: np.ndarray,
    w_arr: np.ndarray,
    region_arr: np.ndarray,
) -> None:
    global KMAP, MASK, L_ARR, B_ARR, D_ARR, W_ARR, REGION_ARR
    KMAP, MASK = kmap, mask
    L_ARR, B_ARR, D_ARR, W_ARR, REGION_ARR = l_arr, b_arr, d_arr, w_arr, region_arr


# ===========================================================================
# Stacking core
# ===========================================================================
def stack_index_range(start: int, stop: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Accumulate weighted kappa sums over objects ``[start, stop)``.

    Returns ``(sum_wk, sum_w, n_used)`` flattened to GRID_SIZE**2, where
    ``n_used`` counts objects with at least one unmasked grid pixel.
    """
    sum_wk = np.zeros(GRID_SIZE**2, dtype=np.float64)
    sum_w = np.zeros(GRID_SIZE**2, dtype=np.float64)
    n_used = 0

    for lo in range(start, stop, STACK_BATCH):
        hi = min(lo + STACK_BATCH, stop)
        l_gal = L_ARR[lo:hi, None]
        b_gal = B_ARR[lo:hi, None]
        d_gal = D_ARR[lo:hi, None]
        w_gal = W_ARR[lo:hi, None]

        cosb = np.clip(np.cos(np.radians(b_gal)), 1e-6, None)
        l_grid = l_gal + (OFF_X[None, :] / d_gal) * DEG_PER_RAD / cosb
        b_grid = b_gal + (OFF_Y[None, :] / d_gal) * DEG_PER_RAD

        theta = np.clip(np.radians(90.0 - b_grid), 0.0, np.pi)
        pix = hp.ang2pix(NSIDE, theta, np.radians(l_grid))

        w = w_gal * MASK[pix]
        n_used += int(np.count_nonzero(w.sum(axis=1)))
        sum_wk += (w * KMAP[pix]).sum(axis=0)
        sum_w += w.sum(axis=0)

    return sum_wk, sum_w, n_used


def process_chunk(chunk_meta: tuple[int, int, int]) -> dict[str, object]:
    """Stack one chunk, splitting it by jackknife region.

    Objects are region-sorted before chunking, so a chunk spans only one or
    two regions and the returned accumulator stays a few hundred kilobytes
    regardless of how many regions the tessellation has.
    """
    chunk_id, start, stop = chunk_meta
    t0 = time.time()

    labels = REGION_ARR[start:stop]
    # Sorted input => region blocks are contiguous, so boundaries are enough.
    edges = np.flatnonzero(np.diff(labels)) + 1
    bounds = np.concatenate(([0], edges, [len(labels)]))

    regions, sums_wk, sums_w, n_used, n_objects = [], [], [], [], []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        wk, w, used = stack_index_range(start + int(lo), start + int(hi))
        regions.append(int(labels[lo]))
        sums_wk.append(wk)
        sums_w.append(w)
        n_used.append(used)
        n_objects.append(int(hi - lo))

    return {
        "chunk_id": chunk_id,
        "chunk_start": start,
        "chunk_stop": stop,
        "regions": np.asarray(regions, dtype=np.int32),
        "sum_wk": np.asarray(sums_wk, dtype=np.float64),
        "sum_w": np.asarray(sums_w, dtype=np.float64),
        "n_used": np.asarray(n_used, dtype=np.int64),
        "n_objects": np.asarray(n_objects, dtype=np.int64),
        "elapsed_seconds": time.time() - t0,
    }


# ===========================================================================
# Catalog loading
# ===========================================================================
def load_joint_catalog(
    args: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load every survey region and concatenate into one joint catalog.

    Returns ``(ra, dec, z, weights, region_source)`` where ``region_source``
    records which survey region each object came from (bookkeeping only -- it
    plays no part in the estimator).
    """
    if args.dataset == "BOSS":
        z_min, z_max = 0.4, 0.7
        weight_scheme = "CMASS" if args.catalog_type == "galaxy" else False
    else:
        z_min, z_max = 0.0, 10000.0
        weight_scheme = args.catalog_type == "galaxy"

    fraction = args.fraction if args.fraction < 1.0 else None
    ra, dec, z, weights, source = [], [], [], [], []

    for idx, region in enumerate(args.regions):
        path = resolve_catalog_path(args.data_dir, args.dataset, region, args.catalog_type)
        logger.info("Loading %s %s: %s", region, args.catalog_type, path)
        if args.catalog_type == "random":
            data = load_catalog_lightweight(
                path,
                columns=("RA", "DEC", "Z"),
                fraction=fraction,
                z_min=z_min,
                z_max=z_max,
                seed=args.seed,
            )
            w = np.ones(len(data), dtype=np.float64)
        else:
            data, w = load_catalog(
                path,
                weights=weight_scheme,
                random_fraction=fraction,
                z_min=z_min,
                z_max=z_max,
                seed=args.seed,
            )
        ra.append(np.asarray(data["RA"], dtype=np.float64))
        dec.append(np.asarray(data["DEC"], dtype=np.float64))
        z.append(np.asarray(data["Z"], dtype=np.float64))
        weights.append(np.asarray(w, dtype=np.float64))
        source.append(np.full(len(data), idx, dtype=np.int8))
        logger.info("  %s: %d objects", region, len(data))

    return (
        np.concatenate(ra),
        np.concatenate(dec),
        np.concatenate(z),
        np.concatenate(weights),
        np.concatenate(source),
    )


def load_footprint_tracer(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    """Load the random positions that define the joint jackknife footprint.

    Always a seeded subsample so the tessellation is reproducible across the
    galaxy and random stacking runs that must share it.
    """
    z_min, z_max = (0.4, 0.7) if args.dataset == "BOSS" else (0.0, 10000.0)
    ra, dec = [], []
    for region in args.regions:
        path = resolve_catalog_path(args.data_dir, args.dataset, region, "random")
        data = load_catalog_lightweight(
            path,
            columns=("RA", "DEC", "Z"),
            fraction=args.jk_fraction,
            z_min=z_min,
            z_max=z_max,
            seed=args.jk_seed,
        )
        ra.append(np.asarray(data["RA"], dtype=np.float64))
        dec.append(np.asarray(data["DEC"], dtype=np.float64))
    return np.concatenate(ra), np.concatenate(dec)


def get_regions(args: argparse.Namespace, regions_path: str) -> JackknifeRegions:
    """Load the cached tessellation, or build and cache it."""
    if os.path.exists(regions_path) and not args.rebuild_regions:
        regions = JackknifeRegions.load(regions_path)
        logger.info(
            "Loaded jackknife regions from %s (%d regions, digest=%s)",
            regions_path,
            regions.n_regions,
            regions.digest,
        )
        if regions.nside != args.jk_nside:
            raise ValueError(
                f"Cached regions use nside={regions.nside} but --jk-nside={args.jk_nside}. "
                "Pass --rebuild-regions to rebuild (and re-run every stack that uses them)."
            )
        return regions

    logger.info("Building jackknife tessellation from the joint random footprint ...")
    ra, dec = load_footprint_tracer(args)
    logger.info("Footprint tracer: %d randoms across %s", len(ra), ", ".join(args.regions))
    regions = build_jackknife_regions(ra, dec, nside=args.jk_nside, min_fill=args.jk_min_fill)
    os.makedirs(os.path.dirname(regions_path), exist_ok=True)
    regions.save(regions_path)
    return regions


# ===========================================================================
# Reducer / checkpoint
# ===========================================================================
def init_reducer(n_regions: int, n_chunks: int) -> dict[str, object]:
    return {
        "sum_wk": np.zeros((n_regions, GRID_SIZE**2), dtype=np.float64),
        "sum_w": np.zeros((n_regions, GRID_SIZE**2), dtype=np.float64),
        "n_used": np.zeros(n_regions, dtype=np.int64),
        "n_objects": np.zeros(n_regions, dtype=np.int64),
        "completed_chunks": set(),
        "n_chunks": n_chunks,
    }


def apply_chunk_result(state: dict[str, object], result: dict[str, object]) -> None:
    chunk_id = int(result["chunk_id"])
    if chunk_id in state["completed_chunks"]:
        raise RuntimeError(f"Chunk {chunk_id} applied twice.")
    for k, region in enumerate(result["regions"]):
        state["sum_wk"][region] += result["sum_wk"][k]
        state["sum_w"][region] += result["sum_w"][k]
        state["n_used"][region] += int(result["n_used"][k])
        state["n_objects"][region] += int(result["n_objects"][k])
    state["completed_chunks"].add(chunk_id)


def save_checkpoint_atomic(path: str, state: dict[str, object], args_json: str) -> None:
    tmp = f"{path}.tmp"
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(tmp, "wb") as handle:
        np.savez(
            handle,
            sum_wk=state["sum_wk"],
            sum_w=state["sum_w"],
            n_used=state["n_used"],
            n_objects=state["n_objects"],
            completed_chunks=np.asarray(sorted(state["completed_chunks"]), dtype=np.int64),
            n_chunks=np.int64(state["n_chunks"]),
            args_json=np.array(args_json),
        )
    os.replace(tmp, path)


def load_checkpoint(path: str) -> dict[str, object]:
    with np.load(path, allow_pickle=False) as data:
        return {
            "sum_wk": data["sum_wk"],
            "sum_w": data["sum_w"],
            "n_used": data["n_used"],
            "n_objects": data["n_objects"],
            "completed_chunks": {int(x) for x in data["completed_chunks"].tolist()},
            "n_chunks": int(data["n_chunks"]),
            "args_json": str(data["args_json"]),
        }


def build_chunk_ranges(n_objects: int, chunk_size: int) -> list[tuple[int, int, int]]:
    return [
        (i, start, min(start + chunk_size, n_objects))
        for i, start in enumerate(range(0, n_objects, chunk_size))
    ]


# Rows per slice for the elementwise preprocessing passes.  Both the distance
# integral and the ICRS->Galactic rotation are elementwise, so slicing is
# bitwise identical to a single call -- but done whole on a 41M-row random
# catalog they allocate ~1 GB temporaries each and an 8 GB machine swaps for
# over an hour before stacking even starts.
PREPROCESS_CHUNK = 2_000_000


def comoving_distance_hmpc(z: np.ndarray) -> np.ndarray:
    """Comoving distance in h^-1 Mpc, matching ``preprocess_catalog_galactic``.

    Interpolating D(z) on a grid was tried and rejected.  A grid accurate to
    1e-8 in distance still shifts the stacked map by ~3e-3 relative on a
    1500-galaxy test, because grid points sitting on a HEALPix boundary flip
    to the neighbouring pixel; the effect averages down as 1/sqrt(N) but there
    is no reason to accept it.  The exact call costs ~3.6 min on 44M rows,
    which is noise against a multi-hour stack.
    """
    z = np.asarray(z, dtype=np.float64)
    out = np.empty(len(z), dtype=np.float64)
    for lo in range(0, len(z), PREPROCESS_CHUNK):
        hi = min(lo + PREPROCESS_CHUNK, len(z))
        out[lo:hi] = cosmo.comoving_distance(z[lo:hi]).value * cosmo.h
    return out


def galactic_in_chunks(
    ra: np.ndarray, dec: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """``fast_icrs_to_galactic`` over slices, to bound peak memory.

    The shared implementation builds a (3, N) stack of direction cosines and a
    rotated copy of it -- two ~1 GB allocations at 41M rows.  Slicing gives
    identical results without the spike.
    """
    l_out = np.empty(len(ra), dtype=np.float64)
    b_out = np.empty(len(ra), dtype=np.float64)
    for lo in range(0, len(ra), PREPROCESS_CHUNK):
        hi = min(lo + PREPROCESS_CHUNK, len(ra))
        l_out[lo:hi], b_out[lo:hi] = fast_icrs_to_galactic(ra[lo:hi], dec[lo:hi])
    return l_out, b_out


def mean_map(sum_wk: np.ndarray, sum_w: np.ndarray) -> np.ndarray:
    """Weighted-mean kappa map from summed accumulators."""
    out = np.zeros_like(sum_wk)
    good = sum_w > 0
    out[good] = sum_wk[good] / sum_w[good]
    return out.reshape(GRID_SIZE, GRID_SIZE)


# ===========================================================================
# CLI
# ===========================================================================
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stack single objects over a joint footprint with per-region "
        "jackknife accumulators.",
    )
    parser.add_argument("--dataset", default="BOSS", choices=["BOSS", "eBOSS"])
    parser.add_argument(
        "--regions",
        default="North,South",
        help="Comma-separated survey regions stacked jointly (default: North,South).",
    )
    parser.add_argument("--catalog-type", required=True, choices=["galaxy", "random"])
    parser.add_argument(
        "--fraction",
        type=float,
        default=1.0,
        help="Fraction of the catalog to stack (default: 1.0).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Seed for catalog subsampling when --fraction < 1. Auto-generated "
        "and recorded in the checkpoint if omitted.",
    )
    parser.add_argument(
        "--sigma-crit-weight",
        dest="sigma_crit_weight",
        action="store_true",
        default=True,
        help="Apply 1/Sigma_crit^2 inverse-variance weights (default: on).",
    )
    parser.add_argument(
        "--no-sigma-crit-weight",
        dest="sigma_crit_weight",
        action="store_false",
        help="Disable inverse-variance weighting (unweighted null test).",
    )
    parser.add_argument("--jk-nside", type=int, default=10,
                        help="HEALPix nside of the jackknife tessellation (default: 10, "
                             "34.4 deg^2 cells -> ~275 CMASS regions).")
    parser.add_argument("--jk-min-fill", type=float, default=0.7,
                        help="Minimum cell occupancy, as a fraction of the mean occupied "
                             "cell, to stand as its own region (default: 0.7).")
    parser.add_argument("--jk-fraction", type=float, default=0.10,
                        help="Random subsample fraction used to trace the footprint "
                             "(default: 0.10).")
    parser.add_argument("--jk-seed", type=int, default=20260727,
                        help="Seed for the footprint tracer subsample. Must match across "
                             "the galaxy and random runs (default: 20260727).")
    parser.add_argument("--rebuild-regions", action="store_true",
                        help="Rebuild the tessellation even if a cached one exists.")
    parser.add_argument("--n-processes", type=int, default=None,
                        help="Worker processes (default: 75%% of cores). Use 1 locally.")
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--checkpoint-interval", type=int, default=50)
    parser.add_argument("--resume-checkpoint", action="store_true")
    parser.add_argument("--label", default=None,
                        help="Extra tag inserted into output filenames.")
    parser.add_argument("--output-dir", default="analysis/boss/results")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument(
        "--fwhm-arcmin", type=float, default=FWHM_ARCMIN,
        help="Gaussian FWHM applied to the Planck kappa map, in arcmin "
             f"(default: {FWHM_ARCMIN}). Lower values keep more small-scale "
             "signal but admit more reconstruction noise. Outputs are NOT "
             "auto-tagged -- pass --label to keep products separate.")
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args(argv)
    args.regions = [r.strip() for r in args.regions.split(",") if r.strip()]
    if not args.regions:
        raise ValueError("--regions must name at least one survey region.")
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be positive.")
    if args.n_processes is not None and args.n_processes <= 0:
        raise ValueError("--n-processes must be positive.")
    if not 0.0 < args.fraction <= 1.0:
        raise ValueError("--fraction must lie in (0, 1].")
    return args


def output_paths(args: argparse.Namespace) -> dict[str, str]:
    region_label = "_".join(args.regions)
    tag = "_scw" if args.sigma_crit_weight else ""
    if args.label:
        tag = f"{tag}_{args.label}"
    stem = f"{args.catalog_type}{tag}_{args.dataset}_{region_label}"
    return {
        "region_label": region_label,
        "tag": tag,
        "regions_npz": os.path.join(
            args.output_dir, "jk",
            f"regions_{args.dataset}_{region_label}_nside{args.jk_nside}.npz",
        ),
        "acc_npz": os.path.join(args.output_dir, "jk", f"acc_single_{stem}.npz"),
        "kappa_csv": os.path.join(args.output_dir, f"kappa_single_{stem}.csv"),
        "meta_json": os.path.join(args.output_dir, "jk", f"acc_single_{stem}.meta.json"),
        "checkpoint": os.path.join(
            args.output_dir, "checkpoints", f"single_jk_{stem}.npz"
        ),
    }


def serialize_config(args: argparse.Namespace, regions: JackknifeRegions) -> str:
    payload = {k: v for k, v in vars(args).items()
               if k not in {"overwrite", "resume_checkpoint", "n_processes",
                            "checkpoint_interval", "rebuild_regions"}}
    payload["jk_digest"] = regions.digest
    payload["n_regions"] = regions.n_regions
    payload["grid_size"] = GRID_SIZE
    return json.dumps(payload, sort_keys=True)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    t0 = time.time()
    started_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")

    paths = output_paths(args)
    os.makedirs(os.path.join(args.output_dir, "jk"), exist_ok=True)

    if args.resume_checkpoint:
        if not os.path.exists(paths["checkpoint"]):
            raise FileNotFoundError(f"Checkpoint not found: {paths['checkpoint']}")
    elif os.path.exists(paths["acc_npz"]) and not args.overwrite:
        logger.info("Output exists: %s (use --overwrite). Skipping.", paths["acc_npz"])
        return

    # Resolve the subsampling seed before loading, so a resume reproduces the
    # same subset -- otherwise accumulator state would mix two catalogs.
    if args.fraction < 1.0:
        if args.resume_checkpoint:
            saved = json.loads(load_checkpoint(paths["checkpoint"])["args_json"])
            if saved.get("seed") is None:
                raise ValueError(
                    "Cannot resume: checkpoint has no seed but --fraction < 1.0."
                )
            if args.seed is not None and args.seed != saved["seed"]:
                raise ValueError(
                    f"--seed {args.seed} conflicts with checkpoint seed {saved['seed']}."
                )
            args.seed = saved["seed"]
        elif args.seed is None:
            args.seed = int(np.random.SeedSequence().entropy) % (2**31)
            logger.info("Auto-generated subsampling seed: %d", args.seed)

    logger.info("=" * 62)
    logger.info("stack_single_jk.py -- joint-footprint stack with region accumulators")
    logger.info("=" * 62)
    logger.info("Dataset / regions : %s / %s", args.dataset, ", ".join(args.regions))
    logger.info("Catalog type      : %s", args.catalog_type)
    logger.info("Fraction          : %.4f (seed %s)", args.fraction, args.seed)
    logger.info("Sigma_crit weights: %s", args.sigma_crit_weight)

    jk_regions = get_regions(args, paths["regions_npz"])

    ra, dec, z, weights, source = load_joint_catalog(args)
    logger.info("Joint catalog: %d objects", len(ra))

    if args.sigma_crit_weight:
        scw = sigma_crit_weights(z)
        weights = weights * scw
        logger.info(
            "Applied 1/Sigma_crit^2 weights (relative range %.3f-%.3f)",
            scw.min() / scw.mean(), scw.max() / scw.mean(),
        )

    # Comoving distance and the validity cut come first: dropping rows after
    # region assignment would desynchronize the labels from the coordinates.
    d_arr = comoving_distance_hmpc(z)
    valid = np.isfinite(d_arr) & (d_arr > 0) & np.isfinite(ra) & np.isfinite(dec)
    if not valid.all():
        logger.info("Dropping %d object(s) with invalid position or distance.",
                    int((~valid).sum()))
        ra, dec, weights, d_arr, source = (
            arr[valid] for arr in (ra, dec, weights, d_arr, source)
        )

    labels = jk_regions.assign(ra, dec)
    per_region = np.bincount(labels, minlength=jk_regions.n_regions)
    logger.info(
        "Objects per jackknife region: min=%d median=%d max=%d",
        per_region.min(), int(np.median(per_region)), per_region.max(),
    )
    if (per_region == 0).any():
        logger.warning(
            "%d region(s) hold no objects of this catalog; they contribute no "
            "leave-one-out variance.", int((per_region == 0).sum()),
        )

    # Region-sorting keeps each chunk inside one or two regions, which bounds
    # the accumulator a worker has to ship back to the parent.
    order = np.argsort(labels, kind="stable")
    ra, dec, weights, d_arr, labels = (
        arr[order] for arr in (ra, dec, weights, d_arr, labels)
    )
    l_arr, b_arr = galactic_in_chunks(ra, dec)
    w_arr = weights
    del ra, dec, z, weights, source

    chunk_ranges = build_chunk_ranges(len(l_arr), args.chunk_size)
    if not chunk_ranges:
        raise RuntimeError("No objects left to stack.")

    alm_path, mask_path = resolve_planck_paths(args.data_dir)
    logger.info("Loading Planck map ...")
    kmap, planck_mask, _ = load_kappa_map(
        alm_file=alm_path, mask_file=mask_path, fwhm_arcmin=args.fwhm_arcmin)

    _set_runtime_globals(
        kmap=kmap, mask=planck_mask,
        l_arr=np.ascontiguousarray(l_arr, dtype=np.float64),
        b_arr=np.ascontiguousarray(b_arr, dtype=np.float64),
        d_arr=np.ascontiguousarray(d_arr, dtype=np.float64),
        w_arr=np.ascontiguousarray(w_arr, dtype=np.float64),
        region_arr=np.ascontiguousarray(labels, dtype=np.int32),
    )

    config_json = serialize_config(args, jk_regions)
    state = init_reducer(jk_regions.n_regions, len(chunk_ranges))

    if args.resume_checkpoint:
        ckpt = load_checkpoint(paths["checkpoint"])
        if ckpt["args_json"] != config_json:
            raise ValueError("Checkpoint configuration does not match this run.")
        if ckpt["n_chunks"] != len(chunk_ranges):
            raise ValueError("Checkpoint chunk layout does not match this catalog.")
        state.update({k: ckpt[k] for k in ("sum_wk", "sum_w", "n_used", "n_objects")})
        state["completed_chunks"] = ckpt["completed_chunks"]
        logger.info("Resumed with %d completed chunks.", len(state["completed_chunks"]))

    pending = [c for c in chunk_ranges if c[0] not in state["completed_chunks"]]
    logger.info("Processing %d of %d chunks.", len(pending), len(chunk_ranges))

    def handle(result: dict[str, object]) -> None:
        apply_chunk_result(state, result)
        done = len(state["completed_chunks"])
        if done % max(1, len(chunk_ranges) // 200) == 0 or done == len(chunk_ranges):
            logger.info(
                "  %d/%d chunks | %.1f min elapsed",
                done, len(chunk_ranges), (time.time() - t0) / 60.0,
            )
        if done % args.checkpoint_interval == 0:
            save_checkpoint_atomic(paths["checkpoint"], state, config_json)

    n_processes = args.n_processes or max(1, int(mp.cpu_count() * 0.75))
    if n_processes == 1:
        logger.info("Single-process mode.")
        for chunk in pending:
            handle(process_chunk(chunk))
    else:
        try:
            ctx = mp.get_context("fork")
        except ValueError as exc:
            raise RuntimeError(
                "Multiprocessing needs the 'fork' start method so workers inherit "
                "the kappa map. Use --n-processes 1, or run on Linux."
            ) from exc
        logger.info("Multiprocessing with %d fork workers.", n_processes)
        with ctx.Pool(processes=n_processes) as pool:
            for result in pool.imap_unordered(process_chunk, pending):
                handle(result)

    save_checkpoint_atomic(paths["checkpoint"], state, config_json)

    total_wk = state["sum_wk"].sum(axis=0)
    total_w = state["sum_w"].sum(axis=0)
    kappa = mean_map(total_wk, total_w)

    np.savez(
        paths["acc_npz"],
        sum_wk=state["sum_wk"],
        sum_w=state["sum_w"],
        n_used=state["n_used"],
        n_objects=state["n_objects"],
        seed_pix=jk_regions.seed_pix,
        jk_digest=np.array(jk_regions.digest),
        jk_nside=np.int32(jk_regions.nside),
        grid_size=np.int32(GRID_SIZE),
        box_size_hmpc=np.float64(BOX_SIZE_HMPC),
        args_json=np.array(config_json),
    )
    logger.info("Saved per-region accumulators -> %s", paths["acc_npz"])

    # Radially symmetrized map, matching the stack_single.py convention so the
    # two pipelines can be diffed directly. The unsymmetrized accumulators in
    # the npz are what combine_jackknife.py profiles.
    pd.DataFrame(symmetrize_map(kappa)).to_csv(paths["kappa_csv"], index=True)
    logger.info("Saved symmetrized kappa map -> %s", paths["kappa_csv"])

    runtime = time.time() - t0
    with open(paths["meta_json"], "w", encoding="ascii") as handle:
        json.dump(
            {
                "args": json.loads(config_json),
                "n_regions": jk_regions.n_regions,
                "n_objects_stacked": int(state["n_objects"].sum()),
                "n_objects_unmasked": int(state["n_used"].sum()),
                "n_chunks": len(chunk_ranges),
                "n_processes": n_processes,
                "started_at": started_at,
                "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "runtime_seconds": runtime,
                "accumulator_npz": paths["acc_npz"],
                "kappa_csv": paths["kappa_csv"],
            },
            handle, indent=2, sort_keys=True,
        )
        handle.write("\n")

    logger.info("Done in %.1f s (%.1f min).", runtime, runtime / 60.0)


if __name__ == "__main__":
    main()
