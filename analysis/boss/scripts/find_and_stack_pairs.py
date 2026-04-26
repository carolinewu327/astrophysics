#!/usr/bin/env python
"""Find galaxy/random pairs and stack them without materializing a pair CSV.

This script is intended for large catalogs, especially full random catalogs,
where the ``find_pairs.py -> stack_pairs.py`` workflow becomes too large to
store or shuffle through disk. It keeps the pair-finding criteria identical to
``find_pairs.py`` and applies the same pair-stacking geometry as
``stack_pairs.py``, but accumulates the stacked map directly in memory.
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
    preprocess_catalog_galactic,
    resolve_catalog_path,
    resolve_planck_paths,
    setup_logging,
)
from constants import BOX_SIZE_HMPC
from geometry import angular_separation, reflect_symmetrize_map

logger = logging.getLogger(__name__)

# Pair stacks use a 101x101 grid so the center pixel is well-defined.
GRID_RES: int = 101
HALF_SIZE: float = BOX_SIZE_HMPC / 2.0

# Globals shared by workers. They are populated once in the parent process
# before any worker pool is created.
KMAP: np.ndarray | None = None
MASK: np.ndarray | None = None
NSIDE: int | None = None
X_GRID: np.ndarray | None = None
Y_GRID: np.ndarray | None = None

L_ARR: np.ndarray | None = None
B_ARR: np.ndarray | None = None
D_ARR: np.ndarray | None = None
Z_ARR: np.ndarray | None = None
WEIGHTS_ARR: np.ndarray | None = None
SORTED_INDICES: np.ndarray | None = None
SORTED_D: np.ndarray | None = None

R_PAR_MAX: float | None = None
R_PERP_MIN: float | None = None
R_PERP_MAX: float | None = None


def _set_runtime_globals(
    *,
    kmap: np.ndarray,
    mask: np.ndarray,
    nside: int,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    l_arr: np.ndarray,
    b_arr: np.ndarray,
    d_arr: np.ndarray,
    z_arr: np.ndarray,
    weights_arr: np.ndarray,
    sorted_indices: np.ndarray,
    sorted_d: np.ndarray,
    r_par_max: float,
    r_perp_min: float,
    r_perp_max: float,
) -> None:
    """Populate module globals used by the worker path."""
    global KMAP, MASK, NSIDE, X_GRID, Y_GRID
    global L_ARR, B_ARR, D_ARR, Z_ARR, WEIGHTS_ARR, SORTED_INDICES, SORTED_D
    global R_PAR_MAX, R_PERP_MIN, R_PERP_MAX

    KMAP = kmap
    MASK = mask
    NSIDE = nside
    X_GRID = x_grid
    Y_GRID = y_grid

    L_ARR = l_arr
    B_ARR = b_arr
    D_ARR = d_arr
    Z_ARR = z_arr
    WEIGHTS_ARR = weights_arr
    SORTED_INDICES = sorted_indices
    SORTED_D = sorted_d

    R_PAR_MAX = r_par_max
    R_PERP_MIN = r_perp_min
    R_PERP_MAX = r_perp_max


def _require_globals() -> None:
    """Raise if worker globals are missing."""
    required = {
        "KMAP": KMAP,
        "MASK": MASK,
        "NSIDE": NSIDE,
        "X_GRID": X_GRID,
        "Y_GRID": Y_GRID,
        "L_ARR": L_ARR,
        "B_ARR": B_ARR,
        "D_ARR": D_ARR,
        "Z_ARR": Z_ARR,
        "WEIGHTS_ARR": WEIGHTS_ARR,
        "SORTED_INDICES": SORTED_INDICES,
        "SORTED_D": SORTED_D,
        "R_PAR_MAX": R_PAR_MAX,
        "R_PERP_MIN": R_PERP_MIN,
        "R_PERP_MAX": R_PERP_MAX,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise RuntimeError(f"Worker globals not initialized: {', '.join(missing)}")


def stack_pairs_inline(i_final: np.ndarray, j_final: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Stack valid pairs directly into worker-local accumulators."""
    _require_globals()

    sum_wk = np.zeros((GRID_RES, GRID_RES), dtype=np.float64)
    sum_w = np.zeros((GRID_RES, GRID_RES), dtype=np.float64)

    if len(i_final) == 0:
        return sum_wk, sum_w, 0

    l1_arr = np.deg2rad(L_ARR[i_final])
    b1_arr = np.deg2rad(B_ARR[i_final])
    l2_arr = np.deg2rad(L_ARR[j_final])
    b2_arr = np.deg2rad(B_ARR[j_final])
    dmid_arr = cosmo.comoving_distance((Z_ARR[i_final] + Z_ARR[j_final]) / 2.0).value * cosmo.h
    pair_weights = WEIGHTS_ARR[i_final] * WEIGHTS_ARR[j_final]

    n_skipped = 0

    for idx in range(len(i_final)):
        l1 = l1_arr[idx]
        b1 = b1_arr[idx]
        l2 = l2_arr[idx]
        b2 = b2_arr[idx]

        dl_raw = (l2 - l1 + np.pi) % (2 * np.pi) - np.pi
        if dl_raw < 0:
            l1, l2 = l2, l1
            b1, b2 = b2, b1
            dl_raw = -dl_raw

        lc = l1 + 0.5 * dl_raw
        bc = 0.5 * (b1 + b2)
        dc = dmid_arr[idx]

        dl = dl_raw * np.cos(bc)
        db = b2 - b1
        norm = np.hypot(dl, db)
        if norm == 0 or not np.isfinite(norm):
            n_skipped += 1
            continue

        cos_theta = dl / norm
        sin_theta = db / norm

        dl_cosbc = (cos_theta * (X_GRID / dc)) - (sin_theta * (Y_GRID / dc))
        db_grid = (sin_theta * (X_GRID / dc)) + (cos_theta * (Y_GRID / dc))
        l_grid = lc + dl_cosbc / np.cos(bc)
        b_grid = bc + db_grid

        l_grid_deg_wrapped = np.mod(np.rad2deg(l_grid), 360.0)
        b_grid_deg = np.rad2deg(b_grid)

        theta = np.radians(90.0 - b_grid_deg)
        if np.min(theta) < 0 or np.max(theta) > np.pi:
            n_skipped += 1
            continue

        phi = np.radians(l_grid_deg_wrapped)
        pix = hp.ang2pix(NSIDE, theta.ravel(), phi.ravel())

        kappa_vals = KMAP[pix].reshape(GRID_RES, GRID_RES)
        mask_vals = MASK[pix].reshape(GRID_RES, GRID_RES)

        valid_mask = (mask_vals != 0) & np.isfinite(kappa_vals)
        kappa_vals = np.where(valid_mask, kappa_vals, 0.0)

        weight = pair_weights[idx]
        sum_wk += weight * kappa_vals
        sum_w += weight * valid_mask.astype(np.float64)

    return sum_wk, sum_w, n_skipped


def process_chunk_and_stack(chunk_meta: tuple[int, int, int]) -> dict[str, object]:
    """Find valid pairs in a chunk and stack them immediately."""
    _require_globals()

    chunk_id, chunk_start, chunk_end = chunk_meta
    t0 = time.time()

    i_blocks: list[np.ndarray] = []
    j_blocks: list[np.ndarray] = []

    for idx in range(chunk_start, chunk_end):
        i = SORTED_INDICES[idx]
        dc1 = D_ARR[i]

        left_bound = np.searchsorted(SORTED_D, dc1 - R_PAR_MAX, side="left")
        right_bound = np.searchsorted(SORTED_D, dc1 + R_PAR_MAX, side="right")
        search_start = max(idx + 1, left_bound)
        search_end = right_bound

        if search_end <= search_start:
            continue

        j_candidates = SORTED_INDICES[search_start:search_end]
        dc2_candidates = D_ARR[j_candidates]

        l1_val = L_ARR[i]
        b1_val = B_ARR[i]
        l2_candidates = L_ARR[j_candidates]
        b2_candidates = B_ARR[j_candidates]

        dl = np.abs(l2_candidates - l1_val)
        dl = np.minimum(dl, 360.0 - dl)
        db = np.abs(b2_candidates - b1_val)

        d_avg = (dc1 + dc2_candidates) / 2.0
        min_angle = R_PERP_MIN / d_avg
        max_angle = R_PERP_MAX / d_avg

        cosb1 = np.cos(np.radians(b1_val))
        rough_angle = np.sqrt((np.radians(dl * cosb1)) ** 2 + (np.radians(db)) ** 2)
        angle_mask = (rough_angle >= min_angle * 0.7) & (rough_angle <= max_angle * 2.0)

        if not np.any(angle_mask):
            continue

        j_filtered = j_candidates[angle_mask]
        if len(j_filtered) == 0:
            continue

        i_blocks.append(np.full(len(j_filtered), i, dtype=np.int64))
        j_blocks.append(j_filtered.astype(np.int64, copy=False))

    if not i_blocks:
        elapsed = time.time() - t0
        return {
            "chunk_id": chunk_id,
            "chunk_start": chunk_start,
            "chunk_end": chunk_end,
            "sum_wk": np.zeros((GRID_RES, GRID_RES), dtype=np.float64),
            "sum_w": np.zeros((GRID_RES, GRID_RES), dtype=np.float64),
            "n_pairs": 0,
            "n_skipped": 0,
            "elapsed_seconds": elapsed,
        }

    i_arr = np.concatenate(i_blocks)
    j_arr = np.concatenate(j_blocks)

    theta = angular_separation(L_ARR[i_arr], B_ARR[i_arr], L_ARR[j_arr], B_ARR[j_arr])
    d_avg_arr = (D_ARR[i_arr] + D_ARR[j_arr]) / 2.0
    r_perp_arr = d_avg_arr * theta
    final_mask = (r_perp_arr >= R_PERP_MIN) & (r_perp_arr <= R_PERP_MAX)

    i_final = i_arr[final_mask]
    j_final = j_arr[final_mask]

    sum_wk, sum_w, n_skipped = stack_pairs_inline(i_final, j_final)
    elapsed = time.time() - t0

    return {
        "chunk_id": chunk_id,
        "chunk_start": chunk_start,
        "chunk_end": chunk_end,
        "sum_wk": sum_wk,
        "sum_w": sum_w,
        "n_pairs": int(len(i_final)),
        "n_skipped": int(n_skipped),
        "elapsed_seconds": elapsed,
    }


def _format_number(value: float) -> str:
    """Format floats for default labels without introducing dots."""
    if float(value).is_integer():
        return str(int(value))
    return str(value).replace(".", "p")


def make_output_path(output_dir: str, label: str, dataset: str, region: str) -> str:
    """Return the output CSV path for a stacked pair map."""
    filename = f"kappa_pairs_{label}_{dataset}_{region}.csv"
    return os.path.join(output_dir, filename)


def make_checkpoint_path(output_dir: str, label: str, dataset: str, region: str) -> str:
    """Return the default checkpoint path for a run."""
    filename = f"kappa_pairs_{label}_{dataset}_{region}.npz"
    return os.path.join(output_dir, "checkpoints", filename)


def build_chunk_ranges(n_objects: int, chunk_size: int) -> list[tuple[int, int, int]]:
    """Split the sorted catalog into chunk descriptors."""
    n_chunks = int(np.ceil(n_objects / chunk_size))
    return [
        (chunk_id, chunk_id * chunk_size, min((chunk_id + 1) * chunk_size, n_objects))
        for chunk_id in range(n_chunks)
    ]


def serialize_run_config(args: argparse.Namespace, output_path: str) -> dict[str, object]:
    """Build the subset of run config that must match on resume.

    The seed is included so a resumed run is forced to use the same catalog
    subset as the initial run; otherwise checkpoint accumulator state would
    be silently combined with chunks from a different subsample.
    """
    return {
        "catalog_type": args.catalog_type,
        "chunk_size": int(args.chunk_size),
        "data_dir": args.data_dir,
        "dataset": args.dataset,
        "fraction": float(args.fraction),
        "grid_res": GRID_RES,
        "label": args.label,
        "output_path": output_path,
        "region": args.region,
        "rpar": float(args.rpar),
        "rperp_max": float(args.rperp_max),
        "rperp_min": float(args.rperp_min),
        "seed": int(args.seed) if args.seed is not None else None,
    }


def init_reducer(chunk_ranges: list[tuple[int, int, int]], output_path: str) -> dict[str, object]:
    """Create the reducer state accumulated by the main process."""
    return {
        "total_sum_wk": np.zeros((GRID_RES, GRID_RES), dtype=np.float64),
        "total_sum_w": np.zeros((GRID_RES, GRID_RES), dtype=np.float64),
        "total_pairs": 0,
        "total_skipped": 0,
        "completed_chunks": set(),
        "chunk_ranges": np.asarray(chunk_ranges, dtype=np.int64),
        "output_path": output_path,
    }


def apply_chunk_result(reducer_state: dict[str, object], result: dict[str, object]) -> None:
    """Merge one worker result into the reducer state."""
    reducer_state["total_sum_wk"] += result["sum_wk"]
    reducer_state["total_sum_w"] += result["sum_w"]
    reducer_state["total_pairs"] += int(result["n_pairs"])
    reducer_state["total_skipped"] += int(result["n_skipped"])
    reducer_state["completed_chunks"].add(int(result["chunk_id"]))


def save_checkpoint_atomic(
    checkpoint_path: str,
    reducer_state: dict[str, object],
    args_dict: dict[str, object],
) -> None:
    """Write the current reducer state to an atomic `.npz` checkpoint."""
    checkpoint_dir = os.path.dirname(checkpoint_path)
    if checkpoint_dir:
        os.makedirs(checkpoint_dir, exist_ok=True)

    tmp_path = f"{checkpoint_path}.tmp"
    with open(tmp_path, "wb") as handle:
        np.savez(
            handle,
            total_sum_wk=reducer_state["total_sum_wk"],
            total_sum_w=reducer_state["total_sum_w"],
            total_pairs=np.int64(reducer_state["total_pairs"]),
            total_skipped=np.int64(reducer_state["total_skipped"]),
            completed_chunks=np.asarray(
                sorted(reducer_state["completed_chunks"]),
                dtype=np.int32,
            ),
            chunk_ranges=reducer_state["chunk_ranges"],
            args_json=np.array(json.dumps(args_dict, sort_keys=True)),
            output_path=np.array(reducer_state["output_path"]),
            saved_at=np.array(time.time()),
        )
    os.replace(tmp_path, checkpoint_path)


def load_checkpoint(checkpoint_path: str) -> dict[str, object]:
    """Load a reducer checkpoint from disk."""
    with np.load(checkpoint_path, allow_pickle=False) as data:
        return {
            "total_sum_wk": data["total_sum_wk"],
            "total_sum_w": data["total_sum_w"],
            "total_pairs": int(data["total_pairs"]),
            "total_skipped": int(data["total_skipped"]),
            "completed_chunks": set(int(x) for x in data["completed_chunks"].tolist()),
            "chunk_ranges": data["chunk_ranges"],
            "args_dict": json.loads(str(data["args_json"])),
            "output_path": str(data["output_path"]),
        }


def finalize_map(reducer_state: dict[str, object]) -> np.ndarray:
    """Compute the mean map and apply reflection symmetrization."""
    total_sum_wk = reducer_state["total_sum_wk"]
    total_sum_w = reducer_state["total_sum_w"]
    kappa_mean = np.zeros_like(total_sum_wk)
    nonzero = total_sum_w > 0
    kappa_mean[nonzero] = total_sum_wk[nonzero] / total_sum_w[nonzero]
    return reflect_symmetrize_map(kappa_mean)


def write_metadata_json(path: str, metadata: dict[str, object]) -> None:
    """Write a small JSON sidecar for provenance."""
    with open(path, "w", encoding="ascii") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def load_preprocessed_catalog(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load the requested catalog and return preprocessed arrays."""
    fits_path = resolve_catalog_path(args.data_dir, args.dataset, args.region, args.catalog_type)
    logger.info("Loading catalog from %s", fits_path)

    if args.dataset == "BOSS":
        z_min, z_max = 0.4, 0.7
        weight_scheme = "CMASS" if args.catalog_type == "galaxy" else False
    else:
        z_min, z_max = 0, 10000
        weight_scheme = True if args.catalog_type == "galaxy" else False

    random_fraction = args.fraction if args.fraction < 1.0 else None
    subsample_seed = args.seed if random_fraction is not None else None

    if args.catalog_type == "random":
        data = load_catalog_lightweight(
            fits_path,
            columns=("RA", "DEC", "Z"),
            fraction=random_fraction,
            z_min=z_min,
            z_max=z_max,
            seed=subsample_seed,
        )
        weights = np.ones(len(data), dtype=np.float64)
    else:
        data, weights = load_catalog(
            fits_path,
            weights=weight_scheme,
            random_fraction=random_fraction,
            z_min=z_min,
            z_max=z_max,
            seed=subsample_seed,
        )

    l_arr, b_arr, d_arr, data_filtered, weights_valid = preprocess_catalog_galactic(data, weights)
    z_arr = np.asarray(data_filtered["Z"], dtype=np.float64)
    weights_arr = np.asarray(weights_valid, dtype=np.float64)

    logger.info("Catalog preprocessed: %d objects", len(z_arr))
    return (
        np.asarray(l_arr, dtype=np.float64),
        np.asarray(b_arr, dtype=np.float64),
        np.asarray(d_arr, dtype=np.float64),
        z_arr,
        weights_arr,
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Find valid pairs and stack them directly against the Planck kappa map.",
    )
    parser.add_argument(
        "--dataset",
        type=str,
        default="BOSS",
        choices=["BOSS", "eBOSS"],
        help="Survey dataset (default: BOSS).",
    )
    parser.add_argument(
        "--region",
        type=str,
        default="North",
        help="Catalog region, e.g. North/South for BOSS (default: North).",
    )
    parser.add_argument(
        "--catalog-type",
        type=str,
        default="random",
        choices=["galaxy", "random"],
        help="Type of catalog to process (default: random).",
    )
    parser.add_argument(
        "--fraction",
        type=float,
        default=1.0,
        help="Fraction of the catalog to use (default: 1.0).",
    )
    parser.add_argument(
        "--rpar",
        type=float,
        default=20.0,
        help="Max parallel (line-of-sight) distance in Mpc/h (default: 20).",
    )
    parser.add_argument(
        "--rperp-min",
        type=float,
        default=18.0,
        help="Min perpendicular distance in Mpc/h (default: 18).",
    )
    parser.add_argument(
        "--rperp-max",
        type=float,
        default=22.0,
        help="Max perpendicular distance in Mpc/h (default: 22).",
    )
    parser.add_argument(
        "--n-processes",
        type=int,
        default=None,
        help="Number of processes. Use 1 for a local serial test (default: 75%% of cores).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Objects per chunk (default: 5000).",
    )
    parser.add_argument(
        "--label",
        type=str,
        default=None,
        help="Label used in output filenames. Default: catalog type + rpar.",
    )
    parser.add_argument(
        "--checkpoint-path",
        type=str,
        default=None,
        help="Optional explicit path for the reducer checkpoint (.npz).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help=(
            "Seed for catalog subsampling RNG. Required-but-auto-generated "
            "when --fraction < 1.0; ignored when --fraction == 1.0. "
            "Auto-generated seeds are saved in the checkpoint so resume uses "
            "the same subset. Pass explicitly to reproduce a specific run."
        ),
    )
    parser.add_argument(
        "--checkpoint-interval",
        type=int,
        default=10,
        help="Save a checkpoint every N completed chunks (default: 10).",
    )
    parser.add_argument(
        "--resume-checkpoint",
        action="store_true",
        help="Resume from an existing checkpoint instead of starting fresh.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="analysis/boss/results",
        help="Directory for stacked outputs (default: analysis/boss/results).",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="data",
        help="Root data directory (default: data).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output CSV.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """CLI entry point."""
    args = parse_args(argv)
    setup_logging()

    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be positive.")
    if args.n_processes is not None and args.n_processes <= 0:
        raise ValueError("--n-processes must be positive.")
    if args.checkpoint_interval <= 0:
        raise ValueError("--checkpoint-interval must be positive.")

    if args.label is None:
        args.label = f"{args.catalog_type}_{_format_number(args.rpar)}"

    output_path = make_output_path(args.output_dir, args.label, args.dataset, args.region)
    checkpoint_path = args.checkpoint_path or make_checkpoint_path(
        args.output_dir,
        args.label,
        args.dataset,
        args.region,
    )
    metadata_path = output_path.replace(".csv", ".meta.json")

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    if args.resume_checkpoint:
        if not os.path.exists(checkpoint_path):
            raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")
        if os.path.exists(output_path):
            raise FileExistsError(
                f"Output already exists at {output_path}; remove it before resuming."
            )
    elif os.path.exists(output_path) and not args.overwrite:
        logger.info("Output already exists: %s  (use --overwrite to replace). Skipping.", output_path)
        return

    # Resolve subsampling seed BEFORE loading the catalog.
    # Resume must use the checkpoint's seed so the catalog is identical to
    # the initial run; otherwise accumulator state would combine with chunks
    # from a different subsample (silently wrong output).
    if args.fraction < 1.0:
        if args.resume_checkpoint:
            checkpoint_peek = load_checkpoint(checkpoint_path)
            saved_seed = checkpoint_peek["args_dict"].get("seed")
            if saved_seed is None:
                raise ValueError(
                    "Cannot resume: checkpoint was written without a seed but "
                    "--fraction < 1.0. The original run's catalog subset cannot "
                    "be reproduced. Re-run from scratch instead."
                )
            if args.seed is not None and args.seed != saved_seed:
                raise ValueError(
                    f"--seed {args.seed} conflicts with checkpoint seed {saved_seed}. "
                    "Omit --seed on resume."
                )
            args.seed = saved_seed
            logger.info("Resume: using seed %d from checkpoint.", args.seed)
        elif args.seed is None:
            args.seed = int(np.random.SeedSequence().entropy)
            logger.info("Auto-generated subsampling seed: %d", args.seed)

    n_processes = args.n_processes
    if n_processes is None:
        n_processes = max(1, int(mp.cpu_count() * 0.75))

    logger.info("=" * 60)
    logger.info("find_and_stack_pairs.py - Combined Pair Finder + Stacker")
    logger.info("=" * 60)
    logger.info("Dataset:            %s", args.dataset)
    logger.info("Region:             %s", args.region)
    logger.info("Catalog type:       %s", args.catalog_type)
    logger.info("Fraction:           %s", args.fraction)
    logger.info("Subsample seed:     %s", args.seed if args.seed is not None else "(none)")
    logger.info("r_par_max:          %.3f Mpc/h", args.rpar)
    logger.info("r_perp range:       %.3f to %.3f Mpc/h", args.rperp_min, args.rperp_max)
    logger.info("Processes:          %d", n_processes)
    logger.info("Chunk size:         %d", args.chunk_size)
    logger.info("Output CSV:         %s", output_path)
    logger.info("Checkpoint path:    %s", checkpoint_path)

    t0 = time.time()
    started_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")

    l_arr, b_arr, d_arr, z_arr, weights_arr = load_preprocessed_catalog(args)
    sorted_indices = np.argsort(d_arr)
    sorted_d = d_arr[sorted_indices]
    chunk_ranges = build_chunk_ranges(len(d_arr), args.chunk_size)

    if not chunk_ranges:
        raise RuntimeError("No valid objects found after preprocessing.")

    x_vals = np.linspace(-HALF_SIZE, HALF_SIZE, GRID_RES)
    y_vals = np.linspace(-HALF_SIZE, HALF_SIZE, GRID_RES)
    x_grid, y_grid = np.meshgrid(x_vals, y_vals)

    alm_file, mask_file = resolve_planck_paths(args.data_dir)
    logger.info("Loading Planck map from %s and %s", alm_file, mask_file)
    kmap, mask, nside = load_kappa_map(alm_file=alm_file, mask_file=mask_file)

    _set_runtime_globals(
        kmap=kmap,
        mask=mask,
        nside=nside,
        x_grid=x_grid,
        y_grid=y_grid,
        l_arr=l_arr,
        b_arr=b_arr,
        d_arr=d_arr,
        z_arr=z_arr,
        weights_arr=weights_arr,
        sorted_indices=sorted_indices,
        sorted_d=sorted_d,
        r_par_max=args.rpar,
        r_perp_min=args.rperp_min,
        r_perp_max=args.rperp_max,
    )

    serialized_args = serialize_run_config(args, output_path)
    reducer_state = init_reducer(chunk_ranges, output_path)

    if args.resume_checkpoint:
        checkpoint = load_checkpoint(checkpoint_path)
        if checkpoint["args_dict"] != serialized_args:
            raise ValueError(
                "Checkpoint metadata does not match the current run configuration."
            )
        if not np.array_equal(checkpoint["chunk_ranges"], reducer_state["chunk_ranges"]):
            raise ValueError("Checkpoint chunk layout does not match the current catalog.")

        reducer_state["total_sum_wk"] = checkpoint["total_sum_wk"]
        reducer_state["total_sum_w"] = checkpoint["total_sum_w"]
        reducer_state["total_pairs"] = checkpoint["total_pairs"]
        reducer_state["total_skipped"] = checkpoint["total_skipped"]
        reducer_state["completed_chunks"] = checkpoint["completed_chunks"]
        logger.info(
            "Resumed checkpoint with %d completed chunks.",
            len(reducer_state["completed_chunks"]),
        )

    pending_chunks = [
        chunk_meta
        for chunk_meta in chunk_ranges
        if chunk_meta[0] not in reducer_state["completed_chunks"]
    ]
    total_chunks = len(chunk_ranges)

    logger.info(
        "Processing %d pending chunk(s) out of %d total.",
        len(pending_chunks),
        total_chunks,
    )

    def handle_result(result: dict[str, object]) -> None:
        apply_chunk_result(reducer_state, result)
        completed = len(reducer_state["completed_chunks"])
        elapsed = time.time() - t0
        logger.info(
            "Chunk %d/%d (id=%d [%d:%d]): %d pair(s), %d skipped in %.1fs | total pairs: %d | elapsed: %.1f min",
            completed,
            total_chunks,
            int(result["chunk_id"]),
            int(result["chunk_start"]),
            int(result["chunk_end"]),
            int(result["n_pairs"]),
            int(result["n_skipped"]),
            float(result["elapsed_seconds"]),
            int(reducer_state["total_pairs"]),
            elapsed / 60.0,
        )
        if completed % args.checkpoint_interval == 0:
            save_checkpoint_atomic(checkpoint_path, reducer_state, serialized_args)
            logger.info("Checkpoint saved to %s", checkpoint_path)

    if n_processes == 1:
        logger.info("Running in single-process mode.")
        for chunk_meta in pending_chunks:
            handle_result(process_chunk_and_stack(chunk_meta))
    else:
        try:
            ctx = mp.get_context("fork")
        except ValueError as exc:
            raise RuntimeError(
                "Multiprocessing requires the 'fork' start method. "
                "Use --n-processes 1 locally, or run on Linux for multi-process execution."
            ) from exc

        logger.info("Running with multiprocessing using the 'fork' start method.")
        with ctx.Pool(processes=n_processes) as pool:
            for result in pool.imap_unordered(process_chunk_and_stack, pending_chunks):
                handle_result(result)

    final_map = finalize_map(reducer_state)
    pd.DataFrame(final_map).to_csv(output_path, index=True)

    finished_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    runtime_seconds = time.time() - t0
    metadata = {
        "catalog_type": args.catalog_type,
        "checkpoint_interval": args.checkpoint_interval,
        "checkpoint_path": checkpoint_path,
        "chunk_size": args.chunk_size,
        "dataset": args.dataset,
        "finished_at": finished_at,
        "fraction": args.fraction,
        "grid_res": GRID_RES,
        "label": args.label,
        "n_chunks_completed": len(reducer_state["completed_chunks"]),
        "seed": args.seed,
        "n_chunks_total": total_chunks,
        "n_processes": n_processes,
        "output_csv": output_path,
        "region": args.region,
        "rpar": args.rpar,
        "rperp_max": args.rperp_max,
        "rperp_min": args.rperp_min,
        "runtime_seconds": runtime_seconds,
        "started_at": started_at,
        "total_pairs": int(reducer_state["total_pairs"]),
        "total_skipped": int(reducer_state["total_skipped"]),
    }
    write_metadata_json(metadata_path, metadata)

    logger.info("Saved stacked map to %s", output_path)
    logger.info("Saved metadata to %s", metadata_path)
    logger.info("Completed in %.1f s (%.1f min).", runtime_seconds, runtime_seconds / 60.0)


if __name__ == "__main__":
    main()
