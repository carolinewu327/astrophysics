#!/usr/bin/env python
"""Joint pair stack with per-region jackknife accumulators.

The pair analogue of ``stack_single_jk.py``.  It stacks the North and South
pair catalogs into one joint measurement -- rather than producing two maps for
``plot_results.average_region_maps`` to average 50/50 -- and accumulates
``sum(w*kappa)`` / ``sum(w)`` per jackknife region so leave-one-out estimates
are subtractions rather than re-stacks.

Why this matters for pairs specifically: the filament map is a small residual
between two much larger maps (corrected pairs minus a superposed-singles
control).  Measured with ``weighting_sensitivity.py``, moving from the 50/50
average to a count-weighted one moves the filament bridge excess by -66% at
5 h^-1 Mpc, because modest shifts in the inputs are amplified in the residual.

Pairs are assigned to regions by their **midpoint**, converted from the
catalog's Galactic coordinates back to equatorial so they land on the same
tessellation as the single-galaxy stack.  That sharing is what lets the
filament -- which subtracts a control built from the single stack -- be
jackknifed coherently, deleting the same sky patch from both terms.

Weighting note: pair weights are ``w1 * w2`` only.  1/Sigma_crit^2 weighting is
deliberately *not* applied here, because these galaxy-pair stacks are combined
with the archived random-pair maps, which were stacked without it.  Applying it
to one side of the subtraction and not the other would be inconsistent.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/stack_pairs_jk.py \\
        --pair-catalogs data/paircatalogs/BOSS/galaxy_pairs_BOSS_North_5.0_4.0_6.0hmpc.csv.gz,\\
data/paircatalogs/BOSS/galaxy_pairs_BOSS_South_5.0_4.0_6.0hmpc.csv.gz \\
        --label galaxy_5 --n-processes 6
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

from catalog import load_kappa_map, resolve_planck_paths, setup_logging
from constants import BOX_SIZE_HMPC, FWHM_ARCMIN, NSIDE
from geometry import fast_galactic_to_icrs, reflect_symmetrize_map
from jackknife import JackknifeRegions

logger = logging.getLogger(__name__)

# Pair stacking uses a 101x101 grid so the midpoint pixel is well defined.
GRID_RES: int = 101
HALF_SIZE: float = BOX_SIZE_HMPC / 2.0
STACK_BATCH: int = 64

_X_VALS = np.linspace(-HALF_SIZE, HALF_SIZE, GRID_RES)
X_GRID, Y_GRID = (a.ravel() for a in np.meshgrid(_X_VALS, _X_VALS))

KMAP: np.ndarray | None = None
MASK: np.ndarray | None = None
LC: np.ndarray | None = None
BC: np.ndarray | None = None
COS_T: np.ndarray | None = None
SIN_T: np.ndarray | None = None
DMID: np.ndarray | None = None
WPAIR: np.ndarray | None = None
REGION_ARR: np.ndarray | None = None


def _set_runtime_globals(**kw) -> None:
    globals().update(kw)


def pair_geometry(
    l1: np.ndarray, b1: np.ndarray, l2: np.ndarray, b2: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Midpoint and pair-axis orientation, matching ``stack_pairs.stack_pairs``.

    Inputs are in radians.  Returns ``(lc, bc, cos_theta, sin_theta)`` with the
    galaxy ordering normalized so the pair axis always points toward increasing
    longitude.
    """
    dl_raw = (l2 - l1 + np.pi) % (2 * np.pi) - np.pi
    flip = dl_raw < 0
    # Swapping galaxy 1 and 2 is equivalent to negating the longitude offset
    # and the latitude difference.
    l1 = np.where(flip, l2, l1)
    b1_o = np.where(flip, b2, b1)
    b2_o = np.where(flip, b1, b2)
    dl_raw = np.abs(dl_raw)

    lc = l1 + 0.5 * dl_raw
    bc = 0.5 * (b1_o + b2_o)

    dl = dl_raw * np.cos(bc)
    db = b2_o - b1_o
    norm = np.hypot(dl, db)
    norm = np.where(norm == 0, 1.0, norm)
    return lc, bc, dl / norm, db / norm


def stack_index_range(start: int, stop: int) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Accumulate weighted kappa sums over pairs ``[start, stop)``."""
    sum_wk = np.zeros(GRID_RES**2, dtype=np.float64)
    sum_w = np.zeros(GRID_RES**2, dtype=np.float64)
    n_used = 0
    n_skipped = 0

    for lo in range(start, stop, STACK_BATCH):
        hi = min(lo + STACK_BATCH, stop)
        lc = LC[lo:hi, None]
        bc = BC[lo:hi, None]
        cos_t = COS_T[lo:hi, None]
        sin_t = SIN_T[lo:hi, None]
        dc = DMID[lo:hi, None]
        w = WPAIR[lo:hi, None]

        x_over_d = X_GRID[None, :] / dc
        y_over_d = Y_GRID[None, :] / dc
        dl_cosbc = cos_t * x_over_d - sin_t * y_over_d
        db_grid = sin_t * x_over_d + cos_t * y_over_d

        l_grid = lc + dl_cosbc / np.cos(bc)
        b_grid = bc + db_grid

        # Written as the degree round-trip that stack_pairs.py uses, rather than
        # the algebraically equivalent pi/2 - b_grid, so the two pipelines agree
        # bit-for-bit instead of differing by rounding at HEALPix boundaries.
        theta = np.radians(90.0 - np.degrees(b_grid))
        # Same guard as stack_pairs.py: pairs whose grid runs off the pole are
        # dropped whole, not clipped, so no pair contributes a folded position.
        ok = (theta.min(axis=1) >= 0.0) & (theta.max(axis=1) <= np.pi)
        n_skipped += int((~ok).sum())
        if not ok.any():
            continue

        theta = theta[ok]
        phi = np.radians(np.mod(np.degrees(l_grid[ok]), 360.0))
        w = w[ok]

        pix = hp.ang2pix(NSIDE, theta, phi)
        kappa_vals = KMAP[pix]
        valid = (MASK[pix] != 0) & np.isfinite(kappa_vals)
        kappa_vals = np.where(valid, kappa_vals, 0.0)

        sum_wk += (w * kappa_vals).sum(axis=0)
        sum_w += (w * valid).sum(axis=0)
        n_used += int(ok.sum())

    return sum_wk, sum_w, n_used, n_skipped


def process_chunk(chunk_meta: tuple[int, int, int]) -> dict[str, object]:
    chunk_id, start, stop = chunk_meta
    labels = REGION_ARR[start:stop]
    edges = np.flatnonzero(np.diff(labels)) + 1
    bounds = np.concatenate(([0], edges, [len(labels)]))

    regions, sums_wk, sums_w, n_used, n_pairs, n_skipped = [], [], [], [], [], 0
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        wk, w, used, skipped = stack_index_range(start + int(lo), start + int(hi))
        regions.append(int(labels[lo]))
        sums_wk.append(wk)
        sums_w.append(w)
        n_used.append(used)
        n_pairs.append(int(hi - lo))
        n_skipped += skipped

    return {
        "chunk_id": chunk_id,
        "regions": np.asarray(regions, dtype=np.int32),
        "sum_wk": np.asarray(sums_wk, dtype=np.float64),
        "sum_w": np.asarray(sums_w, dtype=np.float64),
        "n_used": np.asarray(n_used, dtype=np.int64),
        "n_pairs": np.asarray(n_pairs, dtype=np.int64),
        "n_skipped": n_skipped,
    }


def load_pair_catalogs(paths: list[str]) -> pd.DataFrame:
    frames = []
    for path in paths:
        frame = pd.read_csv(path)
        logger.info("Loaded %d pairs from %s", len(frame), path)
        frames.append(frame)
    pairs = pd.concat(frames, ignore_index=True)
    logger.info("Joint pair catalog: %d pairs", len(pairs))
    return pairs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stack pairs jointly with per-region jackknife accumulators.")
    parser.add_argument("--pair-catalogs", required=True,
                        help="Comma-separated pair catalog CSVs (one per survey region).")
    parser.add_argument("--label", required=True,
                        help="Label for output filenames, e.g. galaxy_5.")
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South",
                        help="Survey regions, used to locate the shared tessellation.")
    parser.add_argument("--output-dir", default="analysis/boss/results")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument(
        "--fwhm-arcmin", type=float, default=FWHM_ARCMIN,
        help="Gaussian FWHM applied to the Planck kappa map, in arcmin "
             f"(default: {FWHM_ARCMIN}). Lower values keep more small-scale "
             "signal but admit more reconstruction noise. Outputs are NOT "
             "auto-tagged -- pass --label to keep products separate.")
    parser.add_argument("--jk-nside", type=int, default=10)
    parser.add_argument("--n-processes", type=int, default=None)
    parser.add_argument("--chunk-size", type=int, default=20000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    args.pair_catalogs = [p.strip() for p in args.pair_catalogs.split(",") if p.strip()]
    args.regions = [r.strip() for r in args.regions.split(",") if r.strip()]
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    t0 = time.time()

    region_label = "_".join(args.regions)
    jk_dir = os.path.join(args.output_dir, "jk")
    os.makedirs(jk_dir, exist_ok=True)
    acc_path = os.path.join(jk_dir, f"acc_pairs_{args.label}_{args.dataset}_{region_label}.npz")
    csv_path = os.path.join(
        args.output_dir, f"kappa_pairs_{args.label}_{args.dataset}_{region_label}_joint.csv")
    if os.path.exists(acc_path) and not args.overwrite:
        logger.info("Output exists: %s (use --overwrite). Skipping.", acc_path)
        return

    regions_path = os.path.join(
        jk_dir, f"regions_{args.dataset}_{region_label}_nside{args.jk_nside}.npz")
    if not os.path.exists(regions_path):
        raise FileNotFoundError(
            f"Tessellation not found: {regions_path}. Run stack_single_jk.py first -- "
            "the pair stack must reuse the single stack's regions for the filament "
            "jackknife to be coherent."
        )
    jk_regions = JackknifeRegions.load(regions_path)
    logger.info("Using tessellation %s (%d regions, digest=%s)",
                regions_path, jk_regions.n_regions, jk_regions.digest)

    pairs = load_pair_catalogs(args.pair_catalogs)

    l1 = np.radians(pairs["l1"].to_numpy(dtype=np.float64))
    b1 = np.radians(pairs["b1"].to_numpy(dtype=np.float64))
    l2 = np.radians(pairs["l2"].to_numpy(dtype=np.float64))
    b2 = np.radians(pairs["b2"].to_numpy(dtype=np.float64))
    lc, bc, cos_t, sin_t = pair_geometry(l1, b1, l2, b2)
    dmid = pairs["Dmid"].to_numpy(dtype=np.float64)
    wpair = (pairs["w1"].to_numpy(dtype=np.float64)
             * pairs["w2"].to_numpy(dtype=np.float64))

    ra_mid, dec_mid = fast_galactic_to_icrs(np.degrees(lc), np.degrees(bc))
    labels = jk_regions.assign(ra_mid, dec_mid)
    per_region = np.bincount(labels, minlength=jk_regions.n_regions)
    logger.info("Pairs per jackknife region: min=%d median=%d max=%d; %d region(s) empty",
                per_region.min(), int(np.median(per_region)), per_region.max(),
                int((per_region == 0).sum()))

    order = np.argsort(labels, kind="stable")
    lc, bc, cos_t, sin_t, dmid, wpair, labels = (
        arr[order] for arr in (lc, bc, cos_t, sin_t, dmid, wpair, labels))

    alm_path, mask_path = resolve_planck_paths(args.data_dir)
    logger.info("Loading Planck map ...")
    kmap, planck_mask, _ = load_kappa_map(
        alm_file=alm_path, mask_file=mask_path, fwhm_arcmin=args.fwhm_arcmin)

    _set_runtime_globals(
        KMAP=kmap, MASK=planck_mask,
        LC=np.ascontiguousarray(lc), BC=np.ascontiguousarray(bc),
        COS_T=np.ascontiguousarray(cos_t), SIN_T=np.ascontiguousarray(sin_t),
        DMID=np.ascontiguousarray(dmid), WPAIR=np.ascontiguousarray(wpair),
        REGION_ARR=np.ascontiguousarray(labels, dtype=np.int32),
    )

    chunks = [(i, s, min(s + args.chunk_size, len(labels)))
              for i, s in enumerate(range(0, len(labels), args.chunk_size))]
    sum_wk = np.zeros((jk_regions.n_regions, GRID_RES**2), dtype=np.float64)
    sum_w = np.zeros((jk_regions.n_regions, GRID_RES**2), dtype=np.float64)
    n_used = np.zeros(jk_regions.n_regions, dtype=np.int64)
    n_pairs = np.zeros(jk_regions.n_regions, dtype=np.int64)
    total_skipped = 0
    done = 0

    def handle(result):
        nonlocal total_skipped, done
        for k, reg in enumerate(result["regions"]):
            sum_wk[reg] += result["sum_wk"][k]
            sum_w[reg] += result["sum_w"][k]
            n_used[reg] += int(result["n_used"][k])
            n_pairs[reg] += int(result["n_pairs"][k])
        total_skipped += int(result["n_skipped"])
        done += 1
        if done % max(1, len(chunks) // 20) == 0 or done == len(chunks):
            logger.info("  %d/%d chunks | %.1f min", done, len(chunks),
                        (time.time() - t0) / 60.0)

    n_processes = args.n_processes or max(1, int(mp.cpu_count() * 0.75))
    if n_processes == 1:
        for chunk in chunks:
            handle(process_chunk(chunk))
    else:
        ctx = mp.get_context("fork")
        logger.info("Multiprocessing with %d fork workers.", n_processes)
        with ctx.Pool(processes=n_processes) as pool:
            for result in pool.imap_unordered(process_chunk, chunks):
                handle(result)

    if total_skipped:
        logger.warning("Skipped %d pair(s) whose grid ran off the pole.", total_skipped)

    total_wk, total_w = sum_wk.sum(axis=0), sum_w.sum(axis=0)
    mean = np.zeros_like(total_wk)
    good = total_w > 0
    mean[good] = total_wk[good] / total_w[good]

    np.savez(
        acc_path, sum_wk=sum_wk, sum_w=sum_w, n_used=n_used, n_pairs=n_pairs,
        seed_pix=jk_regions.seed_pix, jk_digest=np.array(jk_regions.digest),
        jk_nside=np.int32(jk_regions.nside), grid_size=np.int32(GRID_RES),
        box_size_hmpc=np.float64(BOX_SIZE_HMPC),
        args_json=np.array(json.dumps(vars(args), sort_keys=True)),
    )
    logger.info("Saved per-region accumulators -> %s", acc_path)

    pd.DataFrame(reflect_symmetrize_map(mean.reshape(GRID_RES, GRID_RES))).to_csv(
        csv_path, index=True)
    logger.info("Saved symmetrized pair map -> %s", csv_path)
    logger.info("Done in %.1f min (%d pairs stacked).",
                (time.time() - t0) / 60.0, int(n_pairs.sum()))


if __name__ == "__main__":
    main()
