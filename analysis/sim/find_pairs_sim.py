#!/usr/bin/env python
"""Find BigMDPL host-halo pairs with RSD applied only to LOS selection."""

from __future__ import annotations

import argparse
import logging
import os
import time

import numpy as np
import pandas as pd

from sim_utils import (
    BIGMDPL_BOX_SIZE_HMPC,
    BIGMDPL_OMEGA_M,
    BIGMDPL_Z_SNAPSHOT,
    ensure_parent,
    minimal_periodic_delta,
    rsd_displacement_hmpc,
    setup_logging,
)


logger = logging.getLogger(__name__)


PAIR_COLUMNS = [
    "id1",
    "id2",
    "x1",
    "y1",
    "z1",
    "z1_redshiftspace",
    "x2",
    "y2",
    "z2",
    "z2_redshiftspace",
    "M1",
    "M2",
    "r_perp",
    "r_parallel",
    "pair_center_x",
    "pair_center_y",
    "cos_theta",
    "sin_theta",
]


def make_pair_frame(
    halos: pd.DataFrame,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    z_rsd: np.ndarray,
    i1: np.ndarray,
    i2: np.ndarray,
    dx: np.ndarray,
    dy: np.ndarray,
    rperp: np.ndarray,
    rpar: np.ndarray,
    box_size_hmpc: float,
) -> pd.DataFrame:
    center_x = (x[i1] + 0.5 * dx) % box_size_hmpc
    center_y = (y[i1] + 0.5 * dy) % box_size_hmpc

    return pd.DataFrame(
        {
            "id1": halos["id"].to_numpy()[i1],
            "id2": halos["id"].to_numpy()[i2],
            "x1": x[i1],
            "y1": y[i1],
            "z1": z[i1],
            "z1_redshiftspace": z_rsd[i1],
            "x2": x[i2],
            "y2": y[i2],
            "z2": z[i2],
            "z2_redshiftspace": z_rsd[i2],
            "M1": halos["Mvir"].to_numpy()[i1],
            "M2": halos["Mvir"].to_numpy()[i2],
            "r_perp": rperp,
            "r_parallel": rpar,
            "pair_center_x": center_x,
            "pair_center_y": center_y,
            "cos_theta": dx / rperp,
            "sin_theta": dy / rperp,
        },
        columns=PAIR_COLUMNS,
    )


def find_pairs(
    halos: pd.DataFrame,
    rperp_min: float,
    rperp_max: float,
    rpar_max: float,
    box_size_hmpc: float,
    z_snapshot: float,
    omega_m: float,
    max_pairs: int | None,
) -> pd.DataFrame:
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError("scipy is required for periodic pair finding") from exc

    x = halos["x"].to_numpy(dtype=np.float64, copy=False) % box_size_hmpc
    y = halos["y"].to_numpy(dtype=np.float64, copy=False) % box_size_hmpc
    z = halos["z"].to_numpy(dtype=np.float64, copy=False) % box_size_hmpc
    vz = halos["vz"].to_numpy(dtype=np.float64, copy=False)
    z_rsd = (z + rsd_displacement_hmpc(vz, z_snapshot, omega_m)) % box_size_hmpc

    xy = np.column_stack([x, y])
    tree = cKDTree(xy, boxsize=box_size_hmpc)
    candidates = tree.query_pairs(r=rperp_max, output_type="ndarray")
    logger.info("Found %d projected candidates with r_perp <= %.3f", len(candidates), rperp_max)

    if len(candidates) == 0:
        return pd.DataFrame(columns=PAIR_COLUMNS)

    i = candidates[:, 0]
    j = candidates[:, 1]
    dx = minimal_periodic_delta(x[i], x[j], box_size_hmpc)
    dy = minimal_periodic_delta(y[i], y[j], box_size_hmpc)
    rperp = np.hypot(dx, dy)
    rpar = np.abs(minimal_periodic_delta(z_rsd[i], z_rsd[j], box_size_hmpc))

    mask = (rperp >= rperp_min) & (rperp <= rperp_max) & (rpar <= rpar_max) & (rperp > 0)
    i = i[mask]
    j = j[mask]
    dx = dx[mask]
    dy = dy[mask]
    rperp = rperp[mask]
    rpar = rpar[mask]

    swap = (dx < 0) | ((dx == 0) & (dy < 0))
    i1 = np.where(swap, j, i)
    i2 = np.where(swap, i, j)
    dx = np.where(swap, -dx, dx)
    dy = np.where(swap, -dy, dy)

    if max_pairs is not None and len(i1) > max_pairs:
        i1 = i1[:max_pairs]
        i2 = i2[:max_pairs]
        dx = dx[:max_pairs]
        dy = dy[:max_pairs]
        rperp = rperp[:max_pairs]
        rpar = rpar[:max_pairs]

    pairs = make_pair_frame(
        halos=halos,
        x=x,
        y=y,
        z=z,
        z_rsd=z_rsd,
        i1=i1,
        i2=i2,
        dx=dx,
        dy=dy,
        rperp=rperp,
        rpar=rpar,
        box_size_hmpc=box_size_hmpc,
    )
    return pairs


def filter_neighbor_chunk(
    halos: pd.DataFrame,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    z_rsd: np.ndarray,
    neighbor_lists: list[np.ndarray],
    chunk_start: int,
    rperp_min: float,
    rperp_max: float,
    rpar_max: float,
    box_size_hmpc: float,
) -> tuple[pd.DataFrame, int]:
    i1_parts = []
    i2_parts = []
    dx_parts = []
    dy_parts = []
    rperp_parts = []
    rpar_parts = []
    projected_candidates = 0

    for local_i, neighbors in enumerate(neighbor_lists):
        i = chunk_start + local_i
        j = np.asarray(neighbors, dtype=np.int64)
        j = j[j > i]
        projected_candidates += len(j)
        if len(j) == 0:
            continue

        dx = minimal_periodic_delta(x[i], x[j], box_size_hmpc)
        dy = minimal_periodic_delta(y[i], y[j], box_size_hmpc)
        rperp = np.hypot(dx, dy)
        rpar = np.abs(minimal_periodic_delta(z_rsd[i], z_rsd[j], box_size_hmpc))

        mask = (rperp >= rperp_min) & (rperp <= rperp_max) & (rpar <= rpar_max) & (rperp > 0)
        if not np.any(mask):
            continue

        j = j[mask]
        dx = dx[mask]
        dy = dy[mask]
        rperp = rperp[mask]
        rpar = rpar[mask]

        swap = (dx < 0) | ((dx == 0) & (dy < 0))
        i_arr = np.full(len(j), i, dtype=np.int64)
        i1 = np.where(swap, j, i_arr)
        i2 = np.where(swap, i_arr, j)
        dx = np.where(swap, -dx, dx)
        dy = np.where(swap, -dy, dy)

        i1_parts.append(i1)
        i2_parts.append(i2)
        dx_parts.append(dx)
        dy_parts.append(dy)
        rperp_parts.append(rperp)
        rpar_parts.append(rpar)

    if not i1_parts:
        return pd.DataFrame(columns=PAIR_COLUMNS), projected_candidates

    i1 = np.concatenate(i1_parts)
    i2 = np.concatenate(i2_parts)
    dx = np.concatenate(dx_parts)
    dy = np.concatenate(dy_parts)
    rperp = np.concatenate(rperp_parts)
    rpar = np.concatenate(rpar_parts)

    return (
        make_pair_frame(
            halos=halos,
            x=x,
            y=y,
            z=z,
            z_rsd=z_rsd,
            i1=i1,
            i2=i2,
            dx=dx,
            dy=dy,
            rperp=rperp,
            rpar=rpar,
            box_size_hmpc=box_size_hmpc,
        ),
        projected_candidates,
    )


def write_pairs_streaming(
    halos: pd.DataFrame,
    output: str,
    rperp_min: float,
    rperp_max: float,
    rpar_max: float,
    box_size_hmpc: float,
    z_snapshot: float,
    omega_m: float,
    max_pairs: int | None,
    candidate_chunk_size: int,
) -> int:
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError("scipy is required for periodic pair finding") from exc

    x = halos["x"].to_numpy(dtype=np.float64, copy=False) % box_size_hmpc
    y = halos["y"].to_numpy(dtype=np.float64, copy=False) % box_size_hmpc
    z = halos["z"].to_numpy(dtype=np.float64, copy=False) % box_size_hmpc
    vz = halos["vz"].to_numpy(dtype=np.float64, copy=False)
    z_rsd = (z + rsd_displacement_hmpc(vz, z_snapshot, omega_m)) % box_size_hmpc

    xy = np.column_stack([x, y])
    tree = cKDTree(xy, boxsize=box_size_hmpc)

    tmp_output = f"{output}.tmp"
    if os.path.exists(tmp_output):
        os.remove(tmp_output)

    wrote_header = False
    total_pairs = 0
    total_candidates = 0
    n_halos = len(halos)

    for chunk_start in range(0, n_halos, candidate_chunk_size):
        chunk_end = min(chunk_start + candidate_chunk_size, n_halos)
        neighbor_lists = tree.query_ball_point(xy[chunk_start:chunk_end], r=rperp_max, return_sorted=False)
        chunk_pairs, chunk_candidates = filter_neighbor_chunk(
            halos=halos,
            x=x,
            y=y,
            z=z,
            z_rsd=z_rsd,
            neighbor_lists=neighbor_lists,
            chunk_start=chunk_start,
            rperp_min=rperp_min,
            rperp_max=rperp_max,
            rpar_max=rpar_max,
            box_size_hmpc=box_size_hmpc,
        )
        total_candidates += chunk_candidates

        if max_pairs is not None:
            remaining = max_pairs - total_pairs
            if remaining <= 0:
                break
            chunk_pairs = chunk_pairs.iloc[:remaining]

        if len(chunk_pairs) > 0:
            chunk_pairs.to_csv(
                tmp_output,
                mode="a" if wrote_header else "w",
                header=not wrote_header,
                index=False,
            )
            wrote_header = True
            total_pairs += len(chunk_pairs)

        logger.info(
            "Processed halos %d-%d/%d: %d projected candidates, %d kept pairs, %d total pairs",
            chunk_start,
            chunk_end,
            n_halos,
            chunk_candidates,
            len(chunk_pairs),
            total_pairs,
        )

    if not wrote_header:
        pd.DataFrame(columns=PAIR_COLUMNS).to_csv(tmp_output, index=False)

    os.replace(tmp_output, output)
    logger.info("Scanned %d projected candidates with r_perp <= %.3f", total_candidates, rperp_max)
    return total_pairs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find periodic BigMDPL halo pairs for kappa stacking.")
    parser.add_argument("--halos", default="analysis/sim/results/halos_mass13.csv")
    parser.add_argument("--output", default="analysis/sim/results/pairs_mass13_rperp5.csv")
    parser.add_argument("--rperp-min", type=float, default=4.0)
    parser.add_argument("--rperp-max", type=float, default=6.0)
    parser.add_argument("--rpar-max", type=float, default=22.0)
    parser.add_argument("--box-size", type=float, default=BIGMDPL_BOX_SIZE_HMPC)
    parser.add_argument("--z-snapshot", type=float, default=BIGMDPL_Z_SNAPSHOT)
    parser.add_argument("--omega-m", type=float, default=BIGMDPL_OMEGA_M)
    parser.add_argument("--max-halos", type=int, default=None, help="Optional halo cap for smoke tests.")
    parser.add_argument("--max-pairs", type=int, default=None, help="Optional output pair cap.")
    parser.add_argument(
        "--candidate-chunk-size",
        type=int,
        default=10000,
        help="Number of halos to query at a time. Use 0 for the old in-memory query_pairs path.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    t0 = time.time()

    if os.path.exists(args.output) and not args.overwrite:
        logger.info("Output exists: %s (use --overwrite)", args.output)
        return

    halos = pd.read_csv(args.halos)
    if args.max_halos is not None:
        halos = halos.iloc[: args.max_halos].copy()
    logger.info("Loaded %d halos", len(halos))
    logger.info(
        "Using BigMDPL RSD convention: Omega_m=%.3f, H0=100 in h^-1 Mpc units, z=%.3f",
        args.omega_m,
        args.z_snapshot,
    )

    ensure_parent(args.output)
    if args.candidate_chunk_size > 0:
        n_pairs = write_pairs_streaming(
            halos=halos,
            output=args.output,
            rperp_min=args.rperp_min,
            rperp_max=args.rperp_max,
            rpar_max=args.rpar_max,
            box_size_hmpc=args.box_size,
            z_snapshot=args.z_snapshot,
            omega_m=args.omega_m,
            max_pairs=args.max_pairs,
            candidate_chunk_size=args.candidate_chunk_size,
        )
    else:
        pairs = find_pairs(
            halos=halos,
            rperp_min=args.rperp_min,
            rperp_max=args.rperp_max,
            rpar_max=args.rpar_max,
            box_size_hmpc=args.box_size,
            z_snapshot=args.z_snapshot,
            omega_m=args.omega_m,
            max_pairs=args.max_pairs,
        )
        pairs.to_csv(args.output, index=False)
        n_pairs = len(pairs)

    logger.info("Saved %d pairs -> %s in %.1f s", n_pairs, args.output, time.time() - t0)


if __name__ == "__main__":
    main()
