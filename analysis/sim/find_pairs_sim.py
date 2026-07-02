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

    center_x = (x[i1] + 0.5 * dx) % box_size_hmpc
    center_y = (y[i1] + 0.5 * dy) % box_size_hmpc

    pairs = pd.DataFrame(
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
    return pairs


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

    ensure_parent(args.output)
    pairs.to_csv(args.output, index=False)
    logger.info("Saved %d pairs -> %s in %.1f s", len(pairs), args.output, time.time() - t0)


if __name__ == "__main__":
    main()
