#!/usr/bin/env python
"""Stack a simulated kappa map around BigMDPL halos."""

from __future__ import annotations

import argparse
import logging
import os
import time

import numpy as np
import pandas as pd
from tqdm import tqdm

from sim_utils import (
    BIGMDPL_BOX_SIZE_HMPC,
    BIGMDPL_DOWNSAMPLE_FACTOR,
    BIGMDPL_PARTICLE_MASS_HMSUN,
    OBS_BOX_SIZE_HMPC,
    OBS_SINGLE_GRID_SIZE,
    SIGMA_C_HMSUN_PER_MPC2,
    apply_periodic_gaussian_smoothing,
    ensure_parent,
    open_kappa_memmap,
    periodic_bilinear_sample,
    radial_profile_from_map,
    save_map_csv,
    save_profile_csv,
    setup_logging,
    single_stack_offsets,
)


logger = logging.getLogger(__name__)


def stack_single_halos(
    halos: pd.DataFrame,
    kappa_map: np.ndarray,
    pixel_size_hmpc: float,
    map_box_size_hmpc: float,
    stack_box_size_hmpc: float,
    grid_size: int,
) -> np.ndarray:
    _, off_x, off_y = single_stack_offsets(stack_box_size_hmpc, grid_size)
    stack_sum = np.zeros((grid_size, grid_size), dtype=np.float64)

    x = halos["x"].to_numpy(dtype=np.float64, copy=False)
    y = halos["y"].to_numpy(dtype=np.float64, copy=False)

    for i in tqdm(range(len(halos)), desc="Stacking single halos"):
        xs = x[i] + off_x
        ys = y[i] + off_y
        stack_sum += periodic_bilinear_sample(
            kappa_map,
            xs,
            ys,
            pixel_size_hmpc=pixel_size_hmpc,
            box_size_hmpc=map_box_size_hmpc,
        )

    if len(halos) == 0:
        return np.zeros((grid_size, grid_size), dtype=np.float32)
    return (stack_sum / len(halos)).astype(np.float32)


def iter_particle_xy_chunks(
    particle_file: str,
    x_col: int,
    y_col: int,
    chunksize: int,
    max_particles: int | None,
):
    rows_seen = 0
    reader = pd.read_csv(
        particle_file,
        sep=r"\s+",
        header=None,
        usecols=[x_col, y_col],
        names=["x", "y"],
        chunksize=chunksize,
        compression="infer",
        engine="c",
    )
    for chunk in reader:
        if max_particles is not None:
            remaining = max_particles - rows_seen
            if remaining <= 0:
                break
            chunk = chunk.iloc[:remaining]
        rows_seen += len(chunk)
        yield chunk
        if max_particles is not None and rows_seen >= max_particles:
            break


def direct_annular_profile(
    halos: pd.DataFrame,
    particle_file: str,
    output: str,
    r_max_hmpc: float,
    n_bins: int,
    box_size_hmpc: float,
    x_col: int,
    y_col: int,
    chunksize: int,
    max_particles: int | None,
    max_halos: int | None,
) -> None:
    """Independent validator from direct particle counts in transverse annuli."""

    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError("scipy is required for the direct-annular validator") from exc

    if max_halos is not None and len(halos) > max_halos:
        halos = halos.sample(n=max_halos, random_state=0).reset_index(drop=True)

    halo_xy = halos[["x", "y"]].to_numpy(dtype=np.float64)
    bins = np.linspace(0.0, r_max_hmpc, n_bins + 1)
    counts = np.zeros(n_bins, dtype=np.float64)
    total_particles = 0
    tree_query_halos = halo_xy

    for chunk_id, chunk in enumerate(
        iter_particle_xy_chunks(particle_file, x_col, y_col, chunksize, max_particles),
        start=1,
    ):
        particle_xy = np.mod(chunk[["x", "y"]].to_numpy(dtype=np.float64), box_size_hmpc)
        total_particles += len(particle_xy)
        particle_tree = cKDTree(particle_xy, boxsize=box_size_hmpc)
        neighbors = particle_tree.query_ball_point(tree_query_halos, r=r_max_hmpc)

        for h_idx, p_idx in enumerate(neighbors):
            if not p_idx:
                continue
            pts = particle_xy[np.asarray(p_idx, dtype=np.int64)]
            dx = (pts[:, 0] - halo_xy[h_idx, 0] + 0.5 * box_size_hmpc) % box_size_hmpc
            dx -= 0.5 * box_size_hmpc
            dy = (pts[:, 1] - halo_xy[h_idx, 1] + 0.5 * box_size_hmpc) % box_size_hmpc
            dy -= 0.5 * box_size_hmpc
            radii = np.hypot(dx, dy)
            counts += np.histogram(radii, bins=bins)[0]

        if chunk_id == 1 or chunk_id % 10 == 0:
            logger.info("Direct-annular profile processed %.3g particles", total_particles)

    annulus_area = np.pi * (bins[1:] ** 2 - bins[:-1] ** 2)
    mean_count_density = total_particles / box_size_hmpc**2
    mean_count_annulus = counts / max(len(halos), 1) / annulus_area
    mp_eff = BIGMDPL_DOWNSAMPLE_FACTOR * BIGMDPL_PARTICLE_MASS_HMSUN
    kappa = mp_eff * (mean_count_annulus - mean_count_density) / SIGMA_C_HMSUN_PER_MPC2
    radius = 0.5 * (bins[:-1] + bins[1:])
    save_profile_csv(
        output,
        radius,
        kappa,
        extra_columns={
            "particle_count": counts,
            "annulus_area_hmpc2": annulus_area,
        },
    )
    logger.info("Saved direct-annular profile -> %s", output)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stack simulated kappa around compact halo catalog.")
    parser.add_argument("--halos", default="analysis/sim/results/halos_mass13.csv")
    parser.add_argument("--kappa-map", default="analysis/sim/results/kappa_map_l0p5.float32")
    parser.add_argument("--output", default="analysis/sim/results/kappa_single_sim_mass13.csv")
    parser.add_argument("--profile-output", default="analysis/sim/results/radial_profile_single_sim_mass13.csv")
    parser.add_argument(
        "--direct-profile-output",
        default="analysis/sim/results/radial_profile_single_direct_annuli_mass13.csv",
    )
    parser.add_argument("--box-size", type=float, default=OBS_BOX_SIZE_HMPC)
    parser.add_argument("--grid-size", type=int, default=OBS_SINGLE_GRID_SIZE)
    parser.add_argument("--smooth", choices=["none", "2arcmin", "4arcmin", "8arcmin"], default="none")
    parser.add_argument("--max-halos", type=int, default=None, help="Optional halo cap for smoke tests.")
    parser.add_argument(
        "--direct-particles",
        default=None,
        help="If set, also compute direct-annular validator from this particle file.",
    )
    parser.add_argument("--direct-rmax", type=float, default=50.0)
    parser.add_argument("--direct-bins", type=int, default=50)
    parser.add_argument("--direct-max-halos", type=int, default=500)
    parser.add_argument("--direct-max-particles", type=int, default=None)
    parser.add_argument("--particle-x-col", type=int, default=0)
    parser.add_argument("--particle-y-col", type=int, default=1)
    parser.add_argument("--particle-chunksize", type=int, default=500_000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    t0 = time.time()

    for path in [args.output, args.profile_output]:
        if os.path.exists(path) and not args.overwrite:
            logger.info("Output exists: %s (use --overwrite)", path)
            return

    halos = pd.read_csv(args.halos)
    if args.max_halos is not None:
        halos = halos.iloc[: args.max_halos].copy()
    logger.info("Loaded %d halos", len(halos))

    kappa_map, info = open_kappa_memmap(args.kappa_map)
    map_for_stack = apply_periodic_gaussian_smoothing(kappa_map, info.pixel_size_hmpc, args.smooth)

    stack = stack_single_halos(
        halos=halos,
        kappa_map=map_for_stack,
        pixel_size_hmpc=info.pixel_size_hmpc,
        map_box_size_hmpc=info.box_size_hmpc,
        stack_box_size_hmpc=args.box_size,
        grid_size=args.grid_size,
    )

    save_map_csv(args.output, stack)
    radius, profile = radial_profile_from_map(stack, args.box_size)
    save_profile_csv(args.profile_output, radius, profile)
    logger.info("Saved single stack -> %s", args.output)
    logger.info("Saved radial profile -> %s", args.profile_output)

    if args.direct_particles:
        if os.path.exists(args.direct_profile_output) and not args.overwrite:
            logger.info("Direct profile exists: %s (use --overwrite)", args.direct_profile_output)
        else:
            direct_annular_profile(
                halos=halos,
                particle_file=args.direct_particles,
                output=args.direct_profile_output,
                r_max_hmpc=args.direct_rmax,
                n_bins=args.direct_bins,
                box_size_hmpc=BIGMDPL_BOX_SIZE_HMPC,
                x_col=args.particle_x_col,
                y_col=args.particle_y_col,
                chunksize=args.particle_chunksize,
                max_particles=args.direct_max_particles,
                max_halos=args.direct_max_halos,
            )

    elapsed = time.time() - t0
    logger.info("Done in %.1f s", elapsed)


if __name__ == "__main__":
    main()

