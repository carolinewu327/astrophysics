#!/usr/bin/env python
"""Project BigMDPL downsampled particles into a 2D kappa memmap."""

from __future__ import annotations

import argparse
import logging
import os
import time

import numpy as np
import pandas as pd

from sim_utils import (
    BIGMDPL_BOX_SIZE_HMPC,
    BIGMDPL_DOWNSAMPLE_FACTOR,
    BIGMDPL_PARTICLE_MASS_HMSUN,
    KappaMapInfo,
    SIGMA_C_HMSUN_PER_MPC2,
    ensure_parent,
    setup_logging,
    write_kappa_metadata,
)


logger = logging.getLogger(__name__)


def iter_particle_xy_chunks(
    particle_file: str,
    x_col: int,
    y_col: int,
    chunksize: int,
    max_particles: int | None = None,
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


def build_counts_memmap(
    particle_file: str,
    counts_path: str,
    n_grid: int,
    pixel_size_hmpc: float,
    box_size_hmpc: float,
    x_col: int,
    y_col: int,
    chunksize: int,
    max_particles: int | None,
) -> tuple[np.memmap, int]:
    ensure_parent(counts_path)
    if os.path.exists(counts_path):
        os.remove(counts_path)

    counts = np.memmap(counts_path, dtype="uint32", mode="w+", shape=(n_grid, n_grid))
    counts[:] = 0
    counts.flush()

    total = 0
    t0 = time.time()

    for chunk_id, chunk in enumerate(
        iter_particle_xy_chunks(particle_file, x_col, y_col, chunksize, max_particles),
        start=1,
    ):
        x = np.mod(chunk["x"].to_numpy(dtype=np.float64, copy=False), box_size_hmpc)
        y = np.mod(chunk["y"].to_numpy(dtype=np.float64, copy=False), box_size_hmpc)
        ix = np.floor(x / pixel_size_hmpc).astype(np.int64)
        iy = np.floor(y / pixel_size_hmpc).astype(np.int64)
        ix = np.clip(ix, 0, n_grid - 1)
        iy = np.clip(iy, 0, n_grid - 1)
        np.add.at(counts, (iy, ix), 1)
        total += len(chunk)

        if chunk_id == 1 or chunk_id % 10 == 0:
            logger.info(
                "Binned %.3g particles into %dx%d grid in %.1f s",
                total,
                n_grid,
                n_grid,
                time.time() - t0,
            )
            counts.flush()

    counts.flush()
    return counts, total


def counts_to_kappa(
    counts: np.memmap,
    output: str,
    total_particles: int,
    pixel_size_hmpc: float,
    box_size_hmpc: float,
    mp_eff_hmsun: float,
    sigma_c: float,
) -> np.memmap:
    ensure_parent(output)
    if os.path.exists(output):
        os.remove(output)

    kappa = np.memmap(output, dtype="float32", mode="w+", shape=counts.shape)
    mean_count_density = total_particles / box_size_hmpc**2

    for row0 in range(0, counts.shape[0], 1024):
        row1 = min(row0 + 1024, counts.shape[0])
        sigma = mp_eff_hmsun * (
            counts[row0:row1].astype(np.float64) / pixel_size_hmpc**2 - mean_count_density
        )
        kappa[row0:row1] = (sigma / sigma_c).astype("float32")

    kappa.flush()
    return kappa


def make_kappa_map(
    particle_file: str,
    output: str,
    pixel_size_hmpc: float,
    box_size_hmpc: float,
    x_col: int,
    y_col: int,
    chunksize: int,
    max_particles: int | None = None,
    keep_counts: bool = False,
) -> None:
    n_grid_float = box_size_hmpc / pixel_size_hmpc
    n_grid = int(round(n_grid_float))
    if not np.isclose(n_grid, n_grid_float):
        raise ValueError("box_size_hmpc must be exactly divisible by pixel_size_hmpc")

    counts_path = f"{output}.counts.uint32"
    logger.info("Building count map: %s", counts_path)
    counts, total_particles = build_counts_memmap(
        particle_file=particle_file,
        counts_path=counts_path,
        n_grid=n_grid,
        pixel_size_hmpc=pixel_size_hmpc,
        box_size_hmpc=box_size_hmpc,
        x_col=x_col,
        y_col=y_col,
        chunksize=chunksize,
        max_particles=max_particles,
    )

    logger.info("Converting counts to float32 kappa map: %s", output)
    mp_eff = BIGMDPL_DOWNSAMPLE_FACTOR * BIGMDPL_PARTICLE_MASS_HMSUN
    counts_to_kappa(
        counts=counts,
        output=output,
        total_particles=total_particles,
        pixel_size_hmpc=pixel_size_hmpc,
        box_size_hmpc=box_size_hmpc,
        mp_eff_hmsun=mp_eff,
        sigma_c=SIGMA_C_HMSUN_PER_MPC2,
    )

    info = KappaMapInfo(
        path=os.path.abspath(output),
        shape=(n_grid, n_grid),
        dtype="float32",
        pixel_size_hmpc=pixel_size_hmpc,
        box_size_hmpc=box_size_hmpc,
        total_particles=total_particles,
    )
    write_kappa_metadata(
        info,
        extra={
            "particle_file": particle_file,
            "x_col": x_col,
            "y_col": y_col,
            "mp_eff_hmsun": mp_eff,
            "sigma_c_hmsun_per_mpc2": SIGMA_C_HMSUN_PER_MPC2,
            "counts_path": counts_path if keep_counts else None,
        },
    )

    del counts
    if not keep_counts and os.path.exists(counts_path):
        os.remove(counts_path)
    logger.info("Finished kappa map with %d particles", total_particles)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Project downsampled BigMDPL particles into a mean-subtracted kappa map.",
    )
    parser.add_argument(
        "--particle-file",
        default="data/cosmosim/dm_particles_snap_030.dat.bz2",
        help="Whitespace particle file. Default assumes bz2 text with x/y in columns 0/1.",
    )
    parser.add_argument(
        "--output",
        default="analysis/sim/results/kappa_map_l0p5.float32",
        help="Raw float32 memmap output path.",
    )
    parser.add_argument(
        "--pixel-size",
        type=float,
        default=0.5,
        help="Map pixel size in h^-1 Mpc. Use 0.1 for the final 25000^2 map.",
    )
    parser.add_argument(
        "--box-size",
        type=float,
        default=BIGMDPL_BOX_SIZE_HMPC,
        help="Simulation box size in h^-1 Mpc.",
    )
    parser.add_argument("--x-col", type=int, default=0, help="Zero-based particle x column.")
    parser.add_argument("--y-col", type=int, default=1, help="Zero-based particle y column.")
    parser.add_argument("--chunksize", type=int, default=1_000_000, help="Particle rows per chunk.")
    parser.add_argument("--max-particles", type=int, default=None, help="Optional smoke-test cap.")
    parser.add_argument(
        "--keep-counts",
        action="store_true",
        help="Keep the intermediate uint32 count memmap next to the kappa output.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    make_kappa_map(
        particle_file=args.particle_file,
        output=args.output,
        pixel_size_hmpc=args.pixel_size,
        box_size_hmpc=args.box_size,
        x_col=args.x_col,
        y_col=args.y_col,
        chunksize=args.chunksize,
        max_particles=args.max_particles,
        keep_counts=args.keep_counts,
    )


if __name__ == "__main__":
    main()

