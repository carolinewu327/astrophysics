#!/usr/bin/env python
"""Derive coarser and pre-smoothed kappa maps from the 0.1 h^-1 Mpc build.

Subcommands
-----------
block-sum
    Exactly derive a coarser kappa map by block-summing the kept uint32
    count map (only integer block factors are exact: 2 -> 0.2 h^-1 Mpc,
    5 -> 0.5 h^-1 Mpc). Uses the same counts->kappa conversion and
    constants as ``make_kappa_map.py``, read from the source sidecar.

smooth
    Precompute a periodically smoothed copy of a kappa map so stack runs
    can load it directly instead of re-smoothing the full map every run.
    Output sidecar records the smoothing level; pass ``--smooth none`` to
    the stack scripts when using a pre-smoothed map, or they will smooth
    twice.

compare
    Report max |diff| between two kappa maps of equal shape (validation,
    e.g. derived 0.5 vs the directly built 0.5 map).

Examples
--------
    python derive_kappa_maps.py block-sum --factor 2 \
        --source-map results/kappa_map_l0p1.float32 \
        --output results/kappa_map_l0p2.float32
    python derive_kappa_maps.py smooth --map results/kappa_map_l0p2.float32 \
        --level 4arcmin
    python derive_kappa_maps.py compare \
        --map-a results/kappa_map_l0p5_derived.float32 \
        --map-b results/kappa_map_l0p5.float32
"""

from __future__ import annotations

import argparse
import json
import logging
import os

import numpy as np

from make_kappa_map import counts_to_kappa
from sim_utils import (
    FWHM_TO_SIGMA,
    KappaMapInfo,
    ensure_parent,
    open_kappa_memmap,
    read_kappa_metadata,
    setup_logging,
    write_kappa_metadata,
    SMOOTHING_FWHM_HMPC,
)


logger = logging.getLogger(__name__)

STRIP_COARSE_ROWS = 256


def block_sum(source_map: str, factor: int, output: str, counts_path: str | None) -> None:
    info = read_kappa_metadata(source_map)
    with open(f"{source_map}.json", "r", encoding="utf-8") as f:
        extra = json.load(f)
    if counts_path is None:
        counts_path = extra.get("counts_path")
    if not counts_path or not os.path.exists(counts_path):
        raise FileNotFoundError(
            f"Count map not found ({counts_path}); rebuild the source map with --keep-counts"
        )

    n_fine = info.shape[0]
    if n_fine % factor != 0:
        raise ValueError(f"factor {factor} does not divide the {n_fine}^2 source grid")
    n_coarse = n_fine // factor
    pixel_coarse = info.pixel_size_hmpc * factor
    logger.info(
        "Block-summing %d^2 counts by %dx%d -> %d^2 (%.3g h^-1 Mpc pixels)",
        n_fine, factor, factor, n_coarse, pixel_coarse,
    )

    counts_fine = np.memmap(counts_path, dtype="uint32", mode="r", shape=(n_fine, n_fine))
    coarse = np.empty((n_coarse, n_coarse), dtype=np.uint32)
    total = 0
    for row0 in range(0, n_coarse, STRIP_COARSE_ROWS):
        row1 = min(row0 + STRIP_COARSE_ROWS, n_coarse)
        strip = np.asarray(counts_fine[row0 * factor : row1 * factor]).astype(np.int64)
        summed = strip.reshape(row1 - row0, factor, n_coarse, factor).sum(axis=(1, 3))
        if summed.max() >= 2**32:
            raise OverflowError("coarse counts exceed uint32")
        coarse[row0:row1] = summed.astype(np.uint32)
        total += int(summed.sum())

    if info.total_particles is not None and total != info.total_particles:
        raise AssertionError(
            f"block-summed total {total} != source total_particles {info.total_particles}"
        )
    logger.info("Coarse counts total %d matches source", total)

    ensure_parent(output)
    counts_to_kappa(
        counts=coarse,
        output=output,
        total_particles=total,
        pixel_size_hmpc=pixel_coarse,
        box_size_hmpc=info.box_size_hmpc,
        mp_eff_hmsun=extra["mp_eff_hmsun"],
        sigma_c=extra["sigma_c_hmsun_per_mpc2"],
    )
    write_kappa_metadata(
        KappaMapInfo(
            path=os.path.abspath(output),
            shape=(n_coarse, n_coarse),
            dtype="float32",
            pixel_size_hmpc=pixel_coarse,
            box_size_hmpc=info.box_size_hmpc,
            total_particles=total,
        ),
        extra={
            "particle_file": extra.get("particle_file"),
            "mp_eff_hmsun": extra["mp_eff_hmsun"],
            "sigma_c_hmsun_per_mpc2": extra["sigma_c_hmsun_per_mpc2"],
            "derived_from_counts": os.path.abspath(counts_path),
            "block_factor": factor,
        },
    )
    logger.info("Wrote %s", output)


def smooth_periodic_chunked(
    kappa: np.ndarray,
    out: np.memmap,
    sigma_pix: float,
    strip_rows: int = 2000,
) -> None:
    """Exact periodic Gaussian smoothing in row strips (memmap-safe).

    Each strip is extended by a wrapped margin larger than the kernel's
    truncation radius, so the strip interiors equal full-map mode="wrap"
    smoothing to float precision without holding the whole map in memory.
    """
    from scipy.ndimage import gaussian_filter

    n = kappa.shape[0]
    margin = int(np.ceil(6.0 * sigma_pix))
    for row0 in range(0, n, strip_rows):
        row1 = min(row0 + strip_rows, n)
        rows = np.arange(row0 - margin, row1 + margin) % n
        ext = np.asarray(kappa[rows])
        sm = gaussian_filter(ext, sigma=sigma_pix, mode=("nearest", "wrap"))
        out[row0:row1] = sm[margin : margin + (row1 - row0)]


def smooth(map_path: str, level: str, output: str | None) -> None:
    if level not in SMOOTHING_FWHM_HMPC or level == "none":
        raise ValueError(f"level must be one of {sorted(set(SMOOTHING_FWHM_HMPC) - {'none'})}")
    kappa, info = open_kappa_memmap(map_path)
    with open(f"{map_path}.json", "r", encoding="utf-8") as f:
        extra = json.load(f)
    if extra.get("smoothing_level"):
        raise ValueError(f"{map_path} is already smoothed ({extra['smoothing_level']})")

    if output is None:
        base = str(map_path)
        stem, ext = os.path.splitext(base)
        output = f"{stem}_s{level}{ext}"

    sigma_pix = SMOOTHING_FWHM_HMPC[level] * FWHM_TO_SIGMA / info.pixel_size_hmpc
    logger.info(
        "Smoothing %s at %s (FWHM %.2f h^-1 Mpc, sigma %.2f px)",
        map_path, level, SMOOTHING_FWHM_HMPC[level], sigma_pix,
    )
    ensure_parent(output)
    if os.path.exists(output):
        os.remove(output)
    out = np.memmap(output, dtype="float32", mode="w+", shape=info.shape)
    smooth_periodic_chunked(kappa, out, sigma_pix)
    out.flush()
    del out
    write_kappa_metadata(
        KappaMapInfo(
            path=os.path.abspath(output),
            shape=info.shape,
            dtype="float32",
            pixel_size_hmpc=info.pixel_size_hmpc,
            box_size_hmpc=info.box_size_hmpc,
            total_particles=info.total_particles,
        ),
        extra={
            "smoothed_from": os.path.abspath(map_path),
            "smoothing_level": level,
            "smoothing_fwhm_hmpc": SMOOTHING_FWHM_HMPC[level],
            "note": "already smoothed; use --smooth none in stack scripts",
        },
    )
    logger.info("Wrote %s", output)


def compare(map_a: str, map_b: str) -> None:
    a, info_a = open_kappa_memmap(map_a)
    b, info_b = open_kappa_memmap(map_b)
    if info_a.shape != info_b.shape:
        raise ValueError(f"shape mismatch: {info_a.shape} vs {info_b.shape}")
    max_abs = 0.0
    n = info_a.shape[0]
    for row0 in range(0, n, 1024):
        row1 = min(row0 + 1024, n)
        diff = np.abs(a[row0:row1].astype(np.float64) - b[row0:row1].astype(np.float64))
        max_abs = max(max_abs, float(diff.max()))
    scale = float(np.abs(np.asarray(a[: min(1024, n)])).max())
    logger.info("max |a-b| = %.3e (typical |kappa| scale %.3e)", max_abs, scale)
    print(f"max_abs_diff {max_abs:.6e}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p_bs = sub.add_parser("block-sum", help="Derive a coarser kappa map from the kept count map.")
    p_bs.add_argument("--source-map", default="analysis/sim/results/kappa_map_l0p1.float32")
    p_bs.add_argument("--counts", default=None, help="Count map path (default: from source sidecar).")
    p_bs.add_argument("--factor", type=int, required=True, help="Integer block factor (2 or 5).")
    p_bs.add_argument("--output", required=True)

    p_sm = sub.add_parser("smooth", help="Precompute a smoothed copy of a kappa map.")
    p_sm.add_argument("--map", required=True)
    p_sm.add_argument("--level", required=True, choices=[k for k in SMOOTHING_FWHM_HMPC if k != "none"])
    p_sm.add_argument("--output", default=None, help="Default: <map stem>_s<level><ext>")

    p_cp = sub.add_parser("compare", help="Max |diff| between two same-shape kappa maps.")
    p_cp.add_argument("--map-a", required=True)
    p_cp.add_argument("--map-b", required=True)

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    if args.command == "block-sum":
        block_sum(args.source_map, args.factor, args.output, args.counts)
    elif args.command == "smooth":
        smooth(args.map, args.level, args.output)
    elif args.command == "compare":
        compare(args.map_a, args.map_b)


if __name__ == "__main__":
    main()
