"""Shared utilities for the BigMDPL simulation kappa pipeline."""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


BIGMDPL_BOX_SIZE_HMPC = 2500.0
BIGMDPL_PARTICLE_MASS_HMSUN = 2.359e10
BIGMDPL_DOWNSAMPLE_FACTOR = 200.0
SIGMA_C_HMSUN_PER_MPC2 = 8.88e14
BIGMDPL_OMEGA_M = 0.307
BIGMDPL_H = 0.677
BIGMDPL_Z_SNAPSHOT = 0.547

OBS_BOX_SIZE_HMPC = 100.0
OBS_SINGLE_GRID_SIZE = 100
OBS_PAIR_GRID_SIZE = 101

SMOOTHING_HMPC = {
    "none": 0.0,
    "2arcmin": 0.83,
    "4arcmin": 1.66,
    "8arcmin": 3.31,
}


@dataclass(frozen=True)
class KappaMapInfo:
    """Metadata needed to reopen a raw memmapped kappa map."""

    path: Path
    shape: tuple[int, int]
    dtype: str
    pixel_size_hmpc: float
    box_size_hmpc: float
    total_particles: int | None = None


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def ensure_parent(path: str | os.PathLike[str]) -> None:
    parent = Path(path).expanduser().resolve().parent
    parent.mkdir(parents=True, exist_ok=True)


def metadata_path(map_path: str | os.PathLike[str]) -> Path:
    return Path(f"{map_path}.json")


def write_kappa_metadata(info: KappaMapInfo, extra: dict | None = None) -> None:
    data = {
        "path": str(info.path),
        "shape": list(info.shape),
        "dtype": info.dtype,
        "pixel_size_hmpc": info.pixel_size_hmpc,
        "box_size_hmpc": info.box_size_hmpc,
        "total_particles": info.total_particles,
    }
    if extra:
        data.update(extra)
    meta_path = metadata_path(info.path)
    ensure_parent(meta_path)
    with meta_path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def read_kappa_metadata(map_path: str | os.PathLike[str]) -> KappaMapInfo:
    meta_path = metadata_path(map_path)
    with meta_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    return KappaMapInfo(
        path=Path(map_path),
        shape=tuple(data["shape"]),
        dtype=data.get("dtype", "float32"),
        pixel_size_hmpc=float(data["pixel_size_hmpc"]),
        box_size_hmpc=float(data.get("box_size_hmpc", BIGMDPL_BOX_SIZE_HMPC)),
        total_particles=data.get("total_particles"),
    )


def open_kappa_memmap(
    map_path: str | os.PathLike[str],
    mode: str = "r",
    shape: tuple[int, int] | None = None,
    dtype: str | np.dtype | None = None,
) -> tuple[np.memmap, KappaMapInfo]:
    """Open a raw kappa memmap, reading the sidecar metadata by default."""

    if shape is None or dtype is None:
        info = read_kappa_metadata(map_path)
        shape = info.shape if shape is None else shape
        dtype = info.dtype if dtype is None else dtype
    else:
        info = KappaMapInfo(
            path=Path(map_path),
            shape=shape,
            dtype=str(np.dtype(dtype)),
            pixel_size_hmpc=BIGMDPL_BOX_SIZE_HMPC / shape[0],
            box_size_hmpc=BIGMDPL_BOX_SIZE_HMPC,
        )

    mmap = np.memmap(map_path, dtype=dtype, mode=mode, shape=shape)
    return mmap, info


def minimal_periodic_delta(a: np.ndarray | float, b: np.ndarray | float, box_size: float) -> np.ndarray:
    """Return the minimum-image displacement from ``a`` to ``b``."""

    return (np.asarray(b) - np.asarray(a) + 0.5 * box_size) % box_size - 0.5 * box_size


def periodic_bilinear_sample(
    kappa_map: np.ndarray,
    x_hmpc: np.ndarray,
    y_hmpc: np.ndarray,
    pixel_size_hmpc: float,
    box_size_hmpc: float,
) -> np.ndarray:
    """Sample a periodic cell-centered map at arbitrary x/y positions.

    The map is indexed as ``map[y_index, x_index]``. Pixel centers are at
    ``(i + 0.5) * pixel_size_hmpc``.
    """

    x = np.mod(x_hmpc, box_size_hmpc)
    y = np.mod(y_hmpc, box_size_hmpc)

    u = x / pixel_size_hmpc - 0.5
    v = y / pixel_size_hmpc - 0.5

    nx = kappa_map.shape[1]
    ny = kappa_map.shape[0]

    i0 = np.floor(u).astype(np.int64) % nx
    j0 = np.floor(v).astype(np.int64) % ny
    i1 = (i0 + 1) % nx
    j1 = (j0 + 1) % ny

    tx = u - np.floor(u)
    ty = v - np.floor(v)

    f00 = kappa_map[j0, i0]
    f10 = kappa_map[j0, i1]
    f01 = kappa_map[j1, i0]
    f11 = kappa_map[j1, i1]

    return (
        (1.0 - tx) * (1.0 - ty) * f00
        + tx * (1.0 - ty) * f10
        + (1.0 - tx) * ty * f01
        + tx * ty * f11
    )


def single_stack_offsets(
    box_size_hmpc: float = OBS_BOX_SIZE_HMPC,
    grid_size: int = OBS_SINGLE_GRID_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cell = box_size_hmpc / grid_size
    vals = np.linspace(
        -0.5 * box_size_hmpc + 0.5 * cell,
        0.5 * box_size_hmpc - 0.5 * cell,
        grid_size,
    )
    x_grid, y_grid = np.meshgrid(vals, vals)
    return vals, x_grid, y_grid


def pair_stack_offsets(
    box_size_hmpc: float = OBS_BOX_SIZE_HMPC,
    grid_size: int = OBS_PAIR_GRID_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vals = np.linspace(-0.5 * box_size_hmpc, 0.5 * box_size_hmpc, grid_size)
    x_grid, y_grid = np.meshgrid(vals, vals)
    return vals, x_grid, y_grid


def reflect_symmetrize_map(kappa_map: np.ndarray) -> np.ndarray:
    return 0.25 * (
        kappa_map
        + np.flip(kappa_map, axis=0)
        + np.flip(kappa_map, axis=1)
        + np.flip(np.flip(kappa_map, axis=0), axis=1)
    )


def radial_profile_from_map(
    kappa_map: np.ndarray,
    box_size_hmpc: float = OBS_BOX_SIZE_HMPC,
) -> tuple[np.ndarray, np.ndarray]:
    grid_size = kappa_map.shape[0]
    cell = box_size_hmpc / grid_size
    y, x = np.indices(kappa_map.shape)
    cx = 0.5 * (grid_size - 1)
    cy = 0.5 * (grid_size - 1)
    radius_pix = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    annulus = np.floor(radius_pix).astype(np.int64)
    flat_bin = annulus.ravel()
    flat_val = kappa_map.ravel()
    counts = np.bincount(flat_bin)
    sums = np.bincount(flat_bin, weights=flat_val)
    profile = sums / np.maximum(counts, 1)
    radius = (np.arange(len(profile)) + 0.5) * cell
    return radius, profile


def save_map_csv(path: str | os.PathLike[str], arr: np.ndarray) -> None:
    import pandas as pd

    ensure_parent(path)
    pd.DataFrame(arr).to_csv(path, index=True)


def save_profile_csv(
    path: str | os.PathLike[str],
    radius_hmpc: np.ndarray,
    kappa: np.ndarray,
    extra_columns: dict[str, Iterable[float]] | None = None,
) -> None:
    import pandas as pd

    data = {
        "radius_hmpc": radius_hmpc,
        "kappa": kappa,
    }
    if extra_columns:
        data.update(extra_columns)
    ensure_parent(path)
    pd.DataFrame(data).to_csv(path, index=False)


def hubble_bigmdpl_hmpc(z: float, omega_m: float = BIGMDPL_OMEGA_M) -> float:
    """Return H_RSD(z) in km s^-1 (h^-1 Mpc)^-1.

    This uses the common simulation-coordinate shortcut H0 = 100, not
    100 h, so (1 + z) v / H_RSD(z) is directly in h^-1 Mpc.
    """

    return 100.0 * np.sqrt(omega_m * (1.0 + z) ** 3 + 1.0 - omega_m)


def rsd_displacement_hmpc(
    vz_km_s: np.ndarray,
    z_snapshot: float = BIGMDPL_Z_SNAPSHOT,
    omega_m: float = BIGMDPL_OMEGA_M,
) -> np.ndarray:
    """Convert line-of-sight velocity to RSD displacement in h^-1 Mpc."""

    return (1.0 + z_snapshot) * vz_km_s / hubble_bigmdpl_hmpc(z_snapshot, omega_m)


def apply_periodic_gaussian_smoothing(
    kappa_map: np.ndarray,
    pixel_size_hmpc: float,
    smoothing: str,
) -> np.ndarray:
    if smoothing == "none":
        return kappa_map
    if smoothing not in SMOOTHING_HMPC:
        raise ValueError(f"Unknown smoothing option {smoothing!r}")

    try:
        from scipy.ndimage import gaussian_filter
    except ImportError as exc:
        raise RuntimeError("scipy is required for smoothing") from exc

    sigma_pix = SMOOTHING_HMPC[smoothing] / pixel_size_hmpc
    return gaussian_filter(kappa_map, sigma=sigma_pix, mode="wrap")


def output_name(prefix: str, mass_label: str, suffix: str = ".csv") -> str:
    return f"{prefix}_{mass_label}{suffix}"
