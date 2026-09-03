"""Shared utilities for the BigMDPL simulation kappa pipeline."""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from geometry import assert_unit_pixel_scale, stack_axis


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

# FWHM equivalents (D_l x theta at z=0.547) of the observational smoothing.
# The Planck map is smoothed with hp.smoothalm(fwhm=...), so these are FWHM,
# not Gaussian sigma; convert by FWHM_TO_SIGMA before use as a filter sigma.
SMOOTHING_FWHM_HMPC = {
    "none": 0.0,
    "2arcmin": 0.83,
    "4arcmin": 1.66,
    "8arcmin": 3.31,
}
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


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
    vals, _ = stack_axis(box_size_hmpc, grid_size)
    x_grid, y_grid = np.meshgrid(vals, vals)
    return vals, x_grid, y_grid


def pair_stack_offsets(
    box_size_hmpc: float = OBS_BOX_SIZE_HMPC,
    grid_size: int = OBS_PAIR_GRID_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vals, _ = stack_axis(box_size_hmpc, grid_size)
    x_grid, y_grid = np.meshgrid(vals, vals)
    return vals, x_grid, y_grid


def symmetrized_radial_interpolator(single_map: np.ndarray, validate: bool = True):
    """1D radius -> kappa function from a radially symmetrized single stack.

    The symmetrized map is a function of distance from its symmetrization
    physical center ((N-1)/2 in pixel-index coordinates), so a radial profile
    represents it exactly — and
    extends to the map corners (r ~ 70 px), unlike 2D interpolation, which is
    limited to the -N/2..N/2-1 axis range of the even grid and picks up
    center/edge asymmetries. Radii beyond the corners return 0.
    """
    if validate:
        validate_radially_symmetrized_map(single_map)
    # Radius below is measured in PIXELS and the returned function is then
    # evaluated at physical h^-1 Mpc coordinates, so this is only valid when
    # one pixel is one h^-1 Mpc.  Assert it rather than assume it.
    assert_unit_pixel_scale(single_map.shape[0])
    n = single_map.shape[0]
    center = 0.5 * (n - 1)
    yy, xx = np.indices((n, n))
    r = np.hypot(xx - center, yy - center).ravel()
    v = single_map.ravel()
    bins = np.arange(0.0, r.max() + 1.0, 0.5)
    idx = np.clip(np.digitize(r, bins) - 1, 0, len(bins) - 2)
    sums = np.bincount(idx, weights=v, minlength=len(bins) - 1)
    counts = np.bincount(idx, minlength=len(bins) - 1)
    good = counts > 0
    centers = 0.5 * (bins[:-1] + bins[1:])[good]
    means = sums[good] / counts[good]

    def profile(radius: np.ndarray) -> np.ndarray:
        return np.interp(radius, centers, means, left=means[0], right=0.0)

    return profile


def make_two_halo_template(
    single_map: np.ndarray,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    rperp_center: float,
    validate: bool = True,
) -> np.ndarray:
    """Two copies of the symmetrized single-stack profile at +/- rperp/2.

    Grids are in the single map's physical units (1 h^-1 Mpc per pixel).
    """
    profile = symmetrized_radial_interpolator(single_map, validate=validate)
    offset = 0.5 * rperp_center
    return profile(np.hypot(x_grid + offset, y_grid)) + profile(np.hypot(x_grid - offset, y_grid))


def reflect_symmetrize_map(kappa_map: np.ndarray) -> np.ndarray:
    return 0.25 * (
        kappa_map
        + np.flip(kappa_map, axis=0)
        + np.flip(kappa_map, axis=1)
        + np.flip(np.flip(kappa_map, axis=0), axis=1)
    )


def radial_symmetrize_map(kappa_map: np.ndarray, pwr: float = 2 / 3) -> np.ndarray:
    """Average a square map in radial bins about its physical center.

    For an even grid the origin lies between the four central pixels, at
    ``(N - 1) / 2`` in pixel-index coordinates.  For an odd grid the same
    expression selects the central pixel.  Using ``N // 2`` on an even grid
    introduces a half-pixel shift and makes a supposedly radial map asymmetric
    under a physical reflection.
    """

    if kappa_map.ndim != 2 or kappa_map.shape[0] != kappa_map.shape[1]:
        raise ValueError(f"Expected a square 2D map, got shape {kappa_map.shape}.")
    if not np.isfinite(kappa_map).all():
        raise ValueError("Cannot symmetrize a map containing non-finite values.")
    y, x = np.indices(kappa_map.shape)
    cy = 0.5 * (kappa_map.shape[0] - 1)
    cx = 0.5 * (kappa_map.shape[1] - 1)
    radial_bin = (((x - cx) ** 2 + (y - cy) ** 2) ** pwr).astype(np.int64)
    flat_bin = radial_bin.ravel()
    flat_val = kappa_map.ravel()
    sums = np.bincount(flat_bin, weights=flat_val)
    counts = np.bincount(flat_bin)
    radial_avg = np.zeros_like(sums, dtype=np.float64)
    np.divide(sums, counts, out=radial_avg, where=counts > 0)
    return radial_avg[radial_bin].astype(kappa_map.dtype, copy=False)


def validate_radially_symmetrized_map(single_map: np.ndarray) -> None:
    """Reject a single stack that is not centered on the physical map origin.

    A radially symmetrized map is symmetric under both flips, because index
    ``i -> N-1-i`` preserves ``|i - (N-1)/2|``.  Archived even-grid stacks built
    with the old ``N // 2`` convention are symmetric about pixel 50 instead, and
    re-symmetrizing them cannot recover values already averaged into the wrong
    radial bins -- hence a hard failure rather than a silent repair.

    The message names both possible causes -- never symmetrized, or symmetrized
    about the wrong center -- and reports the measured asymmetry, but does not
    guess between them.  Magnitude cannot separate the cases: the archived
    mis-centered simulation single is 20.5% asymmetric and the mis-centered
    product in the regression test is 18.1%, both squarely in the range a raw
    stack occupies.  An earlier draft of this check branched on a 10% threshold
    and would have mislabeled every real instance.  A caller passing a raw map
    on purpose should use ``validate=False`` rather than read the message.

    ``lib/geometry.py`` carries a byte-identical copy: the simulation tree is
    run with ``PYTHONPATH=analysis/sim`` and cannot import ``lib``.
    ``tests/test_sim_centering.py`` asserts the two agree; change both together.
    """
    if single_map.ndim != 2 or single_map.shape[0] != single_map.shape[1]:
        raise ValueError(f"Expected a square 2D single stack, got {single_map.shape}.")
    if not np.isfinite(single_map).all():
        raise ValueError("Single stack contains non-finite values.")
    scale = max(float(np.max(np.abs(single_map))), 1.0e-12)
    atol = max(1.0e-12, 1.0e-8 * scale)
    residual = max(
        float(np.max(np.abs(single_map - np.flip(single_map, axis=0)))),
        float(np.max(np.abs(single_map - np.flip(single_map, axis=1)))),
    )
    if not (
        np.allclose(single_map, np.flip(single_map, axis=0), rtol=1.0e-7, atol=atol)
        and np.allclose(single_map, np.flip(single_map, axis=1), rtol=1.0e-7, atol=atol)
    ):
        raise ValueError(
            f"Single stack is not reflection-symmetric about its physical "
            f"center (asymmetry {residual / scale:.1%} of peak). Either it was "
            "never radially symmetrized, or it is an archived even-grid product "
            "built with the old N//2 centering -- re-symmetrizing that cannot "
            "recover the lost radial bins, so regenerate it with "
            "stack_single_sim.py. Magnitude does not separate the two cases "
            "(both run ~20% of peak); check how the file was produced. If you "
            "are passing a raw map on purpose, call with validate=False."
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
    if smoothing not in SMOOTHING_FWHM_HMPC:
        raise ValueError(f"Unknown smoothing option {smoothing!r}")

    try:
        from scipy.ndimage import gaussian_filter
    except ImportError as exc:
        raise RuntimeError("scipy is required for smoothing") from exc

    sigma_pix = SMOOTHING_FWHM_HMPC[smoothing] * FWHM_TO_SIGMA / pixel_size_hmpc
    return gaussian_filter(kappa_map, sigma=sigma_pix, mode="wrap")


def output_name(prefix: str, mass_label: str, suffix: str = ".csv") -> str:
    return f"{prefix}_{mass_label}{suffix}"
