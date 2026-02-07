"""Pure math/analytics for the gravitational lensing pipeline.

Coordinate transforms, angular separation, map symmetrization, and radial
profiles.  No file I/O -- every function is a pure computation suitable for
unit testing.
"""

from typing import Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Coordinate transforms
# ---------------------------------------------------------------------------
def fast_icrs_to_galactic(
    ra_deg: np.ndarray, dec_deg: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert ICRS coordinates (RA, Dec in degrees) to Galactic (l, b in degrees).

    Uses the standard rotation matrix -- no astropy dependency.

    Parameters
    ----------
    ra_deg : float or array-like
        Right Ascension in degrees.
    dec_deg : float or array-like
        Declination in degrees.

    Returns
    -------
    l_deg, b_deg : tuple of arrays
        Galactic longitude and latitude in degrees.
    """
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    # Rotation matrix from ICRS to Galactic
    R = np.array(
        [
            [-0.0548755604, -0.8734370902, -0.4838350155],
            [0.4941094279, -0.4448296300, 0.7469822445],
            [-0.8676661490, -0.1980763734, 0.4559837762],
        ]
    )

    xg, yg, zg = np.dot(R, np.array([x, y, z]))

    b_rad = np.arcsin(zg)
    l_rad = np.arctan2(yg, xg)

    l_deg = np.degrees(l_rad) % 360
    b_deg = np.degrees(b_rad)

    return l_deg, b_deg


def angular_separation(
    l1_deg: np.ndarray,
    b1_deg: np.ndarray,
    l2_deg: np.ndarray,
    b2_deg: np.ndarray,
) -> np.ndarray:
    """Great-circle angular separation between points in Galactic coordinates.

    Parameters
    ----------
    l1_deg, b1_deg : float or array-like
        Galactic longitude and latitude of the first point(s), in degrees.
    l2_deg, b2_deg : float or array-like
        Galactic longitude and latitude of the second point(s), in degrees.

    Returns
    -------
    theta : float or ndarray
        Angular separation in **radians**.
    """
    l1 = np.radians(l1_deg)
    b1 = np.radians(b1_deg)
    l2 = np.radians(l2_deg)
    b2 = np.radians(b2_deg)

    cos_theta = (
        np.cos(b1) * np.cos(l1) * np.cos(b2) * np.cos(l2)
        + np.cos(b1) * np.sin(l1) * np.cos(b2) * np.sin(l2)
        + np.sin(b1) * np.sin(b2)
    )
    return np.arccos(np.clip(cos_theta, -1.0, 1.0))


# ---------------------------------------------------------------------------
# Map symmetrization
# ---------------------------------------------------------------------------
def symmetrize_map(kappa_map: np.ndarray, pwr: float = 2 / 3) -> np.ndarray:
    """Symmetrize a 2D map by averaging in radial bins.

    Parameters
    ----------
    kappa_map : ndarray
        2-D map to symmetrize.
    pwr : float
        Controls radial binning.  1/2 = digitize in radius, 1 = radius
        squared, 2/3 (default) = radius^(4/3).
    """
    grid_size = kappa_map.shape[0]
    y, x = np.indices(kappa_map.shape)
    cx, cy = grid_size // 2, grid_size // 2
    r = (((x - cx) ** 2 + (y - cy) ** 2) ** pwr).astype(int)
    r_flat = r.ravel()
    kappa_flat = kappa_map.ravel()
    kappa_avg = np.bincount(r_flat, weights=kappa_flat) / np.bincount(r_flat)
    sym_map = kappa_avg[r]
    return sym_map


def reflect_symmetrize_map(kappa_map: np.ndarray) -> np.ndarray:
    """Apply reflection symmetry: average (+x,+y), (-x,+y), (+x,-y), (-x,-y).

    Preserves asymmetry along the pair axis (X) while enhancing the signal
    perpendicular to it.  Assumes *kappa_map* is square with odd dimensions
    so the center pixel is well-defined.
    """
    sym_map = np.copy(kappa_map)
    n = kappa_map.shape[0]
    center = n // 2

    for i in range(center + 1):
        for j in range(center + 1):
            di, dj = i - center, j - center

            coords = [
                (center + di, center + dj),
                (center - di, center + dj),
                (center + di, center - dj),
                (center - di, center - dj),
            ]

            vals = [kappa_map[x, y] for x, y in coords]
            avg_val = np.mean(vals)

            for x, y in coords:
                sym_map[x, y] = avg_val

    return sym_map


# ---------------------------------------------------------------------------
# Radial profile
# ---------------------------------------------------------------------------
def radial_profile(
    arr: np.ndarray,
    sigma: Optional[np.ndarray] = None,
    zoom: int = 70,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Extract a 1D radial profile from a 2D map.

    Returns (profile, errors).  *errors* is ``None`` when *sigma* is not
    provided.
    """
    grid_size = arr.shape[0]
    y, x = np.indices(arr.shape)
    r = np.sqrt((x - grid_size / 2) ** 2 + (y - grid_size / 2) ** 2).astype(int)
    flat = arr.ravel()
    N = np.bincount(r.ravel())
    S = np.bincount(r.ravel(), weights=flat)
    prof = S / N

    err = None
    if sigma is not None:
        flat_s = sigma.ravel()
        S2 = np.bincount(r.ravel(), weights=flat_s**2)
        err = np.sqrt(S2) / N

    # Trim to requested zoom range
    prof = prof[:zoom]
    if err is not None:
        err = err[:zoom]

    return prof, err
