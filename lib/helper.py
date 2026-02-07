"""Shared utility code for the astrophysics gravitational lensing pipeline."""

import logging
from typing import Optional, Tuple, Union

import healpy as hp
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from astropy.io import fits

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GRID_SIZE: int = 100
BOX_SIZE_HMPC: float = 100.0
FWHM_ARCMIN: float = 8.0
NSIDE: int = 2048

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
def setup_logging(level: int = logging.INFO) -> None:
    """Configure root logger with timestamped format."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# ---------------------------------------------------------------------------
# Catalog I/O
# ---------------------------------------------------------------------------
def load_catalog(
    path: str,
    weights: Union[bool, str] = True,
    random_fraction: Optional[float] = None,
    z_min: float = 0,
    z_max: float = 10000,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load a catalog from a FITS file.
    """
    with fits.open(path) as hd:
        cat = hd[1].data
    cat = cat[
        (cat["Z"] > z_min)
        & (cat["Z"] < z_max)
        & np.isfinite(cat["RA"])
        & np.isfinite(cat["DEC"])
    ]
    if random_fraction:
        cat = cat[
            np.random.choice(len(cat), int(random_fraction * len(cat)), replace=False)
        ]

    if weights == "CMASS":
        w = (
            cat["WEIGHT_SEEING"]
            * cat["WEIGHT_STAR"]
            * (cat["WEIGHT_NOZ"] + cat["WEIGHT_CP"] - 1)
        )
    elif weights:
        w = cat["WEIGHT_NOZ"] * cat["WEIGHT_SYSTOT"]
    else:
        w = np.ones(len(cat))

    logger.info("Loaded %d galaxies from %s", len(cat), path)
    return cat, w


# ---------------------------------------------------------------------------
# Coordinate transforms
# ---------------------------------------------------------------------------
def fast_icrs_to_galactic(
    ra_deg: np.ndarray, dec_deg: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert ICRS coordinates (RA, Dec in degrees) to Galactic (l, b in degrees)
    without using astropy. Based on standard transformation matrix.

    Parameters:
        ra_deg : float or array-like
            Right Ascension in degrees
        dec_deg : float or array-like
            Declination in degrees

    Returns:
        l_deg, b_deg : tuple of arrays
            Galactic longitude and latitude in degrees
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


# ---------------------------------------------------------------------------
# Catalog preprocessing
# ---------------------------------------------------------------------------
def preprocess_catalog_galactic(
    data: np.ndarray, weights: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Convert a catalog with RA, DEC, Z to Galactic coords + comoving distance."""
    z = data["Z"]
    ra = data["RA"]
    dec = data["DEC"]
    D_mpc = cosmo.comoving_distance(z).value  # Mpc
    D_hmpc = D_mpc * cosmo.h  # h^-1 Mpc
    valid = (D_hmpc > 0) & np.isfinite(D_hmpc)

    ra_valid = ra[valid]
    dec_valid = dec[valid]
    D_valid = D_hmpc[valid]

    if weights is not None:
        assert len(weights) == len(data), "Weights must match the length of the data."
        weights_valid = weights[valid]
    else:
        weights_valid = None

    l, b = fast_icrs_to_galactic(ra_valid, dec_valid)

    logger.info(
        "Preprocessed catalog: %d / %d galaxies valid", np.sum(valid), len(data)
    )
    return l, b, D_valid, data[valid], weights_valid


# ---------------------------------------------------------------------------
# Planck kappa map loading
# ---------------------------------------------------------------------------
def load_kappa_map(
    alm_file: str = "data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits",
    mask_file: str = "data/planck/COM_Lensing_4096_R3.00/mask.fits",
    nside: int = NSIDE,
    fwhm_arcmin: float = FWHM_ARCMIN,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """Load Planck ALM data, smooth, and return (kmap, mask, nside).

    Reads the complex ALM from REAL+IMAG columns, applies Gaussian smoothing
    at *fwhm_arcmin* resolution, and synthesises the HEALPix convergence map.
    """
    logger.info("Loading ALM from %s", alm_file)
    alm_data = fits.open(alm_file)[1].data
    alm = alm_data["REAL"] + 1j * alm_data["IMAG"]
    lmax = hp.Alm.getlmax(len(alm))

    fwhm_rad = np.radians(fwhm_arcmin / 60.0)
    kmap = hp.alm2map(hp.smoothalm(alm, fwhm=fwhm_rad), nside=nside, lmax=lmax)

    logger.info("Loading mask from %s", mask_file)
    mask = hp.read_map(mask_file)

    logger.info(
        "Kappa map ready: nside=%d, lmax=%d, fwhm=%.1f arcmin", nside, lmax, fwhm_arcmin
    )
    return kmap, mask, nside


# ---------------------------------------------------------------------------
# Map symmetrization
# ---------------------------------------------------------------------------
def symmetrize_map(kappa_map: np.ndarray, pwr: float = 2 / 3) -> np.ndarray:
    """
    Symmetrize a 2D map by averaging in radial bins.
    :param kappa_map:
    :param pwr: 1/2 means digitize in radius, 1 means digitize in radius squared
    2/3 means digitize in radius^(4/3)
    :return:
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
    provided.  Does **not** call ``plt.show()`` -- callers handle plotting.
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
