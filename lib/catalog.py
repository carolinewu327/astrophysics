"""Data I/O for the gravitational lensing pipeline.

Catalog loading, path resolution, preprocessing, and Planck convergence
map loading.
"""

import logging
import os
from typing import Optional, Tuple, Union

import healpy as hp
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from astropy.io import fits

from constants import FWHM_ARCMIN, NSIDE
from geometry import fast_icrs_to_galactic

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
# Path resolution
# ---------------------------------------------------------------------------
def resolve_catalog_path(
    data_dir: str, dataset: str, region: str, catalog_type: str
) -> str:
    """Return the FITS file path for the requested catalog.

    Parameters
    ----------
    data_dir : str
        Root data directory.
    dataset : str
        "BOSS" or "eBOSS".
    region : str
        Region identifier (e.g. "North"/"South" for BOSS, "NGC"/"SGC" for eBOSS).
    catalog_type : str
        "galaxy" or "random".

    Returns
    -------
    str
        Path to the FITS file.
    """
    if dataset == "BOSS":
        if catalog_type == "galaxy":
            fname = f"galaxy_DR12v5_CMASS_{region}.fits.gz"
        else:
            fname = f"random0_DR12v5_CMASS_{region}.fits.gz"
        return os.path.join(data_dir, "BOSS", fname)
    elif dataset == "eBOSS":
        if catalog_type == "galaxy":
            fname = f"eBOSS_LRG_clustering_data-{region}-vDR16.fits"
        else:
            fname = f"eBOSS_LRG_clustering_random-{region}-vDR16.fits"
        return os.path.join(data_dir, "eBOSS", fname)
    else:
        raise ValueError(f"Unknown dataset: {dataset}")


def resolve_planck_paths(data_dir: str) -> Tuple[str, str]:
    """Return (alm_path, mask_path) for the Planck lensing data."""
    base = os.path.join(data_dir, "planck", "COM_Lensing_4096_R3.00")
    alm_path = os.path.join(base, "MV", "dat_klm.fits")
    mask_path = os.path.join(base, "mask.fits")
    return alm_path, mask_path


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
    """Load a catalog from a FITS file.

    Parameters
    ----------
    path : str
        Path to the FITS file.
    weights : bool or str
        ``"CMASS"`` for BOSS CMASS weighting, ``True`` for standard eBOSS
        weighting, ``False`` for uniform weights.
    random_fraction : float or None
        If set, randomly subsample to this fraction of the catalog.
    z_min, z_max : float
        Redshift cuts.

    Returns
    -------
    cat, w : tuple of ndarray
        Filtered catalog and per-object weights.
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


def load_catalog_lightweight(
    path: str,
    columns: Tuple[str, ...] = ("RA", "DEC", "Z"),
    fraction: Optional[float] = None,
    z_min: float = 0,
    z_max: float = 10000,
) -> np.ndarray:
    """Load only selected columns from a FITS catalog, memory-efficiently.

    Reads the full row count but only the requested columns, then applies
    redshift cuts and optional subsampling.  Much lighter than load_catalog
    for large random catalogues where only positions are needed.

    Parameters
    ----------
    path : str
        Path to the FITS file.
    columns : tuple of str
        Column names to read.
    fraction : float or None
        If set, randomly subsample to this fraction *before* reading data.
    z_min, z_max : float
        Redshift cuts (applied only if "Z" is in *columns*).

    Returns
    -------
    np.ndarray
        Structured array with the requested columns.
    """
    with fits.open(path, memmap=True) as hd:
        nrows = hd[1].header["NAXIS2"]

        # Subsample row indices first to avoid reading all data
        if fraction and fraction < 1.0:
            n_keep = int(fraction * nrows)
            idx = np.sort(np.random.choice(nrows, n_keep, replace=False))
        else:
            idx = np.arange(nrows)

        # Read only the requested columns for the selected rows
        col_data = {col: hd[1].data.field(col)[idx] for col in columns}

    # Build a structured array
    dtype = [(col, col_data[col].dtype) for col in columns]
    cat = np.empty(len(idx), dtype=dtype)
    for col in columns:
        cat[col] = col_data[col]

    # Redshift + finite-position cuts
    valid = np.ones(len(cat), dtype=bool)
    if "Z" in columns:
        valid &= (cat["Z"] > z_min) & (cat["Z"] < z_max)
    if "RA" in columns:
        valid &= np.isfinite(cat["RA"])
    if "DEC" in columns:
        valid &= np.isfinite(cat["DEC"])
    cat = cat[valid]

    logger.info("Loaded %d rows (%d columns) from %s", len(cat), len(columns), path)
    return cat


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
