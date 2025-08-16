import numpy as np
from astropy.cosmology import Planck18 as cosmo
from astropy.io import fits

# h = cosmo.h
grid_size = 100


# --- Catalog loader ---
def load_catalog(path, weights=True, random_fraction=None):
    with fits.open(path) as hd:
        cat = hd[1].data
    cat = cat[(cat['Z'] > 0) & np.isfinite(cat['RA']) & np.isfinite(cat['DEC'])]
    if random_fraction:
        cat = cat[np.random.choice(len(cat), int(random_fraction * len(cat)), replace=False)]
    w = np.ones(len(cat)) if not weights else cat['WEIGHT_NOZ'] * cat['WEIGHT_SYSTOT']
    return cat, w


def fast_icrs_to_galactic(ra_deg, dec_deg):
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
    # Convert to radians
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)

    # Convert spherical to Cartesian
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    # Rotation matrix from ICRS to Galactic
    R = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [ 0.4941094279, -0.4448296300,  0.7469822445],
        [-0.8676661490, -0.1980763734,  0.4559837762]
    ])

    # Apply rotation
    xg, yg, zg = np.dot(R, np.array([x, y, z]))

    # Convert back to spherical
    b_rad = np.arcsin(zg)
    l_rad = np.arctan2(yg, xg)

    # Convert to degrees
    l_deg = np.degrees(l_rad) % 360
    b_deg = np.degrees(b_rad)

    return l_deg, b_deg


# Preprocessing function for a catalog with RA, DEC, and redshift
def preprocess_catalog_galactic(data):
    z = data['Z']
    ra = data['RA']
    dec = data['DEC']
    D_mpc = cosmo.comoving_distance(z).value    # Mpc
    D_hmpc = D_mpc * cosmo.h                    # convert to h⁻¹ Mpc
    valid = (D_hmpc > 0) & np.isfinite(D_hmpc)

    ra_valid = ra[valid]
    dec_valid = dec[valid]
    D_valid = D_hmpc[valid]

    l, b = fast_icrs_to_galactic(ra_valid, dec_valid)

    return l, b, D_valid, data[valid]


# --- Symmetrize map ---
def symmetrize_map(kappa_map):
    y, x = np.indices(kappa_map.shape)
    cx, cy = grid_size // 2, grid_size // 2
    r = np.sqrt((x - cx)**2 + (y - cy)**2).astype(int)
    r_flat = r.ravel()
    kappa_flat = kappa_map.ravel()
    kappa_avg = np.bincount(r_flat, weights=kappa_flat) / np.bincount(r_flat)
    sym_map = kappa_avg[r].reshape(kappa_map.shape)
    return sym_map