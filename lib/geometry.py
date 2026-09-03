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


def fast_galactic_to_icrs(
    l_deg: np.ndarray, b_deg: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert Galactic (l, b in degrees) back to ICRS (RA, Dec in degrees).

    The exact inverse of :func:`fast_icrs_to_galactic` -- the rotation matrix
    is orthogonal, so the inverse is its transpose.

    Needed because pair catalogs are stored in Galactic coordinates while
    jackknife regions are defined on equatorial positions.  Assigning pairs to
    the *same* tessellation as the single-galaxy stack is what makes a
    filament-level jackknife coherent: the filament map subtracts a control
    built from the single stack, so both terms must have the same patch of sky
    removed in each leave-one-out realization.
    """
    l = np.radians(l_deg)
    b = np.radians(b_deg)

    xg = np.cos(b) * np.cos(l)
    yg = np.cos(b) * np.sin(l)
    zg = np.sin(b)

    R = np.array(
        [
            [-0.0548755604, -0.8734370902, -0.4838350155],
            [0.4941094279, -0.4448296300, 0.7469822445],
            [-0.8676661490, -0.1980763734, 0.4559837762],
        ]
    )

    x, y, z = np.dot(R.T, np.array([xg, yg, zg]))

    dec_deg = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    ra_deg = np.degrees(np.arctan2(y, x)) % 360

    return ra_deg, dec_deg


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
# Stacking grid geometry
# ---------------------------------------------------------------------------
def stack_axis(box_size_hmpc: float, grid_size: int) -> Tuple[np.ndarray, float]:
    """Physical coordinates of a square stacking grid, and its pixel spacing.

    Two conventions, picked by parity, so that both give a spacing of exactly
    1 h^-1 Mpc on the 100 h^-1 Mpc box this pipeline uses:

        even N   cell-centred: -L/2 + c/2 ... +L/2 - c/2 with c = L/N.
                 N = 100 -> -49.5 ... +49.5, spacing 1.0.  No pixel sits on
                 zero, so the origin falls between the four central pixels.
        odd  N   endpoint-inclusive: -L/2 ... +L/2.
                 N = 101 -> -50 ... +50, spacing 1.0.  The centre pixel sits
                 exactly on the origin.

    Applying the even rule to an odd grid is what went wrong before: it gives
    c = 100/101 = 0.990, so a 101-pixel single stack came out on a scale 1%
    smaller than the 101-pixel pair stack it was subtracted from.  The two
    conventions were duplicated across sim_utils, stack_single_jk and
    combine_jackknife, which is why the mismatch survived.  Everything that
    needs a stacking axis must come through here.

    Several downstream routines -- notably
    :func:`symmetrized_radial_interpolator`, which measures radius in pixel
    units and is then handed coordinates in h^-1 Mpc -- assume the spacing is
    exactly 1.  :func:`assert_unit_pixel_scale` makes that check explicit.
    """
    if grid_size < 2:
        raise ValueError(f"grid_size must be at least 2, got {grid_size}")
    half = 0.5 * float(box_size_hmpc)
    if grid_size % 2 == 1:
        axis = np.linspace(-half, half, grid_size)
    else:
        cell = float(box_size_hmpc) / grid_size
        axis = np.linspace(-half + 0.5 * cell, half - 0.5 * cell, grid_size)
    return axis, float(axis[1] - axis[0])


def assert_unit_pixel_scale(grid_size: int, box_size_hmpc: float = 100.0,
                            tol: float = 1e-9) -> float:
    """Fail loudly when a grid's pixel is not 1 h^-1 Mpc.

    Radial interpolation over a stacked map indexes pixels and then evaluates
    the result at physical coordinates, which is only valid at unit spacing.
    """
    _, spacing = stack_axis(box_size_hmpc, grid_size)
    if abs(spacing - 1.0) > tol:
        raise ValueError(
            f"grid {grid_size} on a {box_size_hmpc:g} h^-1 Mpc box has pixel "
            f"spacing {spacing:.6f}, not 1.0; radial interpolation over this "
            "map would be off by that factor.")
    return spacing


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
    if kappa_map.ndim != 2:
        raise ValueError("kappa_map must be two-dimensional.")
    y, x = np.indices(kappa_map.shape)
    # The standard single stack has an even 100x100 grid whose origin lies
    # between the four central pixels.  Using N//2 here shifts the profile by
    # half a pixel in each direction and consequently mis-centres both copies
    # of the superposed-singles control.
    cx = 0.5 * (kappa_map.shape[1] - 1)
    cy = 0.5 * (kappa_map.shape[0] - 1)
    r = (((x - cx) ** 2 + (y - cy) ** 2) ** pwr).astype(int)
    r_flat = r.ravel()
    kappa_flat = kappa_map.ravel()
    sums = np.bincount(r_flat, weights=kappa_flat)
    counts = np.bincount(r_flat)
    kappa_avg = np.zeros_like(sums, dtype=np.float64)
    np.divide(sums, counts, out=kappa_avg, where=counts > 0)
    sym_map = kappa_avg[r]
    return sym_map


def reflect_symmetrize_map(kappa_map: np.ndarray) -> np.ndarray:
    """Apply reflection symmetry: average (+x,+y), (-x,+y), (+x,-y), (-x,-y).

    Averaging all four quadrants imposes symmetry in *both* axes, so it does
    NOT preserve asymmetry along the pair axis -- any difference between the
    two galaxies of a pair is averaged away.  That is intended: the pairs are
    ordered arbitrarily, so no physical asymmetry along X survives stacking
    anyway, and folding gains a factor of 2 in each direction.  The
    consequence for downstream code is that columns at +X and -X hold the
    same numbers, so a fit or profile must use one half only -- using both is
    exact double counting with no added information.

    Assumes *kappa_map* is square with odd dimensions so the center pixel is
    well-defined.
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
# Two-halo control template
# ---------------------------------------------------------------------------
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

    ``analysis/sim/sim_utils.py`` carries a byte-identical copy: the simulation
    tree is run with ``PYTHONPATH=analysis/sim`` and cannot import ``lib``.
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


def symmetrized_radial_interpolator(single_map: np.ndarray, validate: bool = True):
    """Build a 1D radius -> kappa function from a radially symmetrized stack.

    After :func:`symmetrize_map`, the map is *by construction* a function of
    distance from its physical center (``(N - 1) / 2`` in pixel-index
    coordinates), so a radial profile represents it exactly.  Radii beyond
    the map corners return 0.

    Radius is in pixels, which equals h^-1 Mpc on the standard 100-pixel /
    100 h^-1 Mpc stacking grid.

    *validate* guards against archived mis-centered stacks.  Pass ``False`` only
    when handing in an unsymmetrized map on purpose -- the radial average taken
    here is then the symmetrization, which is exactly what
    ``combine_filament_jackknife.py --no-symmetrize-single`` is testing.
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
    values = single_map.ravel()

    bins = np.arange(0.0, r.max() + 1.0, 0.5)
    idx = np.clip(np.digitize(r, bins) - 1, 0, len(bins) - 2)
    sums = np.bincount(idx, weights=values, minlength=len(bins) - 1)
    counts = np.bincount(idx, minlength=len(bins) - 1)
    good = counts > 0
    centers = (0.5 * (bins[:-1] + bins[1:]))[good]
    means = sums[good] / counts[good]

    def profile(radius: np.ndarray) -> np.ndarray:
        return np.interp(radius, centers, means, left=means[0], right=0.0)

    return profile


def two_halo_template(
    single_map: np.ndarray,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    separation_hmpc: float,
    validate: bool = True,
) -> np.ndarray:
    """Superposed-singles control: two copies of the single stack at +/- sep/2.

    Evaluates the radial profile of the (already symmetrized) single map at the
    *target* grid's own coordinates, so the control is built natively on
    whatever grid it will be subtracted from.

    This replaces shifting a pixel array with ``np.roll``/``ndimage.shift`` and
    then trimming to a common shape.  That route carried a half-pixel error: the
    pair grid is 101 pixels with centers on integers while the single grid is
    100 pixels with centers on half-integers, and trimming one row does not
    align them.  The residual showed up as false structure at the halo peaks,
    where the map gradient is steepest -- worst at small separation, where half
    a pixel is a large fraction of the pair separation.  The simulation side
    has used this profile-based construction since July 2026; the observational
    chain did not.
    """
    profile = symmetrized_radial_interpolator(single_map, validate=validate)
    offset = 0.5 * separation_hmpc
    return profile(np.hypot(x_grid + offset, y_grid)) + profile(
        np.hypot(x_grid - offset, y_grid)
    )


# ---------------------------------------------------------------------------
# Pair-stack band geometry
# ---------------------------------------------------------------------------
# Collapsing a 2-D pair stack to a 1-D profile means choosing which rows to
# average.  These bands are fixed in h^-1 Mpc and identical at every pair
# separation.  The pair grid is 101 px over 100 h^-1 Mpc, so rows are 1 h^-1 Mpc
# tall and centred on integer Y; the cuts sit halfway between rows, making the
# selection exact and immune to float noise in the axis:
#
#   central     the three rows centred at Y = -1, 0, +1   (|Y| <= 1.5)
#   off-centre  the rows centred at |Y| = 2 ... 10        (1.5 <= |Y| <= 10.5)
#
# An earlier version scaled these with the pair separation (0.15/0.45/0.85 times
# sep), which measured each separation with a different ruler: at sep 5 the
# central band collapsed to the single row Y = 0, and at sep 20 it smeared 7 rows
# together.  Fixed physical bands put every separation on the same footing.
#
# X still scales with the separation -- the bridge is the region *between* the
# two galaxies, so its width has to follow them.
#
# This lives in lib/ rather than in either analysis tree because both the BOSS
# and simulation pipelines score with it, and a duplicated copy is exactly how
# the two definitions drifted apart before.
CENTRAL_HALF_Y_HMPC = 1.5
OFF_LO_Y_HMPC, OFF_HI_Y_HMPC = 1.5, 10.5
BRIDGE_HALF_X_FRAC = 0.35


def band_masks(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Row masks for the central and off-centre bands.

    The two definitions share the boundary value 1.5, which looks like double
    counting until you check the grid: rows sit on integer Y, so no row has
    |Y| = 1.5 and the masks are disjoint -- 3 rows central, 18 rows off-centre.
    The off-centre band is therefore 6x better averaged per column, and will
    carry correspondingly more weight in any fit that treats the two jointly.
    """
    central = np.abs(axis) <= CENTRAL_HALF_Y_HMPC
    far = (np.abs(axis) >= OFF_LO_Y_HMPC) & (np.abs(axis) <= OFF_HI_Y_HMPC)
    if not central.any() or not far.any():
        raise ValueError("Empty band: the grid is too coarse for the fixed Y bands.")
    if (central & far).any():
        raise ValueError("Central and off-centre bands overlap on this grid.")
    return central, far


def band_profiles(arr: np.ndarray, axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Central-band and off-centre-band means *separately*, along the pair axis.

    Returns ``(central, off_centre)``, one value per column of *arr*.  This is
    the un-differenced form: :func:`band_profile` returns their difference,
    which cancels any map-wide additive constant, while keeping them apart
    preserves that information for a joint fit to decide what to do with.
    """
    central, far = band_masks(axis)
    return arr[central, :].mean(axis=0), arr[far, :].mean(axis=0)


def band_row_weights(axis: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Fractional overlap of each row with the physical shell ``lo <= |Y| <= hi``.

    :func:`band_masks` can only take whole rows, so it needs the band edges to
    fall between row centres.  That holds on the 101-pixel pair grid (centres on
    integers) and fails on the 100-pixel single grid (centres on half-integers,
    so |Y| = 1.5 *is* a row centre and the two bands overlap -- band_masks
    raises).  Weighting the boundary rows by how much of them lies inside the
    shell gives the same physical band on either grid:

        100-grid central   rows +/-0.5 (w = 1), +/-1.5 (w = 0.5)   -> width 3
        101-grid central   rows 0, +/-1        (w = 1)             -> width 3

    The overlap is computed in *signed* Y and then folded.  Doing it in |Y|
    halves the row straddling Y = 0 and silently returns a width of 2.5 on the
    101 grid.
    """
    axis = np.asarray(axis, dtype=float)
    cell = float(np.diff(axis).mean())
    edge_lo, edge_hi = axis - 0.5 * cell, axis + 0.5 * cell

    def overlap(a: float, b: float) -> np.ndarray:
        return np.clip(np.minimum(edge_hi, b) - np.maximum(edge_lo, a), 0.0, None)

    return (overlap(lo, hi) + overlap(-hi, -lo)) / cell


def band_profiles_fractional(arr: np.ndarray,
                             axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """:func:`band_profiles` with fractional boundary rows; works on any grid.

    Returns ``(central, off_centre)``.  On the 101-pixel grid this reproduces
    :func:`band_profiles` exactly, because every weight there is 0 or 1.
    """
    w_cen = band_row_weights(axis, 0.0, CENTRAL_HALF_Y_HMPC)
    w_off = band_row_weights(axis, OFF_LO_Y_HMPC, OFF_HI_Y_HMPC)
    if w_cen.sum() <= 0 or w_off.sum() <= 0:
        raise ValueError("Empty band: the grid is too coarse for the fixed Y bands.")
    return ((arr * w_cen[:, None]).sum(axis=0) / w_cen.sum(),
            (arr * w_off[:, None]).sum(axis=0) / w_off.sum())


def band_profile(arr: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Central-band mean minus off-centre-band mean, along the pair axis.

    *axis* is the physical coordinate of the grid in h^-1 Mpc, shared by both
    dimensions of the square map *arr*.  Returns one value per column of *arr*.
    """
    central, far = band_profiles(arr, axis)
    return central - far


def band_response_zones(x: np.ndarray, sep: float,
                        x_bin: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exclusive bridge / halo / outer masks over binned X centres.

    Assigned bridge-first so the three are disjoint and exhaust the range:

        bridge  X <= BRIDGE_HALF_X_FRAC * sep      the filament statistic's window
        halo    |X - sep/2| <= max(x_bin/2, 1.5)   the galaxies themselves
        outer   everything beyond

    A zone can come out empty.  At sep = 5 with 4 h^-1 Mpc bins the bridge
    (X <= 1.75) and the halo (X = 2.5) fall in the same first bin, so the
    binning cannot separate them at all -- which is worth knowing before reading
    any bridge interpretation into a fit at that resolution.

    Lives here rather than in either script because the fit and its diagnostic
    plot must not drift apart on where the halo is.
    """
    halo_half = max(x_bin / 2.0, 1.5)
    bridge = x <= BRIDGE_HALF_X_FRAC * sep
    halo = (~bridge) & (np.abs(x - sep / 2.0) <= halo_half)
    return bridge, halo, ~bridge & ~halo


def x_bin_edges(axis: np.ndarray, x_max: float, x_bin: float) -> np.ndarray:
    """Bin index for each column, or -1 for columns outside [0, x_max].

    Only X >= 0 is kept: the maps are reflection-symmetrized, so the profile is
    even in X and the negative half carries no independent information.

    Bin i spans [i*x_bin, (i+1)*x_bin), and X = x_max gets its own bin rather
    than being folded into the previous one.  That matters for the *unbinned*
    case: at x_bin = 1 on a 1 h^-1 Mpc grid, X = 0..20 must give 21 components,
    one per native column, so the vector is 42 long.  An earlier version used
    ceil(x_max/x_bin) bins, which merged X = 19 and X = 20 and silently returned
    40 -- i.e. it applied a residual binning to a configuration that was
    supposed to have none.

    The cost is that a coarse binning leaves a short final bin: at x_bin = 4,
    x_max = 40 the last bin holds only X = 40.  Bin centres are reported as the
    mean X of the columns actually in each bin, so that stays visible.
    """
    if x_bin <= 0 or x_max <= 0:
        raise ValueError(f"Need positive x_max and x_bin, got {x_max}, {x_bin}.")
    n_bins = int(np.floor(x_max / x_bin)) + 1
    idx = np.full(axis.shape, -1, dtype=int)
    inside = (axis >= 0) & (axis <= x_max)
    idx[inside] = np.clip((axis[inside] // x_bin).astype(int), 0, n_bins - 1)
    return idx


def kappa_band(arr: np.ndarray, axis: np.ndarray, x_max: float = 40.0,
               x_bin: float = 4.0) -> Tuple[np.ndarray, np.ndarray]:
    """Compress a pair map to the joint two-band profile vector.

    Returns ``(x_centres, vector)`` where *vector* is the concatenation

        [ kappa_central(X_0..X_n) , kappa_off(X_0..X_n) ]

    of length ``2 * n_bins``.  This is the compression that makes the
    covariance tractable: the full 101x101 map has ~2600 independent pixels,
    far too many to estimate a covariance for from 287 jackknife samples, while
    this vector is a few tens of numbers.

    Binning to roughly the beam scale (8 arcmin ~ 3.3 h^-1 Mpc, hence the 4
    h^-1 Mpc default) is not just a convenience: finer bins add no independent
    information, only more covariance to estimate from the same 287 samples.
    """
    central, off = band_profiles(arr, axis)
    idx = x_bin_edges(axis, x_max, x_bin)
    n_bins = idx.max() + 1
    if n_bins < 1:
        raise ValueError(f"No columns inside X = 0..{x_max}.")

    counts = np.bincount(idx[idx >= 0], minlength=n_bins).astype(float)
    if (counts == 0).any():
        raise ValueError(f"Empty X bin at x_bin={x_bin}: grid too coarse.")

    def binned(vec: np.ndarray) -> np.ndarray:
        return np.bincount(idx[idx >= 0], weights=vec[idx >= 0], minlength=n_bins) / counts

    x_centres = binned(axis)
    return x_centres, np.concatenate([binned(central), binned(off)])


def bridge_excess(arr: np.ndarray, axis: np.ndarray, rperp_center: float) -> float:
    """Bridge statistic: :func:`band_profile` averaged over the bridge window.

    Defined as the windowed mean of the profile rather than as an independent
    2-D box mean.  The two are algebraically identical here -- the window is a
    rectangle and the row counts do not vary with x -- but deriving one from the
    other means the plotted curve and the quoted number cannot drift apart.
    """
    window = np.abs(axis) <= BRIDGE_HALF_X_FRAC * rperp_center
    if not window.any():
        raise ValueError(f"Empty bridge window at rperp_center={rperp_center}.")
    return float(band_profile(arr, axis)[window].mean())


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
