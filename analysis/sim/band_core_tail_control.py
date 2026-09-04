#!/usr/bin/env python
"""Zheng's band-space core + tail control, central band only.

His instruction, taken literally: do not build the control from the circular
radial profile.  Take the *symmetrized single kappa map*, average the same
Y-rows the pair band profile uses, and decompose that 1-D band profile.

    b(x) = core(x) + tail(x),    core(x) = b(x)                    |x| <= r_c
                                        = C exp(-x^2 / 2 sigma^2)  |x| >  r_c

with value and slope matched at r_c, so sigma^2 = -r_c b(r_c) / b'(r_c) is
measured rather than chosen.  The control is then

    T(x) = A_c [core(u-) + core(u+)] + A_t [tail(s_t u-) + tail(s_t u+)] + k
    u±   = |x ± a/2|

averaged over the measured r_perp distribution, not the nominal separation.

Three things this file does differently from ``core_tail_control.py``, each
because the earlier choice was wrong:

  fractional band rows   the single stack is 100 px (rows on half-integers) and
                         the pair stack is 101 px (rows on integers), so
                         ``band_masks`` cannot even be evaluated on the single
                         grid -- |Y| = 1.5 is a row centre and the bands
                         overlap.  ``band_profiles_fractional`` weights the
                         boundary rows and gives the same physical band, width
                         3 and 18, on either grid.  No re-stack needed, and no
                         pixel array is shifted, so the July half-pixel bug
                         cannot come back.

  a free additive constant  the central band is not a small positive profile:
                         it crosses zero near |x| = 20 and settles on a
                         -6.3e-4 shelf, the periodic box's zero-mean
                         constraint.  The old central-minus-off statistic
                         cancelled that to ~2e-5.  Fitting the central band
                         ALONE over all |x| > a/2 does not, so without a free
                         constant A_c and A_t would absorb an offset six times
                         the bridge signal.  The ``full_noconst`` scheme is
                         kept to show exactly that.

  x >= 0 only            the pair maps are reflection-symmetrized, so the two
                         halves are the same numbers.  Fitting both is exact
                         double counting.

Support: a shifted, stretched template needs the single band profile at
s(|x| + a/2).  Beyond the single's own axis it does not exist, and it is NOT
extended -- splicing in the radial profile would mix the old and new
constructions inside one template.  The fit and every reported number are
restricted to the domain supported for *every* s in the grid, so the s-scan
compares models on identical points.  The excluded range is reported.

Usage
-----
    PYTHONPATH=lib:analysis/sim python analysis/sim/band_core_tail_control.py
    PYTHONPATH=lib:analysis/sim python analysis/sim/band_core_tail_control.py --null
    PYTHONPATH=lib:analysis/sim python analysis/sim/band_core_tail_control.py --figure
    PYTHONPATH=lib:analysis/sim python analysis/sim/band_core_tail_control.py \
        --validate-grid analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import lsq_linear

from geometry import band_profiles, band_profiles_fractional, reflect_symmetrize_map
from deconvolve_pair_profile import axis_for, load_map
from stretched_single_control import (HOD, NULL_BANDS, PHYSICAL, SIM,
                                      jackknife, maps_from_blocks)

# The single stack on the *pair* grid.  The 100-pixel stack can be read with
# band_profiles_fractional and gives the same band width, but not the same band
# *sampling*: its rows sit at |Y| = 0.5, 1.5 against the 101-grid's 0, 1, so the
# even grid weights the band toward larger |Y| (<Y^2> = 0.92 vs 0.67) and, on a
# profile that falls with |Y|, comes out 2-4% low at every X.  Measured against
# this stack that bias moves A_c by 0.02-0.04 -- most of the apparent 3-6%
# excess of the paired halo over the field single -- and the halo residual at
# a = 20 by ~4 sigma.  The bridge moves by less than 1 sigma at a = 10 and 20.
# So the fractional extractor is the right mechanism for reading an even grid,
# but it is not a substitute for stacking on the grid you will subtract from.
SINGLE = "kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv"
SINGLE_EVEN = "kappa_single_sim_hodnmatch_8arcmin_centered.csv"

logger = logging.getLogger(__name__)

FWHM_HMPC = 3.31                                              # 8 arcmin at z = 0.547
SIGMA_BEAM = FWHM_HMPC / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # 1.406
S_GRID = np.round(np.arange(0.70, 1.3001, 0.005), 4)
HALO_HALF = 1.5                                               # halo window half-width
CONST_SCALE = 1e-4    # constant column value, so its coefficient is O(1) like
                      # the amplitudes.  A column of ones against templates of
                      # order 1e-3 makes cond(M'M) ~ 1e8 for purely numerical
                      # reasons and hides whether the fit is really degenerate.
BRIDGE_FRAC = 0.25                                            # bridge half-width / a
R_CORE_DEFAULT = 3.0


# --------------------------------------------------------------------------
# the single, in band space
# --------------------------------------------------------------------------
def single_band_profile(single_map):
    """Central-band profile of the symmetrized single stack, folded to x >= 0.

    Returns (x, b) with x >= 0.  Folding is exact for a radially symmetrized
    map and halves the interpolation noise on the outermost points.
    """
    axis = axis_for(single_map.shape[0])
    cen, _ = band_profiles_fractional(single_map, axis)
    x = np.abs(axis)
    order = np.argsort(x)
    x, cen = x[order], cen[order]
    ux, inv = np.unique(np.round(x, 6), return_inverse=True)
    b = np.bincount(inv, weights=cen) / np.bincount(inv)
    return ux, b


def band_spline(x, b):
    """b(|x|), NaN past the measured range -- never zero-filled, never extended.

    The spline is built on the *even extension* of the folded profile, over
    [-x_max, +x_max].  Without it the 100-pixel grid has no sample inside
    |x| = 0.5, so every template argument below that returns NaN -- and the
    binding argument is |X| - a/2, which is exactly zero at the halo centre.
    That punched a hole in the fit and in the halo diagnostic at |X| = a/2 for
    every separation.  Mirroring is exact for a symmetrized map (the profile is
    even by construction), so x = 0 is bracketed by data and interpolated, not
    extrapolated.  The 101-pixel grid has a sample at x = 0 and does not need
    this, which is one more consequence of the single/pair grid offset.
    """
    xs, idx = np.unique(np.concatenate([-x[::-1], x]), return_index=True)
    bs = np.concatenate([b[::-1], b])[idx]
    spl = CubicSpline(xs, bs, extrapolate=False)
    x_max = float(x[-1])

    def f(q):
        return spl(np.abs(np.asarray(q, dtype=float)))

    return f, x_max


def decompose_band(f, r_c, h=0.05):
    """Measured band profile inside r_c, slope-matched Gaussian outside.

    Same algebra as the radial version -- matching value and slope of
    C exp(-x^2/2 sigma^2) to b at x = r_c gives sigma^2 = -r_c b / b'.  The
    tail is identically zero inside r_c, so A_t and s_t cannot touch the peak.
    """
    r_c = float(r_c)
    bc = float(f(r_c))
    dbc = float((f(r_c + h) - f(r_c - h)) / (2.0 * h))
    if not np.isfinite(bc) or not np.isfinite(dbc) or dbc >= 0:
        raise ValueError(f"band profile is not falling at r_c = {r_c}")
    sigma = float(np.sqrt(-r_c * bc / dbc))
    C = bc * np.exp(r_c ** 2 / (2.0 * sigma ** 2))

    def core(q):
        q = np.abs(np.asarray(q, dtype=float))
        inner = np.asarray(f(q), dtype=float)
        outer = C * np.exp(-q ** 2 / (2.0 * sigma ** 2))
        return np.where(q <= r_c, inner, outer)

    def tail(q):
        return np.asarray(f(q), dtype=float) - core(q)

    return core, tail, sigma, C


# --------------------------------------------------------------------------
# templates
# --------------------------------------------------------------------------
def superpose(g, x, seps, wts, s=1.0):
    """sum of two copies of g at +/- a/2, averaged over the r_perp mix."""
    out = np.zeros_like(np.asarray(x, dtype=float))
    for a, w in zip(seps, wts):
        for sign in (-1.0, +1.0):
            out = out + w * np.asarray(g(s * np.abs(x + sign * 0.5 * a)), float)
    return out


def build_bank(core, tail, x, seps, wts):
    """core template (s = 1) and the tail template at every s in the grid."""
    return (superpose(core, x, seps, wts, 1.0),
            {float(s): superpose(tail, x, seps, wts, float(s)) for s in S_GRID})


def support_mask(x, seps, x_max, s_max):
    """Bins where every template in the s-scan is defined, x >= 0 only.

    The binding argument is s(|x| + a/2) for the largest a in the mix and the
    largest s in the grid.  Using one mask for all s keeps the s-scan an
    apples-to-apples comparison of models rather than of point counts.
    """
    a_max = float(np.max(seps))
    return (x >= 0) & (s_max * (np.abs(x) + 0.5 * a_max) <= x_max)


# --------------------------------------------------------------------------
# fitting
# --------------------------------------------------------------------------
def solve_bounded(M, d, with_const):
    """Least squares with A_c, A_t >= 0 and the constant free-signed.

    A negative amplitude improves chi^2 while destroying the interpretation the
    decomposition exists for, so the amplitudes are floored at zero.  The
    constant must stay free-signed: the shelf it absorbs is negative.
    """
    ok = np.isfinite(M).all(axis=1) & np.isfinite(d)
    n_par = M.shape[1]
    if ok.sum() < n_par + 1:
        return np.full(n_par, np.nan), np.nan, True
    M, d = M[ok], d[ok]
    G = M.T @ M
    scale = np.sqrt(np.diag(G))
    cond = float(np.linalg.cond(G / np.outer(scale, scale)))   # scale-free
    try:
        free = np.linalg.solve(G, M.T @ d)
    except np.linalg.LinAlgError:
        return np.full(n_par, np.nan), cond, True
    neg = bool((free[:2] < 0).any())
    if not neg:
        return free, cond, neg
    lo = np.array([0.0, 0.0] + ([-np.inf] if with_const else []))
    hi = np.full(n_par, np.inf)
    return lsq_linear(M, d, bounds=(lo, hi)).x, cond, neg


def design(core_t, tail_t, mask, with_const):
    cols = [core_t[mask], tail_t[mask]]
    if with_const:
        cols.append(np.full(int(mask.sum()), CONST_SCALE))
    return np.column_stack(cols)


def fit_control(p_cen, x, core_t, tail_bank, fit_mask, with_const, s_free=True):
    """Scan s_t, solve A_c, A_t (+ constant) linearly at each s, keep the best."""
    grid = list(tail_bank) if s_free else [1.0]
    d = p_cen[fit_mask]
    best = None
    for s in grid:
        M = design(core_t, tail_bank[s], fit_mask, with_const)
        p, cond, neg = solve_bounded(M, d, with_const)
        if not np.isfinite(p).all():
            continue
        ok = np.isfinite(M).all(axis=1) & np.isfinite(d)
        r = float(np.sqrt(np.mean((d[ok] - M[ok] @ p) ** 2)))
        if best is None or r < best[0]:
            best = (r, p, s, cond, neg)
    if best is None:
        return None
    r, p, s, cond, neg = best
    return dict(A_c=float(p[0]), A_t=float(p[1]),
                const=float(p[2]) * CONST_SCALE if with_const else 0.0,
                s_t=float(s), fit_rms=r, cond=cond, neg=neg)


def control_curve(fit, core_t, tail_bank):
    return fit["A_c"] * core_t + fit["A_t"] * tail_bank[fit["s_t"]] + fit["const"]


# --------------------------------------------------------------------------
# windows and scoring
# --------------------------------------------------------------------------
HOLD_FRAC = 0.25   # outermost fraction of the supported exterior, never fitted


def hold_edge(x, sep, sup):
    """Where the never-fitted outer holdout starts, in u."""
    u = np.abs(x) - 0.5 * sep
    uu = u[sup & (u > 0.0)]
    if uu.size == 0:
        return np.inf
    return float(uu.min() + (1.0 - HOLD_FRAC) * (uu.max() - uu.min()))


def fit_windows(x, sep, sup):
    """Primary window and the robustness variants, all inside the support.

    ``full`` is the advisor-requested fit: every supported bin with |X| > a/2.
    ``full_hold`` fits the inner three quarters of that range so the outermost
    quarter is a genuine prediction -- necessary because a free constant makes
    the *mean* residual over the fitted region identically zero, so an in-fit
    off-halo mean is not a test of anything.
    """
    u = np.abs(x) - 0.5 * sep
    ue = hold_edge(x, sep, sup)
    return {
        "full":         sup & (u > 0.0),            # Zheng's instruction
        "full_noconst": sup & (u > 0.0),            # same window, no constant
        "full_hold":    sup & (u > 0.0) & (u <= ue),
        "buf2":         sup & (u > 2.0),
        "buf4":         sup & (u > 4.0),
        "win415":       sup & (u > 4.0) & (u <= 15.0),
    }


def score(diff, x, sep, sup):
    """Bridge, halo, off-halo and holdout residuals -- separate, never summed.

    ``outer`` is the mean off-halo residual, and it is close to zero *by
    construction* for any scheme carrying a free constant: least squares sets
    the mean residual over the fitted region to zero, and the fitted region and
    the off-halo region nearly coincide.  ``outer_rms`` and ``hold`` are the
    off-halo numbers that actually carry information.
    """
    u = np.abs(x) - 0.5 * sep
    ue = hold_edge(x, sep, sup)
    bridge = sup & (np.abs(x) <= BRIDGE_FRAC * sep)
    halo = sup & (np.abs(u) <= HALO_HALF)
    outer = sup & (u > HALO_HALF)
    hold = sup & (u > ue)

    def safe(w, square=False):
        if not w.any():
            return np.nan
        v = diff[w]
        v = v[np.isfinite(v)]
        if v.size == 0:
            return np.nan
        return float(np.sqrt(np.mean(v ** 2))) if square else float(np.mean(v))

    return dict(bridge=safe(bridge), halo=safe(halo), outer=safe(outer),
                outer_rms=safe(outer, True), hold=safe(hold),
                hold_rms=safe(hold, True))


# --------------------------------------------------------------------------
# samples
# --------------------------------------------------------------------------
def load_sample(blocks, catalog, space, cut):
    full, loo, n_pairs, sep = maps_from_blocks(blocks)
    full = reflect_symmetrize_map(full.astype(np.float64))
    loo = [reflect_symmetrize_map(m.astype(np.float64)) for m in loo]
    col = "r_parallel_real" if space == "real" else "r_parallel_rsd"
    df = pd.read_csv(catalog, usecols=["r_perp", col])
    lo = -np.inf if cut[0] is None else cut[0]
    sel = df["r_perp"][(df[col] >= lo) & (df[col] <= cut[1])].to_numpy()
    counts, edges = np.histogram(sel, bins=21)
    keep = counts > 0
    seps = (0.5 * (edges[:-1] + edges[1:]))[keep]
    wts = counts[keep] / counts[keep].sum()
    return full, loo, int(n_pairs), sep, seps, wts


def jobs_for(null):
    if null:
        return [(b, HOD / f"nonphys_rperp10_rpar{b}_blocks.npz",
                 HOD / "pairs_hodnmatch50_rperp10_real120.csv", "real",
                 tuple(float(v) for v in b.split("_"))) for b in NULL_BANDS]
    return [(k, HOD / v[0], HOD / v[1], "rsd", (None, 10.0))
            for k, v in PHYSICAL.items()]


def pair_central_band(pair_map, axis):
    return band_profiles(pair_map, axis)[0]


# --------------------------------------------------------------------------
# one (sample, r_c) analysis
# --------------------------------------------------------------------------
def analyse(sample, bfun, x_max_single, full, loo, axis, sep, seps, wts,
            r_c, rows):
    core, tail, sigma, C = decompose_band(bfun, r_c)
    core_t, tail_bank = build_bank(core, tail, axis, seps, wts)
    s_max = float(max(tail_bank))
    sup = support_mask(axis, seps, x_max_single, s_max)
    x_sup = float(np.max(axis[sup])) if sup.any() else np.nan
    n_excl = int(((axis >= 0) & ~sup).sum())
    windows = fit_windows(axis, sep, sup)

    p_full = pair_central_band(full, axis)
    p_loo = [pair_central_band(m, axis) for m in loo]
    logger.info("%s r_c=%.1f: sigma=%.3f (beam sigma %.3f), support |x| <= %.0f, "
                "%d bins dropped at the edge", sample, r_c, sigma, SIGMA_BEAM,
                x_sup, n_excl)

    # --- the fitted controls ------------------------------------------------
    out_curves = {}
    for name, mask in windows.items():
        with_const = name != "full_noconst"
        F = fit_control(p_full, axis, core_t, tail_bank, mask, with_const)
        if F is None:
            logger.warning("  %-12s no fit", name)
            continue
        dF = p_full - control_curve(F, core_t, tail_bank)
        sF = score(dF, axis, sep, sup)
        L = [fit_control(p, axis, core_t, tail_bank, mask, with_const)
             for p in p_loo]
        L = [v for v in L if v is not None]
        dL = [score(p - control_curve(v, core_t, tail_bank), axis, sep, sup)
              for p, v in zip(p_loo, L)]
        get = lambda k: np.array([v[k] for v in L], dtype=float)
        gets = lambda k: np.array([v[k] for v in dL], dtype=float)
        rows.append(dict(
            sample=sample, separation=sep, r_c=r_c, scheme=name,
            sigma=sigma, n_fit=int(mask.sum()), x_support=x_sup,
            n_edge_excluded=n_excl, cond=F["cond"], fit_rms=F["fit_rms"] * 1e4,
            A_c=F["A_c"], A_c_err=jackknife(get("A_c")),
            A_t=F["A_t"], A_t_err=jackknife(get("A_t")),
            s_t=F["s_t"], s_t_err=jackknife(get("s_t")),
            const=F["const"] * 1e4, const_err=jackknife(get("const")) * 1e4,
            bridge=sF["bridge"] * 1e4, bridge_err=jackknife(gets("bridge")) * 1e4,
            halo=sF["halo"] * 1e4, halo_err=jackknife(gets("halo")) * 1e4,
            outer=sF["outer"] * 1e4, outer_err=jackknife(gets("outer")) * 1e4,
            outer_rms=sF["outer_rms"] * 1e4,
            outer_rms_err=jackknife(gets("outer_rms")) * 1e4,
            hold=sF["hold"] * 1e4, hold_err=jackknife(gets("hold")) * 1e4,
            hold_rms=sF["hold_rms"] * 1e4,
            hold_rms_err=jackknife(gets("hold_rms")) * 1e4,
            hold_edge_u=hold_edge(axis, sep, sup),
            n_neg=int(sum(bool(v["neg"]) for v in L)),
            n_s_edge=int(np.sum((get("s_t") <= S_GRID[0]) | (get("s_t") >= S_GRID[-1]))),
        ))
        r = rows[-1]
        logger.info("  %-12s A_c=%.3f+/-%.3f A_t=%.3f+/-%.3f s_t=%.3f+/-%.3f "
                    "k=%+.3f+/-%.3f  cond %6.1f fit %.3f | bridge %+7.3f+/-%5.3f "
                    " halo %+7.3f+/-%5.3f  outer rms %.3f  hold %+7.3f+/-%5.3f",
                    name, r["A_c"], r["A_c_err"], r["A_t"], r["A_t_err"],
                    r["s_t"], r["s_t_err"], r["const"], r["const_err"],
                    r["cond"], r["fit_rms"], r["bridge"], r["bridge_err"],
                    r["halo"], r["halo_err"], r["outer_rms"],
                    r["hold"], r["hold_err"])
        if name in ("full", "full_noconst"):
            out_curves[name] = (F, dF)

    # --- the plain superposition, no fitted parameters -----------------------
    plain = superpose(bfun, axis, seps, wts, 1.0)
    dP = p_full - plain
    sP = score(dP, axis, sep, sup)
    dPl = [score(p - plain, axis, sep, sup) for p in p_loo]
    g = lambda k: np.array([v[k] for v in dPl], dtype=float)
    rows.append(dict(
        sample=sample, separation=sep, r_c=r_c, scheme="plain",
        sigma=sigma, n_fit=0, x_support=x_sup, n_edge_excluded=n_excl,
        cond=np.nan, fit_rms=np.nan, A_c=1.0, A_c_err=0.0, A_t=1.0, A_t_err=0.0,
        s_t=1.0, s_t_err=0.0, const=0.0, const_err=0.0,
        bridge=sP["bridge"] * 1e4, bridge_err=jackknife(g("bridge")) * 1e4,
        halo=sP["halo"] * 1e4, halo_err=jackknife(g("halo")) * 1e4,
        outer=sP["outer"] * 1e4, outer_err=jackknife(g("outer")) * 1e4,
        outer_rms=sP["outer_rms"] * 1e4,
        outer_rms_err=jackknife(g("outer_rms")) * 1e4,
        hold=sP["hold"] * 1e4, hold_err=jackknife(g("hold")) * 1e4,
        hold_rms=sP["hold_rms"] * 1e4,
        hold_rms_err=jackknife(g("hold_rms")) * 1e4,
        hold_edge_u=hold_edge(axis, sep, sup),
        n_neg=0, n_s_edge=0))
    r = rows[-1]
    logger.info("  %-12s no free parameters                                  "
                "                    | bridge %+7.3f+/-%5.3f  halo %+7.3f+/-%5.3f"
                "  outer rms %.3f  hold %+7.3f+/-%5.3f", "plain",
                r["bridge"], r["bridge_err"], r["halo"], r["halo_err"],
                r["outer_rms"], r["hold"], r["hold_err"])
    out_curves["plain"] = (None, dP)
    out_curves["meta"] = dict(p_full=p_full, plain=plain, core_t=core_t,
                              tail_bank=tail_bank, sup=sup, sigma=sigma,
                              windows=windows,
                              hold_edge_u=hold_edge(axis, sep, sup))
    return out_curves


# --------------------------------------------------------------------------
# validation against a real 101-pixel stack
# --------------------------------------------------------------------------
def validate_grid(single_100, path_101):
    m101 = load_map(path_101)
    x0, b0 = single_band_profile(single_100)
    x1, b1 = single_band_profile(m101)
    f0, _ = band_spline(x0, b0)
    both = x1[(x1 <= x0[-1]) & (x1 >= x0[0])]
    a, b = np.asarray(f0(both), float), np.interp(both, x1, b1)
    rel = np.abs(a - b) / np.abs(b)
    inner = both <= 15.0
    logger.info("grid validation vs %s", Path(path_101).name)
    logger.info("  100-px fractional + spline vs native 101-px central band:")
    logger.info("  |x| <= 15 : max %.2e abs, %.3f%% median rel",
                np.abs(a - b)[inner].max(), 100 * np.median(rel[inner]))
    logger.info("  all bins  : max %.2e abs, %.3f%% median rel",
                np.abs(a - b).max(), 100 * np.median(rel))
    for q in (0., 1., 2., 3., 5., 10., 20., 30., 40.):
        i = int(np.argmin(np.abs(both - q)))
        logger.info("    x=%5.1f  100frac %+.5e   native101 %+.5e   ratio %.4f",
                    both[i], a[i], b[i], a[i] / b[i])
    return pd.DataFrame(dict(x=both, from_100_fractional=a, native_101=b))


# --------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------
def make_figure(store, out_dir, r_c):
    """Two separate figures, one per control, as Zheng asked for them.

    Each is 1 x 3 (r_perp = 5, 10, 20) and shows the measured pair central
    band, the control central band, and their difference.  The difference is
    ~1e-4 against a peak of ~85e-4, so it gets its own strip under each panel
    rather than being an invisible line on the same axes.  Everything is
    mirrored about X = 0: the pair maps are reflection-symmetrized, so the
    negative half is the same numbers and the fit uses X >= 0 only.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from deconvolve_pair_profile import GREEN, INK, INK_MUTED, tidy
    BLUE, ORANGE, PURPLE = "#2a78d6", "#eb6834", "#7d5ba6"

    def mirror(x, y):
        keep = x > 0
        return (np.concatenate([-x[keep][::-1], x[keep]]),
                np.concatenate([y[keep][::-1], y[keep]]))

    keys = [k for k in ("5", "10", "20") if k in store]

    FIGS = {
        "fitted": dict(
            stem="band_core_tail_fitted_central",
            title=("Core + tail control, fitted outside the galaxies "
                   f"(matched core, $r_c$ = {r_c:g})\n"
                   "Dotted purple is the fit as specified: $A_c$, $A_t$, $s_t$ "
                   "free on all supported $|X| > a/2$.  Dashed blue adds a free "
                   "additive constant.\n"
                   "Grey = between the galaxies (never fitted).  Dotted green = "
                   "galaxy positions, green strip = halo window, orange = "
                   "never-fitted outer holdout."),
        ),
        "plain": dict(
            stem="band_superposed_central",
            title=("Plain superposition of singles — no fitted parameters\n"
                   "The single's own central-band profile added at $\\pm a/2$, "
                   "subtracted from the measured pair band as it stands.\n"
                   "Grey = between the galaxies, dotted green = galaxy "
                   "positions, green strip = halo window."),
        ),
    }

    for which, cfg in FIGS.items():
        fig = plt.figure(figsize=(4.9 * len(keys), 4.9), constrained_layout=True)
        gs = GridSpec(2, len(keys), figure=fig, height_ratios=[3.0, 1.3])

        for col, key in enumerate(keys):
            cur = store[key]
            meta = cur["meta"]
            axis, sep = meta["axis"], meta["sep"]
            sup, p = meta["sup"], meta["p_full"]
            x_sup = float(np.max(axis[sup]))
            u_hold = meta["hold_edge_u"]

            ax = fig.add_subplot(gs[0, col])
            rx = fig.add_subplot(gs[1, col], sharex=ax)

            if which == "fitted":
                F, diff = cur["full"]
                F0, diff0 = cur["full_noconst"]
                ctl = control_curve(F, meta["core_t"], meta["tail_bank"])
                ctl0 = control_curve(F0, meta["core_t"], meta["tail_bank"])
                sub = (rf"$A_c$={F0['A_c']:.3f}  $A_t$={F0['A_t']:.3f}  "
                       rf"$s_t$={F0['s_t']:.3f}" + "\n"
                       rf"+const:  $A_c$={F['A_c']:.3f}  $A_t$={F['A_t']:.3f}  "
                       rf"$s_t$={F['s_t']:.3f}  $k$={F['const']*1e4:+.2f}"
                       rf"$\times10^{{-4}}$")
            else:
                diff = cur["plain"][1]
                ctl = meta["plain"]
                ctl0 = diff0 = None
                sub = "no free parameters"

            xm, pm = mirror(axis, p)
            _, cm = mirror(axis, ctl)
            _, dm = mirror(axis, diff)

            ax.plot(xm, pm * 1e4, lw=2.2, color=INK, label="measured pair",
                    zorder=4)
            if ctl0 is not None:
                _, c0m = mirror(axis, ctl0)
                _, d0m = mirror(axis, diff0)
                ax.plot(xm, c0m * 1e4, lw=1.7, color=PURPLE, ls=(0, (1, 1.6)),
                        label=r"control: $A_c, A_t, s_t$", zorder=6)
                rx.plot(xm, d0m * 1e4, lw=1.6, color=PURPLE, ls=(0, (1, 1.6)),
                        zorder=5)
            ax.plot(xm, cm * 1e4, lw=1.9, color=BLUE, ls=(0, (5, 2)),
                    label=("control: $+$ free constant" if ctl0 is not None
                           else "control (superposed singles)"), zorder=5)
            rx.plot(xm, dm * 1e4, lw=2.0, color=ORANGE, zorder=4)
            rx.axhline(0, color=INK, lw=0.9)

            for a_ in (ax, rx):
                a_.axvspan(-0.5 * sep, 0.5 * sep, color=INK_MUTED,
                           alpha=0.055, lw=0)
                for xx in (-0.5 * sep, 0.5 * sep):
                    a_.axvline(xx, color=GREEN, lw=1.2, ls=":", zorder=2)
                if which == "fitted":
                    for sgn in (-1, 1):
                        lo, hi = sgn * (0.5 * sep + u_hold), sgn * x_sup
                        a_.axvspan(min(lo, hi), max(lo, hi), color=ORANGE,
                                   alpha=0.07, lw=0)
                a_.set_xlim(-x_sup, x_sup)
            for sgn in (-1, 1):
                rx.axvspan(sgn * (0.5 * sep - HALO_HALF),
                           sgn * (0.5 * sep + HALO_HALF),
                           color=GREEN, alpha=0.10, lw=0)

            ax.set_ylabel(r"central band $\kappa$  [$10^{-4}$]")
            rx.set_ylabel(r"difference  [$10^{-4}$]", fontsize=8.5)
            rx.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
            ax.tick_params(labelbottom=False)
            ax.set_title(rf"$r_\perp$ = {key} $h^{{-1}}$Mpc", loc="left",
                         color=INK, fontsize=11.5,
                         pad=26 if which == "fitted" else 14)
            # above the axes, never over the data
            ax.text(0.0, 1.012, sub, transform=ax.transAxes, fontsize=8.2,
                    color=INK_MUTED, va="bottom", linespacing=1.45)
            if col == 0:
                ax.legend(loc="upper right", fontsize=8.2, frameon=False)
            tidy(ax); tidy(rx)

        fig.suptitle(cfg["title"], fontsize=10.5, color=INK, ha="left", x=0.01)
        base = Path(out_dir) / cfg["stem"]
        base.parent.mkdir(parents=True, exist_ok=True)
        for e in ("png", "pdf"):
            fig.savefig(f"{base}.{e}", bbox_inches="tight", dpi=150)
        plt.close(fig)
        logger.info("Saved %s.png / .pdf", base)


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Band-space core+tail control.")
    ap.add_argument("--single", default=SINGLE)
    ap.add_argument("--r-core-scan", default="2,3,4")
    ap.add_argument("--null", action="store_true")
    ap.add_argument("--figure", action="store_true",
                    help="also write the 2x3 central-band figure (r_c = 3)")
    ap.add_argument("--figure-r-core", type=float, default=R_CORE_DEFAULT)
    ap.add_argument("--validate-grid", default="",
                    help="path to a native 101-px single stack to check the "
                         "fractional extractor against")
    ap.add_argument("--output-dir", default="output/plots")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    single = load_map(SIM / args.single)
    x_s, b_s = single_band_profile(single)
    bfun, x_max_single = band_spline(x_s, b_s)
    logger.info("single %s: grid %d, central band on |x| = %.1f..%.1f, "
                "b(0) = %.4e, far shelf b(%.1f) = %.4e",
                args.single, single.shape[0], x_s[0], x_max_single,
                b_s[0], x_s[-1], b_s[-1])

    if args.validate_grid:
        tab = validate_grid(single, args.validate_grid)
        p = SIM / "band_extractor_grid_validation.csv"
        tab.to_csv(p, index=False)
        logger.info("Saved -> %s", p)

    axis = axis_for(101)
    rcs = [float(v) for v in args.r_core_scan.split(",")]
    rows, store = [], {}
    for label, blocks, catalog, space, cut in jobs_for(args.null):
        if not Path(blocks).exists():
            logger.warning("missing %s -- skipped", blocks)
            continue
        full, loo, n_pairs, sep, seps, wts = load_sample(blocks, catalog, space, cut)
        logger.info("=== %s: %d pairs, a = %.1f, %d r_perp bins ===",
                    label, n_pairs, sep, len(seps))
        for r_c in rcs:
            cur = analyse(label, bfun, x_max_single, full, loo, axis, sep,
                          seps, wts, r_c, rows)
            if not args.null and r_c == args.figure_r_core:
                cur["meta"].update(axis=axis, sep=sep)
                store[label] = cur

    stem = "band_core_tail_null" if args.null else "band_core_tail_physical"
    out = SIM / f"{stem}.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    logger.info("Saved -> %s", out)

    if args.figure and store:
        make_figure(store, args.output_dir, args.figure_r_core)


if __name__ == "__main__":
    main()
