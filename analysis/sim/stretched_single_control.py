#!/usr/bin/env python
"""Zheng's two-parameter stretched single-galaxy control.

The one-parameter scaled control (A only) cannot match the halo peak and the
outskirts at the same time -- one amplitude absorbs a normalization error, not
a shape error.  Zheng's proposal adds a stretch:

    b(X, Y) = A * [ f(sqrt(s^2 (X - a/2)^2 + Y^2))
                  + f(sqrt(s^2 (X + a/2)^2 + Y^2)) ]

with (A, s) fitted to the measured profile outside the pair, |X| > a/2.
s < 1 broadens the profile along X, s > 1 narrows it; the halo centres stay
at +/- a/2 either way.

Three points of construction that are easy to get wrong:

  2-D first, bands second   The stretch has to be applied inside the radius,
                            on the 2-D template, and the bands extracted from
                            that.  Evaluating an already-Y-averaged band
                            profile at sX gets the Y dependence wrong: the
                            band average and the stretch do not commute.

  average over the real     The pairs occupy a bin (4-6, 9-11, 18-22), and the
  r_perp distribution       stack places each pair at *its own* separation.  A
                            control built at the nominal a alone is sharper
                            along X than the stack it is subtracted from, so a
                            free s would partly just absorb the bin width.

  spline, not lookup        With a stretch the template is evaluated at radii
                            between the tabulated bins, so the radial profile
                            is splined.  Beyond the measured range it returns
                            NaN, never 0 -- zero-fill silently deletes the
                            outer template (this bit us once already).

Variants
--------
joint    (A, s) fitted together on |X| > a/2 + buffer.  This is Zheng's
         proposal as written, and the primary result.
split    A from the halo-peak windows, s from the outskirts, iterated to
         convergence.  A diagnostic for whether the joint fit is degenerate,
         not a rival estimator -- fitting A at the peaks forces the halo
         residual small, so it stops being an independent check.
radial   stretch applied to the full radius rather than to X alone, so the
         halo stays isotropic.  Not Zheng's formula; included because an
         X-only stretch deforms the template along the same axis as the tidal
         signal we are trying to measure.

Usage
-----
    PYTHONPATH=lib:analysis/sim python analysis/sim/stretched_single_control.py
    PYTHONPATH=lib:analysis/sim python analysis/sim/stretched_single_control.py --null
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

from geometry import band_masks, band_profiles, reflect_symmetrize_map
from deconvolve_pair_profile import axis_for, load_map
from sim_utils import validate_radially_symmetrized_map

logger = logging.getLogger(__name__)

SIM = Path("analysis/sim/results")
HOD = SIM / "hod_pairs"
SINGLE = "kappa_single_sim_hodnmatch_8arcmin_centered.csv"

# Fit geometry, all h^-1 Mpc.  Fixed here, in simulation, before BOSS.
BUFFER = 4.0            # start of the exterior fit, past a/2
OUTER = 25.0            # end of the exterior fit
HOLDOUT = (25.0, 40.0)  # never fitted; used only to score the prediction
HALO_HALF = 1.5         # half-width of the halo window, for the split variant
S_GRID = np.round(np.arange(0.70, 1.3001, 0.005), 4)

PHYSICAL = {"5": ("stack_rperp5_rpar10_blocks.npz", "pairs_hodnmatch50_rperp5_rsd50.csv"),
            "10": ("stack_rperp10_blocks.npz", "pairs_hodnmatch50_rperp10_rsd50.csv"),
            "20": ("stack_rperp20_blocks.npz", "pairs_hodnmatch50_rperp20_rsd50.csv")}
NULL_BANDS = ("40_60", "60_90", "90_120")


# --------------------------------------------------------------------------
# radial profile
# --------------------------------------------------------------------------
def radial_spline(single_map: np.ndarray, validate: bool = True):
    """Cubic-spline f(r) from the symmetrized single stack; NaN past the data.

    Returns NaN rather than 0 outside the measured range so that a template
    which reaches past the map is caught by an explicit NaN check instead of
    quietly losing its outskirts.
    """
    if validate:
        validate_radially_symmetrized_map(single_map)
    n = single_map.shape[0]
    c = 0.5 * (n - 1)
    yy, xx = np.indices((n, n))
    r = np.hypot(xx - c, yy - c).ravel()
    v = single_map.ravel()
    bins = np.arange(0.0, r.max() + 1.0, 0.5)
    idx = np.clip(np.digitize(r, bins) - 1, 0, len(bins) - 2)
    sums = np.bincount(idx, weights=v, minlength=len(bins) - 1)
    counts = np.bincount(idx, minlength=len(bins) - 1)
    good = counts > 0
    centers = (0.5 * (bins[:-1] + bins[1:]))[good]
    means = sums[good] / counts[good]
    spl = CubicSpline(centers, means, extrapolate=False)
    r_max = float(centers[-1])

    def f(radius):
        radius = np.asarray(radius, dtype=float)
        out = spl(np.clip(radius, centers[0], r_max))
        return np.where(radius > r_max, np.nan, out)

    return f, r_max


# --------------------------------------------------------------------------
# templates
# --------------------------------------------------------------------------
def rperp_weights(catalog: Path, rpar_max: float = 10.0, n_bins: int = 21):
    """Histogram of the actual pair separations entering the stack."""
    cols = ["r_perp", "r_parallel_rsd"]
    df = pd.read_csv(catalog, usecols=cols)
    sel = df["r_perp"][df["r_parallel_rsd"] <= rpar_max].to_numpy()
    counts, edges = np.histogram(sel, bins=n_bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    keep = counts > 0
    w = counts[keep].astype(float)
    return centers[keep], w / w.sum(), float(sel.mean()), len(sel)


def band_template(f, axis, seps, weights, s, mode):
    """(central, off-centre) band profiles of the stretched two-halo template.

    Built on the 2-D (X, Y) grid and then band-averaged, and averaged over the
    real distribution of pair separations.
    """
    cen_rows, far_rows = band_masks(axis)
    X, Y = np.meshgrid(axis, axis)
    acc = np.zeros_like(X)
    for a, w in zip(seps, weights):
        for sign in (-1.0, +1.0):
            dx = X + sign * 0.5 * a
            r = np.hypot(s * dx, Y) if mode == "x" else s * np.hypot(dx, Y)
            acc += w * f(r)
    return acc[cen_rows].mean(axis=0), acc[far_rows].mean(axis=0)


def template_bank(f, axis, seps, weights, mode):
    """Templates do not depend on the data, so build them once per (s, mode)."""
    return {s: band_template(f, axis, seps, weights, s, mode) for s in S_GRID}


# --------------------------------------------------------------------------
# fitting
# --------------------------------------------------------------------------
def _solve_A(model, data):
    m = np.isfinite(model) & np.isfinite(data)
    if m.sum() < 3 or not np.any(model[m]):
        return np.nan
    return float(model[m] @ data[m] / (model[m] @ model[m]))


def _rms(model, data, A):
    resid = data - A * model
    resid = resid[np.isfinite(resid)]
    return float(np.sqrt(np.mean(resid ** 2))) if resid.size else np.nan


def exterior_mask(axis, sep, lo, hi):
    return (axis >= 0.5 * sep + lo) & (axis <= hi)


def halo_mask(axis, sep):
    return np.abs(np.abs(axis) - 0.5 * sep) <= HALO_HALF


def fit(p_cen, p_off, bank, axis, sep, variant, buffer=BUFFER, outer=OUTER,
        s_values=None):
    """Return (A, s, ext_rms, holdout_rms, at_edge)."""
    ext = exterior_mask(axis, sep, buffer, outer)
    hold = (axis > HOLDOUT[0]) & (axis <= HOLDOUT[1])
    d_ext = np.concatenate([p_cen[ext], p_off[ext]])
    best = None
    for s in (S_GRID if s_values is None else s_values):
        t_cen, t_off = bank[s]
        m_ext = np.concatenate([t_cen[ext], t_off[ext]])
        if variant == "split":
            hm = halo_mask(axis, sep)
            A = _solve_A(np.concatenate([t_cen[hm], t_off[hm]]),
                         np.concatenate([p_cen[hm], p_off[hm]]))
        else:
            A = _solve_A(m_ext, d_ext)
        if not np.isfinite(A):
            continue
        score = _rms(m_ext, d_ext, A)
        if best is None or score < best[2]:
            best = (A, s, score, t_cen, t_off)
    if best is None:
        return np.nan, np.nan, np.nan, np.nan, False
    A, s, ext_rms, t_cen, t_off = best
    m_h = np.concatenate([t_cen[hold], t_off[hold]])
    d_h = np.concatenate([p_cen[hold], p_off[hold]])
    grid = S_GRID if s_values is None else np.asarray(s_values)
    edge = bool(len(grid) > 1 and (s <= grid[0] or s >= grid[-1]))
    return A, s, ext_rms, _rms(m_h, d_h, A), edge


def residual_stats(p_cen, p_off, bank, axis, sep, A, s):
    t_cen, t_off = bank[s]
    diff = (p_cen - A * t_cen) - (p_off - A * t_off)
    bridge = np.abs(axis) <= 0.25 * sep
    return (float(np.nanmean(diff[bridge])), float(np.nanmean(diff[halo_mask(axis, sep)])), diff)


def evaluate(pair_map, bank, axis, sep, variant, **kw):
    p_cen, p_off = band_profiles(pair_map, axis)
    A, s, ext_rms, hold_rms, edge = fit(p_cen, p_off, bank, axis, sep, variant, **kw)
    if not np.isfinite(A):
        return dict(A=np.nan, s=np.nan, ext_rms=np.nan, hold_rms=np.nan,
                    bridge=np.nan, halo=np.nan,
                    diff=np.full_like(axis, np.nan), at_edge=edge)
    bridge, halo, diff = residual_stats(p_cen, p_off, bank, axis, sep, A, s)
    return dict(A=A, s=s, ext_rms=ext_rms * 1e4, hold_rms=hold_rms * 1e4,
                bridge=bridge * 1e4, halo=halo * 1e4, diff=diff * 1e4,
                at_edge=edge)


def evaluate_fixed(pair_map, bank, axis, sep, A, s):
    """Apply a control with parameters carried in, not refitted."""
    p_cen, p_off = band_profiles(pair_map, axis)
    bridge, halo, _ = residual_stats(p_cen, p_off, bank, axis, sep, A, s)
    return dict(bridge=bridge * 1e4, halo=halo * 1e4)


def jackknife(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    k = len(v)
    return float(np.sqrt((k - 1) / k * np.sum((v - v.mean()) ** 2))) if k > 1 else np.nan


def maps_from_blocks(path):
    d = np.load(path, allow_pickle=True)
    sums, counts = d["sums"], d["counts"]
    tot_s, tot_c = sums.sum(axis=0), counts.sum()
    full = tot_s / tot_c
    loo = np.array([(tot_s - sums[b]) / (tot_c - counts[b]) for b in range(len(counts))])
    return full, loo, int(d["n_pairs"]), float(d["rperp_center"])


# --------------------------------------------------------------------------
# variants
# --------------------------------------------------------------------------
# name -> (stretch mode, fit style, s values; None = free over S_GRID)
VARIANTS = {
    # "_nominal" builds the template at the single nominal separation instead
    # of averaging over the bin.  Comparing joint against joint_nominal is the
    # test of whether a free s is doing physics or just absorbing bin width.
    "unscaled_nominal": ("x_nominal", "joint", np.array([1.0])),
    "unscaled":         ("x", "joint", np.array([1.0])),
    "aonly":            ("x", "joint", np.array([1.0])),
    "joint_nominal":    ("x_nominal", "joint", None),
    "joint":            ("x", "joint", None),
    "split":            ("x", "split", None),
    "radial":           ("radial", "joint", None),
}


def run_one(pair_map, banks, axis, sep, name, **kw):
    mode, style, sv = VARIANTS[name]
    if name.startswith("unscaled"):
        p_cen, p_off = band_profiles(pair_map, axis)
        b, h, d = residual_stats(p_cen, p_off, banks[mode], axis, sep, 1.0, 1.0)
        return dict(A=1.0, s=1.0, ext_rms=np.nan, hold_rms=np.nan,
                    bridge=b * 1e4, halo=h * 1e4, diff=d * 1e4, at_edge=False)
    return evaluate(pair_map, banks[mode], axis, sep, style, s_values=sv, **kw)


def analyse(full_map, loo_maps, banks, axis, sep, label, rows, sens_rows):
    for name in VARIANTS:
        mode, style, sv = VARIANTS[name]
        full = run_one(full_map, banks, axis, sep, name)

        # refit inside every realization -- the honest error
        loo = [run_one(m, banks, axis, sep, name) for m in loo_maps]
        A_k = np.array([r["A"] for r in loo])
        s_k = np.array([r["s"] for r in loo])
        b_k = np.array([r["bridge"] for r in loo])
        h_k = np.array([r["halo"] for r in loo])

        # same maps, parameters carried in from the full sample -- the error
        # you would report if you forgot to refit
        fixed = [evaluate_fixed(m, banks[mode], axis, sep, full["A"], full["s"])
                 for m in loo_maps] if np.isfinite(full["A"]) else []
        b_fix = np.array([r["bridge"] for r in fixed]) if fixed else np.array([np.nan])

        # A bridge mean of zero can hide equal-and-opposite structure, so
        # score the worst single bin inside |X| <= sep as well as the mean.
        prof = np.array([r["diff"] for r in loo])
        per_bin = np.array([jackknife(prof[:, i]) for i in range(len(axis))])
        with np.errstate(divide="ignore", invalid="ignore"):
            sig_all = np.abs(full["diff"]) / per_bin

        def worst(mask):
            v = sig_all[mask]
            return float(np.nanmax(v)) if np.isfinite(v).any() else np.nan

        # |X| <= a spans the halo peaks too, so it mostly restates the halo
        # residual; the bridge-only version is the one the null cares about.
        struct = worst(np.abs(axis) <= sep)
        struct_bridge = worst(np.abs(axis) <= 0.25 * sep)

        ok = np.isfinite(A_k) & np.isfinite(s_k)
        corr = (float(np.corrcoef(A_k[ok], s_k[ok])[0, 1])
                if ok.sum() > 2 and np.std(s_k[ok]) > 0 else np.nan)

        rows.append(dict(
            sample=label, separation=sep, variant=name,
            A=full["A"], s=full["s"],
            A_scatter=jackknife(A_k), s_scatter=jackknife(s_k), A_s_corr=corr,
            ext_rms=full["ext_rms"], hold_rms=full["hold_rms"],
            bridge=full["bridge"], bridge_err=jackknife(b_k),
            bridge_err_fixed=jackknife(b_fix),
            halo=full["halo"], halo_err=jackknife(h_k),
            struct_max_sigma=struct, struct_bridge_sigma=struct_bridge,
            at_edge=full["at_edge"], n_edge_loo=int(sum(r["at_edge"] for r in loo)),
        ))
        logger.info(
            "%-14s %-8s A=%6.3f+/-%.3f  s=%6.3f+/-%.3f  r(A,s)=%+.2f  "
            "ext=%6.3f hold=%6.3f  bridge=%+7.3f+/-%5.3f  halo=%+7.3f+/-%5.3f",
            label, name, full["A"], jackknife(A_k), full["s"], jackknife(s_k),
            corr, full["ext_rms"], full["hold_rms"], full["bridge"],
            jackknife(b_k), full["halo"], jackknife(h_k))
        logger.info("%-14s %-8s   worst bin: bridge %.2f sigma, |X|<=a %.2f sigma",
                    label, name, struct_bridge, struct)

        if sens_rows is not None and not name.startswith("unscaled"):
            for buf in (2.0, 4.0, 6.0):
                for outer in (20.0, 25.0, 30.0):
                    r = run_one(full_map, banks, axis, sep, name,
                                buffer=buf, outer=outer)
                    sens_rows.append(dict(sample=label, separation=sep, variant=name,
                                          buffer=buf, outer=outer,
                                          **{k: v for k, v in r.items() if k != "diff"}))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Zheng's stretched single control.")
    ap.add_argument("--null", action="store_true",
                    help="Stage 2: the non-physical-pair null instead of Stage 1.")
    ap.add_argument("--single", default=SINGLE)
    ap.add_argument("--out-tag", default="")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    f, r_max = radial_spline(load_map(SIM / args.single))
    logger.info("Radial profile splined to r = %.1f h^-1 Mpc", r_max)
    axis = axis_for(101)
    rows, sens_rows = [], []

    jobs = ([(b, HOD / f"nonphys_rperp10_rpar{b}_blocks.npz",
              HOD / "pairs_hodnmatch50_rperp10_real120.csv", "real",
              tuple(float(v) for v in b.split("_"))) for b in NULL_BANDS]
            if args.null else
            [(k, HOD / v[0], HOD / v[1], "rsd", (None, 10.0))
             for k, v in PHYSICAL.items()])

    for label, blocks, catalog, space, cut in jobs:
        if not blocks.exists():
            logger.warning("missing %s -- skipped", blocks)
            continue
        full, loo, n_pairs, sep = maps_from_blocks(blocks)
        full = reflect_symmetrize_map(full.astype(np.float64))
        loo = [reflect_symmetrize_map(m.astype(np.float64)) for m in loo]

        col = "r_parallel_real" if space == "real" else "r_parallel_rsd"
        df = pd.read_csv(catalog, usecols=["r_perp", col])
        sel = df["r_perp"][(df[col] >= (cut[0] if cut[0] is not None else -np.inf))
                           & (df[col] <= cut[1])].to_numpy()
        counts, edges = np.histogram(sel, bins=21)
        centers = 0.5 * (edges[:-1] + edges[1:])
        keep = counts > 0
        seps, wts = centers[keep], counts[keep] / counts[keep].sum()
        logger.info("%s: %d pairs, nominal a=%.0f, actual <a>=%.3f, spread %.3f",
                    label, n_pairs, sep, float(sel.mean()), float(sel.std()))

        banks = {m: template_bank(f, axis, seps, wts, m) for m in ("x", "radial")}
        banks["x_nominal"] = template_bank(f, axis, np.array([sep]),
                                           np.array([1.0]), "x")
        analyse(full, loo, banks, axis, sep, label, rows, sens_rows)

    stem = "stretched_control_null" if args.null else "stretched_control_physical"
    stem += args.out_tag
    pd.DataFrame(rows).to_csv(SIM / f"{stem}.csv", index=False)
    logger.info("Saved -> %s", SIM / f"{stem}.csv")
    if sens_rows:
        pd.DataFrame(sens_rows).to_csv(SIM / f"{stem}_sensitivity.csv", index=False)
        logger.info("Saved -> %s", SIM / f"{stem}_sensitivity.csv")


if __name__ == "__main__":
    main()
