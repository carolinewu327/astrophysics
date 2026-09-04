#!/usr/bin/env python
"""Recover the single-halo profile from a pair stack by deconvolving the companion.

Zheng's proposal (2026-08-13).  Centre on one halo and look *outward*, away from
the companion.  What you measure there is the halo itself plus the companion a
separation away:

    g(x) = f(x) + f(x + a)

Both galaxies share the ensemble-average profile f, so the companion's
contamination at x is the same number as the halo's own signal at x + a -- which
is another measurement.  Rearranging and substituting repeatedly,

    f(x) = g(x) - f(x + a)
         = g(x) - g(x + a) + f(x + 2a)
         = g(x) - g(x + a) + g(x + 2a) - g(x + 3a) + ...

i.e.  f(x) = sum_{k=0}^{N} (-1)^k g(x + k a)  +  (-1)^{N+1} f(x + (N+1) a),

and truncating at order N means asserting the remainder is negligible.  Note the
alternation runs over *every* multiple of a; Zheng's email has g(x + 2a) in the
second term, which appears to be an indexing slip -- at a = 10, x = 5 the two
readings differ by 31%.

With f in hand the control is an honest superposition, because f is a single
halo rather than something that already contains the pair:

    b(X) = f(X - a/2) + f(X + a/2)

Three things this file is careful about:

**Never zero-fill a missing term.**  If g(x + k a) falls off the map it is
unknown, not zero.  Zero-filling silently lowers the effective order at large x
and stamps a boundary artifact into f.  Here every order carries its validity
range x <= x_max - N a, and everything outside is NaN.

**Both bands.**  The kappa_band fit uses the central and off-centre profiles
jointly, so reconstructing only the central band cannot feed it.  The recursion
is linear and applies to each band separately.

**Both exterior sides.**  g is averaged over the +X and -X exteriors.  For
reflection-symmetrized pair stacks this is a no-op (the two agree to ~1e-10,
which is CSV round-off), but it removes an arbitrary choice of side and is
correct if the input is ever unsymmetrized.

What the interior residual means, stated plainly: the exterior fixes f, and
evenness f(-x) = f(x) carries it inside.  So the "filament" is *the part of the
interior that two even profiles inferred from the exterior cannot explain*.  The
exterior agreement is NOT independent evidence for that model -- summing the
series fully makes f(X-a/2) + f(X+a/2) = g(X-a/2) identically, so agreement
there tests truncation, boundaries and interpolation, not the physics.  Only a
sample with no bridge (large-LOS non-physical pairs) tests the estimator.

Usage
-----
    PYTHONPATH=lib:analysis/sim python analysis/sim/deconvolve_pair_profile.py
    PYTHONPATH=lib:analysis/sim python analysis/sim/deconvolve_pair_profile.py \\
        --separations 10 --orders 0,1,2,3
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from geometry import stack_axis, band_masks, band_profiles
from sim_utils import setup_logging, symmetrized_radial_interpolator

logger = logging.getLogger(__name__)

OBS_BOX_HMPC = 100.0

PAIR_STACK = {
    # At separation 5 the endpoint work and the current analysis both use
    # r_par <= 10; the archived r_par <= 5 stack is a different sample.
    "5": "stack_rperp5_rpar10_matched",
    "10": "stack_rperp10_matched",
    "20": "stack_rperp20_matched",
}

INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e3e2df"
BLUE, ORANGE, GREEN, PURPLE = "#2a78d6", "#eb6834", "#1f9d55", "#7d5ba6"

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 7.5,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def load_map(path) -> np.ndarray:
    return pd.read_csv(path, index_col=0).to_numpy(dtype=float)


def axis_for(n: int) -> np.ndarray:
    """Stacking axis for an n-pixel map -- see geometry.stack_axis."""
    return stack_axis(OBS_BOX_HMPC, n)[0]


def tidy(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


# ---------------------------------------------------------------------------
# The recursion
# ---------------------------------------------------------------------------

def exterior_profile(band: np.ndarray, axis: np.ndarray, sep: float,
                     x: np.ndarray) -> np.ndarray:
    """g(x): the band profile centred on a halo, looking outward, x >= 0.

    Averaged over both exteriors -- +a/2 + x and -a/2 - x -- so the result does
    not depend on which halo is called "the right one".
    """
    right = np.interp(0.5 * sep + x, axis, band, left=np.nan, right=np.nan)
    left = np.interp(-0.5 * sep - x, axis, band, left=np.nan, right=np.nan)
    return 0.5 * (right + left)


def deconvolve(g: np.ndarray, x: np.ndarray, sep: float, order: int,
               x_max: float) -> np.ndarray:
    """f(x) = sum_{k=0}^{order} (-1)^k g(x + k*sep), NaN where unavailable.

    A term whose argument leaves the measured range is *unknown*, so the whole
    estimate at that x is NaN rather than being quietly truncated to a lower
    order.  That keeps the effective order constant across x.
    """
    out = np.zeros_like(x, dtype=float)
    for k in range(order + 1):
        shifted = np.interp(x + k * sep, x, g, left=np.nan, right=np.nan)
        out = out + ((-1) ** k) * shifted
    out[x + order * sep > x_max] = np.nan
    return out


def even_extend(f: np.ndarray, x: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Evaluate f at arbitrary (signed) positions using f(-x) = f(x).

    This is the load-bearing assumption: the exterior-inferred profile is
    asserted to hold on the interior side of each halo too.
    """
    ok = np.isfinite(f)
    if not ok.any():
        return np.full_like(query, np.nan, dtype=float)
    return np.interp(np.abs(query), x[ok], f[ok], left=np.nan, right=np.nan)


def band_from_radial(profile, axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Central and off-centre band profiles of a circular map, on `axis`.

    Lets a radially defined reference (the field single) be compared with the
    reconstructed f on exactly the same footing as the pair stack's bands.
    """
    central_rows, far_rows = band_masks(axis)
    gx, gy = np.meshgrid(axis, axis)
    vals = profile(np.hypot(gx, gy))
    return vals[central_rows, :].mean(axis=0), vals[far_rows, :].mean(axis=0)


# ---------------------------------------------------------------------------

def run_separation(key: str, sim: Path, single_name: str, orders: list[int],
                   out_dir: Path) -> list[dict]:
    sep = float(key)
    pair = load_map(sim / "hod_pairs" / f"{PAIR_STACK[key]}.csv")
    axis = axis_for(pair.shape[0])
    p_cen, p_off = band_profiles(pair, axis)

    step = float(np.diff(axis).mean())
    x_max = 0.5 * OBS_BOX_HMPC - 0.5 * sep          # furthest measured exterior
    x = np.arange(0.0, x_max + 0.5 * step, step)

    g_cen = exterior_profile(p_cen, axis, sep, x)
    g_off = exterior_profile(p_off, axis, sep, x)

    field = load_map(sim / single_name)
    f_ref_cen, f_ref_off = band_from_radial(
        symmetrized_radial_interpolator(field), axis)
    ref_cen = np.interp(x, axis[axis >= 0], f_ref_cen[axis >= 0])

    interior = np.abs(axis) <= 0.5 * sep
    exterior = (np.abs(axis) > 0.5 * sep) & (np.abs(axis) <= 1.6 * sep)

    rows, curves = [], {}
    logger.info("=== separation %g  (x measured to %.1f, step %.1f) ===",
                sep, x_max, step)
    for N in orders:
        reach = x_max - N * sep
        if reach < sep:
            logger.info("  N=%d: needs f out to x=%g but only %g is available "
                        "-- skipped", N, sep, max(reach, 0.0))
            continue

        f_cen = deconvolve(g_cen, x, sep, N, x_max)
        f_off = deconvolve(g_off, x, sep, N, x_max)

        ctrl_cen = (even_extend(f_cen, x, axis - 0.5 * sep)
                    + even_extend(f_cen, x, axis + 0.5 * sep))
        ctrl_off = (even_extend(f_off, x, axis - 0.5 * sep)
                    + even_extend(f_off, x, axis + 0.5 * sep))

        res_cen = p_cen - ctrl_cen
        res_off = p_off - ctrl_off
        band_res = res_cen - res_off

        # The exterior check needs f out to |X| + a/2, so at high order it can
        # go entirely undefined.  Report the coverage rather than letting a
        # silent NaN read as "no residual".
        ext_ok = exterior & np.isfinite(res_cen)
        cover = ext_ok.sum() / max(exterior.sum(), 1)
        ext = np.abs(res_cen[ext_ok]).mean() if ext_ok.any() else np.nan

        bridge = np.abs(axis) <= 0.25 * sep
        rows.append({
            "separation_hmpc": sep, "order": N, "x_reach": reach,
            "f_peak": f_cen[0] * 1e4,
            "f_peak_over_field": f_cen[0] / ref_cen[0],
            "exterior_abs_residual": ext * 1e4,
            "exterior_coverage": cover,
            "interior_central": np.nanmean(res_cen[interior]) * 1e4,
            "bridge_central": np.nanmean(res_cen[bridge]) * 1e4,
            "bridge_band_diff": np.nanmean(band_res[bridge]) * 1e4,
        })
        curves[N] = (f_cen, f_off, res_cen, band_res)
        logger.info("  N=%d  reach x<=%5.1f   f(0)=%7.3fe-4 (%.3f x field)   "
                    "|ext resid|=%7.4fe-4 (cover %3.0f%%)   bridge(central)=%+7.3fe-4"
                    "   bridge(c-off)=%+7.3fe-4",
                    N, reach, f_cen[0] * 1e4, f_cen[0] / ref_cen[0], ext * 1e4,
                    100 * cover, np.nanmean(res_cen[bridge]) * 1e4,
                    np.nanmean(band_res[bridge]) * 1e4)

    # An alternating series with decreasing terms brackets its own limit:
    # consecutive partial sums straddle it.  So the midpoint of orders N and
    # N+1 is a far better estimate than either, and half their gap is an honest
    # truncation uncertainty -- which matters here because the raw orders swing
    # by tens of 1e-4 while the answer sits near zero.
    for lo, hi in zip(rows[:-1], rows[1:]):
        if hi["order"] != lo["order"] + 1:
            continue
        for field_name in ("bridge_central", "bridge_band_diff", "f_peak"):
            mid = 0.5 * (lo[field_name] + hi[field_name])
            half = 0.5 * abs(hi[field_name] - lo[field_name])
            hi[f"{field_name}_paired"] = mid
            hi[f"{field_name}_trunc_err"] = half
        logger.info("  N=%d&%d paired:  f(0)=%7.3f+/-%.3f   bridge(central)="
                    "%+7.3f+/-%.3f   bridge(c-off)=%+7.3f+/-%.3f  [e-4]",
                    lo["order"], hi["order"],
                    hi["f_peak_paired"], hi["f_peak_trunc_err"],
                    hi["bridge_central_paired"], hi["bridge_central_trunc_err"],
                    hi["bridge_band_diff_paired"], hi["bridge_band_diff_trunc_err"])

    if curves:
        make_figure(key, sep, axis, x, g_cen, ref_cen, curves, p_cen, out_dir)
    return rows


def make_figure(key, sep, axis, x, g_cen, ref_cen, curves, p_cen, out_dir):
    orders = sorted(curves)
    shades = plt.cm.viridis(np.linspace(0.15, 0.85, len(orders)))

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3), constrained_layout=True)

    ax0 = axes[0]
    ax0.plot(x, g_cen * 1e4, lw=2.2, color=INK, label="measured $g(x)$ (N=0)")
    for c, N in zip(shades, orders):
        if N == 0:
            continue
        ax0.plot(x, curves[N][0] * 1e4, lw=1.8, color=c, label=f"$f(x)$, N={N}")
    ax0.plot(x, ref_cen * 1e4, lw=1.6, color=ORANGE, ls="--",
             label="field single (reference)")
    ax0.set_yscale("log")
    ax0.set_xlim(0, min(3.0 * sep, x.max()))
    ax0.set_xlabel(r"$x$ from the halo, outward  [$h^{-1}$ Mpc]")
    ax0.set_ylabel(r"central-band $\kappa$  [$10^{-4}$]")
    ax0.set_title("(a) Recovered single-halo profile", loc="left", color=INK, pad=5)
    ax0.legend(loc="upper right")
    tidy(ax0)

    ax1 = axes[1]
    ax1.plot(axis, p_cen * 1e4, lw=2.2, color=INK, label="measured pair profile")
    for c, N in zip(shades, orders):
        ctrl = p_cen - curves[N][2]
        ax1.plot(axis, ctrl * 1e4, lw=1.6, color=c, ls="--", label=f"control, N={N}")
    for xx in (-0.5 * sep, 0.5 * sep):
        ax1.axvline(xx, color=GREEN, lw=1.0, ls="--", alpha=0.8)
    ax1.set_xlim(-2.2 * sep, 2.2 * sep)
    ax1.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
    ax1.set_ylabel(r"central-band $\kappa$  [$10^{-4}$]")
    ax1.set_title("(b) Reconstructed control vs the data", loc="left", color=INK, pad=5)
    ax1.legend(loc="upper right")
    tidy(ax1)

    # Plot the statistic we actually quote -- central MINUS off-centre band --
    # as the solid curve, with the central band alone dashed behind it.  These
    # differ a lot (at r_perp = 20, N=1: +1.65 vs +3.13), and showing only the
    # central band while quoting the difference invites exactly the mismatch
    # that was queried.
    ax2 = axes[2]
    for c, N in zip(shades, orders):
        ax2.plot(axis, curves[N][3] * 1e4, lw=2.0, color=c,
                 label=f"N={N}  (central $-$ off)")
        ax2.plot(axis, curves[N][2] * 1e4, lw=1.1, color=c, ls="--", alpha=0.55)
    ax2.axhline(0, color=INK, lw=1.0)
    for xx in (-0.5 * sep, 0.5 * sep):
        ax2.axvline(xx, color=GREEN, lw=1.0, ls="--", alpha=0.8)
    ax2.axvspan(-0.25 * sep, 0.25 * sep, color=INK_MUTED, alpha=0.07, lw=0)
    ax2.set_xlim(-1.8 * sep, 1.8 * sep)
    ax2.set_xlabel(r"$X$ along the pair axis  [$h^{-1}$ Mpc]")
    ax2.set_ylabel(r"residual $\kappa$  [$10^{-4}$]")
    ax2.set_title("(c) Interior residual — solid = central $-$ off band (quoted); "
                  "dashed = central only", loc="left", color=INK, pad=5)
    ax2.legend(loc="upper right")
    tidy(ax2)

    fig.suptitle(
        f"Companion-deconvolved control — simulation, $r_\\perp$ = {key} "
        f"$h^{{-1}}$Mpc\n"
        r"$f(x)=g(x)-g(x{+}a)+g(x{+}2a)-\cdots$, then control "
        r"$=f(X{-}a/2)+f(X{+}a/2)$.  Green marks the halos; shading marks the bridge."
        "\nOutside the halos the control is forced toward the data by the same "
        "recursion, so agreement there tests truncation, not the halo model.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    base = out_dir / f"deconvolved_control_sep{key}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("  saved %s.png / .pdf", base)


def deconvolve_maps(key: str, sim: Path, order: int, out_dir: Path) -> dict:
    """Row-wise extension of the recursion, to produce maps rather than profiles.

    The recursion g(x) = f(x) + f(x+a) is a statement about translation along X,
    exact only on the pair axis.  Applying it independently to each Y row is the
    natural 2-D extension and gives directly comparable heatmaps, but it is NOT
    a radial deconvolution: nothing forces the recovered f to be circular.

    That is worth having rather than assuming away -- a correct reconstruction
    should come out close to circular on its own, so the departure is a free
    diagnostic.  A genuinely radial inversion is a harder problem and is not
    attempted here.
    """
    sep = float(key)
    pair = load_map(sim / "hod_pairs" / f"{PAIR_STACK[key]}.csv")
    axis = axis_for(pair.shape[0])
    step = float(np.diff(axis).mean())
    x_max = 0.5 * OBS_BOX_HMPC - 0.5 * sep
    x = np.arange(0.0, x_max + 0.5 * step, step)

    f_map = np.full_like(pair, np.nan, dtype=float)
    ctrl = np.full_like(pair, np.nan, dtype=float)
    for j in range(pair.shape[0]):
        row = pair[j, :]
        g_row = exterior_profile(row, axis, sep, x)
        f_row = deconvolve(g_row, x, sep, order, x_max)
        f_map[j, :] = even_extend(f_row, x, axis)
        ctrl[j, :] = (even_extend(f_row, x, axis - 0.5 * sep)
                      + even_extend(f_row, x, axis + 0.5 * sep))

    fil = pair - ctrl
    for arr, stem in ((f_map, "kappa_deconvolved_single_sim"),
                      (ctrl, "kappa_deconvolved_control_sim"),
                      (fil, "kappa_deconvolved_filament_sim")):
        pd.DataFrame(arr, index=pd.Index(axis, name="y"), columns=axis).to_csv(
            sim / f"{stem}_rperp{key}_N{order}.csv")
    return {"pair": pair, "f": f_map, "control": ctrl, "filament": fil,
            "axis": axis, "sep": sep, "order": order}


def make_map_figure(d: dict, key: str, out_dir: Path) -> None:
    from matplotlib.colors import LogNorm, TwoSlopeNorm

    axis, sep = d["axis"], d["sep"]
    z = min(2.5 * sep, 40.0)
    ext = (axis[0], axis[-1], axis[0], axis[-1])
    gx, gy = np.meshgrid(axis, axis)
    win = (np.abs(gx) <= z) & (np.abs(gy) <= z)

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.4), constrained_layout=True)
    panels = ((d["f"], r"(a) deconvolved single $f$", False),
              (d["control"], r"(b) control $f(X{-}a/2)+f(X{+}a/2)$", False),
              (d["filament"], "(c) filament = data $-$ control", True))

    for col, (arr, ttl, div) in enumerate(panels):
        ax = axes[0, col]
        vals = arr[win & np.isfinite(arr)]
        if div:
            s = float(np.nanpercentile(np.abs(vals), 99.5)) * 1e4
            im = ax.imshow(arr * 1e4, origin="lower", extent=ext, cmap="RdBu_r",
                           norm=TwoSlopeNorm(vmin=-s, vcenter=0.0, vmax=s),
                           interpolation="nearest")
        else:
            lo = max(float(np.nanpercentile(vals, 2)), 1e-7)
            im = ax.imshow(arr * 1e4, origin="lower", extent=ext, cmap="magma_r",
                           norm=LogNorm(vmin=lo * 1e4, vmax=float(vals.max()) * 1e4),
                           interpolation="nearest")
        for xx in (-0.5 * sep, 0.5 * sep):
            ax.plot(xx, 0.0, "+", color=GREEN, ms=10, mew=1.3)
        ax.set_xlim(-z, z); ax.set_ylim(-z, z)
        ax.set_xlabel(r"$X$  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")
        ax.set_title(ttl, loc="left", color=INK, pad=5)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(r"$\kappa$  [$10^{-4}$]"); cb.ax.tick_params(labelsize=7)

        # Ring width matches the pixel pitch: narrower rings leave annuli with
        # no pixel centres at all, which shows up as fake gaps in the curve.
        # And a ring that runs outside the defined region is not a radial
        # average -- it samples only the surviving wedge near X = 0 -- so stop
        # where coverage falls below 90% rather than plotting a wedge mean.
        axp = axes[1, col]
        rr = np.hypot(gx, gy).ravel()
        v = arr.ravel()
        pitch = float(np.diff(axis).mean())
        edges = np.arange(0.0, min(30.0, z) + pitch, pitch)
        idx = np.digitize(rr, edges) - 1
        prof, cover = [], []
        for i in range(len(edges) - 1):
            sel = idx == i
            n = int(sel.sum())
            ok = np.isfinite(v[sel])
            cover.append(ok.mean() if n else 0.0)
            prof.append(v[sel][ok].mean() if ok.any() else np.nan)
        prof = np.array(prof)
        prof[np.array(cover) < 0.9] = np.nan
        r = 0.5 * (edges[:-1] + edges[1:])
        axp.plot(r, prof * 1e4, lw=2, color=PURPLE)
        good = np.isfinite(prof)
        if good.any():
            axp.set_xlim(0, r[good].max() * 1.05)
        if div:
            axp.axhline(0, color=INK, lw=0.9)
        elif np.nanmin(prof) > 0:
            axp.set_yscale("log")
        axp.axvline(0.5 * sep, color=GREEN, lw=1.0, ls="--", alpha=0.8)
        axp.set_xlabel(r"$r$ from map centre  [$h^{-1}$ Mpc]")
        axp.set_ylabel(r"$\kappa$  [$10^{-4}$]")
        axp.set_title(f"({'def'[col]}) radial profile", loc="left", color=INK, pad=5)
        tidy(axp)

    # How circular did f come out?  Nothing in the row-wise recursion enforces it.
    inner = np.hypot(gx, gy) <= 1.5 * sep
    fv = d["f"][inner & np.isfinite(d["f"])]
    fig.suptitle(
        f"Companion-deconvolved control, maps — simulation, $r_\\perp$ = {key} "
        f"$h^{{-1}}$Mpc, order N = {d['order']}\n"
        "The recursion is applied to each $Y$ row independently, so this is a "
        "translation-based extension, not a radial deconvolution — nothing forces "
        "$f$ to be circular.\n"
        "Blank margins are where a series term leaves the measured map; those are "
        "unknown, not zero. Radial profiles stop where a ring is less than 90% "
        "inside the defined region.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    base = out_dir / f"deconvolved_maps_sep{key}"
    for e in ("png", "pdf"):
        fig.savefig(f"{base}.{e}", bbox_inches="tight")
    plt.close(fig)
    logger.info("  saved %s.png / .pdf  (f over r<1.5a: %.3f-%.3f e-4)",
                base, fv.min() * 1e4, fv.max() * 1e4)


def main(argv=None):
    p = argparse.ArgumentParser(description="Deconvolve the companion from a pair stack.")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--single", default="kappa_single_sim_hodnmatch_8arcmin_centered.csv")
    p.add_argument("--separations", default="5,10,20")
    p.add_argument("--orders", default="0,1,2,3,4,5,6,7,8")
    p.add_argument("--output-dir", default="output/plots")
    p.add_argument("--maps", action="store_true",
                   help="Also run the row-wise 2-D extension and save maps.")
    args = p.parse_args(argv)
    setup_logging()

    sim = Path(args.sim_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    orders = [int(o) for o in args.orders.split(",")]

    # Best usable order per separation: the highest whose control still covers
    # the interior.  Set by the map size, not by convergence.
    map_order = {"5": 8, "10": 3, "20": 1}

    rows = []
    for key in [s.strip() for s in args.separations.split(",")]:
        rows.extend(run_separation(key, sim, args.single, orders, out_dir))
        if args.maps:
            d = deconvolve_maps(key, sim, map_order.get(key, 1), out_dir)
            make_map_figure(d, key, out_dir)

    table = pd.DataFrame(rows)
    path = sim / "deconvolved_control_summary.csv"
    table.to_csv(path, index=False)
    logger.info("\n%s", table.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    logger.info("Saved -> %s", path)


if __name__ == "__main__":
    main()
