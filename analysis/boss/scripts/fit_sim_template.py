#!/usr/bin/env python
"""Fit the simulation filament map to BOSS as a template with one free amplitude.

"Lensing is low" -- the measured lensing amplitude around BOSS galaxies runs
below a matched LCDM+HOD prediction -- means a direct comparison of bridge
numbers conflates two questions: is the pattern there (shape), and is the
normalisation right (amplitude).  Fitting

    chi^2(A) = sum_ij (obs_ij - A * sim_ij)^2 / sigma_ij^2

separates them.  A absorbs the normalisation and is itself a measurement of the
amplitude deficit for pairs; the residual tests the shape.  Because A enters
linearly the solution is closed-form, no minimiser:

    A_hat  = sum(obs*sim/sigma^2) / sum(sim^2/sigma^2)
    dchi2  = chi^2(A=0) - chi^2(A_hat) = A_hat^2 * sum(sim^2/sigma^2)

Three fits are reported per separation:

  full_grid        every valid pixel -- the literal calculation requested
  unique_quadrant  one quadrant only -- see below
  bridge_quadrant  the bridge band within that quadrant -- see below

**Why the quadrant matters.**  Both maps are reflection-symmetrized, so the
101x101 grid carries only 51x51 independent values replicated four times.
Summing chi^2 over the full grid counts each one ~4 times, inflating dchi2 by
~4x and the significance by ~2x.  A_hat is unaffected -- it is a weighted mean
and both the map and its sigma carry the same symmetry -- but dchi2 is not.

**Why dchi2 is still not the answer.**  Even on one quadrant the pixels are
correlated: the 8 arcmin beam is ~3.3 h^-1 Mpc against 1 h^-1 Mpc pixels, so
neighbours are not independent and the diagonal sigma understates the error on
A.  This is the same failure mode as assuming independent radial bins, which
understated profile errors by 2.2-2.6x.  The trustworthy error is the jackknife:
refit A on each of the 287 leave-one-out maps (with the *full-sample* weights
held fixed, so no spurious variance is injected) and take the delete-one
scatter.  That absorbs both the symmetry duplication and the beam correlations.

    Report A/sigma_A_JK as the result; dchi2 is the requested cross-check.

**Why a bridge-restricted fit.**  A_hat weights by the template's own amplitude,
and the largest values in the simulation filament map sit at the halo positions.
The control subtraction is meant to remove those but demonstrably does not
fully -- BOSS peaks near 20e-4 against the simulation's 9e-4, an open item.  So
a whole-map A may be measuring that halo-peak mismatch rather than the bridge.
Fitting a second A over the bridge band alone tests this: agreement means the
whole-map number is safe, divergence means it is not measuring the filament.

**Naming.**  dchi2 measures agreement with the simulation template as a whole,
which contains halo peaks and bridge as well as the quadrupole.  It is a
matched-pattern detection, not a quadrupole detection per se.

**Known limitation.**  The template carries its own noise, which biases A low
(regression dilution).  Quantifying it needs pixel-level simulation leave-one-out
maps, which are not persisted -- jackknife_pair_stack.py keeps only the scalar
statistics.  Given the simulation's ~20 sigma bridge detection this is a
secondary systematic, recorded rather than corrected.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/fit_sim_template.py
"""

from __future__ import annotations

import argparse
import logging
import os

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from catalog import setup_logging
from geometry import BRIDGE_HALF_X_FRAC, OFF_HI_Y_HMPC, band_profile
from jackknife import jackknife_error
from plot_obs_vs_sim_filament import sim_filament
from plot_separation_summary import build_terms
from weighting_sensitivity import SEPARATIONS

logger = logging.getLogger(__name__)

INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e3e2df"
BLUE = "#2a78d6"

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def unique_quadrant_mask(shape: tuple[int, int]) -> np.ndarray:
    """One quadrant of a reflection-symmetrized square grid, centre included.

    reflect_symmetrize_map averages (+i,+j), (-i,+j), (+i,-j), (-i,-j) about the
    centre pixel, so for an odd grid the independent values are exactly the
    i >= centre, j >= centre block.
    """
    n, m = shape
    if n != m or n % 2 == 0:
        raise ValueError(f"Expected an odd square grid, got {shape}.")
    c = n // 2
    mask = np.zeros(shape, dtype=bool)
    mask[c:, c:] = True
    return mask


def bridge_mask(axis: np.ndarray, sep: float) -> np.ndarray:
    """The band the filament statistic lives in: bridge in x, off-centre in y.

    Excludes the halo peaks at x = +/- sep/2, since BRIDGE_HALF_X_FRAC = 0.35
    is inside 0.5.
    """
    x, y = np.meshgrid(axis, axis)
    return (np.abs(x) <= BRIDGE_HALF_X_FRAC * sep) & (np.abs(y) <= OFF_HI_Y_HMPC)


def fit_amplitude_with_offset(obs: np.ndarray, sim: np.ndarray, inv_var: np.ndarray,
                              mask: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Fit obs = A*sim + c, returning (A, c, the 2x2 normal matrix).

    Diagnostic for a specific confusion.  The bridge statistic is a *difference*
    -- central band minus off-centre band -- so any spatially flat offset cancels
    out of it.  The template fit is not: it works on absolute kappa, so a
    residual monopole from imperfect control subtraction would show up as
    amplitude.  If A survives with a free constant, the fit is responding to the
    template's shape; if A collapses, it was responding to a level.
    """
    o, s, w = obs[mask], sim[mask], inv_var[mask]
    sww, sws, ss = float(np.sum(w)), float(np.sum(w * s)), float(np.sum(w * s * s))
    m = np.array([[ss, sws], [sws, sww]])
    rhs = np.array([float(np.sum(w * o * s)), float(np.sum(w * o))])
    if abs(np.linalg.det(m)) < 1e-30:
        return np.nan, np.nan, m
    a, c = np.linalg.solve(m, rhs)
    return float(a), float(c), m


def fit_amplitude(obs: np.ndarray, sim: np.ndarray, inv_var: np.ndarray,
                  mask: np.ndarray) -> dict:
    """Closed-form single-amplitude fit over *mask*. A is unconstrained."""
    o, s, w = obs[mask], sim[mask], inv_var[mask]
    denom = float(np.sum(s * s * w))
    if denom <= 0:
        raise ValueError("Template carries no weight inside the mask.")
    a = float(np.sum(o * s * w) / denom)
    chi2_0 = float(np.sum(o * o * w))
    chi2_min = float(np.sum((o - a * s) ** 2 * w))
    return {
        "A": a,
        "sigma_A_diag": float(1.0 / np.sqrt(denom)),
        "chi2_A0": chi2_0,
        "chi2_min": chi2_min,
        "delta_chi2": chi2_0 - chi2_min,
        "sqrt_delta_chi2": float(np.sqrt(max(chi2_0 - chi2_min, 0.0))),
        "n_pixels": int(mask.sum()),
        "_denom": denom,
    }


def refit_loo(loo_maps: np.ndarray, shape: tuple[int, int], sim: np.ndarray,
              inv_var: np.ndarray, mask: np.ndarray, denom: float) -> np.ndarray:
    """A for each leave-one-out map, with full-sample weights held fixed.

    Recomputing the weights per realization would let the error map fluctuate
    with the map it weights and inject variance that is not really there.
    """
    s, w = sim[mask], inv_var[mask]
    sw = s * w
    return np.array([float(np.sum(m.reshape(shape)[mask] * sw) / denom)
                     for m in loo_maps])


def fit_profile(obs: np.ndarray, sim: np.ndarray, loo_maps: np.ndarray,
                axis: np.ndarray) -> dict:
    """Fit the template to the 1-D band-differenced profile rather than to pixels.

    The 2-D pixel fit has a structural problem: it weights by the template's own
    amplitude, and the filament map's positive feature is a narrow band along the
    pair axis while the negative quadrupole lobes fill most of the area.  So
    80-98% of the pixel fit's weight lands on the perpendicular deficit, and A
    ends up measuring that rather than the bridge.

    Fitting the band-differenced profile instead -- central rows minus off-centre
    rows, as a function of x -- fixes this without shrinking the region:

      * it is targeted at the filament, since the band difference is what
        isolates the quadrupole contrast in the first place;
      * it is immune to a map-wide additive offset by construction, because a
        difference cancels any constant.  That is what made the small-region
        pixel fits swing by factors of 2-3 depending on whether a constant was
        free;
      * it still spans the whole x range, so it is not a small-region fit.

    Only x >= 0 is used: reflection symmetrization makes the profile even in x,
    so the full range would count every value twice.
    """
    p_obs = band_profile(obs, axis)
    p_sim = band_profile(sim, axis)
    p_loo = np.array([band_profile(m.reshape(obs.shape), axis) for m in loo_maps])
    _, sigma = jackknife_error(p_loo)

    half = (axis >= 0) & np.isfinite(sigma) & (sigma > 0)
    w = np.zeros_like(sigma)
    np.divide(1.0, sigma ** 2, out=w, where=half)

    o, s, ww = p_obs[half], p_sim[half], w[half]
    denom = float(np.sum(s * s * ww))
    a = float(np.sum(o * s * ww) / denom)
    a_loo = np.array([float(np.sum(p[half] * s * ww) / denom) for p in p_loo])
    _, err = jackknife_error(a_loo[:, None])

    chi2_0 = float(np.sum(o * o * ww))
    chi2_min = float(np.sum((o - a * s) ** 2 * ww))
    return {
        "A": a,
        "sigma_A_diag": float(1.0 / np.sqrt(denom)),
        "sigma_A_jk": float(err[0]),
        "chi2_A0": chi2_0,
        "chi2_min": chi2_min,
        "delta_chi2": chi2_0 - chi2_min,
        "sqrt_delta_chi2": float(np.sqrt(max(chi2_0 - chi2_min, 0.0))),
        "n_pixels": int(half.sum()),
    }


def validate(obs: np.ndarray, sim: np.ndarray, sigma: np.ndarray, sep: float) -> np.ndarray:
    """Shape/axis agreement and a finite, positive-sigma pixel mask."""
    if obs.shape != sim.shape:
        raise ValueError(
            f"Grid mismatch at sep={sep}: BOSS {obs.shape} vs simulation {sim.shape}. "
            "Both should be the 101-pixel pair grid over a 100 h^-1 Mpc box.")
    if sigma.shape != obs.shape:
        raise ValueError(f"Error map shape {sigma.shape} does not match {obs.shape}.")
    valid = np.isfinite(obs) & np.isfinite(sim) & np.isfinite(sigma) & (sigma > 0)
    dropped = valid.size - int(valid.sum())
    if dropped:
        logger.warning("  sep=%g: dropped %d of %d pixels (non-finite or sigma<=0)",
                       sep, dropped, valid.size)
    else:
        logger.info("  sep=%g: all %d pixels valid", sep, valid.size)
    return valid


def make_sheet(obs, sim, a_hat, axis, sep, out_base, loo_maps, inv_var, quad):
    """The three requested maps on one shared scale, plus a residual profile."""
    model = a_hat * sim
    resid = obs - model

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.0), constrained_layout=True)
    # Display zooms to the signal region even though the fit spans the whole map.
    # Scaling the colour bar to the full 100 h^-1 Mpc box lets BOSS's outer noise
    # set the range and washes the template out to near-blank.
    zoom = 2.5 * sep
    inner = np.abs(axis) <= zoom
    box = np.ix_(inner, inner)
    scale = max(np.nanpercentile(np.abs(obs[box]), 99.5),
                np.nanpercentile(np.abs(model[box]), 99.5), 1e-12)
    norm = TwoSlopeNorm(vmin=-scale, vcenter=0.0, vmax=scale)

    panels = ((axes[0, 0], obs, r"(a) BOSS filament map  $M_{\rm obs}$"),
              (axes[0, 1], model, rf"(b) Scaled simulation  $\hat{{A}}\,M_{{\rm sim}}$"
                                  rf"   ($\hat{{A}}$ = {a_hat:.2f})"),
              (axes[1, 0], resid, r"(c) Residual  $M_{\rm obs} - \hat{A}\,M_{\rm sim}$"))
    for ax, arr, title in panels:
        im = ax.imshow(arr, origin="lower",
                       extent=(axis[0], axis[-1], axis[0], axis[-1]),
                       cmap="RdBu_r", norm=norm, interpolation="nearest", aspect="equal")
        ax.scatter([-sep / 2, sep / 2], [0, 0], s=34, marker="x", c=INK,
                   linewidths=1.3, zorder=5)
        ax.set_xlim(-zoom, zoom)
        ax.set_ylim(-zoom, zoom)
        ax.set_title(title, loc="left", color=INK, pad=5)
        ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(r"$\kappa$")
        cb.ax.tick_params(labelsize=7)

    # Secondary diagnostic: does anything survive in the residual?
    ax = axes[1, 1]
    shape = obs.shape
    for arr, color, label in ((obs, BLUE, "BOSS"),
                              (resid, "#eb6834", "residual")):
        prof = band_profile(arr, axis) * 1e4
        ax.plot(axis, prof, lw=2, color=color, label=label, zorder=3)
    loo_resid = np.array([band_profile(m.reshape(shape) - model, axis) for m in loo_maps])
    _, err = jackknife_error(loo_resid)
    rprof = band_profile(resid, axis) * 1e4
    ax.fill_between(axis, rprof - err * 1e4, rprof + err * 1e4, color="#eb6834",
                    alpha=0.15, lw=0, zorder=1)
    ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6)
    for xx in (-sep / 2, sep / 2):
        ax.axvline(xx, color=INK_MUTED, lw=0.8, ls="--", alpha=0.7)
    ax.set_xlim(-2.5 * sep, 2.5 * sep)
    ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
    ax.set_ylabel(r"central − off-centre band,  $\kappa$  [$10^{-4}$]")
    ax.set_title("(d) Residual profile  (secondary diagnostic)", loc="left",
                 color=INK, pad=5)
    ax.legend(loc="upper right")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.suptitle(
        f"Simulation template fit — separation {sep:g} $h^{{-1}}$ Mpc\n"
        f"Maps (a)-(c) share one symmetric colour scale, zoomed to $\\pm${zoom:g} "
        f"$h^{{-1}}$Mpc for legibility.  The fit itself spans the whole 100 "
        f"$h^{{-1}}$Mpc map, as requested; $\\hat{{A}}$ = {a_hat:.3f}.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description="Single-amplitude simulation template fit.")
    p.add_argument("--dataset", default="BOSS")
    p.add_argument("--regions", default="North,South")
    p.add_argument("--separations", default="5,10,20")
    p.add_argument("--single-tag", default="_scw")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin.csv")
    p.add_argument("--output-dir", default="output/plots")
    args = p.parse_args(argv)
    regions = [r.strip() for r in args.regions.split(",")]
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    rows = []
    for key in [s.strip() for s in args.separations.split(",")]:
        sep = SEPARATIONS[key]["center"]
        logger.info("=== separation %g h^-1 Mpc ===", sep)

        t = build_terms(args.results_dir, args.dataset, regions, key,
                        args.single_tag, args.single_random_tag)
        obs, axis, shape = t["filament"], t["axis"], t["filament"].shape
        loo_maps = t["loo"]["filament"]
        sim, _ = sim_filament(args.sim_dir, key, args.sim_single)

        _, sigma = jackknife_error(loo_maps)
        sigma = sigma.reshape(shape)
        valid = validate(obs, sim, sigma, sep)

        inv_var = np.zeros_like(sigma)
        np.divide(1.0, sigma ** 2, out=inv_var, where=valid)

        quad = unique_quadrant_mask(shape)
        masks = {
            "full_grid": valid,
            "unique_quadrant": valid & quad,
            "bridge_quadrant": valid & quad & bridge_mask(axis, sep),
        }

        a_for_plot = None
        for name, mask in masks.items():
            fit = fit_amplitude(obs, sim, inv_var, mask)
            a_loo = refit_loo(loo_maps, shape, sim, inv_var, mask, fit.pop("_denom"))
            _, err = jackknife_error(a_loo[:, None])
            sigma_jk = float(err[0])
            if name == "full_grid":
                a_for_plot = fit["A"]

            # Does A survive a free additive constant?
            a_off, c_off, _ = fit_amplitude_with_offset(obs, sim, inv_var, mask)
            a_off_loo = np.array([
                fit_amplitude_with_offset(m.reshape(shape), sim, inv_var, mask)[0]
                for m in loo_maps])
            _, off_err = jackknife_error(a_off_loo[:, None])

            rows.append({"separation_hmpc": sep, "region": name, **fit,
                         "sigma_A_jk": sigma_jk,
                         "A_over_sigma_jk": fit["A"] / sigma_jk if sigma_jk > 0 else np.nan,
                         "jk_over_diag_error": sigma_jk / fit["sigma_A_diag"],
                         "A_with_offset": a_off,
                         "offset_c": c_off,
                         "sigma_A_offset_jk": float(off_err[0])})
            logger.info(
                "  %-16s n=%5d  A = %+.3f +/- %.3f (jk) / %.3f (diag)   "
                "A/sig_jk = %5.2f   sqrt(dchi2) = %6.2f",
                name, fit["n_pixels"], fit["A"], sigma_jk, fit["sigma_A_diag"],
                fit["A"] / sigma_jk if sigma_jk > 0 else np.nan, fit["sqrt_delta_chi2"])

        # Third estimator: the 1-D band-differenced profile (see fit_profile).
        pf = fit_profile(obs, sim, loo_maps, axis)
        rows.append({"separation_hmpc": sep, "region": "profile_1d", **pf,
                     "A_over_sigma_jk": pf["A"] / pf["sigma_A_jk"] if pf["sigma_A_jk"] > 0
                     else np.nan,
                     "jk_over_diag_error": pf["sigma_A_jk"] / pf["sigma_A_diag"]})
        logger.info(
            "  %-16s n=%5d  A = %+.3f +/- %.3f (jk) / %.3f (diag)   "
            "A/sig_jk = %5.2f   sqrt(dchi2) = %6.2f",
            "profile_1d", pf["n_pixels"], pf["A"], pf["sigma_A_jk"], pf["sigma_A_diag"],
            pf["A"] / pf["sigma_A_jk"] if pf["sigma_A_jk"] > 0 else np.nan,
            pf["sqrt_delta_chi2"])

        # Residual chi^2/dof on the unique quadrant. Diagnostic only: one
        # quadrant removes the symmetry duplication but not the beam-induced
        # correlation between neighbouring pixels, so this is not a calibrated
        # goodness-of-fit test.
        qm = masks["unique_quadrant"]
        resid = obs - a_for_plot * sim
        dof = int(qm.sum()) - 1
        chi2_resid = float(np.sum(resid[qm] ** 2 * inv_var[qm]))
        logger.info("  residual chi2/dof (quadrant, diagnostic only) = %.3f  (dof=%d)",
                    chi2_resid / dof, dof)
        rows[-1]["residual_chi2_per_dof_quadrant"] = chi2_resid / dof

        make_sheet(obs, sim, a_for_plot, axis, sep,
                   os.path.join(args.output_dir, f"sim_template_fit_sep{key}"),
                   loo_maps, inv_var, qm)

    out = pd.DataFrame(rows)
    path = os.path.join(args.results_dir, f"sim_template_fit_{args.dataset}.csv")
    out.to_csv(path, index=False)
    logger.info("\nSaved -> %s", path)

    logger.info("\nAmplitude fit summary:")
    logger.info("  sep  region             A      sig_jk   A/sig_jk   sqrt(dchi2)  jk/diag")
    for _, r in out.iterrows():
        logger.info("  %3.0f  %-16s %+6.3f  %7.3f  %8.2f  %11.2f  %7.2f",
                    r.separation_hmpc, r.region, r.A, r.sigma_A_jk,
                    r.A_over_sigma_jk, r.sqrt_delta_chi2, r.jk_over_diag_error)


if __name__ == "__main__":
    main()
