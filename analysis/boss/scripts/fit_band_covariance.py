#!/usr/bin/env python
"""Template fit to the compressed two-band profile, with a full covariance.

The pixel-level template fit could not quote an honest error.  Its chi^2 used
diagonal per-pixel errors, but the 8 arcmin Planck beam is ~3.3 h^-1 Mpc at
z = 0.55 against 1 h^-1 Mpc pixels, so neighbouring pixels move together and the
Delta chi^2 = 1 error came out 14-20x smaller than a delete-one jackknife of the
same amplitude.  Estimating the true covariance was not an option either: the
101x101 map carries ~2600 independent values and only 287 jackknife samples
exist to constrain their covariance.

The fix is to compress before fitting.  Reduce each map to the mean central-band
and mean off-centre-band profiles as a function of X -- a few tens of numbers --
and the covariance becomes small enough to estimate from 287 samples and invert.
The chi^2 then accounts for the correlations, so its Delta chi^2 = 1 error is
the right one to quote.  Maps stay in the paper for illustration; the number
comes from the compressed vector.

    kappa_band = [ kappa_central(X_0..X_n) , kappa_off(X_0..X_n) ],   X >= 0

The two bands are kept *separate* rather than differenced.  Differencing cancels
any map-wide additive constant, which is convenient but throws away information;
keeping them apart lets the covariance decide, and the nested-offset test below
recovers the immunity to a monopole explicitly when that is what is wanted.

Estimator, with B = C^-1:

    A_hat  = (s^T B d) / (s^T B s)
    sigma_A = (s^T B s)^-1/2                        [Delta chi^2 = 1]
    chi2(0) = d^T B d
    dchi2   = A_hat^2 (s^T B s) = (A_hat / sigma_A)^2

**On the precision matrix.**  Hartlap's (N-d-2)/(N-1) factor is derived for
covariances built from *independent* Gaussian realizations.  Delete-one
jackknife samples share 285 of their 287 regions with each other, so C is not
Wishart with N = 287 degrees of freedom and the factor has no derivation here.
Rather than argue the point, the default binning sidesteps it: at 4 h^-1 Mpc
bins the data vector is d = 22, where the factor would be 0.920 and would move
sigma_A by 4.3% -- smaller than the effects being measured, whether or not it
applies.  It is reported as a clearly-labelled sensitivity column, never folded
into the headline.  The default B is unscaled.

**On the curvature/jackknife cross-check.**  With B held fixed,

    A_k - A_bar = s^T B (d_k - d_bar) / (s^T B s)

so the delete-one variance of A_k is s^T B C B s / (s^T B s)^2, which at
B = C^-1 collapses to 1 / (s^T C^-1 s) = sigma_A^2 *identically* -- C is defined
from those same d_k.  So the agreement checked below is an algebraic unit test
that catches coding errors (a dropped (N-1)/N, the wrong centring, mismatched
masks); it is not independent statistical validation, and it only holds for the
unscaled B.

**Known limitations, disclosed rather than fixed.**
  * The random-pair map is held fixed rather than jackknifed, so C omits its
    noise.  Measured floor ~0.5e-4, below the galaxy-side error.
  * The simulation template carries its own noise, which biases A low by
    regression dilution.  Expected to be subdominant but currently unquantified:
    a high scalar bridge S/N does not establish that every mode of the two-band
    template is noise-free.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts python \\
        analysis/boss/scripts/fit_band_covariance.py --separations 5,10,20
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
from scipy.linalg import cho_factor, cho_solve

from catalog import setup_logging
from geometry import band_response_zones, kappa_band
from plot_obs_vs_sim_filament import sim_filament
from plot_separation_summary import build_terms
from weighting_sensitivity import SEPARATIONS

logger = logging.getLogger(__name__)

INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e3e2df"
BLUE, ORANGE = "#2a78d6", "#eb6834"

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


# ---------------------------------------------------------------------------
# Covariance
# ---------------------------------------------------------------------------
def jackknife_covariance(loo: np.ndarray) -> np.ndarray:
    """Delete-one covariance of a stack of leave-one-out vectors.

    Centred on the *mean of the leave-one-out vectors*, not on the full-sample
    measurement.  The stack is a ratio estimator -- sum(w*kappa)/sum(w) -- so
    the two are not identical, and the delete-one formula is defined about the
    LOO mean.  Zheng's note allows either; this is the stricter reading.
    """
    n = loo.shape[0]
    delta = loo - loo.mean(axis=0)
    return (n - 1) / n * (delta.T @ delta)


def precision_diagnostics(cov: np.ndarray, n_samples: int) -> dict:
    """Is this covariance safe to invert with the samples available?"""
    eig = np.linalg.eigvalsh(cov)
    d = cov.shape[0]
    hartlap = (n_samples - d - 2) / (n_samples - 1)
    return {
        "dim": d,
        "n_samples": n_samples,
        "rank": int(np.linalg.matrix_rank(cov)),
        "eig_min": float(eig.min()),
        "eig_max": float(eig.max()),
        "cond": float(eig.max() / eig.min()) if eig.min() > 0 else np.inf,
        "hartlap_alpha": float(hartlap),
        "hartlap_sigma_inflation": float(1.0 / np.sqrt(hartlap)) if hartlap > 0 else np.nan,
    }


class Precision:
    """B = C^-1 applied by Cholesky solve; C is never explicitly inverted."""

    def __init__(self, cov: np.ndarray):
        eig_min = float(np.linalg.eigvalsh(cov).min())
        if eig_min <= 0:
            raise ValueError(
                f"Covariance is not positive definite (min eigenvalue {eig_min:.3e}). "
                "Coarsen the X binning or shorten the fit range.")
        self._cho = cho_factor(cov, lower=True)

    def apply(self, vec: np.ndarray) -> np.ndarray:
        return cho_solve(self._cho, vec)

    def dot(self, a: np.ndarray, b: np.ndarray) -> float:
        """a^T B b"""
        return float(a @ self.apply(b))


# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------
def fit_amplitude(d: np.ndarray, s: np.ndarray, B: Precision) -> dict:
    """Single free amplitude, closed form."""
    sbs = B.dot(s, s)
    sbd = B.dot(s, d)
    dbd = B.dot(d, d)
    a = sbd / sbs
    dchi2 = a * a * sbs
    return {
        "A": float(a),
        "sigma_A": float(1.0 / np.sqrt(sbs)),
        "chi2_A0": float(dbd),
        "chi2_min": float(dbd - dchi2),
        "delta_chi2": float(dchi2),
        "significance": float(a * np.sqrt(sbs)),
        "_sbs": float(sbs),
    }


def fit_amplitude_with_offset(d: np.ndarray, s: np.ndarray, B: Precision) -> dict:
    """Nested-model test against a free monopole.

    A residual monopole from imperfect control subtraction adds equally to both
    bands, so the nuisance is a single constant across the whole vector.  The
    *null* must also carry it, otherwise Delta chi^2 would credit A with the
    improvement that comes from fitting the offset:

        null  d = c*1              (c free)
        alt   d = A*s + c*1        (A, c free)
    """
    one = np.ones_like(d)
    dbd = B.dot(d, d)

    obo, obd = B.dot(one, one), B.dot(one, d)
    chi2_null = dbd - obd * obd / obo

    m = np.array([[B.dot(s, s), B.dot(s, one)],
                  [B.dot(one, s), obo]])
    rhs = np.array([B.dot(s, d), obd])
    theta = np.linalg.solve(m, rhs)
    chi2_alt = dbd - float(rhs @ theta)
    cov_theta = np.linalg.inv(m)
    return {
        "A_offset": float(theta[0]),
        "sigma_A_offset": float(np.sqrt(cov_theta[0, 0])),
        "offset_c": float(theta[1]),
        "chi2_null_offset_only": float(chi2_null),
        "chi2_min_offset": float(chi2_alt),
        "delta_chi2_nested": float(chi2_null - chi2_alt),
        "significance_offset": float(theta[0] / np.sqrt(cov_theta[0, 0])),
    }


def response_decomposition(x: np.ndarray, s: np.ndarray, B: Precision, sbs: float,
                           sep: float, x_bin: float) -> tuple[dict, np.ndarray]:
    """Where does the fitted amplitude actually come from?

    Writing u = B s / (s^T B s), the estimator is the linear form A = u . d, and
    u . s = 1 identically.  So the per-element products u_i s_i partition the
    fit's response into fractions summing to one, and can be grouped to ask the
    question that deflated the two previous estimators: is this measuring the
    bridge, or the halo peaks the control subtraction was supposed to remove?

    Zones come from geometry.band_response_zones, shared with the diagnostic
    plot so the two cannot disagree on where the halo is.
    """
    u = B.apply(s) / sbs
    contrib = u * s                      # sums to exactly 1
    n = len(x)
    central, off = contrib[:n], contrib[n:]

    is_bridge, is_halo, is_outer = band_response_zones(x, sep, x_bin)

    per_bin = central + off
    out = {
        "frac_central_band": float(central.sum()),
        "frac_off_band": float(off.sum()),
        "frac_bridge": float(per_bin[is_bridge].sum()),
        "frac_halo": float(per_bin[is_halo].sum()),
        "frac_outer": float(per_bin[is_outer].sum()),
        "n_bins_bridge": int(is_bridge.sum()),
        "n_bins_halo": int(is_halo.sum()),
        "n_bins_outer": int(is_outer.sum()),
        "bridge_halo_resolved": bool(is_bridge.sum() and is_halo.sum()),
    }
    return out, per_bin


def loo_refit_scatter(loo_vectors: np.ndarray, s: np.ndarray, B: Precision,
                      sbs: float) -> float:
    """Delete-one scatter of A with the covariance (hence the weights) held fixed.

    Algebraically equal to sigma_A for unscaled B -- see the module docstring.
    Computed anyway as a unit test on the covariance construction.
    """
    bs = B.apply(s)
    a_k = (loo_vectors @ bs) / sbs
    n = len(a_k)
    return float(np.sqrt((n - 1) / n * np.sum((a_k - a_k.mean()) ** 2)))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def compress_all(obs: np.ndarray, loo_maps: np.ndarray, sim: np.ndarray,
                 axis: np.ndarray, x_max: float, x_bin: float):
    """Compress the observation, its 287 LOO realizations, and the template.

    Every profile is binned *before* the covariance is built.  Binning is linear
    so this is equivalent to binning C afterwards, but it keeps the matrix at
    its final size throughout.
    """
    shape = obs.shape
    x, d = kappa_band(obs, axis, x_max, x_bin)
    _, s = kappa_band(sim, axis, x_max, x_bin)
    loo = np.array([kappa_band(m.reshape(shape), axis, x_max, x_bin)[1]
                    for m in loo_maps])
    return x, d, s, loo


def run_fit(obs, loo_maps, sim, axis, x_max, x_bin, sep):
    x, d, s, loo = compress_all(obs, loo_maps, sim, axis, x_max, x_bin)
    cov = jackknife_covariance(loo)
    diag = precision_diagnostics(cov, loo.shape[0])
    B = Precision(cov)

    fit = fit_amplitude(d, s, B)
    sbs = fit.pop("_sbs")
    fit.update(fit_amplitude_with_offset(d, s, B))
    zones, per_bin = response_decomposition(x, s, B, sbs, sep, x_bin)
    fit.update(zones)

    sigma_loo = loo_refit_scatter(loo, s, B, sbs)
    fit["sigma_A_loo_refit"] = sigma_loo
    fit["unit_test_ratio"] = sigma_loo / fit["sigma_A"]
    fit["sigma_A_hartlap"] = fit["sigma_A"] * diag["hartlap_sigma_inflation"]
    fit["dof"] = diag["dim"] - 1
    fit["chi2_min_per_dof"] = fit["chi2_min"] / fit["dof"]
    fit.update(diag)
    fit["x_max"] = x_max
    fit["x_bin"] = x_bin
    return fit, (x, d, s, cov, per_bin)


def make_plot(x, d, s, cov, a_hat, sep, x_max, x_bin, out_base):
    """Central-band and off-band profiles vs X, with errors and the best fit."""
    n = len(x)
    err = np.sqrt(np.diag(cov))
    # Shared y: both panels are kappa in the same units, and the contrast
    # between the bands is the signal. Letting each autoscale would make a
    # small off-band excursion look like the large central one.
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), constrained_layout=True,
                             sharex=True, sharey=True)

    panels = ((axes[0], slice(0, n), BLUE, "o", "central band  ($|Y| \\leq 1.5$)"),
              (axes[1], slice(n, 2 * n), ORANGE, "s",
               "off-centre band  ($1.5 \\leq |Y| \\leq 10.5$)"))
    for ax, sl, color, marker, label in panels:
        ax.errorbar(x, d[sl] * 1e4, yerr=err[sl] * 1e4, fmt=marker, ms=5,
                    color=color, lw=0, elinewidth=1.2, capsize=2.5,
                    label="BOSS", zorder=3)
        ax.plot(x, a_hat * s[sl] * 1e4, lw=2, color=INK, alpha=0.75,
                label=rf"$\hat{{A}}\,\times$ simulation  ($\hat{{A}}$ = {a_hat:.2f})",
                zorder=2)
        ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6)
        ax.axvline(sep / 2, color=INK_MUTED, lw=0.8, ls="--", alpha=0.7)
        ax.annotate("halo", (sep / 2, ax.get_ylim()[1]), xytext=(3, -10),
                    textcoords="offset points", color=INK_MUTED, fontsize=7.5)
        ax.set_title(label, loc="left", color=INK, pad=5)
        ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
        ax.set_ylabel(r"$\kappa$  [$10^{-4}$]")
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.legend(loc="upper right")

    fig.suptitle(
        f"Compressed band fit — separation {sep:g} $h^{{-1}}$ Mpc\n"
        f"Data vector is the two profiles jointly ({2 * n} numbers, {x_bin:g} "
        f"$h^{{-1}}$Mpc bins to $X$ = {x_max:g}); error bars are the square root of the "
        f"jackknife covariance diagonal, but the fit uses the full matrix.",
        fontsize=10.5, color=INK, ha="left", x=0.01)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
    plt.close(fig)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Compressed two-band template fit.")
    p.add_argument("--dataset", default="BOSS")
    p.add_argument("--regions", default="North,South")
    p.add_argument("--separations", default="5,10,20")
    p.add_argument("--single-tag", default="_scw")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    p.add_argument("--x-max", type=float, default=40.0)
    p.add_argument("--x-bin", type=float, default=4.0,
                   help="Default 4 h^-1 Mpc ~ the 8 arcmin beam scale.")
    p.add_argument("--scan-x-bins", default="1,2,4")
    p.add_argument("--scan-x-max", default="10,20,40")
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin.csv")
    p.add_argument("--output-dir", default="output/plots")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    regions = [r.strip() for r in args.regions.split(",")]
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    primary, scan = [], []
    for key in [s.strip() for s in args.separations.split(",")]:
        sep = SEPARATIONS[key]["center"]
        logger.info("=== separation %g h^-1 Mpc ===", sep)

        t = build_terms(args.results_dir, args.dataset, regions, key,
                        args.single_tag, args.single_random_tag)
        obs, axis = t["filament"], t["axis"]
        loo_maps = t["loo"]["filament"]
        sim, _ = sim_filament(args.sim_dir, key, args.sim_single)
        if sim.shape != obs.shape:
            raise ValueError(f"Grid mismatch at sep={sep}: {obs.shape} vs {sim.shape}.")

        fit, (x, d, s, cov, per_bin) = run_fit(obs, loo_maps, sim, axis,
                                               args.x_max, args.x_bin, sep)
        primary.append({"separation_hmpc": sep, **fit})
        pd.DataFrame({"x_hmpc": x, "response_fraction": per_bin}).to_csv(
            os.path.join(args.results_dir,
                         f"band_fit_response_sep{key}_{args.dataset}.csv"), index=False)

        logger.info("  primary: %g h^-1 Mpc bins to X=%g -> d=%d",
                    args.x_bin, args.x_max, fit["dim"])
        logger.info("    A = %+.3f +/- %.3f   (%.2f sigma)   dchi2 = %.2f",
                    fit["A"], fit["sigma_A"], fit["significance"], fit["delta_chi2"])
        logger.info("    chi2(A=0) = %.2f   chi2_min/dof = %.3f (dof=%d)",
                    fit["chi2_A0"], fit["chi2_min_per_dof"], fit["dof"])
        logger.info("    with free offset: A = %+.3f +/- %.3f (%.2f sigma), "
                    "nested dchi2 = %.2f, c = %+.3fe-4",
                    fit["A_offset"], fit["sigma_A_offset"], fit["significance_offset"],
                    fit["delta_chi2_nested"], fit["offset_c"] * 1e4)
        logger.info("    covariance: rank %d/%d, cond = %.3e, eig_min = %.3e",
                    fit["rank"], fit["dim"], fit["cond"], fit["eig_min"])
        logger.info("    unit test  sigma_loo/sigma_curv = %.6f  (must be 1)",
                    fit["unit_test_ratio"])
        logger.info("    Hartlap sensitivity (approximate, not applied): "
                    "alpha = %.3f -> sigma_A = %.3f",
                    fit["hartlap_alpha"], fit["sigma_A_hartlap"])
        logger.info("    response: central band %+.1f%% / off band %+.1f%%",
                    fit["frac_central_band"] * 100, fit["frac_off_band"] * 100)
        logger.info("    response: bridge %+.1f%% (%d bins) / halo %+.1f%% (%d bins) "
                    "/ outer %+.1f%% (%d bins)",
                    fit["frac_bridge"] * 100, fit["n_bins_bridge"],
                    fit["frac_halo"] * 100, fit["n_bins_halo"],
                    fit["frac_outer"] * 100, fit["n_bins_outer"])
        if not fit["bridge_halo_resolved"]:
            logger.warning("    bridge and halo are NOT separated at this binning "
                           "(sep=%g, x_bin=%g): do not read A as a bridge amplitude.",
                           sep, args.x_bin)
        if abs(fit["unit_test_ratio"] - 1.0) > 1e-6:
            logger.error("    UNIT TEST FAILED at sep=%g: ratio %.6f", sep,
                         fit["unit_test_ratio"])

        make_plot(x, d, s, cov, fit["A"], sep, args.x_max, args.x_bin,
                  os.path.join(args.output_dir, f"band_fit_sep{key}"))

        for xb in [float(v) for v in args.scan_x_bins.split(",")]:
            for xm in [float(v) for v in args.scan_x_max.split(",")]:
                try:
                    f, _ = run_fit(obs, loo_maps, sim, axis, xm, xb, sep)
                except ValueError as exc:
                    logger.warning("  scan x_bin=%g x_max=%g skipped: %s", xb, xm, exc)
                    continue
                scan.append({"separation_hmpc": sep, **f})

    R = args.results_dir
    pd.DataFrame(primary).to_csv(
        os.path.join(R, f"band_fit_{args.dataset}.csv"), index=False)
    pd.DataFrame(scan).to_csv(
        os.path.join(R, f"band_fit_scan_{args.dataset}.csv"), index=False)

    logger.info("\nPrimary fit (%g h^-1 Mpc bins, X = 0-%g):", args.x_bin, args.x_max)
    logger.info("  sep    d      A     sigma_A   sigma   chi2(0)   dchi2   chi2/dof   cond")
    for r in primary:
        logger.info("  %3.0f  %3d  %+6.3f  %8.3f  %6.2f  %8.2f  %6.2f  %9.3f  %8.1e",
                    r["separation_hmpc"], r["dim"], r["A"], r["sigma_A"],
                    r["significance"], r["chi2_A0"], r["delta_chi2"],
                    r["chi2_min_per_dof"], r["cond"])

    logger.info("\nStability scan (A +/- sigma_A, significance, response split):")
    logger.info("  sep  x_bin  x_max   d      A     sigma_A   sigma   bridge%%   halo%%"
                "   outer%%  cond")
    for r in scan:
        logger.info("  %3.0f  %5.0f  %5.0f  %3d  %+6.3f  %8.3f  %6.2f  %7.1f %7.1f "
                    "%7.1f  %7.1e",
                    r["separation_hmpc"], r["x_bin"], r["x_max"], r["dim"],
                    r["A"], r["sigma_A"], r["significance"],
                    r["frac_bridge"] * 100, r["frac_halo"] * 100,
                    r["frac_outer"] * 100, r["cond"])
    logger.info("\nSaved -> %s", os.path.join(R, f"band_fit_{args.dataset}.csv"))


if __name__ == "__main__":
    main()
