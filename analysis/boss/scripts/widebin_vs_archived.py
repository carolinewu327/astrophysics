#!/usr/bin/env python
"""The paper's bridge study on the archived narrow bins and the wide bins, side by side.

Recomputes every number the draft (paper/styau/3_methods_results_draft.tex,
section "Bridge convergence") quotes, with the same functions, for

    archived   r_perp 4-6 (r_par<=5), 9-11 and 18-22 (r_par<=10); pairs
               unweighted by 1/Sigma_crit^2; nominal single-separation control
    wide       r_perp 4-6, 6-15, 15-25 (all r_par<=5), plus the 6-10 / 10-15
               split; every stack weighted by 1/Sigma_crit^2; control averaged
               over the pairs' actual separations; full sub-bin random run

Per sample:
  * bridge excess +/- jackknife error (mean of B(X) over |X| <= 0.35 R), and
    the mock prediction of the same statistic with the observed-minus-mock
    difference in units of the quadrature error;
  * chi^2 of the 22-element compressed profile against zero (full covariance);
  * amplitude of the mock residual profile fitted to the data over the full
    profile (with the share of the response from bins around the galaxies,
    and with a free constant) and over the bridge window only;
  * the random-pair map's own bridge statistic (noise floor);
  * an inverse-variance pool of the three main bins, with its error both
    neglecting and including the covariance between bins (the bins share the
    287 jackknife regions, so the joint covariance comes from the same
    leave-one-out values), and the minimum-variance combination.

The archived column must reproduce the draft; it is recomputed here rather than
copied so that both columns come from one code path.  Nothing in the archive is
written: outputs go to widebin_rpar5/comparison/.

Usage
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python \\
        analysis/boss/scripts/widebin_vs_archived.py
"""

from __future__ import annotations

import argparse
import json
import logging
import os

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist

from catalog import setup_logging
from combine_filament_jackknife import (BRIDGE_CENTERS, grouping_filament, load_grouping_pairs,
                                        load_random_subbins, loo_and_total, pooled_random_map)
from combine_jackknife import load_accumulator
from fit_band_covariance import Precision, run_fit
from geometry import (BRIDGE_HALF_X_FRAC, band_profile, bridge_excess, separation_averaged_template,
                      separation_mix)
from jackknife import jackknife_error
from plot_obs_vs_sim_filament import sim_filament
from plot_separation_summary import build_terms
from weighting_sensitivity import SEPARATIONS, axis_for_map, load_map

logger = logging.getLogger(__name__)

ARCHIVED = [  # (label, SEPARATIONS key, bin text, r_par)
    ("R=5", "5", "4-6", 5.0), ("R=10", "10", "9-11", 10.0), ("R=20", "20", "18-22", 10.0)]
WIDE = [  # (label, grouping, bin text)
    ("R=5", "g4_6", "4-6"), ("R=10", "g6_15", "6-15"), ("R=20", "g15_25", "15-25"),
    ("split", "g6_10", "6-10"), ("split", "g10_15", "10-15")]
X_MAX, X_BIN = 40.0, 4.0


def bridge_with_error(obs: np.ndarray, loo: np.ndarray, axis: np.ndarray, center: float):
    value = bridge_excess(obs, axis, center)
    loo_vals = np.array([bridge_excess(m.reshape(obs.shape), axis, center) for m in loo])
    return value, jackknife_error(loo_vals)[1], loo_vals


def offset_error(d, s, cov):
    """1-sigma error on the free constant of the amplitude-plus-constant fit."""
    B, one = Precision(cov), np.ones_like(d)
    m = np.array([[B.dot(s, s), B.dot(s, one)], [B.dot(one, s), B.dot(one, one)]])
    return float(np.sqrt(np.linalg.inv(m)[1, 1]))


def profile_with_error(obs, loo, axis):
    prof = band_profile(obs, axis)
    loo_p = np.array([band_profile(m.reshape(obs.shape), axis) for m in loo])
    return prof, jackknife_error(loo_p)[1]


def score(label, bin_text, rpar, center, obs, loo, sim, sim_bridge, n_pairs, mean_r, rand_floor, axis):
    """Every draft number for one sample."""
    b, be, b_loo = bridge_with_error(obs, loo, axis, center)
    sb, sbe = sim_bridge
    full, (_, d_vec, s_vec, cov, _) = run_fit(obs, loo, sim, axis, X_MAX, X_BIN, center)
    c_err = offset_error(d_vec, s_vec, cov)
    win, _ = run_fit(obs, loo, sim, axis, BRIDGE_HALF_X_FRAC * center, X_BIN, center)
    return {
        "sample": label, "bin": bin_text, "rpar_max": rpar, "center": center,
        "n_pairs": n_pairs, "mean_rperp": mean_r,
        "bridge": b, "bridge_err": be, "bridge_snr": b / be,
        "mock_bridge": sb, "mock_err": sbe, "obs_minus_mock_sigma": (b - sb) / np.hypot(be, sbe),
        "chi2_profile": full["chi2_A0"], "chi2_dof": full["dim"],
        "chi2_p": float(chi2_dist.sf(full["chi2_A0"], full["dim"])),
        # Same vector and covariance after the best-fit mock profile is removed:
        # chi2_fit = chi2_profile - A_full_sig^2, with one fewer degree of freedom.
        "chi2_fit": full["chi2_min"], "chi2_fit_dof": full["dof"],
        "chi2_fit_p": float(chi2_dist.sf(full["chi2_min"], full["dof"])),
        "A_full": full["A"], "A_full_err": full["sigma_A"], "A_full_sig": full["significance"],
        "A_full_halo_frac": full["frac_halo"], "A_full_bridge_frac": full["frac_bridge"],
        "A_offset": full["A_offset"], "A_offset_err": full["sigma_A_offset"], "offset_c": full["offset_c"],
        "offset_c_err": c_err, "offset_c_sig": full["offset_c"] / c_err,
        "A_window": win["A"], "A_window_err": win["sigma_A"], "A_window_sig": win["significance"],
        "random_floor": rand_floor,
        "_loo_bridge": b_loo,   # leave-one-out bridge values, for the joint covariance; not written
    }


def archived_samples(results_dir, sim_dir, sim_single, single_tag, single_random_tag):
    rows, profiles = [], {}
    for label, key, bin_text, rpar in ARCHIVED:
        t = build_terms(results_dir, "BOSS", ["North", "South"], key, single_tag, single_random_tag)
        center, axis = SEPARATIONS[key]["center"], t["axis"]
        obs, loo = t["filament"], t["loo"]["filament"]
        sim, meta = sim_filament(sim_dir, key, sim_single)
        mb = meta["stats"]["residual_bridge_excess_fixedband_kappa"]
        acc = load_accumulator(os.path.join(
            results_dir, "jk", f"acc_pairs_galaxy_{SEPARATIONS[key]['galaxy_label']}_BOSS_North_South.npz"))
        floor = pd.read_csv(os.path.join(results_dir, "filament_jackknife_BOSS_North_South.csv"))
        floor = float(floor.loc[floor.separation_hmpc == center, "random_pair_noise_floor"].iloc[0])
        rows.append({"version": "archived", **score(label, bin_text, rpar, center, obs, loo, sim,
                                                   (mb["value"], mb["jackknife_error"]),
                                                   int(acc["n_pairs"].sum()), np.nan, floor, axis)})
        profiles[("archived", label)] = (axis, *profile_with_error(obs, loo, axis), band_profile(sim, axis), center)
        logger.info("archived %s done", label)
    return rows, profiles


def wide_samples(results_dir, wide_dir, sim_dir, sim_single, single_tag, single_random_tag, random_label):
    gal = load_accumulator(os.path.join(results_dir, "jk", f"acc_single_galaxy{single_tag}_BOSS_North_South.npz"))
    rnd = load_accumulator(os.path.join(results_dir, "jk",
                                        f"acc_single_random{single_random_tag}_BOSS_North_South.npz"))
    gt, gl = loo_and_total(gal)
    rt, rl = loo_and_total(rnd)
    st, sl, g = gt - rt, gl - rl, gal["grid_size"]
    sim_single_map = load_map(os.path.join(sim_dir, sim_single))
    ax = np.linspace(-50, 50, 101)
    X, Y = np.meshgrid(ax, ax)

    rows, profiles = [], {}
    for label, name, bin_text in WIDE:
        lo, hi = (float(v) for v in name[1:].split("_"))
        center = BRIDGE_CENTERS[name]
        pair = load_grouping_pairs(os.path.join(
            wide_dir, "jk", f"acc_pairs_galaxy_{name}_rpar5_scw_BOSS_North_South.npz"), lo, hi)
        rand = pooled_random_map(load_random_subbins(wide_dir, random_label, "BOSS",
                                                     ["North", "South"], lo, hi))
        pt, pl = loo_and_total(pair)
        fine = (pair["fine_w"], pair["fine_wr"], pair["fine_cov"])

        # Pair maps reflection-symmetrized, as build_terms does for the archived
        # column (grouping_filament then reflects the averaged control too).
        def fil(p, s, mix):
            return grouping_filament(p, s, rand, mix, center, "separation_averaged", g, symmetrize=True)[0]

        obs = fil(pt, st, separation_mix(*fine))
        loo = np.array([fil(pl[k], sl[k], separation_mix(*fine, exclude_region=k)).ravel()
                        for k in range(sl.shape[0])])
        axis = axis_for_map(obs)

        base = os.path.join(sim_dir, "hod_pairs", "widebin_rpar5", f"stack_{name}_rpar5_avgctrl")
        blocks = np.load(base + "_blocks.npz")
        seps, wts, _ = separation_mix(blocks["fine_w"], blocks["fine_wr"])
        sim = load_map(base + "_map.csv") - separation_averaged_template(sim_single_map, X, Y, seps, wts)
        with open(base + ".json", encoding="utf-8") as fh:
            mb = json.load(fh)["stats"]["residual_bridge_excess_fixedband_kappa"]

        mean_r = float(pair["fine_wr"].sum() / pair["fine_w"].sum())
        rows.append({"version": "wide", **score(label, bin_text, 5.0, center, obs, loo, sim,
                                               (mb["value"], mb["jackknife_error"]),
                                               int(pair["n_pairs"].sum()), mean_r,
                                               bridge_excess(rand, axis_for_map(rand), center), axis)})
        profiles[("wide", bin_text)] = (axis, *profile_with_error(obs, loo, axis), band_profile(sim, axis), center)
        logger.info("wide %s done", bin_text)
    return rows, profiles


def pooled(df: pd.DataFrame, version: str) -> dict:
    """Combine the three main bins.

    pooled_bridge uses inverse-variance weights.  pooled_err treats the bins as
    independent; pooled_err_joint is the error of that same weighted mean with
    the jackknife covariance between bins.  The gls_* columns are the
    minimum-variance combination under that covariance.
    """
    main = df[(df.version == version) & (df["sample"] != "split")]
    w = (1.0 / main.bridge_err ** 2).to_numpy()
    b = main.bridge.to_numpy()
    value = float((w * b).sum() / w.sum())
    err = float(1.0 / np.sqrt(w.sum()))

    loo = np.column_stack(list(main["_loo_bridge"]))      # (regions, bins)
    k = loo.shape[0]
    dev = loo - loo.mean(axis=0)
    cov = (k - 1) / k * dev.T @ dev
    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    err_joint = float(np.sqrt(w @ cov @ w) / w.sum())
    cinv = np.linalg.inv(cov)
    gls_err = float(1.0 / np.sqrt(cinv.sum()))
    gls = float((cinv @ b).sum() * gls_err ** 2)
    return {"version": version, "pooled_bridge": value, "pooled_err": err, "pooled_snr": value / err,
            "pooled_err_joint": err_joint, "pooled_snr_joint": value / err_joint,
            "corr_5_10": corr[0, 1], "corr_5_20": corr[0, 2], "corr_10_20": corr[1, 2],
            "gls_bridge": gls, "gls_err": gls_err, "gls_snr": gls / gls_err}


def plot_profiles(profiles, out_base):
    pairs = [(("archived", "R=5"), ("wide", "4-6")), (("archived", "R=10"), ("wide", "6-15")),
             (("archived", "R=20"), ("wide", "15-25"))]
    titles = {"archived": "archived bin", "wide": "wide bin"}
    fig, axes = plt.subplots(3, 2, figsize=(9.0, 9.0), sharey="row", constrained_layout=True)
    for r, row in enumerate(pairs):
        for c, key in enumerate(row):
            ax_ = axes[r, c]
            axis, prof, err, sim_prof, center = profiles[key]
            m = (axis >= 0) & (axis <= 25)
            ax_.axvspan(0, BRIDGE_HALF_X_FRAC * center, color="#e3e2df", zorder=0)
            ax_.errorbar(axis[m], prof[m] * 1e4, yerr=err[m] * 1e4, fmt="o", ms=3.5, color="#2a78d6",
                         elinewidth=1, capsize=1.5, label="BOSS")
            ax_.plot(axis[m], sim_prof[m] * 1e4, color="#eb6834", lw=2, label="mock")
            ax_.axhline(0, color="#8a8985", lw=0.8)
            if key[0] == "wide":
                lo, hi = (float(v) for v in key[1].split("-"))
                ax_.axvspan(lo / 2, hi / 2, color="#2a78d6", alpha=0.06, zorder=0)
            else:
                ax_.axvline(center / 2, color="#3a9b5c", ls=":", lw=1)
            bin_text = key[1] if key[0] == "wide" else {"R=5": "4-6", "R=10": "9-11", "R=20": "18-22"}[key[1]]
            ax_.set_title(f"{titles[key[0]]} r$_\\perp$ = {bin_text}", loc="left", fontsize=9.5)
            ax_.set_xlabel(r"$X$ along pair axis [$h^{-1}$Mpc]")
            if c == 0:
                ax_.set_ylabel(r"$B(X)$ [$10^{-4}$]")
            for side in ("top", "right"):
                ax_.spines[side].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--results-dir", default="analysis/boss/results")
    p.add_argument("--wide-dir", default="analysis/boss/results/widebin_rpar5")
    p.add_argument("--sim-dir", default="analysis/sim/results")
    p.add_argument("--sim-single", default="kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv")
    p.add_argument("--random-label", default="random_w3_25_rpar5_scw_frac100")
    p.add_argument("--single-tag", default="_scw")
    p.add_argument("--single-random-tag", default="_scw_frac100")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    setup_logging()

    out_dir = os.path.join(args.wide_dir, "comparison")
    outputs = [os.path.join(out_dir, f) for f in
               ("bridge_study_archived_vs_wide.csv", "bridge_study_pooled.csv",
                "band_profiles_archived_vs_wide.png", "band_profiles_archived_vs_wide.pdf")]
    existing = [o for o in outputs if os.path.exists(o)]
    if existing and not args.overwrite:
        raise FileExistsError("Outputs exist (use --overwrite): " + ", ".join(existing))
    os.makedirs(out_dir, exist_ok=True)

    rows_a, prof_a = archived_samples(args.results_dir, args.sim_dir, args.sim_single,
                                      args.single_tag, args.single_random_tag)
    rows_w, prof_w = wide_samples(args.results_dir, args.wide_dir, args.sim_dir, args.sim_single,
                                  args.single_tag, args.single_random_tag, args.random_label)
    df = pd.DataFrame(rows_a + rows_w)
    df.drop(columns="_loo_bridge").to_csv(outputs[0], index=False)
    pd.DataFrame([pooled(df, "archived"), pooled(df, "wide")]).to_csv(outputs[1], index=False)
    plot_profiles({**prof_a, **prof_w}, os.path.join(out_dir, "band_profiles_archived_vs_wide"))
    logger.info("Saved -> %s", out_dir)


if __name__ == "__main__":
    main()
