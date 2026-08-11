#!/usr/bin/env python
"""Validation figure: map-sampled vs direct-annular single-halo radial profile.

Compares the single-halo kappa radial profile obtained by sampling the
simulated kappa map (``stack_single_sim.py``) against the independent
direct-annular validator, which counts particles in periodic transverse
annuli around halos and never touches the kappa map. Agreement in the
signal-dominated inner region validates the map-construction and
map-sampling path used by all downstream stacks.

Usage
-----
    python plot_validation_single_profile.py
    python plot_validation_single_profile.py --n-halos 5000 \
        --direct-profile analysis/sim/results/radial_profile_single_direct_annuli_mass13_5000halos.csv
"""

from __future__ import annotations

import argparse
import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sim_utils import (
    BIGMDPL_DOWNSAMPLE_FACTOR,
    BIGMDPL_PARTICLE_MASS_HMSUN,
    SIGMA_C_HMSUN_PER_MPC2,
    ensure_parent,
    setup_logging,
)


logger = logging.getLogger(__name__)


def poisson_sigma_kappa(
    particle_count: np.ndarray,
    annulus_area_hmpc2: np.ndarray,
    n_halos: int,
) -> np.ndarray:
    """Poisson (shot-noise) uncertainty of the direct-annular kappa.

    This is a lower bound on the true uncertainty: at large radii the
    scatter between halos is dominated by correlated large-scale
    structure, which Poisson counting does not capture.
    """

    mp_eff = BIGMDPL_DOWNSAMPLE_FACTOR * BIGMDPL_PARTICLE_MASS_HMSUN
    return (
        mp_eff
        * np.sqrt(np.maximum(particle_count, 1.0))
        / (n_halos * annulus_area_hmpc2)
        / SIGMA_C_HMSUN_PER_MPC2
    )


def make_figure(
    map_profile: pd.DataFrame,
    direct_profile: pd.DataFrame,
    n_halos: int,
    output: str,
) -> None:
    merged = pd.merge(
        map_profile,
        direct_profile,
        on="radius_hmpc",
        suffixes=("_map", "_direct"),
    )
    sigma = poisson_sigma_kappa(
        merged["particle_count"].to_numpy(dtype=np.float64),
        merged["annulus_area_hmpc2"].to_numpy(dtype=np.float64),
        n_halos,
    )
    pull = (merged["kappa_direct"] - merged["kappa_map"]) / sigma

    fig, (ax_prof, ax_pull) = plt.subplots(
        2,
        1,
        figsize=(7.0, 7.5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )

    ax_prof.errorbar(
        merged["radius_hmpc"],
        merged["kappa_direct"],
        yerr=sigma,
        fmt="o",
        ms=4,
        capsize=2,
        color="tab:orange",
        label=f"Direct annular counts ({n_halos} halos, Poisson errors)",
    )
    ax_prof.plot(
        map_profile["radius_hmpc"],
        map_profile["kappa"],
        "-",
        color="tab:blue",
        lw=2,
        label="Map-sampled stack (all halos)",
    )
    ax_prof.axhline(0.0, color="gray", lw=0.8, ls=":")
    ax_prof.set_yscale("symlog", linthresh=1e-4)
    ax_prof.set_ylabel(r"$\kappa$")
    ax_prof.legend(loc="upper right", fontsize=9)
    ax_prof.set_title("Single-halo profile validation: map sampling vs direct particle counts")

    ax_pull.axhspan(-2, 2, color="gray", alpha=0.15)
    ax_pull.axhspan(-1, 1, color="gray", alpha=0.25)
    ax_pull.axhline(0.0, color="gray", lw=0.8)
    ax_pull.plot(merged["radius_hmpc"], pull, "o", ms=4, color="tab:orange")
    ax_pull.set_xlabel(r"$r_p$ [$h^{-1}$ Mpc]")
    ax_pull.set_ylabel(r"(direct $-$ map) / $\sigma_{\rm Poisson}$")
    ax_pull.set_ylim(-6, 6)

    ensure_parent(output)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved validation figure -> %s", output)

    inner = merged[merged["radius_hmpc"] <= 10]
    outer = merged[merged["radius_hmpc"] > 30]
    logger.info(
        "Inner r<=10: map mean=%.3e, direct mean=%.3e (ratio %.2f); outer r>30 pull RMS=%.1f",
        inner["kappa_map"].mean(),
        inner["kappa_direct"].mean(),
        inner["kappa_direct"].mean() / inner["kappa_map"].mean(),
        np.sqrt((pull[merged["radius_hmpc"] > 30] ** 2).mean()) if len(outer) else float("nan"),
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot map-sampled vs direct-annular single-halo profile validation figure.",
    )
    parser.add_argument(
        "--map-profile",
        default="analysis/sim/results/radial_profile_single_sim_mass13.csv",
    )
    parser.add_argument(
        "--direct-profile",
        default="analysis/sim/results/radial_profile_single_direct_annuli_mass13.csv",
    )
    parser.add_argument(
        "--output",
        default="analysis/sim/results/validation_single_profile_map_vs_direct.png",
    )
    parser.add_argument(
        "--n-halos",
        type=int,
        default=500,
        help="Halo count used by the direct-annular run (not stored in its CSV).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    map_profile = pd.read_csv(args.map_profile)
    direct_profile = pd.read_csv(args.direct_profile)
    make_figure(map_profile, direct_profile, args.n_halos, args.output)


if __name__ == "__main__":
    main()
