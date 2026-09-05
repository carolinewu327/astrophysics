#!/usr/bin/env python
"""Build a central-galaxy mock catalog from the compact host-halo catalog.

Implements the HOD procedure from Z. Zheng's 2026-07-14 notes: each host halo
hosts a central galaxy with probability

    <Ncen(Mh)> = 0.5 * [1 + erf( (ln Mh - ln Mcut) / (sqrt(2) * sigma) )]

using the White et al. (2011) CMASS best-fit parameters
log10(Mcut / h^-1 Msun) = 13.04 and sigma = 0.94 (width in ln M). The paper's
mock overrides Mcut on the command line; see the note on LOG10_MCUT_DEFAULT. A seeded
uniform draw per halo decides occupation; occupied halos are written with the
halo's position and peculiar velocity, in the same CSV schema as
``halos_mass13.csv`` so all downstream scripts consume the mock unchanged.

Known limitation (recorded in the sidecar): the input compact catalog is
floored at log10 Mvir >= 12.5, where <Ncen> ~ 0.093, so centrals in lower-mass
halos (~16% of the full expectation) cannot be realized. Satellites are
deliberately ignored at this stage (CMASS satellite fraction ~10%).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.special import erf

from sim_utils import ensure_parent, setup_logging


logger = logging.getLogger(__name__)

# White et al. (2011) published log10 Mcut = 13.04 for CMASS, and that is the
# default here so the code matches the citation. It is NOT the value the paper's
# mock uses. Applied to Rockstar Mvir on BigMDPL with the log10 Mvir >= 12.5
# floor, 13.04 gives nbar = 4.21e-4 (h/Mpc)^3, 40% above CMASS -- their fit was
# calibrated on a different halo definition and mass function.
#
# The paper's catalog is galaxies_hod_nmatch.csv, built 2026-07-20 with
# --log10-mcut 13.246. That threshold matches the CMASS *central* number
# density, not the total: this HOD places no satellites, so the target is
# (1 - f_sat) * 3e-4 = 2.7e-4 for f_sat = 0.10, and 13.246 solves it over
# hosts_minmass12p5.bin (realized nbar = 2.7027e-4). Solving for the full 3e-4
# would instead give 13.20. The tuning was done by hand and is not reproduced by
# any script here; the per-run record is the .meta.json sidecar written beside
# each output catalog.
LOG10_MCUT_DEFAULT = 13.04
SIGMA_LN_DEFAULT = 0.94
CMASS_NBAR_HMPC3 = 3.0e-4  # approximate CMASS number density, (h/Mpc)^3
OUTPUT_COLUMNS = ["id", "Mvir", "x", "y", "z", "vx", "vy", "vz"]


def mean_ncen(log10_mvir: np.ndarray, log10_mcut: float, sigma_ln: float) -> np.ndarray:
    """<Ncen(Mh)> per White et al. (2011) eq. 12; sigma is the width in ln M."""
    x = (log10_mvir - log10_mcut) * np.log(10.0) / (np.sqrt(2.0) * sigma_ln)
    return 0.5 * (1.0 + erf(x))


def load_compact_halos(path: str) -> tuple[np.ndarray, dict]:
    with open(f"{path}.json", "r", encoding="utf-8") as f:
        info = json.load(f)
    dtype = np.dtype([tuple(field) for field in info["dtype"]])
    arr = np.fromfile(path, dtype=dtype)
    if len(arr) != info["n_rows"]:
        raise AssertionError(f"{path}: row count {len(arr)} != sidecar {info['n_rows']}")
    return arr, info


def build_catalog(
    halos: np.ndarray,
    log10_mcut: float,
    sigma_ln: float,
    seed: int,
    fraction: float,
    fraction_seed: int,
) -> tuple[pd.DataFrame, dict]:
    log10_mvir = np.log10(halos["Mvir"].astype(np.float64))
    ncen = mean_ncen(log10_mvir, log10_mcut, sigma_ln)
    expected = float(ncen.sum())

    rng = np.random.default_rng(seed)
    occupied = rng.uniform(size=len(halos)) < ncen
    galaxies = halos[occupied]
    logger.info(
        "Occupied %d of %d halos (expected %.0f, Poisson sigma %.0f)",
        occupied.sum(), len(halos), expected, np.sqrt(expected),
    )

    stats = {
        "n_halos_in": int(len(halos)),
        "expected_centrals": expected,
        "realized_centrals": int(occupied.sum()),
    }

    if fraction < 1.0:
        keep = np.random.default_rng(fraction_seed).uniform(size=len(galaxies)) < fraction
        galaxies = galaxies[keep]
        stats["downsampled_to"] = int(len(galaxies))
        logger.info("Downsampled to %d galaxies (fraction %.3f)", len(galaxies), fraction)

    frame = pd.DataFrame({name: galaxies[name] for name in OUTPUT_COLUMNS})
    return frame, stats


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the HOD central-galaxy mock catalog.")
    parser.add_argument("--halo-catalog", default="analysis/sim/results/hosts_minmass12p5.bin")
    parser.add_argument("--output", default="analysis/sim/results/galaxies_hod_cen.csv")
    parser.add_argument("--log10-mcut", type=float, default=LOG10_MCUT_DEFAULT)
    parser.add_argument("--sigma-ln", type=float, default=SIGMA_LN_DEFAULT,
                        help="Occupation width in ln M (White et al. 2011 convention).")
    parser.add_argument("--seed", type=int, default=20260716,
                        help="Seed for the occupation draw (recorded in the sidecar).")
    parser.add_argument("--fraction", type=float, default=1.0,
                        help="Optional random downsampling fraction applied after occupation.")
    parser.add_argument("--fraction-seed", type=int, default=1,
                        help="Seed for the downsampling draw.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()

    if os.path.exists(args.output) and not args.overwrite:
        logger.info("Output exists: %s (use --overwrite)", args.output)
        return

    halos, halo_info = load_compact_halos(args.halo_catalog)
    frame, stats = build_catalog(
        halos, args.log10_mcut, args.sigma_ln, args.seed, args.fraction, args.fraction_seed
    )

    ensure_parent(args.output)
    frame.to_csv(args.output, index=False)

    box = 2500.0
    nbar = len(frame) / box**3
    meta = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "hod": {
            "model": "central-only erf (White et al. 2011 eq. 12)",
            "log10_mcut": args.log10_mcut,
            "sigma_ln": args.sigma_ln,
        },
        "seed": args.seed,
        "fraction": args.fraction,
        "fraction_seed": args.fraction_seed,
        "input_catalog": os.path.abspath(args.halo_catalog),
        "input_selection": halo_info.get("selection"),
        "columns": OUTPUT_COLUMNS,
        "stats": stats,
        "number_density_hmpc3": nbar,
        "cmass_reference_nbar_hmpc3": CMASS_NBAR_HMPC3,
        "ncen_at_input_floor": float(
            mean_ncen(np.array([12.5]), args.log10_mcut, args.sigma_ln)[0]
        ),
        "known_limitations": [
            "input halos floored at log10 Mvir >= 12.5; centrals in lower-mass halos "
            "are not realized (see ncen_at_input_floor for how much occupation the "
            "floor truncates for these HOD parameters)",
            "central galaxies only; CMASS satellite fraction ~10% is ignored",
            "halo masses are Rockstar Mvir; published HOD fits may use other "
            "definitions (e.g. M200b in Yuan et al. 2021)",
        ],
    }
    with open(f"{args.output}.meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)
        f.write("\n")

    logger.info(
        "Wrote %d mock galaxies -> %s (nbar %.3e (h/Mpc)^3 vs CMASS ~%.1e)",
        len(frame), args.output, nbar, CMASS_NBAR_HMPC3,
    )


if __name__ == "__main__":
    main()
