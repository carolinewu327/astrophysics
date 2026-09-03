#!/usr/bin/env python
"""Non-physical pair null test: is the superposed-singles control unbiased?

Zheng's point 4.  The filament measurement is

    filament = corrected pair stack - control,

where the control is the single-halo stack superposed at +/- r_perp/2.  That
control is a *model* for what two unassociated galaxies at this projected
separation would produce.  If it is biased, the filament residual is biased --
and the filament is a small difference between two much larger maps, so a
modest error in the control is a large error in the answer.  (On BOSS, a 4.4e-4
shift in the control moved the filament bridge excess from 7.5e-4 to 11.9e-4.)

This tests the control empirically.  Pairs close in projection but far apart
along the line of sight are not physically associated, so they contain no
filament.  Stack them and subtract the same template: the residual bridge
excess should be **zero**.  Whatever it is instead is control bias.

Running a ladder of increasing line-of-sight separation -- Zheng's "try the
test with larger line-of-sight separation" -- shows how fast the residual
approaches zero as physical association is switched off, which separates
genuine control bias from residual real correlation at modest LOS separation.

The simulation's kappa map is a projection through the entire 2500 h^-1 Mpc
box, so LOS-separated pairs still land in the same projection; the test is not
limited by slab depth.

Usage
-----
    PYTHONPATH=../../lib python nonphysical_pair_test.py \\
        --pairs results/hod_pairs/pairs_hodnmatch50_rperp10_real120.csv \\
        --rperp-center 10 \\
        --single results/kappa_single_sim_hodnmatch_8arcmin.csv
"""

from __future__ import annotations

import argparse
import json
import logging
from argparse import Namespace
from pathlib import Path

import pandas as pd

from jackknife_pair_stack import run
from sim_utils import setup_logging

logger = logging.getLogger(__name__)

DEFAULT_BANDS = "0:20,20:40,40:60,60:90,90:120"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Residual bridge excess vs line-of-sight separation.")
    parser.add_argument("--pairs", required=True,
                        help="Deep real-space pair catalog (large --rpar-max).")
    parser.add_argument("--kappa-map",
                        default="results/kappa_map_l0p1_s8arcmin.float32")
    parser.add_argument("--single",
                        default="results/kappa_single_sim_hodnmatch_8arcmin.csv")
    parser.add_argument("--rperp-center", type=float, required=True)
    parser.add_argument("--bands", default=DEFAULT_BANDS,
                        help="Comma-separated lo:hi LOS bands in h^-1 Mpc, applied to "
                             "the true (real-space) separation.")
    parser.add_argument("--blocks-per-side", type=int, default=5)
    parser.add_argument("--output-dir", default="results/hod_pairs")
    parser.add_argument("--label", default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    setup_logging()
    label = args.label or f"rperp{int(args.rperp_center)}"
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for band in args.bands.split(","):
        lo, hi = (float(v) for v in band.split(":"))
        stem = out_dir / f"nonphys_{label}_rpar{int(lo)}_{int(hi)}"
        stack_args = Namespace(
            pairs=args.pairs, kappa_map=args.kappa_map, single=args.single,
            rperp_center=args.rperp_center, rpar_min=lo, rpar_max=hi,
            rpar_space="real", blocks_per_side=args.blocks_per_side,
            grid_size=101, box_size=100.0, max_pairs=None, seed=0,
            rpar_half_open=True,
            stack_output=f"{stem}_map.csv",
            # The per-block accumulators are the point of re-running this: the
            # scalar stats below apply only to the built-in field-single
            # template, so the recursive estimator would otherwise have no
            # uncertainty of its own.
            blocks_output=f"{stem}_blocks.npz",
        )
        logger.info("=== LOS band %.0f-%.0f h^-1 Mpc (real space) ===", lo, hi)
        result = run(stack_args)

        res = result["stats"]["residual_bridge_excess_kappa"]
        raw = result["stats"]["raw_bridge_excess_kappa"]
        rows.append({
            "rpar_lo_hmpc": lo, "rpar_hi_hmpc": hi,
            "n_pairs": result["n_pairs"],
            "raw_bridge_excess": raw["value"], "raw_bridge_err": raw["jackknife_error"],
            "residual_bridge_excess": res["value"],
            "residual_bridge_err": res["jackknife_error"],
            "residual_significance": (
                res["value"] / res["jackknife_error"]
                if res["jackknife_error"] > 0 else float("nan")),
        })
        with (out_dir / f"nonphys_{label}_rpar{int(lo)}_{int(hi)}.json").open(
                "w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write("\n")

    out = pd.DataFrame(rows)
    path = out_dir / f"nonphysical_pair_test_{label}.csv"
    out.to_csv(path, index=False)
    logger.info("Saved -> %s", path)

    logger.info(
        "\nResidual bridge excess vs LOS separation (kappa x 1e4).\n"
        "Physically associated pairs sit in the first band; the last bands are\n"
        "non-physical. A control with no bias drives the residual to zero there.\n")
    logger.info("  LOS band      pairs     raw      residual      err    sigma")
    for _, r in out.iterrows():
        logger.info("  %3.0f-%3.0f  %9d  %7.3f  %10.3f  %7.3f  %7.1f",
                    r.rpar_lo_hmpc, r.rpar_hi_hmpc, r.n_pairs,
                    r.raw_bridge_excess * 1e4, r.residual_bridge_excess * 1e4,
                    r.residual_bridge_err * 1e4, r.residual_significance)


if __name__ == "__main__":
    main()
