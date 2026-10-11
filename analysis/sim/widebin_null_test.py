#!/usr/bin/env python
"""Non-physical pair null test for a wide r_perp bin with the separation-averaged control.

The wide-bin counterpart of nonphysical_pair_test.py.  Mock pairs close on the
sky but far apart along the line of sight (real space) are not physically
associated, so they hold no filament: pair stack minus control should leave a
bridge excess of zero.  Whatever it leaves instead is control bias.

Three steps, each skipped if its output already exists:

  1. draw a random fraction of the mock galaxies (default 25 per cent, seeded);
  2. find their pairs in the r_perp bin out to 120 h^-1 Mpc along the line of
     sight in real space (find_pairs_sim.py --rpar-space real);
  3. stack each line-of-sight band with jackknife_pair_stack.py
     --control separation_averaged and collect the fixed-band bridge excess.

A quarter of the galaxies is one sixteenth of the pairs.  That is deliberate:
the stacker reads the whole pair file into memory (the full 6-15 catalogue is
~22M pairs, ~6 GB), and the mock error is set by the structure in the simulated
map rather than the pair count, so the full sample would add little.

Writes to results/hod_pairs/widebin_rpar5/nulltest/:
    nonphysical_pair_test_<bin>.csv           the table (one row per band)
    nonphys_<bin>_rpar<lo>_<hi>{.json,_map.csv,_blocks.npz}   per band
    halos_*_frac*.csv, pairs_*_real120.csv    intermediate, large, not in git

Usage (from the repository root, about 20 minutes per bin)
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python \\
        analysis/sim/widebin_null_test.py --bin g6_15
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

import find_pairs_sim
import jackknife_pair_stack
from combine_filament_jackknife import BRIDGE_CENTERS
from sim_utils import setup_logging

logger = logging.getLogger(__name__)

DEFAULT_BANDS = "0:20,20:40,40:60,60:90,90:120"
RPERP_OUTER = 25.0   # the only closed upper edge; every other bin is [lo, hi)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bin", required=True, help="Grouping name, e.g. g6_15.")
    p.add_argument("--halos", default="analysis/sim/results/galaxies_hod_nmatch_frac50.csv")
    p.add_argument("--halo-fraction", type=float, default=0.25)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--rpar-max", type=float, default=120.0)
    p.add_argument("--bands", default=DEFAULT_BANDS,
                   help="Comma-separated lo:hi real-space line-of-sight bands in h^-1 Mpc.")
    p.add_argument("--kappa-map", default="analysis/sim/results/kappa_map_l0p1_s8arcmin.float32")
    p.add_argument("--single",
                   default="analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv")
    p.add_argument("--output-dir", default="analysis/sim/results/hod_pairs/widebin_rpar5/nulltest")
    p.add_argument("--overwrite", action="store_true", help="Replace existing band stacks and the table.")
    args = p.parse_args(argv)
    setup_logging()

    lo, hi = (float(v) for v in args.bin[1:].split("_"))
    center = BRIDGE_CENTERS[args.bin]
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table = out / f"nonphysical_pair_test_{args.bin}.csv"
    if table.exists() and not args.overwrite:
        raise FileExistsError(f"{table} exists (use --overwrite)")

    pct = round(100 * args.halo_fraction)
    halos = out / f"halos_hodnmatch50_frac{pct}_seed{args.seed}.csv"
    if not halos.exists():
        full = pd.read_csv(args.halos)
        full.sample(frac=args.halo_fraction, random_state=args.seed).to_csv(halos, index=False)
        logger.info("Drew %d%% of %d mock galaxies -> %s", pct, len(full), halos)

    pairs = out / f"pairs_hodnmatch50_frac{pct}_seed{args.seed}_rperp{lo:g}_{hi:g}_real{args.rpar_max:g}.csv"
    if not pairs.exists():
        find_pairs_sim.main(["--halos", str(halos), "--output", str(pairs),
                             "--rperp-min", f"{lo:g}", "--rperp-max", f"{hi:g}",
                             "--rpar-max", f"{args.rpar_max:g}", "--rpar-space", "real"])

    rows = []
    for band in args.bands.split(","):
        b_lo, b_hi = (float(v) for v in band.split(":"))
        stem = out / f"nonphys_{args.bin}_rpar{b_lo:g}_{b_hi:g}"
        logger.info("=== %s, line-of-sight band %g-%g h^-1 Mpc (real space) ===", args.bin, b_lo, b_hi)
        jackknife_pair_stack.main(
            ["--pairs", str(pairs), "--kappa-map", args.kappa_map, "--single", args.single,
             "--rperp-min", f"{lo:g}", "--rperp-max", f"{hi:g}", "--rperp-center", f"{center:g}",
             "--rpar-min", f"{b_lo:g}", "--rpar-max", f"{b_hi:g}", "--rpar-half-open",
             "--rpar-space", "real", "--control", "separation_averaged",
             "--output", f"{stem}.json", "--stack-output", f"{stem}_map.csv",
             "--blocks-output", f"{stem}_blocks.npz"]
            + ([] if hi == RPERP_OUTER else ["--rperp-half-open"])
            + (["--overwrite"] if args.overwrite else []))
        with open(f"{stem}.json", encoding="utf-8") as fh:
            result = json.load(fh)
        raw = result["stats"]["raw_bridge_excess_fixedband_kappa"]
        res = result["stats"]["residual_bridge_excess_fixedband_kappa"]
        rows.append({"grouping": args.bin, "rpar_lo_hmpc": b_lo, "rpar_hi_hmpc": b_hi,
                     "n_pairs": result["n_pairs"], "mean_rperp_hmpc": result["mean_rperp_hmpc"],
                     "raw_bridge_excess": raw["value"], "raw_bridge_err": raw["jackknife_error"],
                     "residual_bridge_excess": res["value"], "residual_bridge_err": res["jackknife_error"],
                     "residual_significance": res["value"] / res["jackknife_error"],
                     "halo_fraction": args.halo_fraction, "seed": args.seed, "control": result["control"]})

    df = pd.DataFrame(rows)
    df.to_csv(table, index=False)
    logger.info("\n%s bridge excess after control (kappa x 1e4):", args.bin)
    for r in rows:
        logger.info("  %3g-%-3g  %8d pairs   raw %6.2f   residual %6.2f +/- %.2f  (%+.1f sigma)",
                    r["rpar_lo_hmpc"], r["rpar_hi_hmpc"], r["n_pairs"], r["raw_bridge_excess"] * 1e4,
                    r["residual_bridge_excess"] * 1e4, r["residual_bridge_err"] * 1e4,
                    r["residual_significance"])
    logger.info("Saved %s", table)


if __name__ == "__main__":
    main()
