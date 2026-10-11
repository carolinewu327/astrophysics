#!/usr/bin/env python
"""How clean is the r_par <= 5 pair selection in the wide bins?  (mock only)

Pairs are selected on the redshift-space line-of-sight separation, which
includes peculiar velocities.  Using the mock's true positions, for each wide
bin and each candidate cut C this reports

    completeness         of the truly associated pairs (true |r_par| <= 20
                         h^-1 Mpc), the fraction with redshift-space r_par <= C
    interloper_fraction  of the pairs selected with redshift-space r_par <= C,
                         the fraction that is not truly associated

with the definitions and the table code of los_contamination_table.py.

No new pair catalogue is needed.  The deep real-space catalogue the null test
builds (widebin_null_test.py: true |r_par| <= 120 for a seeded quarter of the
mock galaxies) carries both separations for every pair, so the truth sample and
every redshift-space selection are subsets of it.  That holds as long as no
selected pair sits near the 120 h^-1 Mpc edge; the script checks this and
records the largest true separation among the selected pairs.

Writes results/hod_pairs/widebin_rpar5/los_purity_wide.csv.

Usage (from the repository root; seconds once the null test has been run)
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python analysis/sim/widebin_los_purity.py
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from los_contamination_table import build_table
from sim_utils import setup_logging

logger = logging.getLogger(__name__)

BINS = ["g4_6", "g6_15", "g15_25"]
RPERP_OUTER = 25.0   # the only closed upper edge; every other bin is [lo, hi)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nulltest-dir", default="analysis/sim/results/hod_pairs/widebin_rpar5/nulltest")
    p.add_argument("--pairs-template",
                   default="pairs_hodnmatch50_frac25_seed12345_rperp{lo:g}_{hi:g}_real120.csv")
    p.add_argument("--deep-edge", type=float, default=120.0,
                   help="True line-of-sight extent of the deep catalogue.")
    p.add_argument("--cuts", nargs="+", type=float, default=[5.0, 10.0, 20.0, 30.0])
    p.add_argument("--assoc-rpar", type=float, default=20.0,
                   help="A pair is truly associated if its true |r_par| is at most this.")
    p.add_argument("--output", default="analysis/sim/results/hod_pairs/widebin_rpar5/los_purity_wide.csv")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    setup_logging()

    out = Path(args.output)
    if out.exists() and not args.overwrite:
        raise FileExistsError(f"{out} exists (use --overwrite)")

    tables = []
    for name in BINS:
        lo, hi = (float(v) for v in name[1:].split("_"))
        path = Path(args.nulltest_dir) / args.pairs_template.format(lo=lo, hi=hi)
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing. It is built by: python analysis/sim/widebin_null_test.py --bin {name}")
        df = pd.read_csv(path, usecols=["r_perp", "r_parallel_real", "r_parallel_rsd"])
        keep = (df.r_perp >= lo) & ((df.r_perp <= hi) if hi == RPERP_OUTER else (df.r_perp < hi))
        df = df[keep]
        if (df.r_parallel_real < 0).any() or (df.r_parallel_rsd < 0).any():
            raise ValueError(f"{path}: separations are expected as absolute values")

        selected = df[df.r_parallel_rsd <= max(args.cuts)]
        deepest = float(selected.r_parallel_real.max())
        if deepest > args.deep_edge - 20.0:
            raise ValueError(
                f"{name}: a pair selected at r_par <= {max(args.cuts):g} has true separation {deepest:.1f}, "
                f"within 20 of the catalogue edge ({args.deep_edge:g}); the catalogue is not deep enough.")

        truth = df[df.r_parallel_real <= args.assoc_rpar]
        table = build_table(truth, df, args.cuts, args.assoc_rpar, assoc_3d_factor=2.0)
        table.insert(0, "grouping", name)
        table["n_truth_pairs"] = len(truth)
        table["max_true_los_of_selected_hmpc"] = deepest
        tables.append(table)
        for r in table.itertuples():
            logger.info("%-7s r_par <= %-3g  %8d pairs   completeness %5.1f%%   interlopers %5.1f%%   "
                        "median true separation %.1f", name, r.rpar_cut_hmpc, r.pair_count,
                        100 * r.completeness, 100 * r.interloper_fraction, r.true_los_p50_hmpc)

    pd.concat(tables, ignore_index=True).to_csv(out, index=False)
    logger.info("Saved %s", out)


if __name__ == "__main__":
    main()
