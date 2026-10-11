#!/usr/bin/env python
"""How does the mock bridge signal change with the line-of-sight cut, in the wide bins?

The data use redshift-space r_par <= 5 in every bin.  widebin_los_purity.py
shows that cut is pure but keeps under half of the associated pairs, and that
r_par <= 10 would nearly double the pairs with under 1 per cent interlopers.
This measures what the wider cuts do to the signal: for each wide bin and each
cut C it stacks the mock pairs with redshift-space r_par <= C, subtracts the
separation-averaged control, and records the bridge excess S(C) and the pair
count N(C), each also as a ratio to C = 5.

The ratio S(C)/S(5) gets a paired jackknife error: the cuts are nested and
share the same 25 blocks, so most of the block-to-block scatter cancels.

What this can and cannot say.  The mock has no CMB reconstruction noise, so its
own S/N is not the data's.  The table gives the signal change and the pair
count; `fom_sqrtN` = S(C)/S(5) x sqrt(N(C)/N(5)) is the gain in S/N if the
error fell as 1/sqrt(N).  That is an upper bound: pairs added by deepening the
cut sit on the same sky as the ones already there and share their Planck
noise.

No new pair catalogue: the cuts are subsets of the deep real-space catalogue
the null test builds (widebin_null_test.py; seeded quarter of the mock
galaxies), which widebin_los_purity.py checks is deep enough.

Writes results/hod_pairs/widebin_rpar5/los_signal/:
    los_signal_wide.csv                                the table
    los_<bin>_rsd<C>{.json,_map.csv,_blocks.npz}       per bin and cut (not in git)

Usage (from the repository root, about 40 minutes)
-----
    PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python analysis/sim/widebin_los_signal.py
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import jackknife_pair_stack
from combine_filament_jackknife import BRIDGE_CENTERS
from geometry import bridge_excess, reflect_symmetrize_map, separation_averaged_template, separation_mix
from sim_utils import setup_logging
from summarize_sim_sensitivity import load_map

logger = logging.getLogger(__name__)

BINS = ["g4_6", "g6_15", "g15_25"]
RPERP_OUTER = 25.0   # the only closed upper edge; every other bin is [lo, hi)
REFERENCE_CUT = 5.0


def loo_bridge(blocks, single, center):
    """Residual bridge excess for the full stack and each leave-one-block-out stack."""
    sums, counts = blocks["sums"], blocks["counts"]
    n = sums.shape[-1]
    axis = np.linspace(-(n - 1) / 2, (n - 1) / 2, n)
    gx, gy = np.meshgrid(axis, axis)

    def score(stack, exclude=None):
        kwargs = {} if exclude is None else {"exclude_region": exclude}
        seps, wts, _ = separation_mix(blocks["fine_w"], blocks["fine_wr"], **kwargs)
        residual = stack - separation_averaged_template(single, gx, gy, seps, wts)
        return bridge_excess(reflect_symmetrize_map(residual.astype(np.float32)).astype(np.float64),
                             axis, center)

    full = score(sums.sum(0) / counts.sum())
    keep = [b for b in range(len(counts)) if counts.sum() - counts[b] > 0]
    loo = np.array([score((sums.sum(0) - sums[b]) / (counts.sum() - counts[b]), b) for b in keep])
    return full, loo


def jk_err(loo):
    k = len(loo)
    return float(np.sqrt((k - 1) / k * np.sum((loo - loo.mean()) ** 2)))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bins", default=",".join(BINS))
    p.add_argument("--cuts", default="5,10,20", help="Redshift-space r_par cuts in h^-1 Mpc.")
    p.add_argument("--nulltest-dir", default="analysis/sim/results/hod_pairs/widebin_rpar5/nulltest")
    p.add_argument("--pairs-template",
                   default="pairs_hodnmatch50_frac25_seed12345_rperp{lo:g}_{hi:g}_real120.csv")
    p.add_argument("--kappa-map", default="analysis/sim/results/kappa_map_l0p1_s8arcmin.float32")
    p.add_argument("--single",
                   default="analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv")
    p.add_argument("--output-dir", default="analysis/sim/results/hod_pairs/widebin_rpar5/los_signal")
    p.add_argument("--overwrite", action="store_true", help="Replace existing stacks and the table.")
    args = p.parse_args(argv)
    setup_logging()

    cuts = [float(c) for c in args.cuts.split(",")]
    if REFERENCE_CUT not in cuts:
        raise ValueError(f"--cuts must include the reference cut {REFERENCE_CUT:g}")
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table = out / "los_signal_wide.csv"
    if table.exists() and not args.overwrite:
        raise FileExistsError(f"{table} exists (use --overwrite)")
    single = load_map(Path(args.single))

    rows = []
    for name in args.bins.split(","):
        lo, hi = (float(v) for v in name[1:].split("_"))
        center = BRIDGE_CENTERS[name]
        pairs = Path(args.nulltest_dir) / args.pairs_template.format(lo=lo, hi=hi)
        if not pairs.exists():
            raise FileNotFoundError(
                f"{pairs} is missing. It is built by: python analysis/sim/widebin_null_test.py --bin {name}")
        scored = {}
        for cut in cuts:
            stem = out / f"los_{name}_rsd{cut:g}"
            if not Path(f"{stem}.json").exists() or args.overwrite:
                logger.info("=== %s, redshift-space r_par <= %g ===", name, cut)
                jackknife_pair_stack.main(
                    ["--pairs", str(pairs), "--kappa-map", args.kappa_map, "--single", args.single,
                     "--rperp-min", f"{lo:g}", "--rperp-max", f"{hi:g}", "--rperp-center", f"{center:g}",
                     "--rpar-max", f"{cut:g}", "--rpar-space", "redshift",
                     "--control", "separation_averaged",
                     "--output", f"{stem}.json", "--stack-output", f"{stem}_map.csv",
                     "--blocks-output", f"{stem}_blocks.npz"]
                    + ([] if hi == RPERP_OUTER else ["--rperp-half-open"])
                    + (["--overwrite"] if args.overwrite else []))
            with open(f"{stem}.json", encoding="utf-8") as fh:
                result = json.load(fh)
            full, loo = loo_bridge(np.load(f"{stem}_blocks.npz"), single, center)
            stat = result["stats"]["residual_bridge_excess_fixedband_kappa"]
            if not np.isclose(full, stat["value"], rtol=1e-4, atol=1e-9):
                raise ValueError(f"{stem}: rebuilt bridge {full:.6e} != stacker's {stat['value']:.6e}")
            scored[cut] = (result, full, loo)

        ref_result, ref_full, ref_loo = scored[REFERENCE_CUT]
        for cut in cuts:
            result, full, loo = scored[cut]
            n_ratio = result["n_pairs"] / ref_result["n_pairs"]
            same_blocks = len(loo) == len(ref_loo)
            ratio = full / ref_full
            ratio_err = jk_err(loo / ref_loo) if same_blocks and cut != REFERENCE_CUT else 0.0
            rows.append({"grouping": name, "rpar_cut_hmpc": cut, "n_pairs": result["n_pairs"],
                         "mean_rperp_hmpc": result["mean_rperp_hmpc"],
                         "bridge": full, "bridge_err": jk_err(loo),
                         "n_ratio_to_cut5": n_ratio, "signal_ratio_to_cut5": ratio,
                         "signal_ratio_err": ratio_err, "fom_sqrtN": ratio * np.sqrt(n_ratio),
                         "fom_sqrtN_err": ratio_err * np.sqrt(n_ratio)})

    df = pd.DataFrame(rows)
    df.to_csv(table, index=False)
    logger.info("\nMock bridge excess versus line-of-sight cut (kappa x 1e4):")
    for r in df.itertuples():
        logger.info("  %-7s r_par <= %-3g %8d pairs (x%.2f)   S = %5.2f +/- %.2f   S/S(5) = %.2f +/- %.2f"
                    "   S/S(5) x sqrt(N/N5) = %.2f", r.grouping, r.rpar_cut_hmpc, r.n_pairs,
                    r.n_ratio_to_cut5, r.bridge * 1e4, r.bridge_err * 1e4, r.signal_ratio_to_cut5,
                    r.signal_ratio_err, r.fom_sqrtN)
    logger.info("Saved %s", table)


if __name__ == "__main__":
    main()
