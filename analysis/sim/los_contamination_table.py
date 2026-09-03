#!/usr/bin/env python
"""Phase 3 LOS-cut decision table: completeness and interloper fractions.

Inputs are two pair catalogs from ``find_pairs_sim.py`` over the same halo
sample and transverse bin:

- a *truth* catalog found with ``--rpar-space real --rpar-max <assoc-rpar>``:
  every pair satisfying the predefined real-space association criterion,
  regardless of where RSD scatters it;
- a *deep RSD* catalog found with ``--rpar-space redshift`` and a cut well
  above the largest candidate (e.g. 50), from which nested candidate-cut
  subsets are taken.

For each candidate cut C the table reports:

- ``completeness``: fraction of truth pairs with ``r_parallel_rsd <= C``;
- ``interloper_fraction``: fraction of RSD-selected pairs (``r_parallel_rsd
  <= C``) failing the association criterion in real space;
- the same two quantities under the robustness criterion (true 3D
  separation <= ``--assoc-3d-factor`` x r_perp);
- pair counts and percentiles of the true LOS separation distribution.

The decision rule (preregistered): choose the smallest C with completeness
>= the pre-agreed target and interloper_fraction <= the pre-agreed ceiling.
This script computes the table only; it does not choose.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from sim_utils import ensure_parent, setup_logging


logger = logging.getLogger(__name__)


def load_catalog(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=["r_perp", "r_parallel_real", "r_parallel_rsd"])
    logger.info("%s catalog: %d pairs from %s", label, len(df), path)
    return df


def classify_associated(
    df: pd.DataFrame,
    assoc_rpar: float,
    assoc_3d_factor: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Primary and robustness association flags for each pair."""
    primary = df["r_parallel_real"].to_numpy() <= assoc_rpar
    r3d = np.hypot(df["r_perp"].to_numpy(), df["r_parallel_real"].to_numpy())
    robustness = r3d <= assoc_3d_factor * df["r_perp"].to_numpy()
    return primary, robustness


def build_table(
    truth: pd.DataFrame,
    rsd: pd.DataFrame,
    cuts: list[float],
    assoc_rpar: float,
    assoc_3d_factor: float,
) -> pd.DataFrame:
    truth_rsd_sep = truth["r_parallel_rsd"].to_numpy()
    n_truth = len(truth)
    truth_primary, truth_robust = classify_associated(truth, assoc_rpar, assoc_3d_factor)
    # The truth catalog IS the primary-criterion sample by construction;
    # under the robustness criterion only a subset of it qualifies.
    n_truth_robust = int(truth_robust.sum())

    rsd_primary, rsd_robust = classify_associated(rsd, assoc_rpar, assoc_3d_factor)
    rsd_sep = rsd["r_parallel_rsd"].to_numpy()

    rows = []
    for cut in cuts:
        in_cut = rsd_sep <= cut
        n_sel = int(in_cut.sum())
        truth_in = truth_rsd_sep <= cut
        row = {
            "rpar_cut_hmpc": cut,
            "pair_count": n_sel,
            "completeness": float(truth_in.mean()) if n_truth else np.nan,
            "interloper_fraction": float((~rsd_primary[in_cut]).mean()) if n_sel else np.nan,
            "completeness_3d": (
                float(truth_in[truth_robust].mean()) if n_truth_robust else np.nan
            ),
            "interloper_fraction_3d": float((~rsd_robust[in_cut]).mean()) if n_sel else np.nan,
        }
        if n_sel:
            true_los = rsd["r_parallel_real"].to_numpy()[in_cut]
            for q in (50, 90, 99):
                row[f"true_los_p{q}_hmpc"] = float(np.percentile(true_los, q))
        rows.append(row)
    return pd.DataFrame(rows)


def boundary_checks(truth: pd.DataFrame, rsd_max_available: float, cuts: list[float]) -> dict:
    """How much truth lies beyond the largest candidate cut / the deep catalog edge."""
    sep = truth["r_parallel_rsd"].to_numpy()
    largest_cut = max(cuts)
    return {
        "truth_pairs": int(len(truth)),
        "truth_beyond_largest_cut": int((sep > largest_cut).sum()),
        "truth_beyond_largest_cut_fraction": float((sep > largest_cut).mean()),
        "truth_within_5hmpc_of_deep_edge": int((sep > rsd_max_available - 5.0).sum()),
        "deep_catalog_rpar_edge_hmpc": rsd_max_available,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Phase 3 LOS-cut completeness/contamination table.")
    parser.add_argument("--truth-pairs", type=Path, required=True,
                        help="Catalog from --rpar-space real --rpar-max <assoc-rpar>.")
    parser.add_argument("--rsd-pairs", type=Path, required=True,
                        help="Deep catalog from --rpar-space redshift (cut >> largest candidate).")
    parser.add_argument("--cuts", nargs="+", type=float, default=[5.0, 10.0, 20.0, 30.0])
    parser.add_argument("--assoc-rpar", type=float, default=20.0,
                        help="Primary association criterion: true |r_parallel| <= this (h^-1 Mpc). "
                             "Must match the truth catalog's selection cut.")
    parser.add_argument("--assoc-3d-factor", type=float, default=2.0,
                        help="Robustness criterion: true 3D separation <= factor x r_perp.")
    parser.add_argument("--output", type=Path,
                        default=Path("analysis/sim/results/los_sweep/los_contamination_table.csv"))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()

    truth = load_catalog(args.truth_pairs, "truth")
    rsd = load_catalog(args.rsd_pairs, "deep RSD")

    max_truth_real = float(truth["r_parallel_real"].max())
    if max_truth_real > args.assoc_rpar + 1e-6:
        raise ValueError(
            f"truth catalog contains r_parallel_real up to {max_truth_real:.2f} > "
            f"--assoc-rpar {args.assoc_rpar}; it was not built with the association criterion"
        )

    table = build_table(truth, rsd, sorted(args.cuts), args.assoc_rpar, args.assoc_3d_factor)
    checks = boundary_checks(truth, float(rsd["r_parallel_rsd"].max()), args.cuts)

    ensure_parent(args.output)
    table.to_csv(args.output, index=False)
    meta = {
        "truth_pairs": str(args.truth_pairs),
        "rsd_pairs": str(args.rsd_pairs),
        "assoc_rpar_hmpc": args.assoc_rpar,
        "assoc_3d_factor": args.assoc_3d_factor,
        "cuts_hmpc": sorted(args.cuts),
        "boundary_checks": checks,
    }
    with open(f"{args.output}.meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)
        f.write("\n")

    logger.info("Boundary checks: %s", checks)
    print(table.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    logger.info("Wrote %s", args.output)


if __name__ == "__main__":
    main()
