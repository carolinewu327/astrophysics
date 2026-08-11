#!/usr/bin/env python
"""Verify the compact host-halo binary against the mass13 reference catalog.

Safety gate before deleting the raw Rockstar hlist: applying the mass13
selection to the compact ``.bin`` catalog must reproduce ``halos_mass13.csv``
(which was extracted directly from the raw hlist) row for row. Positions,
velocities, and Mvir in the binary are float32, so values are compared at
float32 tolerance rather than exactly.

Usage
-----
    python verify_compact_halos.py
    python verify_compact_halos.py --compact results/hosts_minmass12.bin \
        --reference results/halos_mass13.csv
"""

from __future__ import annotations

import argparse
import json
import logging
import sys

import numpy as np
import pandas as pd

from sim_utils import setup_logging


logger = logging.getLogger(__name__)

FLOAT_COLUMNS = ["Mvir", "x", "y", "z", "vx", "vy", "vz"]


def load_compact(path: str) -> tuple[np.ndarray, dict]:
    with open(f"{path}.json", "r", encoding="utf-8") as f:
        info = json.load(f)
    dtype = np.dtype([tuple(field) for field in info["dtype"]])
    arr = np.fromfile(path, dtype=dtype)
    if len(arr) != info["n_rows"]:
        raise AssertionError(
            f"Row count mismatch: file has {len(arr)}, sidecar says {info['n_rows']}"
        )
    return arr, info


def verify(compact_path: str, reference_path: str, log10_mvir: float, half_width: float) -> bool:
    arr, info = load_compact(compact_path)
    logger.info("Compact catalog: %d rows, selection %s", len(arr), info["selection"])

    ref = pd.read_csv(reference_path)
    logger.info("Reference catalog: %d rows", len(ref))

    ok = True

    # --- Global invariants of the compact catalog ---
    if not (arr["upid"] == -1).all():
        logger.error("FAIL: compact catalog contains non-host rows (upid != -1)")
        ok = False
    floor = info["selection"].get("min_log10_mvir")
    if floor is not None and np.log10(arr["Mvir"].min()) < floor - 1e-6:
        logger.error("FAIL: compact catalog violates its own mass floor")
        ok = False
    if len(np.unique(arr["id"])) != len(arr):
        logger.error("FAIL: duplicate halo ids in compact catalog")
        ok = False

    # --- Reproduce the mass13 selection from the compact catalog ---
    # The reference cut ran on full-precision text values; the binary holds
    # float32 Mvir. Halos within float32 rounding of the bin edges can land
    # on either side, so compare the interiors and allow edge stragglers.
    log_m = np.log10(arr["Mvir"].astype(np.float64))
    sel = arr[np.abs(log_m - log10_mvir) < half_width]
    logger.info("Compact mass13 selection: %d rows (reference %d)", len(sel), len(ref))

    sel_ids = set(sel["id"].tolist())
    ref_ids = set(ref["id"].astype(np.int64).tolist())
    only_ref = ref_ids - sel_ids
    only_sel = sel_ids - ref_ids
    edge_tol = 4e-7  # float32 relative rounding of Mvir in log10
    for label, ids in [("reference-only", only_ref), ("compact-only", only_sel)]:
        if not ids:
            continue
        if label == "reference-only":
            masses = ref.set_index("id").loc[sorted(ids), "Mvir"].to_numpy()
        else:
            lookup = {int(i): m for i, m in zip(sel["id"], sel["Mvir"])}
            masses = np.array([lookup[i] for i in sorted(ids)])
        edge_dist = np.abs(np.abs(np.log10(masses) - log10_mvir) - half_width)
        if edge_dist.max() > edge_tol:
            logger.error(
                "FAIL: %d %s ids are not bin-edge rounding cases (max edge distance %.2e)",
                len(ids),
                label,
                edge_dist.max(),
            )
            ok = False
        else:
            logger.info(
                "%d %s ids are all within float32 rounding of the bin edge (max %.2e) -- acceptable",
                len(ids),
                label,
                edge_dist.max(),
            )

    # --- Value agreement on the common ids ---
    common = np.array(sorted(sel_ids & ref_ids), dtype=np.int64)
    ref_common = ref.set_index("id").loc[common]
    sel_df = pd.DataFrame({name: sel[name] for name in ["id"] + FLOAT_COLUMNS}).set_index("id").loc[common]
    for col in FLOAT_COLUMNS:
        ref_vals = ref_common[col].to_numpy(dtype=np.float64)
        bin_vals = sel_df[col].to_numpy(dtype=np.float64)
        # float32 storage: relative tolerance 1e-6, absolute floor for velocities near 0
        if not np.allclose(bin_vals, ref_vals, rtol=2e-6, atol=1e-3):
            worst = np.abs(bin_vals - ref_vals).max()
            logger.error("FAIL: column %s disagrees beyond float32 tolerance (max abs diff %.3e)", col, worst)
            ok = False
    logger.info("Compared %d common rows across %s", len(common), FLOAT_COLUMNS)

    if ok:
        logger.info("PASS: compact catalog reproduces the mass13 reference.")
    else:
        logger.error("VERIFICATION FAILED -- do NOT delete the raw hlist.")
    return ok


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify compact host-halo binary against mass13 CSV.")
    parser.add_argument("--compact", default="analysis/sim/results/hosts_minmass12.bin")
    parser.add_argument("--reference", default="analysis/sim/results/halos_mass13.csv")
    parser.add_argument("--log10-mvir", type=float, default=13.0)
    parser.add_argument("--half-width", type=float, default=0.05)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    ok = verify(args.compact, args.reference, args.log10_mvir, args.half_width)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
