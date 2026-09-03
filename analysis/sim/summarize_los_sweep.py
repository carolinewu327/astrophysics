#!/usr/bin/env python
"""Summarize the Phase 3 LOS sweep: residual bridge stats per (cut, smoothing).

Reads the manifest written by ``run_los_sweep.py``. For each stack it builds
the two-halo template by placing two copies of the *matched* single-halo
stack (same map/smoothing) at +/- rperp_center/2 on the fixed physical axis,
subtracts it from the pair stack, and reports bridge statistics of the
residual using the same bridge/sideband definitions as
``summarize_sim_sensitivity.py`` (regions scale with the pair separation).

These stacks document the consequence of the LOS choice; the choice itself
comes from the preregistered completeness/interloper rule applied to
``los_contamination_table.py`` output.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from sim_utils import make_two_halo_template, setup_logging
from summarize_sim_sensitivity import count_pairs, load_map, map_stats


logger = logging.getLogger(__name__)


def make_fixed_template(single: np.ndarray, pair_axis: np.ndarray, rperp_center: float) -> np.ndarray:
    # radial-profile construction: exact for the symmetrized single map,
    # symmetric by construction, and covers radii out to the map corners
    x_grid, y_grid = np.meshgrid(pair_axis, pair_axis)
    return make_two_halo_template(single, x_grid, y_grid, rperp_center)


def summarize(manifest_path: Path) -> pd.DataFrame:
    with manifest_path.open("r", encoding="utf-8") as f:
        manifest = json.load(f)
    rperp_center = float(manifest["rperp_center_hmpc"])

    singles: dict[str, np.ndarray] = {}
    rows = []
    for stack_cfg in manifest["stacks"]:
        stack_path = Path(stack_cfg["stack_path"])
        if not stack_path.exists():
            logger.warning("Skipping missing stack %s", stack_path)
            continue
        label = stack_cfg["smoothing_label"]
        if label not in singles:
            singles[label] = load_map(Path(stack_cfg["single_path"]))

        pair_map = load_map(stack_path)
        n = pair_map.shape[0]
        pair_axis = np.linspace(-50.0, 50.0, n)
        template = make_fixed_template(singles[label], pair_axis, rperp_center)
        residual = pair_map - template

        # bridge/sideband regions are defined in units of the pair separation
        norm_axis = pair_axis / rperp_center
        row = {
            "smoothing_label": label,
            "rpar_cut_hmpc": stack_cfg["rpar_cut_hmpc"],
            "pair_count": count_pairs(Path(stack_cfg["pair_path"])),
        }
        row.update(map_stats(pair_map, norm_axis, "raw"))
        row.update(map_stats(residual, norm_axis, "residual"))
        rows.append(row)

    return pd.DataFrame(rows).sort_values(["smoothing_label", "rpar_cut_hmpc"]).reset_index(drop=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize LOS sweep residual bridge statistics.")
    parser.add_argument("--manifest", type=Path,
                        default=Path("analysis/sim/results/los_sweep/manifest.json"))
    parser.add_argument("--output", type=Path,
                        default=Path("analysis/sim/results/los_sweep/los_sweep_summary.csv"))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    stats = summarize(args.manifest)
    stats.to_csv(args.output, index=False)
    print(stats.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    logger.info("Wrote %s", args.output)


if __name__ == "__main__":
    main()
