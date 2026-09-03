#!/usr/bin/env python
"""Phase 3 LOS-cut sweep: nested pair subsets stacked with matched templates.

Takes one deep redshift-space pair catalog (from ``find_pairs_sim.py``,
``--rpar-max`` well above the largest candidate cut), filters it into nested
``r_parallel_rsd <= C`` subsets, and stacks each subset on each configured
kappa map. Single-halo templates are produced once per map, so every pair
stack has a template from the *same* map/smoothing (the matched-template
rule). Fixed-separation stacks only — normalized stacking is a Phase 4
sensitivity axis, not part of the LOS decision.

Smoothing discipline: the configured maps are expected to be pre-smoothed
(``derive_kappa_maps.py smooth``); every stack runs with ``--smooth none``.

Writes ``manifest.json`` for ``summarize_los_sweep.py``. Default is a dry
run that prints commands; pass ``--execute`` to run.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd

from sim_utils import ensure_parent, setup_logging


logger = logging.getLogger(__name__)


def default_maps() -> list[tuple[str, str]]:
    return [
        ("none", "analysis/sim/results/kappa_map_l0p1.float32"),
        ("8arcmin", "analysis/sim/results/kappa_map_l0p1_s8arcmin.float32"),
    ]


def parse_map_specs(specs: list[str] | None) -> list[tuple[str, str]]:
    if not specs:
        return default_maps()
    out = []
    for spec in specs:
        label, _, path = spec.partition("=")
        if not path:
            raise ValueError(f"--map expects label=path, got {spec!r}")
        out.append((label, path))
    return out


def write_subsets(
    deep_pairs: Path,
    cuts: list[float],
    out_dir: Path,
    rperp_label: str,
) -> dict[float, Path]:
    df = pd.read_csv(deep_pairs)
    if "r_parallel_rsd" not in df.columns:
        raise ValueError(
            f"{deep_pairs} lacks r_parallel_rsd; regenerate with the current find_pairs_sim.py"
        )
    edge = float(df["r_parallel_rsd"].max())
    if edge <= max(cuts):
        raise ValueError(
            f"deep catalog edge {edge:.1f} <= largest cut {max(cuts):g}; "
            "the largest subset would be the whole catalog by construction"
        )
    logger.info("Deep catalog: %d pairs, r_parallel_rsd edge %.2f", len(df), edge)

    paths = {}
    for cut in cuts:
        subset = df[df["r_parallel_rsd"] <= cut]
        path = out_dir / f"pairs_{rperp_label}_rpar{cut:g}.csv"
        ensure_parent(path)
        subset.to_csv(path, index=False)
        logger.info("Subset rpar<=%g: %d pairs -> %s", cut, len(subset), path)
        paths[cut] = path
    return paths


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Phase 3 LOS-cut stack sweep.")
    parser.add_argument("--deep-pairs", type=Path, required=True,
                        help="Deep redshift-space pair catalog (cut >> largest candidate).")
    parser.add_argument("--cuts", nargs="+", type=float, default=[5.0, 10.0, 20.0, 30.0])
    parser.add_argument("--rperp-label", default="rperp10")
    parser.add_argument("--rperp-center", type=float, default=10.0)
    parser.add_argument("--map", action="append", dest="maps", metavar="LABEL=PATH",
                        help="Kappa map to stack on (repeatable). Default: none=<0.1 map>, "
                             "8arcmin=<pre-smoothed 0.1 map>. Maps must already be smoothed; "
                             "stacks always run --smooth none.")
    parser.add_argument("--halos", type=Path, default=Path("analysis/sim/results/halos_mass13.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("analysis/sim/results/los_sweep"))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--script-dir", type=Path, default=Path("analysis/sim"))
    parser.add_argument("--max-pairs-per-stack", type=int, default=None,
                        help="Pair cap per stack (timing pilots only).")
    parser.add_argument("--execute", action="store_true",
                        help="Run the commands. Default prints them and writes the manifest.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    maps = parse_map_specs(args.maps)
    cuts = sorted(args.cuts)

    subset_paths = write_subsets(args.deep_pairs, cuts, args.output_dir, args.rperp_label)

    commands: list[list[str]] = []
    singles: dict[str, Path] = {}
    for label, map_path in maps:
        single_out = args.output_dir / f"kappa_single_{label}.csv"
        singles[label] = single_out
        commands.append([
            args.python, str(args.script_dir / "stack_single_sim.py"),
            "--halos", str(args.halos),
            "--kappa-map", map_path,
            "--smooth", "none",
            "--output", str(single_out),
            "--profile-output", str(args.output_dir / f"radial_profile_single_{label}.csv"),
            "--overwrite",
        ])

    stacks: list[dict] = []
    for label, map_path in maps:
        for cut in cuts:
            stack_out = args.output_dir / f"kappa_pairs_{args.rperp_label}_rpar{cut:g}_{label}.csv"
            cmd = [
                args.python, str(args.script_dir / "stack_pairs_sim.py"),
                "--pairs", str(subset_paths[cut]),
                "--kappa-map", map_path,
                "--smooth", "none",
                "--output", str(stack_out),
                "--overwrite",
            ]
            if args.max_pairs_per_stack:
                cmd += ["--max-pairs", str(args.max_pairs_per_stack)]
            commands.append(cmd)
            stacks.append({
                "smoothing_label": label,
                "rpar_cut_hmpc": cut,
                "pair_path": str(subset_paths[cut]),
                "stack_path": str(stack_out),
                "single_path": str(singles[label]),
            })

    manifest = {
        "deep_pairs": str(args.deep_pairs),
        "rperp_label": args.rperp_label,
        "rperp_center_hmpc": args.rperp_center,
        "cuts_hmpc": cuts,
        "maps": {label: path for label, path in maps},
        "max_pairs_per_stack": args.max_pairs_per_stack,
        "stacks": stacks,
    }
    manifest_path = args.output_dir / "manifest.json"
    ensure_parent(manifest_path)
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")
    logger.info("Wrote %s (%d commands)", manifest_path, len(commands))

    for cmd in commands:
        logger.info("CMD: %s", " ".join(cmd))
        if args.execute:
            subprocess.run(cmd, check=True)

    if not args.execute:
        logger.info("Dry run only. Re-run with --execute to launch.")


if __name__ == "__main__":
    main()
