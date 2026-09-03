#!/usr/bin/env python
"""Generate or run BigMDPL simulation sensitivity variants.

The baseline simulation products are expensive, so this script separates
variant definition from execution. By default it writes a manifest and prints
the commands. Add ``--execute`` to run the selected variants.
"""

from __future__ import annotations

import argparse
import json
import logging
import shlex
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


logger = logging.getLogger(__name__)


RPERP_BINS = {
    "rperp5": (5.0, 4.0, 6.0),
    "rperp10": (10.0, 9.0, 11.0),
    "rperp20": (20.0, 19.0, 21.0),
}

SEPARATION_BIN_VARIANTS = {
    "narrow": 0.5,
    "wide": 2.0,
}

MASS_VARIANTS = {
    "mass12p9": (12.9, 0.05),
    "mass13p1": (13.1, 0.05),
    "mass13_wide": (13.0, 0.10),
}


@dataclass
class Variant:
    name: str
    kind: str
    parameters: dict
    single_path: Path
    pairs: dict[str, dict]
    commands: list[list[str]]

    def to_manifest(self) -> dict:
        return {
            "name": self.name,
            "kind": self.kind,
            "parameters": self.parameters,
            "single_path": str(self.single_path),
            "pairs": self.pairs,
            "commands": self.commands,
        }


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def py_cmd(args: argparse.Namespace, script_name: str, *script_args: str) -> list[str]:
    return [args.python, str(args.script_dir / script_name), *script_args]


def pair_outputs(out_dir: Path, mass_label: str, rperp_label: str) -> dict[str, Path]:
    return {
        "pairs": out_dir / f"pairs_{mass_label}_{rperp_label}.csv",
        "fixed": out_dir / f"kappa_pairs_sim_{mass_label}_{rperp_label}.csv",
        "normalized": out_dir / f"kappa_pairs_sim_{mass_label}_{rperp_label}_normalized.csv",
    }


def add_stack_pair_commands(
    commands: list[list[str]],
    args: argparse.Namespace,
    pair_path: Path,
    fixed_path: Path,
    normalized_path: Path,
    smooth: str,
    map_path: str | None = None,
) -> None:
    kappa_map = str(args.kappa_map) if map_path is None else map_path
    if not args.normalized_only:
        commands.append(
            py_cmd(
                args,
                "stack_pairs_sim.py",
                "--pairs",
                str(pair_path),
                "--kappa-map",
                kappa_map,
                "--output",
                str(fixed_path),
                "--smooth",
                smooth,
                "--overwrite",
            )
        )
    commands.append(
        py_cmd(
            args,
            "stack_pairs_sim.py",
            "--pairs",
            str(pair_path),
            "--kappa-map",
            kappa_map,
            "--output",
            str(normalized_path),
            "--normalize-separation",
            "--smooth",
            smooth,
            "--overwrite",
        )
    )


def make_baseline_variant(args: argparse.Namespace) -> Variant:
    pairs = {}
    for label, (center, _, _) in RPERP_BINS.items():
        paths = pair_outputs(args.results_dir, args.mass_label, label)
        pairs[label] = {
            "rperp_center_hmpc": center,
            "pair_path": str(paths["pairs"]),
            "fixed_path": str(paths["fixed"]),
            "normalized_path": str(paths["normalized"]),
        }
    return Variant(
        name="baseline",
        kind="baseline",
        parameters={"mass_label": args.mass_label, "smooth": "none", "rpar_max": 22.0},
        single_path=args.results_dir / f"kappa_single_sim_{args.mass_label}.csv",
        pairs=pairs,
        commands=[],
    )


def smoothed_map_for(args: argparse.Namespace, smooth: str) -> str:
    """Pre-smoothed map for a smoothing level (never smooth the 0.1 map in RAM)."""
    if smooth == "none":
        return str(args.kappa_map)
    stem, ext = str(args.kappa_map).rsplit(".", 1)
    return f"{stem}_s{smooth}.{ext}"


def make_smoothing_variant(args: argparse.Namespace, smooth: str) -> Variant:
    out_dir = args.output_dir / f"smooth_{smooth}"
    mass_label = args.mass_label
    map_path = smoothed_map_for(args, smooth)
    single_path = out_dir / f"kappa_single_sim_{mass_label}.csv"
    profile_path = out_dir / f"radial_profile_single_sim_{mass_label}.csv"
    commands = [
        py_cmd(
            args,
            "stack_single_sim.py",
            "--halos",
            str(args.halos),
            "--kappa-map",
            map_path,
            "--output",
            str(single_path),
            "--profile-output",
            str(profile_path),
            "--smooth",
            "none",
            "--overwrite",
        )
    ]
    pairs = {}
    for label, (center, _, _) in RPERP_BINS.items():
        paths = pair_outputs(out_dir, mass_label, label)
        baseline_pairs = args.results_dir / f"pairs_{mass_label}_{label}.csv"
        add_stack_pair_commands(
            commands, args, baseline_pairs, paths["fixed"], paths["normalized"], "none", map_path
        )
        pairs[label] = {
            "rperp_center_hmpc": center,
            "pair_path": str(baseline_pairs),
            "fixed_path": str(paths["fixed"]),
            "normalized_path": str(paths["normalized"]),
        }
    return Variant(
        name=f"smooth_{smooth}",
        kind="smoothing",
        parameters={"smooth": smooth, "mass_label": mass_label, "kappa_map": map_path},
        single_path=single_path,
        pairs=pairs,
        commands=commands,
    )


def make_los_variant(args: argparse.Namespace, rpar_max: float) -> Variant:
    out_dir = args.output_dir / f"los_rpar{rpar_max:g}"
    mass_label = args.mass_label
    commands: list[list[str]] = []
    pairs = {}
    for label, (center, rmin, rmax) in RPERP_BINS.items():
        paths = pair_outputs(out_dir, mass_label, label)
        commands.append(
            py_cmd(
                args,
                "find_pairs_sim.py",
                "--halos",
                str(args.halos),
                "--output",
                str(paths["pairs"]),
                "--rperp-min",
                str(rmin),
                "--rperp-max",
                str(rmax),
                "--rpar-max",
                str(rpar_max),
                "--overwrite",
            )
        )
        add_stack_pair_commands(commands, args, paths["pairs"], paths["fixed"], paths["normalized"], "none")
        pairs[label] = {
            "rperp_center_hmpc": center,
            "pair_path": str(paths["pairs"]),
            "fixed_path": str(paths["fixed"]),
            "normalized_path": str(paths["normalized"]),
        }
    return Variant(
        name=f"los_rpar{rpar_max:g}",
        kind="line_of_sight",
        parameters={"rpar_max_hmpc": rpar_max, "mass_label": mass_label},
        single_path=args.results_dir / f"kappa_single_sim_{mass_label}.csv",
        pairs=pairs,
        commands=commands,
    )


def make_separation_bin_variant(args: argparse.Namespace, name: str, half_width: float) -> Variant:
    out_dir = args.output_dir / f"sepbin_{name}"
    mass_label = args.mass_label
    commands: list[list[str]] = []
    pairs = {}
    for label, (center, _, _) in RPERP_BINS.items():
        paths = pair_outputs(out_dir, mass_label, label)
        rmin = center - half_width
        rmax = center + half_width
        commands.append(
            py_cmd(
                args,
                "find_pairs_sim.py",
                "--halos",
                str(args.halos),
                "--output",
                str(paths["pairs"]),
                "--rperp-min",
                str(rmin),
                "--rperp-max",
                str(rmax),
                "--rpar-max",
                str(args.baseline_rpar_max),
                "--overwrite",
            )
        )
        add_stack_pair_commands(commands, args, paths["pairs"], paths["fixed"], paths["normalized"], "none")
        pairs[label] = {
            "rperp_center_hmpc": center,
            "pair_path": str(paths["pairs"]),
            "fixed_path": str(paths["fixed"]),
            "normalized_path": str(paths["normalized"]),
            "rperp_min_hmpc": rmin,
            "rperp_max_hmpc": rmax,
        }
    return Variant(
        name=f"sepbin_{name}",
        kind="separation_bin",
        parameters={"half_width_hmpc": half_width, "mass_label": mass_label},
        # matched template: the smooth_none variant's single stack is built on
        # the same (unsmoothed) map these pair stacks use — run smoothing first
        single_path=args.output_dir / "smooth_none" / f"kappa_single_sim_{mass_label}.csv",
        pairs=pairs,
        commands=commands,
    )


def make_mass_variant(args: argparse.Namespace, mass_label: str, log10_mvir: float, half_width: float) -> Variant:
    out_dir = args.output_dir / mass_label
    halos_path = out_dir / f"halos_{mass_label}.csv"
    single_path = out_dir / f"kappa_single_sim_{mass_label}.csv"
    profile_path = out_dir / f"radial_profile_single_sim_{mass_label}.csv"
    commands = [
        py_cmd(
            args,
            "prepare_halos.py",
            "--halo-catalog",
            str(args.rockstar_halos),
            "--output",
            str(halos_path),
            "--log10-mvir",
            str(log10_mvir),
            "--half-width",
            str(half_width),
        ),
        py_cmd(
            args,
            "stack_single_sim.py",
            "--halos",
            str(halos_path),
            "--kappa-map",
            str(args.kappa_map),
            "--output",
            str(single_path),
            "--profile-output",
            str(profile_path),
            "--overwrite",
        ),
    ]
    pairs = {}
    for label, (center, rmin, rmax) in RPERP_BINS.items():
        paths = pair_outputs(out_dir, mass_label, label)
        commands.append(
            py_cmd(
                args,
                "find_pairs_sim.py",
                "--halos",
                str(halos_path),
                "--output",
                str(paths["pairs"]),
                "--rperp-min",
                str(rmin),
                "--rperp-max",
                str(rmax),
                "--rpar-max",
                str(args.baseline_rpar_max),
                "--overwrite",
            )
        )
        add_stack_pair_commands(commands, args, paths["pairs"], paths["fixed"], paths["normalized"], "none")
        pairs[label] = {
            "rperp_center_hmpc": center,
            "pair_path": str(paths["pairs"]),
            "fixed_path": str(paths["fixed"]),
            "normalized_path": str(paths["normalized"]),
        }
    return Variant(
        name=mass_label,
        kind="mass_cut",
        parameters={"log10_mvir": log10_mvir, "half_width_dex": half_width},
        single_path=single_path,
        pairs=pairs,
        commands=commands,
    )


def selected_groups(groups: list[str]) -> set[str]:
    group_set = set(groups)
    if "all" in group_set:
        return {"smoothing", "los", "separation-bin", "mass"}
    return group_set


def build_variants(args: argparse.Namespace) -> list[Variant]:
    groups = selected_groups(args.groups)
    variants = [make_baseline_variant(args)]
    if "smoothing" in groups:
        variants.extend(make_smoothing_variant(args, smooth) for smooth in args.smooth_values)
    if "los" in groups:
        variants.extend(make_los_variant(args, rpar) for rpar in args.rpar_values)
    if "separation-bin" in groups:
        variants.extend(
            make_separation_bin_variant(args, name, half_width)
            for name, half_width in SEPARATION_BIN_VARIANTS.items()
        )
    if "mass" in groups:
        variants.extend(
            make_mass_variant(args, label, log10_mvir, half_width)
            for label, (log10_mvir, half_width) in MASS_VARIANTS.items()
        )
    return variants


def write_manifest(args: argparse.Namespace, variants: list[Variant]) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "execute": args.execute,
        "baseline_results_dir": str(args.results_dir),
        "kappa_map": str(args.kappa_map),
        "variants": [variant.to_manifest() for variant in variants],
    }
    args.manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    logger.info("Wrote manifest: %s", args.manifest)


def run_commands(variants: list[Variant], execute: bool) -> None:
    for variant in variants:
        if not variant.commands:
            continue
        logger.info("Variant %s (%s): %d commands", variant.name, variant.kind, len(variant.commands))
        for command in variant.commands:
            logger.info("%s", shlex.join(command))
            if execute:
                for output_flag in ("--output", "--profile-output"):
                    if output_flag in command:
                        output_path = Path(command[command.index(output_flag) + 1])
                        output_path.parent.mkdir(parents=True, exist_ok=True)
                subprocess.run(command, check=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate or run simulation sensitivity variants.")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--script-dir", type=Path, default=Path("analysis/sim"))
    parser.add_argument("--results-dir", type=Path, default=Path("analysis/sim/results"))
    parser.add_argument("--output-dir", type=Path, default=Path("analysis/sim/results/sensitivity"))
    parser.add_argument("--manifest", type=Path, default=Path("analysis/sim/results/sensitivity/manifest.json"))
    parser.add_argument("--mass-label", default="mass13")
    parser.add_argument("--halos", type=Path, default=Path("analysis/sim/results/halos_mass13.csv"))
    parser.add_argument("--rockstar-halos", type=Path, default=Path("analysis/sim/results/hosts_minmass12p5.bin"))
    parser.add_argument("--kappa-map", type=Path, default=Path("analysis/sim/results/kappa_map_l0p1.float32"),
                        help="Unsmoothed base map; smoothing variants use its pre-smoothed "
                             "siblings (<stem>_s<level>.<ext> from derive_kappa_maps.py).")
    parser.add_argument("--normalized-only", action="store_true",
                        help="Skip fixed-separation pair stacks (halves stack time; the "
                             "sensitivity summary uses normalized stacks only).")
    parser.add_argument(
        "--groups",
        nargs="+",
        default=["smoothing", "los", "separation-bin"],
        choices=["smoothing", "los", "separation-bin", "mass", "all"],
    )
    parser.add_argument("--smooth-values", nargs="+", default=["none", "2arcmin", "4arcmin", "8arcmin"])
    parser.add_argument("--rpar-values", nargs="+", type=float, default=[10.0, 35.0])
    parser.add_argument("--baseline-rpar-max", type=float, default=22.0)
    parser.add_argument("--execute", action="store_true", help="Run commands. Default only writes manifest and prints commands.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    variants = build_variants(args)
    write_manifest(args, variants)
    run_commands(variants, args.execute)
    if not args.execute:
        logger.info("Dry run only. Re-run with --execute to launch these variants.")


if __name__ == "__main__":
    main()
