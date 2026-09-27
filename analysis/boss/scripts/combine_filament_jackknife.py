#!/usr/bin/env python
"""Filament map and bridge excess with jackknife errors.

Chains the joint accumulators into the derived-map sequence that
``plot_results.py`` builds, but does it once per leave-one-out region so the
filament signal gets an uncertainty:

    corrected pairs  = galaxy pairs - random pairs
    corrected single = single galaxies - single randoms
    control          = corrected single superposed at +/- sep/2
    filament         = corrected pairs - control

Every term that carries galaxy shot noise -- the galaxy pairs and both single
stacks -- is deleted region by region *together*, using the shared
tessellation.  That coherence is the point: the control is built from the
single stack, so a jackknife that dropped a region from the pair term but not
from the control would mis-estimate the very cancellation the filament
measurement depends on.

The random-pair map is held fixed rather than jackknifed, because per-region
random-pair accumulators would require re-running the 23-36 hour full-random
pair jobs.  That approximation is quantified rather than assumed: the measured
bridge excess of the frac100 random-pair maps is reported, and since randoms
carry no filament signal its departure from zero is a direct estimate of the
noise this term contributes.  Empirically it is ~0.5e-4, well below the
galaxy-side error.

The frac10 random-pair maps are *not* usable here.  Subsampling the random
catalog to 10% leaves 1% of the pairs -- pair counts go as the square of the
fraction -- so their noise is ~10x the frac100 maps.  Re-measured under the
fixed Y bands (2026-08-06), swapping frac10 for frac100 shifts the per-region
random bridge excess by -8.8 and +24.2e-4 at sep 5 and by -4.8 and -6.3e-4 at
sep 10.  Jointly that is a +5.2e-4 shift in the sep-10 filament against a 5.1e-4
error: a full 1 sigma of pure bookkeeping.  Sep 20 survives the swap (0.2e-4)
only because it has 16x more random pairs than sep 5, so do not generalize from
it -- any new separation needs its own frac100 run.

Wide r_perp groupings (``--groupings``)
----------------------------------------
A grouping such as ``g6_15`` holds pairs from r_perp 6 to 15, stacked by
``stack_pairs_jk.py --rperp-min 6 --rperp-max 15`` with fine-bin separation
summaries.  Two things change relative to the ``--separations`` path:

- **Random pairs** come from one ``find_and_stack_pairs.py --rperp-bin-edges``
  run per region, whose raw per-sub-bin sums are added over the grouping's
  sub-bins *and* over regions before dividing:
  ``sum_r sum_b sum_wk / sum_r sum_b sum_w``, then reflection-symmetrized once.
  That is the same per-pixel pooling the galaxy side uses.  The archived rule
  (per-region maps weighted by pair count) is still computed and its effect on
  the filament is reported as ``random_ns_counts_shift``.  The random map stays
  fixed across jackknife samples, as before.
- **The control** (``--control separation_averaged``, the default here) averages
  the two-halo template over the pairs' actual separations, weighted per pixel
  by each fine bin's coverage (``geometry.separation_averaged_template``).
  Every leave-one-out sample rebuilds it from that sample's separation mix.
  ``--control nominal`` uses the single ``bridge_center`` separation instead;
  the nominal-control filament is always reported as a diagnostic.

``bridge_center`` sets the bridge window (|x| <= 0.35 * center) and the nominal
control's separation.  The defaults below are provisional until the bin
choice is final; override with ``--bridge-centers g6_15=10.5,...``.  The value
used is written to every output row.

Inputs are read from two places: the single-galaxy stacks and everything the
``--separations`` path uses come from ``--results-dir`` (read only); the
grouping pair accumulators and random sub-bin sums come from ``--output-dir``,
which also receives every output.  Existing outputs are never replaced
without ``--overwrite``.

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/combine_filament_jackknife.py \\
        --separations 5,10,20

    W=analysis/boss/results/widebin_rpar5
    PYTHONPATH=lib python analysis/boss/scripts/combine_filament_jackknife.py \\
        --groupings g4_6,g6_15,g15_25 \\
        --random-subbins-label random_w3_25_rpar5_scw_frac100 \\
        --single-tag _scw --single-random-tag _scw_frac100 \\
        --results-dir analysis/boss/results --output-dir $W
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re

import numpy as np
import pandas as pd

from catalog import setup_logging
from combine_jackknife import leave_one_out_maps, load_accumulator, total_map
from geometry import (
    reflect_symmetrize_map,
    separation_averaged_template,
    separation_mix,
    symmetrize_map,
    two_halo_template,
)
from jackknife import jackknife_error
from stack_pairs_jk import CATALOG_CUTS_RE
from weighting_sensitivity import (
    SEPARATIONS,
    axis_for_map,
    bridge_excess,
    build_control_pair_map,
    combine,
    load_map,
    reconcile_shapes,
)

logger = logging.getLogger(__name__)

GALAXY_COUNTS = {"North": 579089.0, "South": 213205.0}

# Provisional bridge centres (h^-1 Mpc) per grouping: Zheng's R labels where
# they exist, the bin midpoint otherwise.  The final values are an open
# decision (plan section 9) to be fixed before reading the results.
BRIDGE_CENTERS = {
    "g4_6": 5.0, "g3_7": 5.0, "g6_10": 8.0, "g10_15": 12.5,
    "g6_15": 10.0, "g7_15": 11.0, "g15_25": 20.0,
}
GROUPING_RE = re.compile(r"^g(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)$")
CONTROL_TOKENS = {"separation_averaged": "avgctrl", "nominal": "nomctrl"}


def loo_and_total(acc: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return (total_map, leave_one_out_maps) as flattened arrays."""
    return (
        total_map(acc["sum_wk"], acc["sum_w"]),
        leave_one_out_maps(acc["sum_wk"], acc["sum_w"]),
    )


def parse_grouping(name: str) -> tuple[float, float]:
    match = GROUPING_RE.match(name)
    if not match:
        raise ValueError(f"Grouping {name!r} is not of the form g{{lo}}_{{hi}}.")
    lo, hi = float(match.group(1)), float(match.group(2))
    if not lo < hi:
        raise ValueError(f"Grouping {name!r} needs lo < hi.")
    return lo, hi


def load_random_subbins(output_dir: str, label: str, dataset: str, regions: list[str],
                        lo: float, hi: float) -> dict[str, dict[str, object]]:
    """Per-region random sums over the sub-bins that make up [lo, hi].

    Refuses a run whose sub-bin edges do not land exactly on lo and hi, or
    that did not complete every chunk (e.g. an unmerged shard).
    """
    out = {}
    for reg in regions:
        path = os.path.join(output_dir, f"kappa_pairs_{label}_{dataset}_{reg}.npz")
        with np.load(path, allow_pickle=False) as data:
            d = {key: data[key] for key in data.files}
        config = json.loads(str(d["config_json"]))
        if config.get("region") != reg:
            raise ValueError(f"{path} was run for region {config.get('region')!r}, not {reg!r}.")
        n_done, n_total = len(d["completed_chunks"]), int(d["n_chunks_total"])
        if n_done != n_total:
            raise ValueError(f"{path} holds {n_done} of {n_total} chunks; merge or finish it first.")
        edges = np.asarray(d["rperp_edges"], dtype=np.float64)
        i_lo = np.flatnonzero(np.isclose(edges, lo, rtol=0, atol=1e-9))
        i_hi = np.flatnonzero(np.isclose(edges, hi, rtol=0, atol=1e-9))
        if i_lo.size != 1 or i_hi.size != 1:
            raise ValueError(
                f"{path}: sub-bin edges {edges.tolist()} do not include both {lo:g} and {hi:g}.")
        bins = slice(int(i_lo[0]), int(i_hi[0]))
        out[reg] = {
            "path": path,
            "S": d["sum_wk"][bins].sum(axis=0),
            "W": d["sum_w"][bins].sum(axis=0),
            "n_pairs": int(d["n_pairs"][bins].sum()),
            "rpar": float(config["rpar"]),
            "sigma_crit_weight": bool(d["sigma_crit_weight"]),
            "closed_top": bool(i_hi[0] == len(edges) - 1),
            "subbin_edges": edges[bins.start:bins.stop + 1].tolist(),
        }
    for key in ("rpar", "sigma_crit_weight"):
        values = {reg: r[key] for reg, r in out.items()}
        if len(set(values.values())) != 1:
            raise ValueError(f"Random sub-bin runs differ in {key}: {values}")
    return out


def ratio_map(S: np.ndarray, W: np.ndarray) -> np.ndarray:
    """S / W with pixels no pair reached set to NaN."""
    out = np.full(S.shape, np.nan)
    good = W > 0
    out[good] = S[good] / W[good]
    return out


def pooled_random_map(rand: dict[str, dict[str, object]]) -> np.ndarray:
    """Primary random-pair map: sums pooled over sub-bins and regions."""
    S = np.sum([r["S"] for r in rand.values()], axis=0)
    W = np.sum([r["W"] for r in rand.values()], axis=0)
    return reflect_symmetrize_map(ratio_map(S, W))


def counts_random_map(rand: dict[str, dict[str, object]]) -> np.ndarray:
    """Archived rule, as a diagnostic: per-region maps weighted by pair count."""
    maps = {reg: reflect_symmetrize_map(ratio_map(r["S"], r["W"])) for reg, r in rand.items()}
    return combine(maps, {reg: float(r["n_pairs"]) for reg, r in rand.items()})


def galaxy_rpar(pair: dict) -> float:
    """r_par cut of the catalogs a pair accumulator was stacked from."""
    rpars = set()
    for path in pair["args"].get("pair_catalogs", []):
        match = CATALOG_CUTS_RE.search(os.path.basename(path))
        if match is None:
            raise ValueError(f"Cannot read the r_par cut from catalog name {path}.")
        rpars.add(float(match.group(1)))
    if len(rpars) != 1:
        raise ValueError(f"Pair accumulator was built from catalogs with r_par {sorted(rpars)}.")
    return rpars.pop()


def single_sigma_crit_weight(acc: dict, path: str) -> bool:
    """1/Sigma_crit^2 setting a single-stack accumulator was built with."""
    value = acc["args"].get("sigma_crit_weight")
    if value is None:
        raise ValueError(f"{path} does not record sigma_crit_weight; cannot check it "
                         "against the pair weighting.")
    return bool(value)


def load_grouping_pairs(path: str, lo: float, hi: float) -> dict:
    """Pair accumulator plus the r_perp cut and fine-bin fields C3 adds."""
    pair = load_accumulator(path)
    with np.load(path, allow_pickle=False) as data:
        missing = [k for k in ("rperp_min", "fine_w") if k not in data.files]
        if missing:
            raise ValueError(f"{path} has no {missing}; rerun stack_pairs_jk.py with --rperp-min/max.")
        for key in ("rperp_min", "rperp_max", "rperp_half_open", "sigma_crit_weight",
                    "fine_edges", "fine_w", "fine_wr", "fine_cov"):
            pair[key] = data[key]
    if not (np.isclose(pair["rperp_min"], lo) and np.isclose(pair["rperp_max"], hi)):
        raise ValueError(
            f"{path} was cut to r_perp {float(pair['rperp_min']):g}-{float(pair['rperp_max']):g}, "
            f"not {lo:g}-{hi:g}.")
    if pair["fine_edges"].size == 0:
        raise ValueError(f"{path} has no fine-bin summaries (--fine-bin-width 0?).")
    return pair


def prepare_single(single_flat: np.ndarray, g: int, symmetrize_single: bool) -> np.ndarray:
    single_map = single_flat.reshape(g, g)
    return symmetrize_map(single_map) if symmetrize_single else single_map


def grouping_filament(pair_flat: np.ndarray, single_flat: np.ndarray, rand_pair: np.ndarray,
                      mix: tuple, center: float, control: str, g: int,
                      symmetrize: bool = False, symmetrize_single: bool = True
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(filament, corrected pairs, control) for one realization of a grouping.

    ``mix`` is ``geometry.separation_mix``'s (separations, weights, coverage)
    for the same realization; it is ignored by the nominal control.
    """
    gp = rand_pair.shape[0]
    pair_map = pair_flat.reshape(gp, gp)
    if symmetrize:
        pair_map = reflect_symmetrize_map(pair_map)
    corrected_pairs = pair_map - rand_pair
    single_map = prepare_single(single_flat, g, symmetrize_single)
    axis = axis_for_map(corrected_pairs)
    gx, gy = np.meshgrid(axis, axis)
    if control == "nominal":
        ctrl = two_halo_template(single_map, gx, gy, center, validate=symmetrize_single)
    else:
        seps, weights, coverage = mix
        ctrl = separation_averaged_template(single_map, gx, gy, seps, weights,
                                            coverage=coverage, validate=symmetrize_single)
    if symmetrize:
        # Per-pixel coverage makes the averaged control slightly asymmetric;
        # reflect it exactly as the pair map was, so the filament map carries
        # no residual that the (symmetric) bridge statistic would hide.
        ctrl = reflect_symmetrize_map(ctrl)
    return corrected_pairs - ctrl, corrected_pairs, ctrl


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Filament bridge excess with jackknife errors.")
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument("--separations", default=None,
                        help="Archived single-separation products (default 5,10,20 "
                             "when --groupings is not given).")
    parser.add_argument("--groupings", default=None,
                        help="Wide r_perp groupings, e.g. g4_6,g6_15,g15_25.")
    parser.add_argument("--random-subbins-label", default=None,
                        help="Label of the find_and_stack_pairs.py --rperp-bin-edges run "
                             "(kappa_pairs_{label}_{dataset}_{region}.npz in --output-dir).")
    parser.add_argument("--galaxy-label-template", default="galaxy_{grouping}_rpar5_scw",
                        help="Pair accumulator label per grouping, read from "
                             "{output-dir}/jk/acc_pairs_{label}_{dataset}_{regions}.npz.")
    parser.add_argument("--control", choices=sorted(CONTROL_TOKENS), default=None,
                        help="Control for --groupings (default separation_averaged). "
                             "The --separations path always uses the nominal control.")
    parser.add_argument("--bridge-centers", default="",
                        help="Override bridge centres, e.g. g6_15=10.5,g15_25=19.")
    parser.add_argument("--tag", default="",
                        help="Extra token in --groupings output names, e.g. pilot.")
    parser.add_argument("--single-tag", default="_scw")
    parser.add_argument("--single-random-tag", default=None,
                        help="Defaults to --single-tag; use e.g. _scw_frac100.")
    parser.add_argument("--results-dir", default="analysis/boss/results",
                        help="Read-only inputs: single stacks and archived pair products.")
    parser.add_argument("--output-dir", default=None,
                        help="Grouping inputs and all outputs (default: --results-dir).")
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace existing outputs.")
    parser.add_argument("--no-symmetrize-single", action="store_true",
                        help="Build the control from the unsymmetrized single map. "
                             "Off by default: the archived pipeline symmetrizes, and "
                             "not doing so leaks azimuthal noise into the filament.")
    parser.add_argument("--symmetrize", action="store_true",
                        help="Reflection-symmetrize each pair realization before "
                             "deriving statistics (matches the archived maps).")
    args = parser.parse_args(argv)
    args.regions = [r.strip() for r in args.regions.split(",")]
    if args.groupings is not None and args.separations is not None:
        parser.error("--groupings and --separations are separate runs; give one.")
    if args.groupings is None:
        args.separations = [s.strip() for s in (args.separations or "5,10,20").split(",")]
        args.groupings = []
        if args.control == "separation_averaged" or args.random_subbins_label:
            parser.error("--control separation_averaged and --random-subbins-label need --groupings.")
    else:
        args.separations = []
        args.groupings = [x.strip() for x in args.groupings.split(",") if x.strip()]
        if not args.random_subbins_label:
            parser.error("--groupings needs --random-subbins-label.")
        if args.control is None:
            args.control = "separation_averaged"
    centers = dict(BRIDGE_CENTERS)
    for item in filter(None, (x.strip() for x in args.bridge_centers.split(","))):
        name, _, value = item.partition("=")
        centers[name.strip()] = float(value)
    args.centers = {}
    for name in args.groupings:
        try:
            parse_grouping(name)
        except ValueError as exc:
            parser.error(str(exc))
        if name not in centers:
            parser.error(f"No bridge centre for {name}; pass --bridge-centers {name}=<value>.")
        args.centers[name] = centers[name]
    if args.single_random_tag is None:
        args.single_random_tag = args.single_tag
    if args.output_dir is None:
        args.output_dir = args.results_dir
    return args


def output_paths(args, region_label: str) -> dict[str, dict[str, str]]:
    """Every file this run writes, keyed by separation/grouping."""
    O = args.output_dir
    paths = {}
    for sep in args.separations:
        label = SEPARATIONS[sep]["galaxy_label"]
        paths[sep] = {
            "map": os.path.join(O, f"kappa_filament_{label}_{args.dataset}_{region_label}_joint.csv"),
            "err": os.path.join(O, f"error_filament_{label}_{args.dataset}_{region_label}_joint.csv"),
        }
    if args.separations:
        paths["summary"] = {"csv": os.path.join(
            O, f"filament_jackknife_{args.dataset}_{region_label}.csv")}
    ctrl = CONTROL_TOKENS.get(args.control, "")
    tag = f"_{args.tag}" if args.tag else ""
    for name in args.groupings:
        label = args.galaxy_label_template.format(grouping=name)
        stem = f"{label}{tag}_{ctrl}_{args.dataset}_{region_label}_joint.csv"
        paths[name] = {
            "map": os.path.join(O, f"kappa_filament_{stem}"),
            "err": os.path.join(O, f"error_filament_{stem}"),
        }
    if args.groupings:
        base = os.path.join(O, f"filament_jackknife_wide{tag}_{ctrl}_{args.dataset}_{region_label}")
        paths["summary"] = {"csv": base + ".csv", "meta": base + ".meta.json"}
    return paths


def run_grouping(args, name: str, gal: dict, single_total: np.ndarray, single_loo: np.ndarray,
                 region_label: str, paths: dict) -> tuple[dict, dict]:
    lo, hi = parse_grouping(name)
    center = args.centers[name]
    label = args.galaxy_label_template.format(grouping=name)
    pair_path = os.path.join(args.output_dir, "jk",
                             f"acc_pairs_{label}_{args.dataset}_{region_label}.npz")
    pair = load_grouping_pairs(pair_path, lo, hi)
    if pair["jk_digest"] != gal["jk_digest"]:
        raise ValueError(
            f"{pair_path} uses tessellation {pair['jk_digest']} but the single stacks use "
            f"{gal['jk_digest']}; the jackknife would delete different sky from the two terms.")

    rand = load_random_subbins(args.output_dir, args.random_subbins_label, args.dataset,
                               args.regions, lo, hi)
    first = next(iter(rand.values()))
    # All four inputs must share one weighting: a pair stack weighted by
    # 1/Sigma_crit^2 against an unweighted single stack is exactly the
    # inconsistency the wide-bin run exists to remove.
    weighting = {
        "galaxy pairs": bool(pair["sigma_crit_weight"]),
        "random pairs": first["sigma_crit_weight"],
        "single galaxies": args.single_weighting["galaxy"],
        "single randoms": args.single_weighting["random"],
    }
    if len(set(weighting.values())) != 1:
        raise ValueError(f"{name}: inputs disagree on sigma_crit_weight: {weighting}")
    rpar = galaxy_rpar(pair)
    if not np.isclose(rpar, first["rpar"]):
        raise ValueError(f"{name}: galaxy pairs use r_par <= {rpar:g}, random pairs {first['rpar']:g}.")
    # Random sub-bins are [e_i, e_i+1) with the top one closed; the galaxy cut
    # must treat the upper edge the same way or one side keeps pairs at r = hi.
    if bool(pair["rperp_half_open"]) == first["closed_top"]:
        raise ValueError(
            f"{name}: galaxy cut is {'[lo, hi)' if pair['rperp_half_open'] else '[lo, hi]'} but "
            f"the random sub-bins end {'closed' if first['closed_top'] else 'open'} at {hi:g}.")

    gp = int(pair["grid_size"])
    rand_pair = pooled_random_map(rand)
    rand_counts = counts_random_map(rand)
    if rand_pair.shape != (gp, gp):
        raise ValueError(f"Random map is {rand_pair.shape}, galaxy pair grid is {gp}x{gp}.")
    n_bad = int(np.isnan(rand_pair).sum())
    if n_bad:
        logger.warning("%s: %d random-map pixels have no pairs (NaN).", name, n_bad)

    g = gal["grid_size"]
    sym_single = not args.no_symmetrize_single
    pair_total, pair_loo = loo_and_total(pair)
    fine = (pair["fine_w"], pair["fine_wr"], pair["fine_cov"])

    def filament(pair_flat, single_flat, mix, control=args.control, rp=rand_pair):
        return grouping_filament(pair_flat, single_flat, rp, mix, center, control, g,
                                 symmetrize=args.symmetrize, symmetrize_single=sym_single)

    mix_total = separation_mix(*fine)
    fil_total, cpairs_total, control_total = filament(pair_total, single_total, mix_total)
    be_total = bridge_excess(fil_total, center)

    # Diagnostics on the full sample: the other control, histogram-only
    # weights instead of per-pixel coverage (plan test T8), and the archived
    # N/S rule for the random map.
    other = "nominal" if args.control == "separation_averaged" else "separation_averaged"
    fil_other, _, _ = filament(pair_total, single_total, mix_total, control=other)
    fil_hist, _, _ = filament(pair_total, single_total, (mix_total[0], mix_total[1], None),
                              control="separation_averaged")
    fil_avg = fil_total if args.control == "separation_averaged" else fil_other
    fil_counts, _, _ = filament(pair_total, single_total, mix_total, rp=rand_counts)

    n_regions = single_loo.shape[0]
    be_loo = np.empty(n_regions)
    fil_loo = np.empty((n_regions, fil_total.size))
    for k in range(n_regions):
        fil_k, _, _ = filament(pair_loo[k], single_loo[k], separation_mix(*fine, exclude_region=k))
        fil_loo[k] = fil_k.ravel()
        be_loo[k] = bridge_excess(fil_k, center)
    _, be_err = jackknife_error(be_loo)
    _, fil_err = jackknife_error(fil_loo)

    tot_w = pair["fine_w"].sum()
    row = {
        "grouping": name,
        "rperp_min": lo,
        "rperp_max": hi,
        "rperp_half_open": bool(pair["rperp_half_open"]),
        "rpar_max": rpar,
        "bridge_center": center,
        "control_mode": args.control,
        "n_pairs": int(pair["n_pairs"].sum()),
        "mean_rperp": float(pair["fine_wr"].sum() / tot_w) if tot_w > 0 else np.nan,
        "corrected_pairs_bridge_excess": bridge_excess(cpairs_total, center),
        "control_bridge_excess": bridge_excess(control_total, center),
        "filament_bridge_excess": be_total,
        "filament_bridge_err": be_err,
        "significance": be_total / be_err if be_err > 0 else np.nan,
        f"filament_bridge_excess_{other}_control": bridge_excess(fil_other, center),
        "control_coverage_shift": bridge_excess(fil_hist, center) - bridge_excess(fil_avg, center),
        "random_pair_noise_floor": bridge_excess(rand_pair, center),
        "random_ns_counts_shift": bridge_excess(fil_counts, center) - be_total,
        "n_random_pairs": int(sum(r["n_pairs"] for r in rand.values())),
        "galaxy_label": label,
        "random_label": args.random_subbins_label,
    }
    pd.DataFrame(fil_total).to_csv(paths[name]["map"])
    pd.DataFrame(fil_err.reshape(fil_total.shape)).to_csv(paths[name]["err"])
    logger.info("%s (r_perp %g-%g, mean %.2f): filament %.3fe-4 +/- %.3fe-4 (%.1f sigma); "
                "%s control %.3fe-4; coverage shift %.4fe-4; N/S count-rule shift %.4fe-4",
                name, lo, hi, row["mean_rperp"], be_total * 1e4, be_err * 1e4, row["significance"],
                other, row[f"filament_bridge_excess_{other}_control"] * 1e4,
                row["control_coverage_shift"] * 1e4, row["random_ns_counts_shift"] * 1e4)
    inputs = {"galaxy_pairs": pair_path, "random_pairs": [r["path"] for r in rand.values()],
              "random_subbin_edges": first["subbin_edges"]}
    return row, inputs


def main(argv=None):
    args = parse_args(argv)
    setup_logging()
    R = args.results_dir
    O = args.output_dir
    jk_dir = os.path.join(R, "jk")
    region_label = "_".join(args.regions)

    paths = output_paths(args, region_label)
    existing = [p for group in paths.values() for p in group.values() if os.path.exists(p)]
    if existing and not args.overwrite:
        raise FileExistsError("Outputs exist (use --overwrite): " + ", ".join(existing))
    os.makedirs(O, exist_ok=True)

    # ---- single-galaxy side: corrected single, per leave-one-out region ----
    gal_single_path = os.path.join(
        jk_dir, f"acc_single_galaxy{args.single_tag}_{args.dataset}_{region_label}.npz")
    rnd_single_path = os.path.join(
        jk_dir, f"acc_single_random{args.single_random_tag}_{args.dataset}_{region_label}.npz")
    gal = load_accumulator(gal_single_path)
    rnd = load_accumulator(rnd_single_path)
    if args.groupings:
        args.single_weighting = {
            "galaxy": single_sigma_crit_weight(gal, gal_single_path),
            "random": single_sigma_crit_weight(rnd, rnd_single_path),
        }
    if gal["jk_digest"] != rnd["jk_digest"]:
        raise ValueError("Single galaxy/random accumulators use different tessellations.")

    gal_total, gal_loo = loo_and_total(gal)
    rnd_total, rnd_loo = loo_and_total(rnd)
    single_total = gal_total - rnd_total
    single_loo = gal_loo - rnd_loo
    n_regions = single_loo.shape[0]
    g = gal["grid_size"]

    if args.groupings:
        rows, inputs = [], {}
        for name in args.groupings:
            row, inputs[name] = run_grouping(args, name, gal, single_total, single_loo,
                                             region_label, paths)
            rows.append(row)
        out = pd.DataFrame(rows)
        out.to_csv(paths["summary"]["csv"], index=False)
        from find_and_stack_pairs import run_provenance
        meta = {
            "args": {k: v for k, v in vars(args).items()},
            "single_accumulators": [gal_single_path, rnd_single_path],
            "inputs": inputs,
            "outputs": paths,
            "provenance": run_provenance(),
        }
        with open(paths["summary"]["meta"], "w", encoding="ascii") as handle:
            json.dump(meta, handle, indent=2, sort_keys=True, default=str)
            handle.write("\n")
        logger.info("Saved -> %s", paths["summary"]["csv"])
        logger.info("\nWide-bin filament measurement (kappa x 1e4), control=%s:", args.control)
        logger.info("  grouping  mean_r   center   filament   error   sigma   rand-floor")
        for _, r in out.iterrows():
            logger.info("  %-8s  %6.2f  %6.1f  %9.2f  %6.2f  %5.1f  %9.2f",
                        r.grouping, r.mean_rperp, r.bridge_center, r.filament_bridge_excess * 1e4,
                        r.filament_bridge_err * 1e4, r.significance,
                        r.random_pair_noise_floor * 1e4)
        return

    rows = []
    for sep in args.separations:
        cfg = SEPARATIONS[sep]
        center = cfg["center"]

        pair_acc_path = os.path.join(
            jk_dir, f"acc_pairs_galaxy_{cfg['galaxy_label']}_{args.dataset}_{region_label}.npz")
        if not os.path.exists(pair_acc_path):
            logger.warning("Missing %s; skipping sep=%s", pair_acc_path, sep)
            continue
        pair = load_accumulator(pair_acc_path)
        if pair["jk_digest"] != gal["jk_digest"]:
            raise ValueError(
                f"Pair accumulator for sep={sep} uses tessellation {pair['jk_digest']} "
                f"but the single stacks use {gal['jk_digest']}; the filament jackknife "
                "would delete different sky from the two terms.")
        gp = pair["grid_size"]
        pair_total, pair_loo = loo_and_total(pair)

        # Random pairs: fixed, count-weighted across regions (no accumulators).
        rand_pair_maps, rand_pair_counts = {}, {}
        for reg in args.regions:
            rand_pair_maps[reg] = load_map(os.path.join(
                R, f"kappa_pairs_random_{cfg['random_label']}_{args.dataset}_{reg}.csv"))
            with open(os.path.join(
                R, f"kappa_pairs_random_{cfg['random_label']}_{args.dataset}_{reg}.meta.json"),
                encoding="utf-8") as handle:
                rand_pair_counts[reg] = float(json.load(handle)["total_pairs"])
        rand_pair = combine(rand_pair_maps, rand_pair_counts)
        rand_pair_floor = bridge_excess(rand_pair, center)

        def filament_from(single_flat: np.ndarray, pair_flat: np.ndarray) -> np.ndarray:
            pair_map = pair_flat.reshape(gp, gp)
            if args.symmetrize:
                pair_map = reflect_symmetrize_map(pair_map)
            cp, rp = reconcile_shapes(pair_map, rand_pair)
            corrected_pairs = cp - rp
            # Radially symmetrize the single map before superposing it: the
            # single-galaxy signal is radially symmetric by construction, so this
            # suppresses azimuthal noise without biasing the control. Skipping it
            # leaks that noise straight into the filament residual.
            single_map = single_flat.reshape(g, g)
            if not args.no_symmetrize_single:
                single_map = symmetrize_map(single_map)
            # Build the control natively on the pair grid via the radial profile,
            # rather than shifting a pixel array and trimming to a common shape.
            # The latter left a half-pixel misalignment (101-pixel pair grid on
            # integers vs 100-pixel single grid on half-integers) that produced
            # false residuals at the halo peaks. See geometry.two_halo_template.
            pair_axis = axis_for_map(corrected_pairs)
            gx, gy = np.meshgrid(pair_axis, pair_axis)
            # The symmetry check exists to catch archived mis-centered stacks.
            # Under --no-symmetrize-single the caller is handing in a raw map on
            # purpose, so the check would reject the very thing being tested.
            control = two_halo_template(single_map, gx, gy, center,
                                        validate=not args.no_symmetrize_single)
            return corrected_pairs - control, corrected_pairs, control

        fil_total, cpairs_total, control_total = filament_from(single_total, pair_total)
        be_total = bridge_excess(fil_total, center)

        be_loo = np.empty(n_regions)
        fil_loo = np.empty((n_regions, fil_total.size))
        for k in range(n_regions):
            fil_k, _, _ = filament_from(single_loo[k], pair_loo[k])
            fil_loo[k] = fil_k.ravel()
            be_loo[k] = bridge_excess(fil_k, center)

        _, be_err = jackknife_error(be_loo)
        _, fil_err = jackknife_error(fil_loo)

        rows.append({
            "separation_hmpc": center,
            "n_pairs": int(pair["n_pairs"].sum()),
            "corrected_pairs_bridge_excess": bridge_excess(cpairs_total, center),
            "control_bridge_excess": bridge_excess(control_total, center),
            "filament_bridge_excess": be_total,
            "filament_bridge_err": be_err,
            "significance": be_total / be_err if be_err > 0 else np.nan,
            "random_pair_noise_floor": rand_pair_floor,
        })

        pd.DataFrame(fil_total).to_csv(paths[sep]["map"])
        pd.DataFrame(fil_err.reshape(fil_total.shape)).to_csv(paths[sep]["err"])
        logger.info("sep=%s: filament bridge excess %.3fe-4 +/- %.3fe-4 (%.1f sigma)",
                    sep, be_total * 1e4, be_err * 1e4,
                    be_total / be_err if be_err > 0 else np.nan)

    if not rows:
        logger.warning("No separations produced results.")
        return

    out = pd.DataFrame(rows)
    path = paths["summary"]["csv"]
    out.to_csv(path, index=False)
    logger.info("Saved -> %s", path)

    logger.info("\nJoint filament measurement (kappa x 1e4):")
    logger.info("  sep   corrected_pairs   control   filament   error   sigma   rand-floor")
    for _, r in out.iterrows():
        logger.info("  %4.0f  %14.2f  %8.2f  %9.2f  %6.2f  %5.1f  %9.2f",
                    r.separation_hmpc, r.corrected_pairs_bridge_excess * 1e4,
                    r.control_bridge_excess * 1e4, r.filament_bridge_excess * 1e4,
                    r.filament_bridge_err * 1e4, r.significance,
                    r.random_pair_noise_floor * 1e4)


if __name__ == "__main__":
    main()
