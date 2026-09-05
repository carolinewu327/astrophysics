#!/usr/bin/env python
"""Pair stack with spatial-jackknife errors on the residual bridge statistics.

Splits the pairs into K = blocks_per_side^2 spatial blocks by pair-center
position, stacks each block once (total cost = one full stack), and combines
leave-one-out means to estimate jackknife uncertainties of the bridge
statistics of (pair stack - matched two-halo template).

The template is built from the supplied single stack (same map/smoothing as
the pair stack: matched-template rule) and held fixed across jackknife
regions -- its own sampling noise is negligible when the single stack uses
the full galaxy sample.

Supports nested LOS cuts via --rpar-max (filters r_parallel_rsd), so one
deep pair catalog serves every candidate cut. Maps are expected pre-smoothed;
stacking is always fixed-separation (normalized stacking is a Phase 4
sensitivity axis).
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from geometry import (BRIDGE_HALF_X_FRAC, CENTRAL_HALF_Y_HMPC, OFF_HI_Y_HMPC,
                      OFF_LO_Y_HMPC, bridge_excess)
from sim_utils import (
    make_two_halo_template,
    open_kappa_memmap,
    reflect_symmetrize_map,
    save_map_csv,
    setup_logging,
)
from stack_pairs_sim import stack_pairs
from summarize_sim_sensitivity import load_map, map_stats


logger = logging.getLogger(__name__)


def jackknife_stats(theta: np.ndarray) -> tuple[float, float]:
    """Jackknife mean and error from leave-one-out estimates."""
    k = len(theta)
    mean = float(np.mean(theta))
    var = (k - 1) / k * float(np.sum((theta - mean) ** 2))
    return mean, float(np.sqrt(var))


def run(args: argparse.Namespace) -> dict:
    pairs = pd.read_csv(args.pairs)
    n_total = len(pairs)
    rpar_col = "r_parallel_real" if args.rpar_space == "real" else "r_parallel_rsd"
    if args.rpar_min is not None or args.rpar_max is not None:
        if rpar_col not in pairs.columns:
            raise ValueError(
                f"Pair catalog lacks {rpar_col!r}; it predates the real/rsd split. "
                "Regenerate it with find_pairs_sim.py.")
        values = pairs[rpar_col].to_numpy()
        keep = np.ones(len(pairs), dtype=bool)
        if args.rpar_min is not None:
            keep &= values >= args.rpar_min
        if args.rpar_max is not None:
            # Half-open [lo, hi) when requested, so adjacent LOS bands cannot
            # both claim a pair sitting exactly on the boundary.
            if getattr(args, "rpar_half_open", False):
                keep &= values < args.rpar_max
            else:
                keep &= values <= args.rpar_max
        pairs = pairs[keep]
    # The r_perp bin is what makes a mock sample comparable to a BOSS sample,
    # and nothing used to check it: the 19-21 mock stack was compared against
    # the 18-22 BOSS bin for weeks because only r_par was asserted.  Declare
    # the bin, filter to it, and record both the declared bounds and what the
    # catalogue actually contains, so a later reader can tell them apart.
    if "r_perp" not in pairs.columns:
        raise ValueError(
            "Pair catalog lacks 'r_perp'; regenerate it with find_pairs_sim.py.")
    rperp_vals = pairs["r_perp"].to_numpy()
    inside = (rperp_vals >= args.rperp_min) & (rperp_vals <= args.rperp_max)
    if not inside.any():
        raise ValueError(
            f"No pairs with {args.rperp_min} <= r_perp <= {args.rperp_max}; the "
            f"catalogue spans {rperp_vals.min():g}-{rperp_vals.max():g}. Wrong "
            "catalogue for this bin.")
    if not inside.all():
        logger.warning("Dropping %d of %d pairs outside r_perp [%g, %g]",
                       int((~inside).sum()), len(pairs), args.rperp_min,
                       args.rperp_max)
        pairs = pairs[inside]
        rperp_vals = rperp_vals[inside]
    rperp_obs = (float(rperp_vals.min()), float(rperp_vals.max()))

    # A catalogue built for a narrower bin cannot be widened by asking for a
    # wider one -- it simply has no pairs out there.  Catch that rather than
    # silently returning the narrow sample under a wide label.
    tol = 0.05 * (args.rperp_max - args.rperp_min)
    if (rperp_obs[0] - args.rperp_min > tol) or (args.rperp_max - rperp_obs[1] > tol):
        raise ValueError(
            f"Declared r_perp bin [{args.rperp_min:g}, {args.rperp_max:g}] but the "
            f"catalogue only spans [{rperp_obs[0]:g}, {rperp_obs[1]:g}]. This is the "
            "19-21 vs 18-22 failure; regenerate the pair catalogue for the "
            "declared bin with find_pairs_sim.py.")

    if args.max_pairs is not None:
        # .iloc[:n] would take the FIRST n rows, which are ordered by the
        # pair-finder and are therefore spatially correlated -- not a random
        # subsample.  Sample explicitly, with a recorded seed.
        if len(pairs) > args.max_pairs:
            pairs = pairs.sample(n=args.max_pairs,
                                 random_state=getattr(args, "seed", 0))
    logger.info("Pairs: %d of %d after cuts (%s in [%s, %s])", len(pairs), n_total,
                rpar_col, args.rpar_min, args.rpar_max)

    kappa_map, info = open_kappa_memmap(args.kappa_map)
    box = info.box_size_hmpc
    cell = box / args.blocks_per_side
    bx = np.clip((pairs["pair_center_x"].to_numpy() // cell).astype(int), 0, args.blocks_per_side - 1)
    by = np.clip((pairs["pair_center_y"].to_numpy() // cell).astype(int), 0, args.blocks_per_side - 1)
    block = bx * args.blocks_per_side + by

    grid = args.grid_size
    sums = []
    counts = []
    for b in range(args.blocks_per_side ** 2):
        sub = pairs[block == b]
        if len(sub) == 0:
            continue
        mean_map = stack_pairs(
            pairs=sub,
            kappa_map=kappa_map,
            pixel_size_hmpc=info.pixel_size_hmpc,
            map_box_size_hmpc=box,
            stack_box_size_hmpc=args.box_size,
            grid_size=grid,
            normalize_separation=False,
            normalized_half_size=2.5,
        ).astype(np.float64)
        sums.append(mean_map * len(sub))
        counts.append(len(sub))
    sums = np.array(sums)
    counts = np.array(counts, dtype=np.float64)
    k = len(counts)
    logger.info("Stacked %d non-empty blocks (%d pairs)", k, int(counts.sum()))

    # stack_single_sim.py already writes a radially symmetrized map.  Do not
    # re-symmetrize it here: doing so used to hide archived products centered
    # on N//2 instead of the physical (N-1)/2 origin.  The template builder
    # validates the saved product and rejects the old convention.
    single = load_map(Path(args.single))
    if single.shape[0] != grid:
        logger.warning(
            "single %s is %d px against a %d px pair grid; its band is sampled "
            "2-4%% low because its rows sit at half-integer Y. Prefer the "
            "matched-grid stack (…_centered_g101.csv).",
            Path(args.single).name, single.shape[0], grid)
    axis = np.linspace(-0.5 * args.box_size, 0.5 * args.box_size, grid)
    x_grid, y_grid = np.meshgrid(axis, axis)
    template = make_two_halo_template(single, x_grid, y_grid, args.rperp_center)
    norm_axis = axis / args.rperp_center

    def fixed_band(arr):
        """Bridge statistic in the BOSS estimator's fixed physical bands.

        ``map_stats`` below is scored on ``norm_axis`` (= axis / r_perp), so its
        Y bands scale with the separation and each separation is measured with a
        different ruler.  That is a different statistic from the one BOSS
        reports, which is why the two were never comparable.  This one uses
        geometry.bridge_excess: fixed |Y| bands, only the X window follows the
        separation.

        Symmetrized first so the number describes exactly the map that
        --stack-output writes.  It makes no numerical difference -- the bands
        and the window are both symmetric, so the fold is linear on this
        statistic (measured: 3e-8 relative) -- but it removes the question.
        """
        return bridge_excess(
            reflect_symmetrize_map(arr.astype(np.float32)).astype(np.float64),
            axis, args.rperp_center)

    full_map = sums.sum(axis=0) / counts.sum()
    full_stats = {}
    full_stats.update(map_stats(full_map, norm_axis, "raw"))
    full_stats.update(map_stats(full_map - template, norm_axis, "residual"))
    full_stats["raw_bridge_excess_fixedband_kappa"] = fixed_band(full_map)
    full_stats["residual_bridge_excess_fixedband_kappa"] = fixed_band(
        full_map - template)

    loo_keys = None
    loo_values = []
    total_sum = sums.sum(axis=0)
    total_count = counts.sum()
    for b in range(k):
        loo_map = (total_sum - sums[b]) / (total_count - counts[b])
        stats = {}
        stats.update(map_stats(loo_map, norm_axis, "raw"))
        stats.update(map_stats(loo_map - template, norm_axis, "residual"))
        # Template held fixed across realizations, as above: it comes from the
        # single stack, which these blocks do not resample.
        stats["raw_bridge_excess_fixedband_kappa"] = fixed_band(loo_map)
        stats["residual_bridge_excess_fixedband_kappa"] = fixed_band(
            loo_map - template)
        if loo_keys is None:
            loo_keys = list(stats)
        loo_values.append([stats[key] for key in loo_keys])
    loo_values = np.array(loo_values)

    result = {
        "pairs": str(args.pairs),
        "kappa_map": str(args.kappa_map),
        "single": str(args.single),
        "rperp_center_hmpc": args.rperp_center,
        "rpar_max_hmpc": args.rpar_max,
        "rpar_min_hmpc": args.rpar_min,
        "rpar_space": args.rpar_space,
        "rperp_min_hmpc": args.rperp_min,
        "rperp_max_hmpc": args.rperp_max,
        "rperp_observed_min_hmpc": rperp_obs[0],
        "rperp_observed_max_hmpc": rperp_obs[1],
        # Two scoring conventions live in "stats" and they are NOT
        # interchangeable.  Keys ending in _fixedband_ use the fixed physical
        # bands below and match the BOSS estimator; every other bridge/side key
        # comes from map_stats on axis / r_perp and so scales its Y bands with
        # the separation.  The scaled ones are kept as a legacy sensitivity
        # axis; production readers should require the _fixedband_ keys.
        "band_definition": {
            "fixedband_keys": {
                "central_half_y_hmpc": CENTRAL_HALF_Y_HMPC,
                "off_lo_hmpc": OFF_LO_Y_HMPC,
                "off_hi_hmpc": OFF_HI_Y_HMPC,
                "bridge_half_x_frac": BRIDGE_HALF_X_FRAC,
                "source": "lib/geometry.py",
            },
            "legacy_scaled_keys": {
                "central_half_y_frac": 0.15,
                "off_lo_frac": 0.45,
                "off_hi_frac": 0.85,
                "bridge_half_x_frac": 0.35,
                "source": "summarize_sim_sensitivity.map_stats on axis/r_perp",
                "note": "not comparable with the BOSS estimator",
            },
        },
        "n_pairs": int(counts.sum()),
        "n_blocks": k,
        "blocks_per_side": args.blocks_per_side,
        "stats": {},
    }
    for j, key in enumerate(loo_keys):
        _, err = jackknife_stats(loo_values[:, j])
        result["stats"][key] = {"value": full_stats[key], "jackknife_error": err}

    if args.stack_output:
        save_map_csv(args.stack_output, reflect_symmetrize_map(full_map.astype(np.float32)))

    # Persist the per-block accumulators, not just the collapsed map.  Any
    # control, band definition or extrapolation can then be recomputed on the
    # full map *and* on every leave-one-out map without re-stacking -- which is
    # what an estimator other than the built-in template needs in order to have
    # a jackknife error at all.  25 blocks x 101 x 101 float64 is ~2 MB.
    if getattr(args, "blocks_output", None):
        np.savez_compressed(
            args.blocks_output,
            sums=sums, counts=counts,
            grid_size=grid, box_size=args.box_size,
            rperp_center=args.rperp_center,
            rpar_min=np.nan if args.rpar_min is None else args.rpar_min,
            rpar_max=np.nan if args.rpar_max is None else args.rpar_max,
            rpar_space=str(args.rpar_space),
            rperp_min=args.rperp_min,
            rperp_max=args.rperp_max,
            n_pairs=int(counts.sum()),
        )
        logger.info("Per-block accumulators -> %s (%d blocks)", args.blocks_output, k)
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pair stack with spatial-jackknife bridge errors.")
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--kappa-map", required=True, help="Pre-smoothed map; no smoothing applied here.")
    parser.add_argument("--single", required=True,
                        help="Matched single stack CSV (same map/smoothing as --kappa-map).")
    parser.add_argument("--rperp-center", type=float, required=True)
    parser.add_argument("--rperp-min", type=float, required=True,
                        help="Lower edge of the transverse bin. Required: the "
                             "bin is what makes a mock sample comparable to a "
                             "BOSS sample, and an unstated one is how a 19-21 "
                             "stack ended up plotted against BOSS's 18-22.")
    parser.add_argument("--rperp-max", type=float, required=True,
                        help="Upper edge of the transverse bin.")
    parser.add_argument("--rpar-max", type=float, default=None,
                        help="Optional upper LOS cut, applied in --rpar-space.")
    parser.add_argument("--rpar-min", type=float, default=None,
                        help="Optional lower LOS cut. Combined with --rpar-max and "
                             "--rpar-space real this selects NON-PHYSICAL pairs: close "
                             "in projection but far apart along the line of sight, so "
                             "they are not physically associated. Their residual bridge "
                             "excess should be zero if the two-halo template is an "
                             "unbiased control.")
    parser.add_argument("--rpar-space", choices=["redshift", "real"], default="redshift",
                        help="Which LOS separation the cuts apply to (default: redshift, "
                             "matching the production/BOSS-like selection).")
    parser.add_argument("--blocks-per-side", type=int, default=5)
    parser.add_argument("--grid-size", type=int, default=101)
    parser.add_argument("--box-size", type=float, default=100.0)
    parser.add_argument("--max-pairs", type=int, default=None, help="Pilot cap.")
    parser.add_argument("--output", required=True, help="Output JSON with stats +/- jackknife errors.")
    parser.add_argument("--stack-output", default=None, help="Optional full stacked-map CSV.")
    parser.add_argument("--blocks-output", default=None,
                        help="Optional .npz of per-block sum maps and counts. Required "
                             "to give any non-built-in estimator a jackknife error.")
    parser.add_argument("--rpar-half-open", action="store_true",
                        help="Treat the LOS cut as [min, max) instead of [min, max].")
    parser.add_argument("--seed", type=int, default=0,
                        help="Seed for --max-pairs subsampling.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    result = run(args)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
    for key, label in (("residual_bridge_excess_fixedband_kappa", "fixed band"),
                       ("residual_bridge_excess_kappa", "scaled band (legacy)")):
        r = result["stats"][key]
        logger.info(
            "residual bridge excess, %-21s = %.4e +/- %.4e (S/N %.2f)",
            label, r["value"], r["jackknife_error"],
            r["value"] / r["jackknife_error"] if r["jackknife_error"] > 0
            else float("nan"))
    logger.info("-> %s", out)


if __name__ == "__main__":
    main()
