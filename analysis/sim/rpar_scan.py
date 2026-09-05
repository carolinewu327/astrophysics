#!/usr/bin/env python
"""How does the filament signal change with the line-of-sight cut?

Zheng's question is whether widening r_par from 5 to 10 h^-1 Mpc for the
r_perp = 5 pairs improves S/N, so that all three separations share one LOS cut.
Widening admits more pairs (noise falls) but the added pairs are less physically
associated (signal dilutes), so the answer is an optimum, not a monotone.

**The simulation cannot answer the S/N question directly.**  It has no CMB
reconstruction noise -- its errors are sample variance over a 2500 h^-1 Mpc box,
while BOSS's error is dominated by Planck noise that is external to the pairs.
Optimising the simulation's own S/N would optimise the wrong quantity.

What the simulation gives is the *signal* curve S(C): how much bridge excess
survives at each cut C.  Combined with the pair count N(C), the figure of merit
for BOSS is

    predicted S/N (C)  proportional to  S(C) * sqrt(N(C))

reported here as a ratio to the current production cut so it is unit-free.

That sqrt(N) is an upper bound, and a generous one.  Deepening the LOS window
adds pairs drawn from the same galaxies at nearly the same sky positions, so
they share Planck noise with the pairs already counted -- far less independent
than pairs added by widening the survey.  BOSS's own scaling already shows this:
across separations the error falls *faster* than 1/sqrt(N) because larger
separations also spread pairs over more sky, a benefit widening r_par does not
get.  Treat the ranking as reliable and the magnitude as optimistic; the real
error bar comes from re-stacking BOSS galaxy pairs at the candidate cut.

Cuts are free: the pair catalogs are built at r_par,rsd <= 50, so every cut is a
subset filter rather than a new pair search.

Scoring uses the shared fixed-Y bands from lib/geometry.py -- the same statistic
as the BOSS chain, which is the point.  Run with lib/ on the path:

    PYTHONPATH=lib:analysis/sim python analysis/sim/rpar_scan.py \\
        --catalog 5=analysis/sim/results/hod_pairs/pairs_hodnmatch50_rperp5_rsd50.csv \\
        --catalog 10=analysis/sim/results/hod_pairs/pairs_hodnmatch50_rperp10_rsd50.csv \\
        --cuts 5,10,15,20,30,50
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from geometry import band_profile, bridge_excess
from sim_utils import (
    make_two_halo_template,
    open_kappa_memmap,
    reflect_symmetrize_map,
    setup_logging,
)
from stack_pairs_sim import stack_pairs
from summarize_sim_sensitivity import load_map

logger = logging.getLogger(__name__)


def jackknife_error(loo: np.ndarray) -> float:
    """Delete-one error with the (K-1)/K factor."""
    k = len(loo)
    return float(np.sqrt((k - 1) / k * np.sum((loo - loo.mean()) ** 2)))


def shell_stacks(pairs: pd.DataFrame, rpar: np.ndarray, cuts: list[float],
                 kappa_map, info, args) -> tuple[np.ndarray, np.ndarray]:
    """Stack every pair exactly once, bucketed by (spatial block, LOS shell).

    The cuts are *nested* -- r_par <= 5 is a subset of r_par <= 10 and so on --
    so stacking each cut independently would re-stack the inner pairs once per
    outer cut.  For six cuts on the r_perp = 5 catalog that is 4.3M pair-stacks
    to extract information carried by 1.26M pairs.

    Instead each pair goes into the disjoint shell (cuts[i-1], cuts[i]] it falls
    in, and any cut's stack is the running sum of the shells inside it.  Exactly
    equivalent, 3.4x less work.

    Returns arrays indexed [block, shell] of summed maps and pair counts.
    """
    cell = info.box_size_hmpc / args.blocks_per_side
    bx = np.clip((pairs["pair_center_x"].to_numpy() // cell).astype(int),
                 0, args.blocks_per_side - 1)
    by = np.clip((pairs["pair_center_y"].to_numpy() // cell).astype(int),
                 0, args.blocks_per_side - 1)
    block = bx * args.blocks_per_side + by

    # np.searchsorted puts each pair in the first shell whose upper edge it is
    # within; anything beyond the last cut is dropped, matching "cut <= C".
    shell = np.searchsorted(np.asarray(cuts), rpar, side="left")

    n_blocks, n_shells = args.blocks_per_side ** 2, len(cuts)
    g = args.grid_size
    sums = np.zeros((n_blocks, n_shells, g, g))
    counts = np.zeros((n_blocks, n_shells))

    for b in range(n_blocks):
        in_block = block == b
        for s in range(n_shells):
            sub = pairs[in_block & (shell == s)]
            if len(sub) == 0:
                continue
            mean_map = stack_pairs(
                pairs=sub,
                kappa_map=kappa_map,
                pixel_size_hmpc=info.pixel_size_hmpc,
                map_box_size_hmpc=info.box_size_hmpc,
                stack_box_size_hmpc=args.box_size,
                grid_size=g,
                normalize_separation=False,
                normalized_half_size=2.5,
            ).astype(np.float64)
            sums[b, s] = mean_map * len(sub)
            counts[b, s] = len(sub)
    return sums, counts


def scan_catalog(rperp_center: float, path: Path, cuts: list[float], args) -> list[dict]:
    pairs_all = pd.read_csv(path)
    if "r_parallel_rsd" not in pairs_all.columns:
        raise ValueError(f"{path} lacks r_parallel_rsd; regenerate with find_pairs_sim.py.")
    rpar = pairs_all["r_parallel_rsd"].to_numpy()
    logger.info("%s: %d pairs, r_par,rsd up to %.1f", path.name, len(pairs_all), rpar.max())

    kappa_map, info = open_kappa_memmap(args.kappa_map)
    axis = np.linspace(-0.5 * args.box_size, 0.5 * args.box_size, args.grid_size)
    x_grid, y_grid = np.meshgrid(axis, axis)
    single = load_map(Path(args.single))
    template = make_two_halo_template(single, x_grid, y_grid, rperp_center)

    # A cut at or beyond the catalog depth is meaningful -- it simply means "no
    # cut", and it is the reference point showing maximum dilution. Only drop a
    # cut that would duplicate an earlier one, i.e. the second and subsequent
    # cuts past the depth. (The catalogs are built at r_par,rsd <= 50 but the
    # largest value sits fractionally below 50.0, so a naive `c <= rpar.max()`
    # guard silently discarded the widest cut.)
    usable, depth_reached = [], False
    for c in sorted(cuts):
        if depth_reached:
            logger.warning("Cut %.0f duplicates the full catalog; dropped.", c)
            continue
        usable.append(c)
        if c >= rpar.max():
            depth_reached = True
    sums, counts = shell_stacks(pairs_all, rpar, usable, kappa_map, info, args)

    # Reflection symmetrization matches the BOSS chain. The band statistics are
    # exactly invariant under it (the masks are themselves symmetric in +/-x and
    # +/-y), so this affects the profile curve only, not the number.
    def residual(sum_map: np.ndarray, n: float) -> np.ndarray:
        return reflect_symmetrize_map((sum_map / n).astype(np.float32)) - template

    rows = []
    for i, cut in enumerate(usable):
        # Cut i is the running sum over shells 0..i.
        inner_sum = sums[:, : i + 1].sum(axis=(0, 1))
        inner_count = counts[:, : i + 1].sum()
        full = residual(inner_sum, inner_count)
        value = bridge_excess(full, axis, rperp_center)

        per_block_sum = sums[:, : i + 1].sum(axis=1)
        per_block_count = counts[:, : i + 1].sum(axis=1)
        loo = np.array([
            bridge_excess(residual(inner_sum - per_block_sum[b],
                                   inner_count - per_block_count[b]),
                          axis, rperp_center)
            for b in range(len(per_block_count)) if per_block_count[b] > 0
        ])
        err = jackknife_error(loo)

        rows.append({
            "rperp_center_hmpc": rperp_center,
            "rpar_cut_hmpc": cut,
            "n_pairs": int(inner_count),
            "bridge_excess": value,
            "bridge_excess_err": err,
            "sim_snr": value / err if err > 0 else np.nan,
            "profile": band_profile(full, axis),
        })
        logger.info("  r_perp %.0f, cut %4.0f: %9d pairs, bridge %+7.3fe-4 +/- %.3fe-4",
                    rperp_center, cut, int(inner_count), value * 1e4, err * 1e4)
    return rows


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Signal versus line-of-sight cut, in simulation.")
    p.add_argument("--catalog", action="append", required=True,
                   help="rperp_center=path/to/pairs.csv (repeatable).")
    p.add_argument("--cuts", default="5,10,15,20,30,50")
    p.add_argument("--kappa-map",
                   default="analysis/sim/results/kappa_map_l0p1_s8arcmin.float32")
    p.add_argument("--single",
                   default="analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv")
    p.add_argument("--reference-cut", type=float, default=None,
                   help="Cut to normalise the figure of merit against. Default: the "
                        "production cut, 5 for r_perp=5 and 10 otherwise.")
    p.add_argument("--blocks-per-side", type=int, default=5)
    p.add_argument("--grid-size", type=int, default=101)
    p.add_argument("--box-size", type=float, default=100.0)
    p.add_argument("--output", default="analysis/sim/results/rpar_scan.csv")
    args = p.parse_args(argv)
    args.cuts = [float(c) for c in args.cuts.split(",")]
    return args


def main(argv=None):
    args = parse_args(argv)
    setup_logging()

    rows = []
    for spec in args.catalog:
        label, _, path = spec.partition("=")
        if not path:
            raise ValueError(f"--catalog expects rperp=path, got {spec!r}")
        rows.extend(scan_catalog(float(label), Path(path), args.cuts, args))

    df = pd.DataFrame(rows)
    profiles = df.pop("profile")

    # Figure of merit, normalised per r_perp against the production cut so the
    # arbitrary overall scale drops out.
    df["fom"] = np.nan
    for rperp, grp in df.groupby("rperp_center_hmpc"):
        ref_cut = args.reference_cut if args.reference_cut is not None else (
            5.0 if rperp == 5.0 else 10.0)
        ref = grp[grp.rpar_cut_hmpc == ref_cut]
        if ref.empty:
            logger.warning("No reference cut %.0f for r_perp %.0f; skipping FoM.", ref_cut, rperp)
            continue
        r = ref.iloc[0]
        fom = (grp.bridge_excess * np.sqrt(grp.n_pairs)) / (r.bridge_excess * np.sqrt(r.n_pairs))
        df.loc[grp.index, "fom"] = fom

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    prof = pd.DataFrame(
        {f"rperp{int(r.rperp_center_hmpc)}_cut{int(r.rpar_cut_hmpc)}": p
         for (_, r), p in zip(df.iterrows(), profiles)},
        index=np.linspace(-0.5 * args.box_size, 0.5 * args.box_size, args.grid_size))
    prof.to_csv(out.with_name(out.stem + "_profiles.csv"))

    logger.info("\nSignal versus LOS cut (bridge excess, kappa x 1e4):")
    logger.info("  r_perp   cut     n_pairs     signal      err   sim S/N   pred. BOSS S/N")
    for rperp, grp in df.groupby("rperp_center_hmpc"):
        for _, r in grp.iterrows():
            logger.info("  %6.0f  %4.0f  %10d  %+8.3f  %7.3f   %7.2f   %13.2f",
                        rperp, r.rpar_cut_hmpc, r.n_pairs, r.bridge_excess * 1e4,
                        r.bridge_excess_err * 1e4, r.sim_snr, r.fom)
        best = grp.loc[grp.fom.idxmax()] if grp.fom.notna().any() else None
        if best is not None:
            logger.info("  -> best cut for r_perp %.0f: %.0f h^-1 Mpc (%.2fx current)",
                        rperp, best.rpar_cut_hmpc, best.fom)
    logger.info("Saved -> %s", out)


if __name__ == "__main__":
    main()
