#!/usr/bin/env python
"""One summary sheet per pair separation: the three maps plus the band-differenced profile.

Requested layout, per separation:
  (a) physical pair stack   -- galaxy pairs minus random pairs
  (b) control pair stack    -- the single-galaxy stack superposed at +/- sep/2
  (c) filament map          -- (a) minus (b)
  (d) profile comparison    -- for each of the three maps, the mean of the central
      rows around the pair axis minus the mean of rows farther away, as a function
      of position along the pair axis.  Panel (d) is the bridge statistic unrolled
      in x: the bridge excess quoted in the tables is the average of the filament
      curve over |x| <= 0.35*sep.

Panel (d) uses the same bands as the bridge statistic, drawn on panel (c) so the
maps and the profile are visibly connected.  The bands are fixed in h^-1 Mpc and
identical at every separation, so the three separations are measured with one
ruler:
  central band      the three rows centred at Y = -1, 0, +1
  off-centre band   the rows centred at |Y| = 2 ... 10

Errors are the delete-one jackknife over the 287 shared HEALPix cells.  Every term
carrying galaxy shot noise -- the pair stack and both single stacks feeding the
control -- has the same sky removed in each realization.  The random-pair map is
held fixed (see combine_filament_jackknife.py).

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/plot_separation_summary.py
"""

from __future__ import annotations

import argparse
import json
import logging
import os

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from catalog import setup_logging
from combine_jackknife import load_accumulator, leave_one_out_maps, total_map
from geometry import (
    BRIDGE_HALF_X_FRAC,
    CENTRAL_HALF_Y_HMPC,
    OFF_HI_Y_HMPC,
    OFF_LO_Y_HMPC,
    band_profile,
    reflect_symmetrize_map,
    symmetrize_map,
    two_halo_template,
)
from jackknife import jackknife_error
from weighting_sensitivity import (
    SEPARATIONS,
    axis_for_map,
    combine,
    load_map,
    reconcile_shapes,
)

logger = logging.getLogger(__name__)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, INK_MUTED = "#0b0b0b", "#52514e", "#8a8985"
GRID = "#e3e2df"

# Band definitions come from lib/geometry.py, which is shared with the
# simulation pipeline; the names below are local aliases for readability.
CENTRAL_HALF_Y = CENTRAL_HALF_Y_HMPC
OFF_LO_Y, OFF_HI_Y = OFF_LO_Y_HMPC, OFF_HI_Y_HMPC
BRIDGE_X_FRAC = BRIDGE_HALF_X_FRAC

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "font.size": 8.5,
    "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white",
})



def boss_rperp_bin(key):
    """The BOSS sample's transverse bin, from its pair-catalog stem.

    The stem is "<rpar_max>_<rperp_min>_<rperp_max>", e.g. "10.0_18.0_22.0".
    """
    _, lo, hi = SEPARATIONS[key]["pair_catalog"].split("_")
    return float(lo), float(hi)


def assert_comparable(key, meta, stem, allow_rperp_mismatch=False):
    """Refuse to plot a mock stack that is not the same measurement as BOSS.

    Three things have to agree before a mock curve can sit beside an observed
    one, and only the first was ever checked:

      r_par cut     already asserted below
      r_perp bin    never asserted, which is how a 19-21 mock stack ended up in
                    a committed figure against BOSS's 18-22
      band scoring  the JSON must carry the fixed physical bands of
                    lib/geometry.py, not the separation-scaled legacy ones

    A JSON written before the fixed-band change has no band_definition block at
    all; that is reported as its own case rather than as a mismatch, because the
    fix is to re-run the stacker, not to edit the file.
    """
    want_lo, want_hi = boss_rperp_bin(key)
    got_lo, got_hi = meta.get("rperp_min_hmpc"), meta.get("rperp_max_hmpc")
    if got_lo is None or got_hi is None:
        raise ValueError(
            f"r_perp = {key}: {stem}.json predates the transverse-bin record. "
            "Re-run jackknife_pair_stack.py with --rperp-min/--rperp-max.")
    if not (np.isclose(want_lo, got_lo) and np.isclose(want_hi, got_hi)):
        msg = (f"r_perp = {key}: BOSS uses the {want_lo:g}-{want_hi:g} h^-1 Mpc "
               f"bin but {stem} was built from {got_lo:g}-{got_hi:g}. These are "
               "different samples; regenerate the pair catalogue for the BOSS "
               "bin with find_pairs_sim.py.")
        if not allow_rperp_mismatch:
            raise ValueError(msg)
        # The opt-out exists for one situation: publishing a provisional figure
        # while the corrected catalogue is still being generated.  It has to be
        # asked for explicitly, it says so at WARNING, and the caption is
        # expected to carry the same caveat -- see paper/README.md.
        logger.warning("PROVISIONAL: %s Proceeding because the caller passed "
                       "allow_rperp_mismatch.", msg)

    band = meta.get("band_definition", {}).get("fixedband_keys")
    if band is None:
        raise ValueError(
            f"r_perp = {key}: {stem}.json has no fixed-band scoring. Its bridge "
            "numbers come from separation-scaled Y bands and are not the "
            "statistic BOSS reports. Re-run jackknife_pair_stack.py.")
    for name, want in (("central_half_y_hmpc", CENTRAL_HALF_Y_HMPC),
                       ("off_lo_hmpc", OFF_LO_Y_HMPC),
                       ("off_hi_hmpc", OFF_HI_Y_HMPC),
                       ("bridge_half_x_frac", BRIDGE_HALF_X_FRAC)):
        if not np.isclose(band.get(name, np.nan), want):
            raise ValueError(
                f"r_perp = {key}: {stem} was scored with {name} = "
                f"{band.get(name)}, but lib/geometry.py now says {want}. The "
                "mock and the observation would be different statistics.")

def build_terms(R: str, dataset: str, regions: list[str], sep_key: str,
                single_tag: str, single_random_tag: str,
                pair_suffix: str = "", symmetrize_pairs: bool = True) -> dict:
    """Total and leave-one-out maps for the three terms at one separation.

    Mirrors combine_filament_jackknife.py; the caller cross-checks the resulting
    bridge excess against that script's CSV so the two cannot drift apart.
    """
    jk = os.path.join(R, "jk")
    region_label = "_".join(regions)
    cfg = SEPARATIONS[sep_key]
    sep = cfg["center"]

    gal = load_accumulator(os.path.join(
        jk, f"acc_single_galaxy{single_tag}_{dataset}_{region_label}.npz"))
    rnd = load_accumulator(os.path.join(
        jk, f"acc_single_random{single_random_tag}_{dataset}_{region_label}.npz"))
    pair = load_accumulator(os.path.join(
        jk, f"acc_pairs_galaxy_{cfg['galaxy_label']}{pair_suffix}_{dataset}_{region_label}.npz"))
    for name, acc in (("random", rnd), ("pair", pair)):
        if acc["jk_digest"] != gal["jk_digest"]:
            raise ValueError(f"{name} accumulator uses a different tessellation.")

    g, gp = gal["grid_size"], pair["grid_size"]
    single_tot = total_map(gal["sum_wk"], gal["sum_w"]) - total_map(rnd["sum_wk"], rnd["sum_w"])
    single_loo = leave_one_out_maps(gal["sum_wk"], gal["sum_w"]) - \
        leave_one_out_maps(rnd["sum_wk"], rnd["sum_w"])
    pair_tot = total_map(pair["sum_wk"], pair["sum_w"])
    pair_loo = leave_one_out_maps(pair["sum_wk"], pair["sum_w"])

    rand_maps, rand_counts = {}, {}
    for reg in regions:
        stem = f"kappa_pairs_random_{cfg['random_label']}_{dataset}_{reg}"
        rand_maps[reg] = load_map(os.path.join(R, f"{stem}.csv"))
        with open(os.path.join(R, f"{stem}.meta.json"), encoding="utf-8") as fh:
            rand_counts[reg] = float(json.load(fh)["total_pairs"])
    rand_pair = combine(rand_maps, rand_counts)

    def terms(single_flat, pair_flat):
        # Reflection-symmetrize the pair stack. The archived random-pair maps are
        # already symmetrized (find_and_stack_pairs.finalize_map), as are the
        # simulation stacks, so leaving the galaxy side raw mixes conventions
        # within a single subtraction. The signal is 4-fold symmetric anyway --
        # which galaxy is "first" is an arbitrary l2 > l1 convention, and there is
        # no preferred side perpendicular to the pair axis -- so this is noise
        # reduction, not cosmetics. Bridge statistics are exactly invariant under
        # it (the bridge and sideband masks are themselves symmetric in +/-x, y).
        pair_map = pair_flat.reshape(gp, gp)
        if symmetrize_pairs:
            pair_map = reflect_symmetrize_map(pair_map)
        cp, rp = reconcile_shapes(pair_map, rand_pair)
        physical = cp - rp
        # Control built natively on the pair grid (see geometry.two_halo_template):
        # shifting a pixel array and trimming left a half-pixel misalignment that
        # showed up as false structure at the halo peaks.
        pair_axis = axis_for_map(physical)
        gx, gy = np.meshgrid(pair_axis, pair_axis)
        control = two_halo_template(symmetrize_map(single_flat.reshape(g, g)), gx, gy, sep)
        return physical, control, physical - control

    phys, ctrl, fil = terms(single_tot, pair_tot)
    n = single_loo.shape[0]
    loo = {"physical": [], "control": [], "filament": []}
    for k in range(n):
        p, c, f = terms(single_loo[k], pair_loo[k])
        loo["physical"].append(p.ravel())
        loo["control"].append(c.ravel())
        loo["filament"].append(f.ravel())

    return {
        "sep": sep, "n_pairs": int(pair["n_pairs"].sum()),
        "physical": phys, "control": ctrl, "filament": fil,
        "loo": {k: np.asarray(v) for k, v in loo.items()},
        "axis": axis_for_map(fil),
    }


def draw_map(ax, arr, axis, sep, title, norm, zoom, show_bands=False):
    im = ax.imshow(arr, origin="lower", extent=(axis[0], axis[-1], axis[0], axis[-1]),
                   cmap="RdBu_r", norm=norm, interpolation="nearest", aspect="equal")
    ax.scatter([-sep / 2, sep / 2], [0, 0], s=34, marker="x", c=INK, linewidths=1.3,
               zorder=5)
    if show_bands:
        for yy in (CENTRAL_HALF_Y, -CENTRAL_HALF_Y):
            ax.axhline(yy, color=INK, lw=0.7, ls="-", alpha=0.55)
        for yy in (OFF_LO_Y, OFF_HI_Y, -OFF_LO_Y, -OFF_HI_Y):
            ax.axhline(yy, color=INK, lw=0.7, ls=":", alpha=0.55)
    ax.set_xlim(-zoom, zoom)
    ax.set_ylim(-zoom, zoom)
    ax.set_title(title, loc="left", color=INK, pad=5)
    ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
    ax.set_ylabel(r"$Y$  [$h^{-1}$ Mpc]")
    return im


def make_sheet(t: dict, out_base: str, bridge_ref: float | None) -> None:
    sep, axis = t["sep"], t["axis"]
    zoom = 2.5 * sep

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.0), constrained_layout=True)

    shared = max(np.nanpercentile(np.abs(t["physical"]), 99.5),
                 np.nanpercentile(np.abs(t["control"]), 99.5))
    norm_pc = TwoSlopeNorm(vmin=-shared, vcenter=0.0, vmax=shared)
    fmax = max(np.nanpercentile(np.abs(t["filament"]), 99.5), 1e-12)
    norm_f = TwoSlopeNorm(vmin=-fmax, vcenter=0.0, vmax=fmax)

    im0 = draw_map(axes[0, 0], t["physical"], axis, sep,
                   "(a) Physical pair stack  (galaxy pairs − random pairs)", norm_pc, zoom)
    im1 = draw_map(axes[0, 1], t["control"], axis, sep,
                   "(b) Control pair stack  (single stack superposed at ±sep/2)",
                   norm_pc, zoom)
    im2 = draw_map(axes[1, 0], t["filament"], axis, sep,
                   "(c) Filament map  (a − b)", norm_f, zoom, show_bands=True)
    for im, ax in ((im0, axes[0, 0]), (im1, axes[0, 1]), (im2, axes[1, 0])):
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(r"$\kappa$")
        cb.ax.tick_params(labelsize=7)

    ax = axes[1, 1]
    curves = [("Physical pairs", "physical", BLUE), ("Control pairs", "control", ORANGE),
              ("Filament (difference)", "filament", AQUA)]
    for label, key, color in curves:
        prof = band_profile(t[key], axis) * 1e4
        loo = np.array([band_profile(m.reshape(t[key].shape), axis)
                        for m in t["loo"][key]])
        _, err = jackknife_error(loo)
        err = err * 1e4
        ax.plot(axis, prof, lw=2, color=color, label=label, zorder=3)
        ax.fill_between(axis, prof - err, prof + err, color=color, alpha=0.15, lw=0,
                        zorder=1)

    ax.axhline(0, color=INK_MUTED, lw=0.8, alpha=0.6)
    for xx in (-sep / 2, sep / 2):
        ax.axvline(xx, color=INK_MUTED, lw=0.8, ls="--", alpha=0.7)
    ax.axvspan(-BRIDGE_X_FRAC * sep, BRIDGE_X_FRAC * sep, color=INK_MUTED, alpha=0.07,
               lw=0, zorder=0)
    ax.set_xlim(-zoom, zoom)
    ax.set_xlabel(r"$X$ along pair axis  [$h^{-1}$ Mpc]")
    ax.set_ylabel(r"central band − off-centre band,  $\kappa$  [$10^{-4}$]")
    ax.set_title("(d) Profile along the pair axis", loc="left", color=INK, pad=5)
    ax.legend(loc="upper right")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.text(0.02, 0.03,
            f"dashed: pair members at $X=\\pm${sep / 2:.1f}\n"
            f"shaded: bridge window $|X|\\leq${BRIDGE_X_FRAC * sep:.1f}\n"
            f"bands: rows $Y=-1,0,1$ vs rows $|Y|=2\\ldots10$",
            transform=ax.transAxes, fontsize=7, color=INK_MUTED, va="bottom")

    note = (f"BOSS CMASS, joint North+South.  $r_\\perp$ = {sep:g} $h^{{-1}}$Mpc, "
            f"{t['n_pairs']:,} pairs.  Bands are jackknife errors over 287 equal-area cells.")
    fig.suptitle(f"BOSS pair stack summary — separation {sep:g} $h^{{-1}}$ Mpc\n" + note,
                 fontsize=10.5, color=INK, ha="left", x=0.01)

    if bridge_ref is not None:
        fil_prof = band_profile(t["filament"], axis)
        window = np.abs(axis) <= BRIDGE_X_FRAC * sep
        recomputed = float(fil_prof[window].mean())
        if not np.isclose(recomputed, bridge_ref, rtol=0.02, atol=1e-7):
            logger.warning(
                "Bridge excess from panel (d) is %.3e but combine_filament_jackknife "
                "reported %.3e. The x-averaged band profile and the 2D box mean differ "
                "slightly because the box weights pixels uniformly in 2D; check if the "
                "gap is large.", recomputed, bridge_ref)
        else:
            logger.info("  panel (d) bridge excess %.3e matches the table (%.3e)",
                        recomputed, bridge_ref)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s.png / .pdf", out_base)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Per-separation summary sheets of the BOSS pair measurement.")
    parser.add_argument("--dataset", default="BOSS")
    parser.add_argument("--regions", default="North,South")
    parser.add_argument("--separations", default="5,10,20")
    parser.add_argument("--single-tag", default="_scw")
    parser.add_argument("--single-random-tag", default="_scw_frac100")
    parser.add_argument("--results-dir", default="analysis/boss/results")
    parser.add_argument("--output-dir", default="output/plots")
    args = parser.parse_args(argv)
    args.regions = [r.strip() for r in args.regions.split(",")]
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)

    ref = {}
    ref_path = os.path.join(
        args.results_dir, f"filament_jackknife_{args.dataset}_{'_'.join(args.regions)}.csv")
    if os.path.exists(ref_path):
        table = pd.read_csv(ref_path)
        ref = dict(zip(table.separation_hmpc, table.filament_bridge_excess))

    for sep_key in (s.strip() for s in args.separations.split(",")):
        sep = SEPARATIONS[sep_key]["center"]
        logger.info("=== separation %g h^-1 Mpc ===", sep)
        t = build_terms(args.results_dir, args.dataset, args.regions, sep_key,
                        args.single_tag, args.single_random_tag)
        make_sheet(t, os.path.join(args.output_dir, f"boss_summary_sep{sep_key}"),
                   ref.get(sep))


if __name__ == "__main__":
    main()
