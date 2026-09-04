#!/usr/bin/env python
"""
plot_footprint.py -- Survey footprint and redshift-distribution figure.

Produces the two panels of fig:footprint_zdist in the paper:

    {output_dir}/footprint.{pdf,png}  Mollweide projection (Galactic coordinates)
                                      of the Planck lensing analysis mask with the
                                      CMASS galaxies of each region overplotted.
    {output_dir}/zdist.{pdf,png}      Weighted redshift distribution of the CMASS
                                      galaxies per region, with the redshift
                                      selection shaded.

It also writes {output_dir}/footprint_stats.txt with the sample statistics
quoted in the Data section: galaxy counts (raw and weighted) per region,
median redshift, Planck mask sky fraction, and the fraction of galaxies that
fall on unmasked Planck pixels.

Usage:
    python plot_footprint.py
    python plot_footprint.py --dataset BOSS --regions North,South --z-min 0.4 --z-max 0.7
"""

import argparse
import logging
import os

import numpy as np
import healpy as hp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from catalog import (  # noqa: E402
    load_catalog, resolve_catalog_path, resolve_planck_paths, setup_logging,
)
from geometry import fast_icrs_to_galactic  # noqa: E402

logger = logging.getLogger(__name__)

REGION_COLORS = {"North": "#1f77b4", "South": "#d62728", "NGC": "#1f77b4", "SGC": "#d62728"}


# ===========================================================================
# Data loading
# ===========================================================================
def load_region(data_dir, dataset, region):
    """Load the full galaxy catalog of one region (no redshift cut).

    Returns a dict with Galactic l, b (deg), redshift z, and weights w.
    """
    path = resolve_catalog_path(data_dir, dataset, region, "galaxy")
    weight_scheme = "CMASS" if dataset == "BOSS" else True
    cat, w = load_catalog(path, weights=weight_scheme)
    l, b = fast_icrs_to_galactic(np.asarray(cat["RA"], float), np.asarray(cat["DEC"], float))
    return {"l": l, "b": b, "z": np.asarray(cat["Z"], float), "w": np.asarray(w, float)}


# ===========================================================================
# Statistics
# ===========================================================================
def compute_stats(regions, mask, z_min, z_max):
    """Per-region and total sample statistics inside the redshift selection."""
    nside_mask = hp.get_nside(mask)
    fsky = float(np.mean(mask > 0))
    stats = {"fsky_planck_mask": fsky, "regions": {}}

    all_z, all_w, all_in = [], [], []
    for name, d in regions.items():
        sel = (d["z"] > z_min) & (d["z"] < z_max)
        z, w = d["z"][sel], d["w"][sel]
        pix = hp.ang2pix(nside_mask, d["l"][sel], d["b"][sel], lonlat=True)
        in_mask = mask[pix] > 0
        stats["regions"][name] = {
            "n_gal": int(sel.sum()),
            "n_gal_weighted": float(w.sum()),
            "median_z": float(np.median(z)),
            "weighted_mean_z": float(np.average(z, weights=w)),
            "frac_on_unmasked_pixel": float(in_mask.mean()),
            "n_gal_on_unmasked_pixel": int(in_mask.sum()),
        }
        all_z.append(z)
        all_w.append(w)
        all_in.append(in_mask)

    z = np.concatenate(all_z)
    w = np.concatenate(all_w)
    in_mask = np.concatenate(all_in)
    stats["total"] = {
        "n_gal": int(len(z)),
        "n_gal_weighted": float(w.sum()),
        "median_z": float(np.median(z)),
        "weighted_mean_z": float(np.average(z, weights=w)),
        "frac_on_unmasked_pixel": float(in_mask.mean()),
        "n_gal_on_unmasked_pixel": int(in_mask.sum()),
    }
    return stats


def write_stats(stats, path, z_min, z_max):
    lines = [
        f"Redshift selection: {z_min} < z < {z_max}",
        f"Planck lensing mask sky fraction: {stats['fsky_planck_mask']:.4f}",
        "",
    ]
    for name, s in list(stats["regions"].items()) + [("Total", stats["total"])]:
        lines += [
            f"[{name}]",
            f"  N galaxies                 : {s['n_gal']:,}",
            f"  Sum of weights             : {s['n_gal_weighted']:,.1f}",
            f"  Median z                   : {s['median_z']:.4f}",
            f"  Weighted mean z            : {s['weighted_mean_z']:.4f}",
            f"  On unmasked Planck pixel   : {s['n_gal_on_unmasked_pixel']:,} "
            f"({100 * s['frac_on_unmasked_pixel']:.2f} %)",
            "",
        ]
    with open(path, "w") as f:
        f.write("\n".join(lines))
    for line in lines:
        logger.info(line)


# ===========================================================================
# Plots
# ===========================================================================
def plot_footprint(regions, mask, z_min, z_max, nside_plot, n_points, out_base):
    """Mollweide map of the Planck mask with galaxies overplotted."""
    mask_dg = hp.ud_grade(mask.astype(float), nside_plot)
    # Two-level colour map: masked (0) light grey, unmasked (1) white.
    cmap = ListedColormap(["#d9d9d9", "#ffffff"])
    # Off-projection pixels are placed below the colour range by healpy;
    # paint them white so the surround of the ellipse is blank.
    cmap.set_under("white")
    cmap.set_bad("white")

    fig = plt.figure(figsize=(9, 5), facecolor="white")
    hp.mollview(
        mask_dg, fig=fig.number, title="", cbar=False, notext=True,
        cmap=cmap, min=0, max=1, badcolor="white", bgcolor="white",
    )
    hp.graticule(dpar=30, dmer=60, color="0.6", alpha=0.5, lw=0.5)

    rng = np.random.default_rng(12345)
    handles = []
    for name, d in regions.items():
        sel = np.flatnonzero((d["z"] > z_min) & (d["z"] < z_max))
        if len(sel) > n_points:
            sel = rng.choice(sel, n_points, replace=False)
        color = REGION_COLORS.get(name, None)
        sc = hp.projscatter(
            d["l"][sel], d["b"][sel], lonlat=True, s=0.4, color=color,
            alpha=0.5, linewidths=0, rasterized=True,
        )
        sc.set_rasterized(True)
        handles.append(Line2D([], [], marker="o", ls="", color=color, ms=5, label=f"CMASS {name}"))

    handles.append(Line2D([], [], marker="s", ls="", color="#d9d9d9", ms=9,
                          label=f"Planck lensing mask ($f_{{\\rm sky}}={np.mean(mask > 0):.2f}$)"))
    ax = plt.gca()
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08),
              ncol=3, frameon=False, fontsize=9, markerscale=1.0)
    # Longitude labels along the equator (astronomical convention: l increases leftwards).
    for lon in (120, 60, 0, -60, -120):
        hp.projtext(lon, -4, f"$l={lon % 360:d}^\\circ$", lonlat=True, fontsize=8,
                    ha="center", va="top", color="0.35")

    for ext in ("pdf", "png"):
        fig.savefig(f"{out_base}.{ext}", dpi=200, bbox_inches="tight", facecolor="white")
        logger.info("Wrote %s.%s", out_base, ext)
    plt.close(fig)


def plot_zdist(regions, z_min, z_max, out_base, dz=0.01, z_plot=(0.3, 0.85)):
    """Weighted redshift histograms per region with the selection shaded."""
    bins = np.arange(z_plot[0], z_plot[1] + dz / 2, dz)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))

    ax.axvspan(z_min, z_max, color="0.92", zorder=0, label=f"${z_min} < z < {z_max}$")

    total_hist = np.zeros(len(bins) - 1)
    for name, d in regions.items():
        color = REGION_COLORS.get(name, None)
        h, _ = np.histogram(d["z"], bins=bins, weights=d["w"])
        total_hist += h
        ax.step(bins[:-1], h, where="post", color=color, lw=1.4, label=f"CMASS {name}")
    ax.step(bins[:-1], total_hist, where="post", color="k", lw=1.0, ls="--", label="Total")

    sel_z = np.concatenate([d["z"][(d["z"] > z_min) & (d["z"] < z_max)] for d in regions.values()])
    zmed = np.median(sel_z)
    ax.axvline(zmed, color="0.3", lw=0.8, ls=":", label=f"median $z={zmed:.3f}$")

    ax.set_xlim(*z_plot)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("redshift $z$")
    ax.set_ylabel(f"weighted number of galaxies per $\\Delta z={dz}$")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()

    for ext in ("pdf", "png"):
        fig.savefig(f"{out_base}.{ext}", dpi=200)
        logger.info("Wrote %s.%s", out_base, ext)
    plt.close(fig)


# ===========================================================================
# CLI
# ===========================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Footprint and redshift-distribution figure for the Data section.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--dataset", default="BOSS", choices=["BOSS", "eBOSS"])
    p.add_argument("--regions", default="North,South",
                   help="Comma-separated regions (BOSS: North,South; eBOSS: NGC,SGC).")
    p.add_argument("--z-min", type=float, default=0.4)
    p.add_argument("--z-max", type=float, default=0.7)
    p.add_argument("--data-dir", default="data")
    p.add_argument("--output-dir", default="output/plots")
    p.add_argument("--nside-plot", type=int, default=256,
                   help="HEALPix nside to which the mask is degraded for plotting.")
    p.add_argument("--n-points", type=int, default=60000,
                   help="Maximum galaxies drawn per region on the footprint map.")
    return p.parse_args()


def main():
    setup_logging()
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    regions = {}
    for region in [r.strip() for r in args.regions.split(",") if r.strip()]:
        logger.info("Loading %s %s galaxies", args.dataset, region)
        regions[region] = load_region(args.data_dir, args.dataset, region)

    _, mask_path = resolve_planck_paths(args.data_dir)
    logger.info("Loading Planck lensing mask from %s", mask_path)
    mask = hp.read_map(mask_path)

    stats = compute_stats(regions, mask, args.z_min, args.z_max)
    write_stats(stats, os.path.join(args.output_dir, "footprint_stats.txt"), args.z_min, args.z_max)

    plot_footprint(regions, mask, args.z_min, args.z_max, args.nside_plot, args.n_points,
                   os.path.join(args.output_dir, "footprint"))
    plot_zdist(regions, args.z_min, args.z_max, os.path.join(args.output_dir, "zdist"))


if __name__ == "__main__":
    main()
