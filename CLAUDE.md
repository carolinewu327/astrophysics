# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an astrophysics research project analyzing gravitational lensing convergence (κ) maps by stacking galaxy positions from BOSS and eBOSS surveys against Planck CMB lensing data. The primary goal is to detect filamentary structures in the cosmic web by analyzing galaxy pairs at various separations.

## Key Data Sources

- **Galaxy Catalogs**: BOSS CMASS and eBOSS LRG catalogs (FITS format)
  - Expected locations: `data/BOSS/galaxy_DR12v5_CMASS_{region}.fits.gz`
  - Expected locations: `data/eBOSS/eBOSS_{catalog}_clustering_data-{region}-vDR16.fits`
  - Regions: "North"/"South" (BOSS) or "NGC"/"SGC" (eBOSS)

- **Planck Lensing Data**: CMB lensing convergence maps
  - ALM file: `data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits`
  - Mask file: `data/planck/COM_Lensing_4096_R3.00/mask.fits`

- **Output Data**:
  - Pair catalogs: `data/paircatalogs/{dataset}/` (e.g. `data/paircatalogs/BOSS/`)
  - Stacked κ maps and derived CSVs: `analysis/boss/results/`
  - Plots: `output/plots/`

## Code Structure

### `lib/` — Shared library modules

- **`lib/jackknife.py`** — Equal-area HEALPix jackknife regions over a joint footprint:
  - `build_jackknife_regions(ra, dec, nside, min_fill)` — Tessellate a footprint traced by randoms. Cells above `min_fill` of the mean occupancy become region cores; partially covered edge cells are **merged into the lightest of their 4 nearest cores** (not dropped), so no object falls outside the tessellation. At `nside=10` (34.4 deg² cells) the joint CMASS footprint gives **287 regions** with ~11% RMS size scatter.
  - `JackknifeRegions` — dataclass with `assign(ra, dec)`, `save`/`load`, and a `digest` that lets the combiner reject accumulators built on different tessellations.

### `analysis/boss/scripts/` — Standalone CLI scripts

All scripts accept `--help` for full argument documentation. Legacy Jupyter notebooks are preserved in `legacy/`.

- **`find_and_stack_pairs.py`** — Combined find-and-stack pipeline that never materializes a pair catalog. Used for full random catalogs where the pair list would be too large to write to disk (~2B pairs / hundreds of GB). Workers stack pairs inline and return per-chunk accumulator grids, which the main process reduces into the final stacked map. Multiprocess mode requires Linux `fork()` for kappa-map sharing via copy-on-write; local serial validation works with `--n-processes 1` on any platform.
- **`stack_single_jk.py`** — **Preferred single-object stacker for anything needing error bars.** Stacks all survey regions as one joint sample and accumulates `sum(w*κ)` / `sum(w)` *per jackknife region* in a single pass, so every leave-one-out estimate is a subtraction rather than a re-stack (287 regions at ~1× the cost of one stack, vs ~287×).
- **`combine_jackknife.py`** — Combine galaxy + random accumulators into the corrected map, per-pixel jackknife errors, and radial profiles with a **full bin-to-bin covariance**. Supports linear, log, or explicit `--bin-edges` binning.
- **`region_split_check.py`** — Diagnostic: re-derives North-only, South-only, unweighted-average, count-weighted, and joint profiles from the same accumulators to isolate the effect of the region-combination rule.

## Analysis Pipeline

The analysis follows this multi-stage workflow, run as CLI scripts:

### 1. Stack single galaxies (`stack_single.py`)
Stacks ALL galaxies from a catalog to create a baseline κ map. Includes jackknife error estimation for galaxy catalogs.

```bash
python stack_single.py --dataset BOSS --region North --catalog-type galaxy
python stack_single.py --dataset BOSS --region North --catalog-type random --fraction 0.10
```

**Output:** `analysis/boss/results/kappa_single_{catalog_type}_{dataset}_{region}.csv`
Galaxy runs also produce: `analysis/boss/results/error_single_{dataset}_{region}.csv`

### 1b. Joint-footprint single stack with jackknife errors (`stack_single_jk.py` + `combine_jackknife.py`)

**This is the default path for the single-galaxy measurement.** North and South are stacked as one joint sample, so there is no per-region averaging step — the relative weighting of the two footprints follows from the objects themselves. Galaxies and randoms are stacked separately but share one tessellation, so the combiner can delete the same patch of sky from both and put the error on the *subtracted* map.

```bash
PYTHONPATH=lib python analysis/boss/scripts/stack_single_jk.py \
    --dataset BOSS --regions North,South --catalog-type galaxy --n-processes 6
PYTHONPATH=lib python analysis/boss/scripts/stack_single_jk.py \
    --dataset BOSS --regions North,South --catalog-type random \
    --fraction 1.0 --label frac100 --n-processes 24 --chunk-size 20000
PYTHONPATH=lib python analysis/boss/scripts/combine_jackknife.py \
    --dataset BOSS --regions North,South --tag _scw --random-tag _scw_frac100
```

**Outputs:** per-region accumulators in `analysis/boss/results/jk/acc_single_*.npz` (the reusable product — any rebinning or hemisphere split follows from these without re-stacking), plus corrected maps, `error_single_*_joint.csv`, `profile_corrected_*.csv`, `profile_cov_corrected_*.csv`, and `profile_bands_corrected_*.csv`.

**Two estimator points that matter:**
- **Never collapse per-pixel errors into a profile error by assuming pixels are independent.** The 8 arcmin Planck beam is ~3.3 h⁻¹ Mpc at z = 0.55, over 3× the 1 h⁻¹ Mpc pixel. Measured on the joint stack, the independent-bin approximation understates band errors by **2.2–2.6×**. `combine_jackknife.py` reports both so the gap stays visible.
- **1/Σ_crit² inverse-variance weighting is the default estimator** (`--no-sigma-crit-weight` for the unweighted null test). Unweighted outputs drop the `_scw` tag and so never overwrite the default products.

**Known cost:** preprocessing the 40.9M-row full random catalog spent ~110 min before stacking began on an 8 GB machine — `fast_icrs_to_galactic` allocates two ~1 GB temporaries for a catalog that size and the machine swapped. Stacking itself was ~60 min on 5 cores. Chunk the coordinate transform if this becomes a bottleneck.

### 2. Find galaxy/random pairs (`find_pairs.py`)
Identifies pairs based on parallel and perpendicular separation criteria using vectorized distance calculations with multiprocessing.

```bash
python find_pairs.py --dataset BOSS --region North --catalog-type galaxy
python find_pairs.py --dataset BOSS --region North --catalog-type random --fraction 0.10
python find_pairs.py --dataset BOSS --region South --catalog-type galaxy --rpar 25 --rperp-min 10 --rperp-max 15
```

The script uses `imap_unordered` for load balancing (chunks complete in arbitrary order) and a vectorized `Dmid` computation. Output rows therefore appear in completion order, not chunk order — sort by full pair coordinates if you need a canonical ordering.

**Output:** `data/paircatalogs/{dataset}/{type}_pairs_{dataset}_{region}_{r_par}_{r_perp_min}_{r_perp_max}hmpc.csv`
Columns: l1, b1, z1, w1, Dc1, ID1, l2, b2, z2, w2, Dc2, ID2, Dmid

### 3. Stack galaxy pairs (`stack_pairs.py`)
Takes a pair catalog and stacks κ values at grid positions oriented along the pair axis. Uses a 101×101 grid (vs 100×100 for single stacking). Applies reflection symmetry.

```bash
python stack_pairs.py \
    --pair-catalog data/paircatalogs/BOSS/galaxy_pairs_BOSS_North_20.0_18.0_22.0hmpc.csv \
    --label galaxy_20 --region North
```

**Output:** `analysis/boss/results/kappa_pairs_{label}_{dataset}_{region}.csv`

### 3b. Combined find-and-stack for full randoms (`find_and_stack_pairs.py`)
For runs where the pair catalog would be too large to materialize on disk (full random catalogs at `--fraction 1.0`), use this script instead of steps 2 + 3. It runs the same pair-finding logic in a worker pool but stacks each pair into per-chunk accumulator grids inline, returning ~hundreds of KB per chunk instead of millions of pair rows. Pairs are never written anywhere.

```bash
python find_and_stack_pairs.py --dataset BOSS --region South --catalog-type random \
    --fraction 1.0 --rpar 5 --rperp-min 4 --rperp-max 6 \
    --label random_5_frac100 \
    --checkpoint-path analysis/boss/results/checkpoints/random_5_frac100_BOSS_South.npz \
    --n-processes 24 --chunk-size 5000
```

**Key behaviors:**
- Loads the kappa map and mask as module globals *before* `Pool()` creation. Multiprocess mode explicitly requests `fork()` (via `mp.get_context("fork")`) so workers can inherit the kappa map via copy-on-write. This is the intended production mode on Linux. Local validation should use `--n-processes 1`; multiprocess macOS may work but is not the production target.
- Atomic checkpoint writes (`os.replace()`); resume via `--resume-checkpoint`. The checkpoint stores the set of completed chunk IDs (not just a "last chunk"), which is required for `imap_unordered`. Runs before 2026-09 never saved a final checkpoint, so their checkpoints are stale and their raw sums survive nowhere; only the symmetrized CSV is valid.
- `--seed` controls catalog subsampling for `--fraction < 1.0`. If omitted, a fresh seed is auto-generated and saved in the checkpoint so resume uses the *same* subset. Resume across `--fraction < 1` runs without a seed would silently mix accumulator state from different subsamples.
- Z is upgraded to float64 inside this script (vs the float32 from BOSS FITS). This produces ~5e-6 differences in the final 101×101 kappa map vs the legacy two-step pipeline (well below physical noise ~1e-3). The new script is more accurate; the old reference is the artifact.
- Sidecar `.meta.json` written next to the output CSV with full run provenance (timing, chunk counts, seed, all CLI args).

**Output:** `analysis/boss/results/kappa_pairs_{label}_{dataset}_{region}.csv` (same naming as `stack_pairs.py`)

### 3c. Wide-bin pair analysis (2026-09)

r⊥ groupings (4–6, 6–15, 15–25 and splits) at r∥ ≤ 5 with 1/Σ_crit² weighting and a separation-averaged control. Design, run sheet, and test record: `notes/wide_bin_implementation_plan.md`. All its products live under `widebin_rpar5/` directories and never overwrite archived results; `--separations` in `combine_filament_jackknife.py` still reproduces the archived 5/10/20 numbers exactly.

### 4. Generate plots and analysis (`plot_results.py`)
Loads all stacked CSV maps, computes derived maps, and generates 8 map plots + 3 profile plots. Works entirely from CSV files (no FITS data needed).

```bash
python plot_results.py --separation 20 --regions North,South
```

**Derived maps computed:**
1. Single galaxies (1) and single randoms (3) → corrected single (5) = (1)−(3)
2. Galaxy pairs (2) and random pairs (4) → corrected pairs (6) = (2)−(4)
3. Corrected single (5) → control pair map (7) via shifted copies at ±(sep/2)
4. Filament map (8) = corrected pairs (6) − control pair (7)

**Output:** Plots in `output/plots/`, derived CSVs in `analysis/boss/results/`

## Coordinate Systems

The analysis uses **Galactic coordinates** (l, b) instead of equatorial (RA, Dec) to avoid coordinate singularities at high declinations.

## Important Implementation Details

### Grid Orientation for Pair Stacking
When stacking galaxy pairs, the coordinate system is oriented such that:
- **X-axis**: Points from galaxy 1 to galaxy 2 (along pair axis)
- **Y-axis**: Perpendicular to pair axis (where filament signal is expected)
- Enforces consistent ordering: galaxies are sorted so l2 > l1 (after wrapping)
- Rotation angle θ computed from Δl and Δb between galaxies

### Grid Size Difference
- Single-galaxy stacking uses **100×100** grid (`lib/constants.py`: `GRID_SIZE = 100`)
- Pair stacking uses **101×101** grid (`stack_pairs.py`: `GRID_RES = 101`)
- `plot_results.py` includes `reconcile_shapes()` to handle this mismatch when subtracting maps

### Symmetrization Methods
Two approaches are used:
1. **Radial symmetry** (`symmetrize_map()` in `lib/geometry.py`): Averages in radial bins about the physical `(N-1)/2` center — between pixels 49/50 for a 100×100 single stack and on pixel 50 for a 101×101 stack. Never use `N//2` for an even grid; archived simulation singles made with that convention must be regenerated, not re-symmetrized.
2. **Reflection symmetry** (`reflect_symmetrize_map()` in `lib/geometry.py`): Averages across Y-axis — used for pair stacks to preserve asymmetry along the pair axis

### Common Issues

**Theta out of bounds**: Pairs near Galactic poles may produce invalid θ values when mapping to HEALPix. These are skipped with a logged warning.

**Memory usage**: Large catalogs (>500k galaxies) with small `r_perp_min` (<5 Mpc/h) can generate millions of pairs. Consider:
- Using `load_catalog_lightweight()` for random catalogs (reads only needed columns via memmap)
- Increasing `r_perp_min` to reduce pair count
- Reducing `--n-processes` to lower memory footprint
- Using `--fraction` to subsample (especially for randoms)

**Coordinate wrapping**: Galactic longitude wraps at 360°. The code handles this in pair stacking, but be careful when manually computing angular separations.

