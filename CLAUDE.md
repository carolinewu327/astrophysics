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

- **`lib/constants.py`** — Pipeline-wide constants:
  - `GRID_SIZE = 100`, `BOX_SIZE_HMPC = 100.0`, `FWHM_ARCMIN = 8.0`, `NSIDE = 2048`

- **`lib/geometry.py`** — Pure math/coordinate transforms (no file I/O):
  - `fast_icrs_to_galactic(ra, dec)` — RA/Dec to Galactic l/b via rotation matrix
  - `angular_separation(l1, b1, l2, b2)` — Great-circle separation (returns radians)
  - `symmetrize_map(kappa_map, pwr)` — Radial symmetry averaging
  - `reflect_symmetrize_map(kappa_map)` — Reflection symmetry across Y-axis (for pair stacking)
  - `radial_profile(arr, sigma, zoom)` — Extract 1D radial profile from 2D map

- **`lib/catalog.py`** — Data I/O, catalog loading, Planck map loading:
  - `setup_logging()` — Configure root logger with timestamps
  - `resolve_catalog_path(data_dir, dataset, region, catalog_type)` — Auto-resolve FITS file paths
  - `resolve_planck_paths(data_dir)` — Return (alm_path, mask_path) for Planck data
  - `load_catalog(path, weights, random_fraction, z_min, z_max)` — Full catalog loader with weighting
  - `load_catalog_lightweight(path, columns, fraction, z_min, z_max)` — Memory-efficient loader for large random catalogs (memmap + column selection)
  - `preprocess_catalog_galactic(data, weights)` — Convert RA/Dec/Z to Galactic coords + comoving distances
  - `load_kappa_map(alm_file, mask_file, nside, fwhm_arcmin)` — Load Planck ALM, smooth, synthesize HEALPix map

### `analysis/boss/scripts/` — Standalone CLI scripts

All scripts accept `--help` for full argument documentation. Legacy Jupyter notebooks are preserved in `legacy/`.

- **`find_pairs.py`** — Find galaxy or random-point pairs by separation criteria
- **`stack_single.py`** — Stack single galaxies/randoms against the Planck κ map (with jackknife errors)
- **`stack_pairs.py`** — Stack galaxy pairs from a pair catalog against the Planck κ map
- **`plot_results.py`** — Load stacked κ CSVs, compute derived maps, generate analysis plots

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

### 2. Find galaxy/random pairs (`find_pairs.py`)
Identifies pairs based on parallel and perpendicular separation criteria using vectorized distance calculations with multiprocessing.

```bash
python find_pairs.py --dataset BOSS --region North --catalog-type galaxy
python find_pairs.py --dataset BOSS --region North --catalog-type random --fraction 0.10
python find_pairs.py --dataset BOSS --region South --catalog-type galaxy --rpar 25 --rperp-min 10 --rperp-max 15
```

**Key parameters:**
- `--rpar`: Maximum parallel (line-of-sight) distance in Mpc/h (default: 20)
- `--rperp-min` / `--rperp-max`: Perpendicular distance range in Mpc/h (default: 18–22)
- `--n-processes`: Parallel workers (default: auto = 75% of cores, capped at 12)
- `--chunk-size`: Galaxies per chunk (default: 10000)

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

**Transformations:**
- `fast_icrs_to_galactic()` in `lib/geometry.py`: Converts RA/Dec to Galactic l/b
- Comoving distance: Uses Planck18 cosmology via `astropy.cosmology`
- Physical units: Distances in h⁻¹ Mpc (comoving distance × h)

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
1. **Radial symmetry** (`symmetrize_map()` in `lib/geometry.py`): Averages in radial bins from center — used for single-galaxy stacks
2. **Reflection symmetry** (`reflect_symmetrize_map()` in `lib/geometry.py`): Averages across Y-axis — used for pair stacks to preserve asymmetry along the pair axis

### Weighting Schemes
- **BOSS CMASS**: `WEIGHT_SEEING × WEIGHT_STAR × (WEIGHT_NOZ + WEIGHT_CP - 1)`
- **eBOSS**: `WEIGHT_NOZ × WEIGHT_SYSTOT`
- Weights account for observational systematics (seeing, stellar density, etc.)
- Random catalogs use uniform weights (w = 1)

### Performance Optimization
The pair-finding algorithm uses:
- Distance-sorted galaxy lists for efficient neighbor search
- Binary search to limit search range based on `r_par_max`
- Two-pass angular filtering: rough cosb-corrected filter, then exact arccos
- Vectorized angular calculations (avoid Python loops)
- Multiprocessing with checkpoint saves every 10 chunks

### Common Issues

**Theta out of bounds**: Pairs near Galactic poles may produce invalid θ values when mapping to HEALPix. These are skipped with a logged warning.

**Memory usage**: Large catalogs (>500k galaxies) with small `r_perp_min` (<5 Mpc/h) can generate millions of pairs. Consider:
- Using `load_catalog_lightweight()` for random catalogs (reads only needed columns via memmap)
- Increasing `r_perp_min` to reduce pair count
- Reducing `--n-processes` to lower memory footprint
- Using `--fraction` to subsample (especially for randoms)

**Coordinate wrapping**: Galactic longitude wraps at 360°. The code handles this in pair stacking, but be careful when manually computing angular separations.

## File Naming Conventions

- Single galaxy maps: `kappa_single_galaxy_{dataset}_{region}.csv`
- Single random maps: `kappa_single_random_{dataset}_{region}.csv`
- Error maps: `error_single_{dataset}_{region}.csv`
- Pair catalogs: `data/paircatalogs/{dataset}/{type}_pairs_{dataset}_{region}_{r_par}_{r_perp_min}_{r_perp_max}hmpc.csv`
- Pair-stacked maps: `kappa_pairs_{label}_{dataset}_{region}.csv`
- Corrected single: `kappa_corrected_single_{dataset}.csv`
- Corrected pairs: `kappa_corrected_pairs_{separation}_{dataset}.csv`
- Control pair: `kappa_control_pair_{separation}_{dataset}.csv`
- Filament: `kappa_filament_{separation}_{dataset}.csv`

## Physical Constants

- Cosmology: Planck18 (from astropy)
- h = 0.6766 (Hubble parameter)
- HEALPix nside: 2048 (resolution for Planck lensing map)
- Smoothing FWHM: 8 arcmin (applied to Planck κ map)
