# Kernel Smoothing for Lensing Analysis

## Overview

The naive approach of increasing `GRID_SIZE` doesn't add information because the fundamental resolution limit is the **Planck kappa map** itself:

- **HEALPix NSIDE=2048** → ~1.7 arcmin/pixel
- **FWHM=8 arcmin** smoothing already applied
- **At z~0.5**: 8 arcmin ≈ **3-4 Mpc/h** physical scale

Simply increasing `GRID_SIZE` oversamples the same blurred data. Instead, we need **kernel-based smoothing** to create continuous fields from discrete grid measurements.

## Solution: Optional Gaussian Kernel Density Estimation (KDE)

Added to `stack_single.py`:

### New Functions

1. **`extract_radial_profile_kde()`**
   - Applies Gaussian kernel smoothing in radial direction
   - Respects Planck resolution via `kernel_width_hmpc` parameter
   - Returns smooth, continuous radial profile
   - Propagates errors correctly through kernel weighting

2. **`extract_radial_profile_binned()`**
   - Traditional azimuthal averaging in discrete radial bins
   - Provided for comparison with KDE approach

### Usage

#### Command-Line Interface

```bash
# Standard stacking (existing behavior, no change)
python stack_single.py --dataset BOSS --region North --catalog-type galaxy

# Enable kernel smoothing (NEW)
python stack_single.py --dataset BOSS --region North --catalog-type galaxy \
    --extract-radial-kde \
    --kernel-width 3.5 \
    --n-radial-bins 500
```

#### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--extract-radial-kde` | False | Enable smooth radial profile extraction |
| `--kernel-width` | 3.0 | Gaussian kernel σ in h⁻¹ Mpc. Should be ≥ Planck resolution (~3-4 Mpc) |
| `--n-radial-bins` | 500 | Number of evaluation points for smooth profile |

### Output Files

When `--extract-radial-kde` is enabled, two additional CSV files are created:

1. **`radial_kde_{catalog_type}_{dataset}_{region}.csv`**
   - Columns: `r_hmpc`, `kappa`, `kappa_err` (if available)
   - Smooth Gaussian KDE profile with `n_radial_bins` points

2. **`radial_binned_{catalog_type}_{dataset}_{region}.csv`**
   - Columns: `r_hmpc`, `kappa`, `kappa_err` (if available)
   - Traditional binned profile (50 bins) for comparison

### Python API

```python
from stack_single import extract_radial_profile_kde, extract_radial_profile_binned

# After stacking
kappa_map, kappa_err, _ = stack_kappa(data, weights, kmap, mask)

# Extract smooth KDE profile
r_kde, kappa_kde, kappa_kde_err = extract_radial_profile_kde(
    kappa_map,
    kappa_err=kappa_err,
    kernel_width_hmpc=3.5,  # Match Planck resolution
    n_radial_bins=500,      # Smooth continuous output
)

# Traditional binned profile for comparison
r_bins, kappa_bins, kappa_bins_err = extract_radial_profile_binned(
    kappa_map,
    kappa_err=kappa_err,
    n_bins=50,
)
```

## Methodology

### Hard Binning (Original)

```
Grid positions → HEALPix lookup → Discrete bins
              ↓
       Each pixel assigned to ONE bin
              ↓
        Jagged, noisy profile
```

### Gaussian KDE (New)

```
Grid positions → HEALPix lookup → Radial coordinates
              ↓
   Each sample contributes to NEARBY radii
         via Gaussian kernel:
              ↓
      w(r) = exp(-Δr²/(2σ²))
              ↓
    Smooth, continuous profile
```

At each evaluation radius `r_eval`, the KDE computes:

```
κ(r_eval) = Σᵢ w(rᵢ, r_eval) × κᵢ / Σᵢ w(rᵢ, r_eval)

where: w(rᵢ, r_eval) = exp(-0.5 × ((rᵢ - r_eval) / σ)²)
```

## Choosing Kernel Width

The kernel width `σ` (sigma) controls the smoothing scale:

| Kernel Width | Effect | Use Case |
|--------------|--------|----------|
| **2-3 Mpc/h** | Minimal smoothing | Clean data, high S/N |
| **3-4 Mpc/h** | Matches Planck resolution | **Recommended** for most cases |
| **5-8 Mpc/h** | Heavy smoothing | Very noisy data, qualitative plots |
| **>8 Mpc/h** | Over-smoothing | Risk losing real features |

**Rule of thumb**: Set `kernel_width ≈ Planck FWHM at your redshift`

For BOSS (z ~ 0.5):
- Planck FWHM = 8 arcmin → **~3.5 Mpc/h**
- Recommended: `--kernel-width 3.5`

## Example: BOSS Galaxy Stack

```bash
# Full pipeline with kernel smoothing
python stack_single.py \
    --dataset BOSS \
    --region North \
    --catalog-type galaxy \
    --extract-radial-kde \
    --kernel-width 3.5 \
    --n-radial-bins 500 \
    --output-dir analysis/boss/results
```

**Output:**
- `kappa_single_galaxy_BOSS_North.csv` — 2D stacked map (50×50)
- `error_single_BOSS_North.csv` — Jackknife errors (50×50)
- `radial_kde_galaxy_BOSS_North.csv` — **Smooth radial profile (500 points)**
- `radial_binned_galaxy_BOSS_North.csv` — Binned radial profile (50 points)

## Plotting Example

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load KDE and binned profiles
kde = pd.read_csv('radial_kde_galaxy_BOSS_North.csv')
binned = pd.read_csv('radial_binned_galaxy_BOSS_North.csv')

plt.figure(figsize=(10, 6))

# Smooth KDE profile
plt.plot(kde['r_hmpc'], kde['kappa'], 'b-', linewidth=2, label='KDE (σ=3.5 Mpc/h)')
if 'kappa_err' in kde.columns:
    plt.fill_between(kde['r_hmpc'],
                     kde['kappa'] - kde['kappa_err'],
                     kde['kappa'] + kde['kappa_err'],
                     alpha=0.3, color='blue')

# Binned profile (for comparison)
plt.plot(binned['r_hmpc'], binned['kappa'], 'o',
         alpha=0.5, markersize=4, label='Binned (50 bins)')

plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('r [h⁻¹ Mpc]')
plt.ylabel('κ(r)')
plt.title('BOSS Galaxy Lensing Signal')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('lensing_profile_kde.png', dpi=150)
```

## Trade-offs

| Method | Pros | Cons |
|--------|------|------|
| **Hard binning** | • Fast<br>• Simple<br>• Traditional | • Jagged<br>• Noisy<br>• Limited by discrete grid |
| **Gaussian KDE** | • Smooth<br>• Continuous<br>• Respects Planck resolution<br>• Better for fitting | • Slightly slower<br>• Requires kernel width choice |

## Technical Details

### Error Propagation

For KDE with weights `wᵢ` and errors `σᵢ`:

```
σ²_smooth(r) = Σᵢ (wᵢ × σᵢ)² / (Σᵢ wᵢ)²
```

This correctly accounts for the Gaussian kernel weighting.

### Computational Complexity

- **Hard binning**: O(N_galaxies × GRID_SIZE²)
- **KDE extraction**: O(GRID_SIZE² × n_radial_bins) — runs AFTER stacking

KDE is a post-processing step, so the stacking time is unchanged.

### Backward Compatibility

The new kernel smoothing is **fully optional**:
- Default behavior (no flags) → unchanged, produces same outputs as before
- Add `--extract-radial-kde` → creates additional KDE profile files
- All existing scripts/notebooks continue to work unchanged

## References

- **Planck 2018 Lensing**: Planck Collaboration (2020), A&A 641, A8
- **Kernel Density Estimation**: Silverman (1986), "Density Estimation"
- **Galaxy-galaxy lensing**: Mandelbaum et al. (2006), MNRAS 368, 715

## See Also

- `stack_single.py` — Main stacking script (updated with KDE)
- `geometry.py` — Contains `symmetrize_map()` for radial averaging
- `constants.py` — Grid parameters (`GRID_SIZE`, `BOX_SIZE_HMPC`)
