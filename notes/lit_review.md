# Literature Review

Summaries of key papers for the filament lensing project, with comparisons to our methodology.

**Papers reviewed:**
1. [Kondo et al. (2019)](#kondo-et-al-2019) -- Filament weak lensing with BOSS+HSC
2. [He et al. (2018)](#he-et-al-2018) -- First CMB lensing x filament detection (density ridges)
3. [de Graaff et al. (2019)](#de-graaff-et-al-2019) -- Missing baryons via tSZ + CMB lensing from filaments
4. [Tanimura et al. (2020)](#tanimura-et-al-2020) -- Filament density & temperature via tSZ + CMB lensing (DisPerSE)

---

## Kondo et al. (2019) -- "Weak Lensing Measurement of Filamentary Structure with the SDSS BOSS and Subaru Hyper Suprime-Cam Data"

**Citation:** Kondo, Miyatake, Shirasaki, Sugiyama & Nishizawa (2019), MNRAS 000, 1-10. [arXiv:1905.08991](https://arxiv.org/abs/1905.08991)
**Local PDF:** `notes/references/1905.08991v1.pdf`

### Overview

This paper reports a **3.9-sigma detection** of filament weak lensing between BOSS CMASS galaxy pairs at z ~ 0.55, using the Subaru Hyper Suprime-Cam (HSC) first-year galaxy shape catalogue. Despite a small overlap area of only ~140 deg^2 between BOSS and HSC, the deep HSC imaging (source galaxy density ~21.8 arcmin^-2) enables a competitive measurement. This is the highest S/N weak lensing detection of filaments between *galaxy-scale* halos at this redshift range.

### Data

**Lens sample (foreground):**
- BOSS CMASS LSS catalogue (spectroscopic redshifts)
- Redshift cut: 0.43 < z < 0.70 (effectively; formal cut 0.14 < z < 0.88 then pair selection restricts range)
- Spatial cut to HSC first-year full-depth full-color footprints: 14,422 CMASS galaxies remain
- Weighting: w_l = w_seeing * w_star * (w_noz + w_cp - 1) -- same BOSS CMASS weight formula we use
- No FKP weight (not needed for lensing, only for clustering measurement)

**Source sample (background):**
- HSC first-year galaxy shape catalogue
- i-band imaging, mean seeing FWHM 0.58", 5-sigma depth i ~ 26
- Magnitude cut i < 24.5 for robust shape calibration
- Weighted source number density: 21.8 arcmin^-2 (much higher than SDSS ~0.5 arcmin^-2)
- Photo-z from MLZ (machine learning self-organising map)
- Shape measurement: re-Gaussianization PSF correction method, calibrated with GalSim image simulations

**Cosmology:** Planck 2015, Omega_m = 0.307, flat LCDM

### Pair Selection

| Parameter | Value |
|-----------|-------|
| Line-of-sight cut | \|Pi\| < 6 h^-1 Mpc |
| Transverse separation | 6 < R < 14 h^-1 Mpc |
| Number of pairs | 70,210 |
| Number of galaxies | 14,422 |
| Mean pair redshift | z ~ 0.55 |

Key points:
- Follows the Clampitt et al. (2016) pair selection procedure
- A galaxy can appear in multiple pairs (pairs > galaxies)
- Pair position/redshift defined as the average of the two galaxy positions/redshifts
- Pair weight = product of individual galaxy weights

**Random pairs:** Created from BOSS random catalogue with the same cuts. Random catalogue divided into subsamples with the same number of points as CMASS galaxies, yielding 96 random pair samples.

**Null test pairs (separated):** Selected with 40 < |Pi| < 60 h^-1 Mpc (same R range). These should have no filament signal since galaxies are far apart along LOS.

### Measurement Method

Uses the **Clampitt et al. (2016) estimator** -- a 4-point cancellation scheme designed to isolate filament signal by removing halo contributions:

1. **Stretch and align:** Each pair is rescaled so the two halos are at fixed positions (separation normalized to unity). The coordinate system is centred on the left halo.

2. **4-point cancellation:** For a point "p1" at position (x, y) between the pair:
   - Adding shear at "p2" (90-degree rotation from p1 w.r.t. left halo) cancels left halo contribution
   - Adding "p3" and "p4" (analogous rotations w.r.t. right halo) cancels right halo contribution
   - Net effect: all spherically symmetric halo signal cancels out, leaving only filament signal

3. **Translational symmetry:** Sum over 10 bins along the pair axis (x-direction)

4. **Line symmetry:** Also sum positions mirrored across the connecting line

5. **Random subtraction:** Subtract signal from random pairs to correct for observational systematics (imperfect PSF correction, shape noise at survey edges)

The estimator output is a 1D function of y (distance from the pair connecting line), with 5 bins (y/R_pair from 0 to 0.5).

The observable is the **excess surface mass density** Delta Sigma (in M_sun h / pc^2), measured from galaxy shapes (tangential shear). This is distinct from our approach using CMB lensing convergence kappa.

### Mock Catalogues and Covariance

- 108 mock realisations from full-sky ray-tracing simulations (Takahashi et al. 2017)
- Mock CMASS galaxies populated via HOD fitting to CMASS clustering + abundance
- HOD fitting done in 5 redshift bins: [0.43,0.44], [0.44,0.51], [0.51,0.57], [0.57,0.64], [0.64,0.70]
- Covariance estimated from scatter across 108 realisations
- **Key finding on covariance composition:**
  - Shape noise contribution: ~40% of total covariance
  - Intrinsic scatter of filament properties + LOS LSS fluctuations: ~60% of total
  - "Our filament lensing measurement is no longer limited by the intrinsic shape noise of source galaxies"
  - This means increasing source density won't help much; need larger survey area

### Results

- **Detection significance:** 3.9 sigma against null hypothesis
- **Consistency with theory:** chi^2 = 8.1 (dof = 5, p = 0.15) -- consistent with mock prediction
- Data points are consistently slightly below theory but chi^2 is still acceptable due to strong bin-to-bin correlation (correlation matrix entries > 0.77)
- 2D convergence map shows diffuse filament between circularly symmetric halos, supporting **"thick" filament model** (from N-body simulations, not a thin string of small halos)

**Null tests (all pass):**

| Test | chi^2/dof | p-value |
|------|-----------|---------|
| Cross-mode | 7.27/5 | 0.20 |
| Separated pair + mode | 2.34/5 | 0.80 |
| Separated pair x mode | 6.05/5 | 0.30 |

### Comparison with Other Studies (from the paper)

| Study | Lens sample | z_lens | Source | Separation cuts | S/N | Notes |
|-------|-------------|--------|--------|-----------------|-----|-------|
| Clampitt+2016 | SDSS LRGs | ~0.25 | SDSS (~0.5/arcmin^2) | 6-14 h^-1 Mpc, Pi<6 | 4.5 sigma | ~8000 deg^2, higher amplitude (more massive halos) |
| Epps & Hudson 2017 | BOSS LOWZ+CMASS | ~0.42 | CFHTLenS | 6-10 h^-1 Mpc, Dz<0.002 | 5 sigma | Used "control" convergence map subtraction |
| He+2018 | CMASS filament catalogue | ~0.55 | Planck CMB kappa | Different method (density ridges) | 5 sigma | Filaments where galaxies reside, not between pairs |
| de Graaff+2019 | CMASS pairs | ~0.55 | Planck CMB kappa | 6-14 h^-1 Mpc, Pi<5 | 1.9 sigma (lensing) | Focus was tSZ; entire CMASS sample, larger area but CMB lensing less sensitive |
| **This paper** | CMASS pairs | ~0.55 | HSC (~21.8/arcmin^2) | 6-14 h^-1 Mpc, Pi<6 | **3.9 sigma** | Small area (140 deg^2) but deep imaging |

Key insight: de Graaff+2019 used the full CMASS sample (much larger area) with Planck CMB lensing and only got 1.9 sigma on the lensing signal. This paper achieves 3.9 sigma with only 140 deg^2 because HSC galaxy shapes have much better lensing sensitivity at z ~ 0.55 than Planck CMB lensing. The galaxy lensing kernel peaks near the lens redshift, while the CMB kernel peaks at z ~ 2.

### Contrast with Our Approach

| Aspect | Kondo+2019 | Our project |
|--------|------------|-------------|
| **Lensing source** | Galaxy shapes (HSC, z_s ~ 1) | CMB lensing convergence (Planck) |
| **Observable** | Excess surface mass density Delta Sigma (from tangential shear) | Convergence kappa (directly from Planck alm) |
| **Halo removal** | Clampitt+2016 4-point cancellation estimator (analytically removes spherical halo contributions from shear field) | Control map subtraction: stack single galaxies, shift by +/-separation/2, subtract from pair-stacked map |
| **Pair selection** | |Pi| < 6 h^-1 Mpc, 6 < R < 14 h^-1 Mpc | r_par < 20 Mpc/h (currently; 5 under discussion), various r_perp |
| **Coordinate stretching** | Yes -- all pairs rescaled to unit separation | No -- fixed physical grid (101x101, 100 Mpc/h box), pairs at different separations produce galaxies at different pixel locations |
| **Symmetrization** | Line symmetry (mirror across pair axis) + translational symmetry (sum along pair axis) -- both built into the estimator | Reflection symmetry across Y-axis (perpendicular to pair axis) |
| **Covariance** | From 108 mock realisations (ray-tracing + HOD) | Jackknife (for single stacks); no formal covariance for pair stacks yet |
| **Random pair handling** | Subtract random pair signal to remove observational systematics | Subtract random-point stacked kappa map |
| **Weighting** | Per lens-source pair weight including Sigma_cr, shape noise, CMASS weights | Equal weight per pair (pair stacking); CMASS weights for single galaxy stacking |
| **Galaxy sample** | CMASS only (14,422 galaxies in HSC overlap) | CMASS North+South (full BOSS footprint, ~500k+ galaxies) |
| **Cosmology** | Planck 2015 (Omega_m = 0.307) | Planck 2018 (via astropy) |

### Key Differences and Implications

1. **Lensing sensitivity at z ~ 0.55:** Galaxy shape lensing (HSC) is much more sensitive than CMB lensing (Planck) for filaments at z ~ 0.55. Kondo+2019 achieves 3.9 sigma with only 140 deg^2, while de Graaff+2019 achieved 1.9 sigma lensing significance with the full CMASS footprint (~10,000 deg^2) using Planck. Our CMB lensing approach trades per-pair sensitivity for full sky coverage. We need many more pairs to compensate.

2. **Halo subtraction method:** The Clampitt estimator is more elegant -- it analytically cancels halo contributions using the 4-point symmetry of tangential shear around spherical halos. Our control-map-subtraction approach is simpler but noisier: errors in the control map (from limited random points, radial profile assumptions) propagate into the filament signal.

3. **Coordinate stretching vs. fixed grid:** Kondo+2019 rescales all pairs to unit separation before stacking, which coherently adds filament signal regardless of pair separation. Our fixed-grid approach means pairs of different separations place their galaxies at different locations in the grid, which is fine if we bin by separation, but we should be aware that within a separation bin (e.g. 18-22 Mpc/h), the filament signal is slightly smeared.

4. **r_par cut:** Kondo+2019 uses |Pi| < 6 h^-1 Mpc, consistent with the literature consensus (Clampitt, Epps & Hudson, de Graaff). This further supports reducing our r_par from 20 to ~5-6 Mpc/h.

5. **Covariance:** Kondo+2019's mock-based covariance reveals that even with deep HSC data, shape noise is only 40% of the error budget. The dominant source is intrinsic filament property scatter + LOS structure fluctuations. For our CMB lensing case, the noise budget is likely even more dominated by Planck map noise (which is large-scale and correlated), making a robust covariance estimate important. Our current jackknife approach for single stacks should be extended to pair stacks.

6. **"Thick" vs "thin" filament:** Kondo+2019's convergence map supports the "thick" filament model (diffuse dark matter distribution), not a "thin" filament (string of small halos). This is the same physical picture our stacking analysis should recover.

### References from this paper relevant to us

- Clampitt et al. (2016) -- foundational estimator method (already in our reference list)
- Epps & Hudson (2017) -- control map approach (already in our reference list)
- de Graaff et al. (2019) -- closest to our approach (CMB lensing + CMASS pairs, already in our reference list)
- He et al. (2018) -- alternative CMB lensing approach using density ridges (different methodology, worth noting)
- Takahashi et al. (2017) -- ray-tracing simulations (useful if we ever want mock-based covariance)

---

## He et al. (2018) -- "The detection of the imprint of filaments on Cosmic Microwave Background lensing"

**Citation:** He, Alam, Ferraro, Chen & Ho (2018), Nature Astronomy 2, 401. [arXiv:1709.02543](https://arxiv.org/abs/1709.02543)
**Local PDF:** `notes/references/1709.02543v2.pdf`

### Overview

This paper reports the **first detection of CMB lensing by filaments**, at **5.0-5.2 sigma** significance. Unlike the pair-stacking approach used by Clampitt, Epps & Hudson, Kondo, and de Graaff, this paper takes a fundamentally different approach: it identifies filaments as **density ridges** in the CMASS galaxy distribution using a topological/statistical algorithm, constructs a 2D filament intensity map, and then **cross-correlates** this map with the Planck CMB lensing convergence map in **Fourier space** (angular power spectrum C_l).

The key physics result is a measurement of the **filament bias** b_f ~ 1.5, quantifying how filaments trace the underlying matter distribution on large scales.

### Data

**Filament catalogue:**
- Derived from the **Cosmic Web Reconstruction** filament catalogue (Chen et al. 2016)
- Based on BOSS DR12 CMASS spectroscopic galaxies
- Redshift range: z = 0.450 to z = 0.700 (effective redshift z_eff = 0.56)
- Universe partitioned into redshift slices of Delta_z = 0.005

**CMB lensing:**
- Planck 2015 lensing convergence map (SMICA foreground cleaning)
- Downgraded to N_side = 512 to match filament map resolution

**Cosmology:** Planck 2013, Omega_m = 0.315, h = 0.673, sigma_8 = 0.829

### Methodology -- Fundamentally Different from Pair-Stacking

#### 1. Filament Finding (Density Ridge Algorithm)

The filament finder works in 2D angular slices:

1. Partition galaxies into redshift bins (Delta_z = 0.005, yielding ~130 bins from z=0.05 to z=0.70)
2. For each bin:
   - Project galaxies onto 2D angular space
   - Reconstruct the galaxy probability density field using kernel density estimation (KDE)
   - Threshold: remove galaxies in low-density regions (stabilises the algorithm)
   - Apply the **ridge-finding algorithm** (Chen et al. 2015) to identify density ridges = filaments
   - Apply galaxy mask to remove filaments outside observed footprint

The Delta_z = 0.005 bin width is chosen as a **bias-variance tradeoff**: too small means too few galaxies (unstable), too large means structures get washed out in projection.

The filament uncertainty is computed via bootstrap: resample galaxies 1000 times, re-run the filament finder each time, compute RMS distance of filament points to their bootstrap counterparts.

#### 2. Filament Intensity Map

Each sky position is assigned a filament intensity:

I(n_hat, z) = (1 / sqrt(2*pi) * sigma_f) * exp(-||n_hat - n_f||^2 / (2 * sigma_f^2))

where n_f is the nearest filament point and sigma_f is the filament uncertainty at that point. This creates a smooth intensity map from the discrete filament catalogue.

The filament intensity overdensity is then:

delta_f(n_hat) = (integral of I over z - I_bar) / I_bar

This is constructed as a HEALPix map at N_side = 512.

#### 3. Cross-Correlation in Fourier Space

The measurement is the **cross angular power spectrum** C_l^{f,kappa}:

C_l^{f,kappa} = (1 / (2l+1) * f_sky) * sum_m (delta_f)_lm * (kappa)_lm*

measured using a pseudo-C_l estimator. This operates entirely in harmonic space -- no real-space stacking, stretching, or grid sampling.

#### 4. Phenomenological Model

The filament-convergence cross-power spectrum is modelled as:

C_l^{f,kappa} = integral over z of [W(z) * phi_f(z) / chi^2 * (1+z) * P_mf(l/chi, z)] dz

where:
- W(z) is the CMB lensing kernel
- phi_f(z) is the filament intensity redshift distribution
- P_mf is the matter-filament power spectrum, modelled as b_f * P_mm with exponential suppression below the filament scale

The suppression models the fact that filaments smooth out density fluctuations below their characteristic length (along the filament) and spacing (perpendicular to the filament).

Two smoothing models are tested:
- **Model 1:** filament length = overall smoothing scale (both parallel and perpendicular)
- **Model 2:** filament length = parallel smoothing; perpendicular scale is a free parameter (best fit: 0.65x filament length)

### Results

| Model | b_f (filament bias) | S/N | chi^2_fit | dof |
|-------|---------------------|-----|-----------|-----|
| Model 1 | 1.68 +/- 0.334 | 5.0 | 25.77 | 15 |
| Model 2 | 1.47 +/- 0.28 | 5.2 | 24.39 | 14 |

Key findings:
- **Filament bias ~ 1.5** -- filaments trace matter with slightly higher bias than unity, but significantly less than CMASS galaxy bias (~2)
- **5 sigma detection** via SNR = sqrt(chi^2_null - chi^2_fit)
- Covariance from jackknife resampling (77 equal-area regions)

**Null tests:**
- Rotated CMB kappa maps (90, 135, 180 degrees): all consistent with null (chi^2/dof ~ 0.75-1.04)
- Masking redMaPPer clusters: less than 4% change in C_l, confirming signal is not dominated by cluster nodes

**Filament-galaxy correlation:**
- Filaments and galaxies are highly correlated on large scales (both trace LSS)
- Correlation coefficient decreases on small scales -- filament map contains information beyond galaxy density field
- This suggests filaments provide **additional cosmological information** beyond standard galaxy clustering

### Filament Properties

- **Filament lengths:** Mean increases with redshift (from ~60 to ~130 Mpc at z=0.45-0.65)
  - This is partly physical (structure growth) and partly selection: CMASS number density decreases with z, so fewer short filaments are detected
- **Length distribution** is non-Gaussian (large difference between mean and median)
- Effective redshift of filament sample: z_eff = 0.56

### Contrast with Our Approach

| Aspect | He+2018 | Our project |
|--------|---------|-------------|
| **Filament identification** | Density ridge algorithm (topological; finds filaments as 1D ridges in the galaxy density field) | Galaxy pair selection (geometric; assumes a filament connects each pair) |
| **What is a "filament"** | Continuous curves tracing density ridges -- filaments where galaxies *reside* | The space *between* two galaxies at a given separation -- filaments connecting pairs |
| **CMB lensing map** | Planck 2015 kappa, N_side=512 | Planck kappa reconstructed from alm, N_side=2048 |
| **Measurement space** | Fourier/harmonic space (angular power spectrum C_l) | Real space (pixel-by-pixel stacking on a 2D grid) |
| **Signal extraction** | Cross-correlation of filament intensity overdensity map with kappa map | Stack kappa at pair positions, subtract control map |
| **Halo contribution** | Masking redMaPPer clusters changes signal by <4%; filament ridges naturally avoid cluster centres | Must explicitly subtract (control map or Clampitt 4-point estimator) |
| **Galaxy sample** | CMASS (same as ours) | CMASS (same) |
| **Redshift range** | z = 0.45-0.70 (z_eff = 0.56) | z = 0.4-0.7 (similar) |
| **Detection significance** | 5.0-5.2 sigma | TBD |
| **Key output** | Filament bias b_f ~ 1.5 | 2D kappa maps showing filament excess between pairs |
| **Covariance** | Jackknife (77 regions) | Jackknife (single stacks); TBD for pairs |

### Key Differences and Implications

1. **Fundamentally different filament definition.** He+2018 finds filaments as density ridges in the galaxy field -- these are filaments where galaxies live. Our pair-stacking approach looks for matter *between* galaxy pairs. As Kondo+2019 notes: "their filament sample is expected to be more massive than ours" since density ridges trace the overdense spines of the cosmic web, while pair-stacking probes the diffuse matter bridges between halos.

2. **Fourier vs. real space.** He+2018 works entirely in harmonic space (C_l), which naturally handles survey geometry, noise properties, and scale-dependent effects. Our real-space stacking is more intuitive and directly produces 2D images, but is harder to model analytically and more susceptible to systematic effects at the map level.

3. **No pair selection needed.** The density ridge approach sidesteps the r_par / r_perp pair selection entirely -- there are no "pairs" in their analysis. This avoids the r_par question that we face, but at the cost of not having a clean geometric picture of individual filaments.

4. **Filament bias vs. filament profile.** He+2018 measures a single number (b_f ~ 1.5) describing how filaments trace matter on large scales. Our approach can potentially measure the **transverse profile** of filaments (kappa as a function of distance from the pair axis), which is more detailed but requires higher S/N.

5. **CMB lensing resolution.** They use N_side=512 (pixel size ~7 arcmin), much lower than our N_side=2048 (~1.7 arcmin). For a cross-correlation measurement on large scales this is fine, but our real-space stacking benefits from the higher resolution to resolve the filament profile.

6. **Same CMASS galaxy sample and redshift range.** The underlying galaxy data is the same, making this a complementary analysis to ours. If we detect a filament signal with our pair-stacking + CMB lensing approach, it would be an independent confirmation using a different methodology on the same underlying cosmic web.

### References from this paper relevant to us

- Chen et al. (2015, 2016) -- Cosmic Web Reconstruction filament catalogue and algorithm
- de Graaff et al. (2019) -- tSZ filament detection (already in our reference list)
- Clampitt et al. (2016), Epps & Hudson (2017) -- pair-stacking predecessors (already in our list)
- Pullen et al. (2016) -- CMB lensing x galaxy velocity cross-correlation methodology

---

## de Graaff et al. (2019) -- "Probing the missing baryons with the Sunyaev-Zel'dovich effect from filaments"

**Citation:** de Graaff, Cai, Heymans & Peacock (2019), A&A 624, A48. [arXiv:1709.10378](https://arxiv.org/abs/1709.10378)
**Local PDF:** `notes/references/1709.10378v3.pdf`

### Overview

This paper is **the closest predecessor to our project**. It uses the same galaxy sample (BOSS CMASS), the same CMB lensing map (Planck), and the same pair-stacking approach. The primary focus is detecting the **thermal Sunyaev-Zel'dovich (tSZ)** signal from filamentary gas (the "missing baryons" / WHIM), but they also measure the **CMB lensing convergence** signal to estimate the total matter density in the filament, enabling a joint constraint on gas density and temperature.

Key results:
- **tSZ filament detection: 2.9 sigma** (between pairs); 3.8 sigma including signal beyond pairs
- **CMB lensing filament detection: 1.9 sigma** (between pairs); 3.1 sigma extended
- Filament gas density: rho_b = (5.5 +/- 2.9) * rho_bar_b
- Filament temperature: T = (2.7 +/- 1.7) x 10^6 K
- Filament baryons account for 11 +/- 7% of total baryon content (28 +/- 12% extended)

### Data

**Galaxy pairs (lens sample):**
- BOSS DR12 CMASS, North + South (full catalogue: ~850,000 galaxies)
- Redshift range: 0.43 < z < 0.75
- Pair selection: transverse separation 6-14 h^-1 Mpc, LOS separation < 5 h^-1 Mpc
- **1,002,334 physical pairs** (mean separation 10.5 h^-1 Mpc, mean angular sep 26.5 arcmin)
- Control sample ("non-physical pairs"): same projected separation but 40-200 h^-1 Mpc LOS separation -> 13,622,456 pairs

**Planck maps:**
- **tSZ:** MILCA Compton y-parameter map (also verified with NILC), FWHM 10 arcmin beam, N_side = 2048, 40% Galactic mask
- **CMB lensing:** Planck kappa map, convolved with Gaussian FWHM 10 arcmin to match y-map resolution, harmonic limit l < 2048

**Cosmology:** flat LCDM, H_0 = 100h km/s/Mpc, h = 0.68, Omega_m = 0.31, Omega_b = 0.049

### Pair Selection Details

| Parameter | Value |
|-----------|-------|
| Transverse separation | 6-14 h^-1 Mpc |
| LOS separation | < 5 h^-1 Mpc |
| Number of physical pairs | 1,002,334 |
| Mean pair separation | 10.5 h^-1 Mpc |
| Control pairs | 13,622,456 (LOS: 40-200 h^-1 Mpc) |

**On the LOS cut (footnote 1, p2):** "Our criterion on line-of-sight separation implicitly assumes that the total redshift differences are purely cosmological; but in reality, peculiar velocities can have an effect -- as can be seen from e.g. Fig. 11 of Jenkins et al. (1998). For filaments of our length that are near to the plane of the sky, the pairwise dispersion in radial velocity is close to 500 km/s, and this has a convolving effect. Thus, some of our pairs will have a true radial separation slightly larger than our limit of 5 h^-1 Mpc, and some pairs with a smaller radial separation will be rejected because the apparent separation in redshift space is above our limit. But our selection still picks pairs that are nearly transverse to the line of sight."

### Methodology

#### 1. Stacking Procedure

For each galaxy pair:
1. **Rotate** the map so the pair aligns with the equator, centred at the origin
2. **Rescale** according to the angular pair separation (all pairs overlap at the same normalised positions)
3. **Project** onto a 2D rectangular grid (nearest-neighbour interpolation)
4. **Mirror** across the vertical axis (reflection symmetry)
5. Masked HEALPix pixels get weight = 0

Working resolution: N_side = 1024 (~3.4 arcmin pixel) for computational efficiency.

#### 2. Halo Subtraction (Isotropic Profile Fitting)

This is the key methodological step and differs from our control-map approach:

1. Extract the **radial profile along the vertical direction** (perpendicular to pair axis), using only data within a 60-degree cone to avoid the filament region
2. Fit a model: F(r_1) = f(r_1) + f(r_2), where f(r) is a 4th-order polynomial x exponential, and r_2 = sqrt(r_1^2 + r_12^2) accounts for the second halo
3. Generate a full 2D model from this isotropic profile
4. **Subtract** the 2D model from the stacked map -> residual = filament signal

Robustness check: introduced a free scaling parameter s for the halo profile width. Best fit gives s = 0.94 (slightly narrower), and the filament signal changes by only ~10%. The filament signal is ~10% of the overlapping isotropic halo signals.

#### 3. Significance Estimation

**Naive estimate:** Treating all pairs as independent gives 5.1 sigma. But this is wrong -- with ~10^6 pairs over 9376 deg^2, and each filament rectangle covering ~0.12 deg^2, each SZ pixel appears in the stack ~13 times.

**Jackknife method:** Split pairs into N_sub = 250 sub-samples (~40 deg^2 each, ~4000 pairs each). Repeat full stacking + halo fitting for each jackknife realisation. Construct covariance from jackknife scatter with Hartlap correction for inverse covariance bias.

Result: **2.9 sigma** for the filament region between the pair (n=12 bins). Including signal beyond the pair: **3.8 sigma**.

### Contamination Assessment

| Source | Contribution | Method |
|--------|-------------|--------|
| Uncorrelated sources (foreground, dust) | Consistent with zero | Stack non-physical pairs (13.6M), draw 500 subsamples |
| Bound gas in correlated haloes | ~20% of filament signal | (a) Stack CMASS galaxy number density map, find filament is 1-2% of peak vs 10% in SZ; (b) Millennium simulation with Y-M relation (slope 5/3), HOD galaxy population |
| CIB leakage into y-map | 2-3 orders of magnitude below signal | Leakage coefficient epsilon_CIB ~ 10^-7 to 10^-6; stacking Planck 857 GHz map gives ~10^-5 K, multiplied by epsilon gives negligible contamination |

Key argument against CIB contamination: the filament-to-halo contrast ratio is **larger** in SZ than in the galaxy density map, so the excess SZ must come from diffuse gas, not from resolved haloes. Since CIB is associated with star formation in haloes, it cannot produce the diffuse excess.

### CMB Lensing Results

They apply the **identical stacking and halo subtraction procedure** to the Planck kappa map (convolved to 10 arcmin FWHM to match the y-map).

- Mean convergence in filament region: kappa_bar = (0.58 +/- 0.31) x 10^-3 -> **1.9 sigma**
- Extended region: kappa_bar = (0.33 +/- 0.11) x 10^-3 -> **3.1 sigma**
- The 2D stacked kappa maps (their Figs 6 & 7) show the same qualitative bridge feature as the tSZ maps

#### Filament density from lensing:
Using a cylindrical Gaussian model (Appendix A):
- Central matter overdensity: delta(z) = (5.5 +/- 2.9) * rho_bar(z)
- Assuming universal baryon fraction in filament -> gas density n_0(z) = (5.5 +/- 2.9) * n_bar_e(z)
- Combined with tSZ measurement -> T_e = (2.7 +/- 1.7) x 10^6 K

#### Baryon budget:
- Filament between pairs: 11 +/- 7% of total baryons
- Extended filament (42 h^-1 Mpc baseline): 28 +/- 12% of total baryons

### Comparison with Tanimura et al. (2019)

Tanimura et al. performed a similar tSZ stacking with SDSS DR12 LRG galaxies at z < 0.4 (262,864 pairs). They found y ~ 1 x 10^-8 at 5.3 sigma. de Graaff finds y ~ 0.6 x 10^-8 at 2.9 sigma with the higher-z CMASS sample. The two studies are independent and complementary in redshift.

### Contrast with Our Approach

| Aspect | de Graaff+2019 | Our project |
|--------|---------------|-------------|
| **Primary observable** | tSZ Compton y-parameter (gas pressure) + CMB lensing kappa (total matter) | CMB lensing kappa only |
| **Planck map** | MILCA y-map (FWHM 10') + kappa map (convolved to 10') | kappa reconstructed from alm at N_side=2048 (~1.7') |
| **Galaxy sample** | BOSS CMASS North+South (full, ~850k galaxies) | BOSS CMASS North+South (same) |
| **Pair selection** | 6-14 h^-1 Mpc transverse, < 5 h^-1 Mpc LOS | Various separations (5, 10, 20), r_par = 20 (under review) |
| **Number of pairs** | 1,002,334 | TBD (634k for South 20 Mpc/h alone) |
| **Coordinate rescaling** | Yes -- all pairs rescaled to same normalised separation | No -- fixed physical grid (101x101, 100 Mpc/h box) |
| **Halo subtraction** | Fit isotropic radial profile from 60-degree vertical cone, subtract 2D model | Subtract control map (shifted single-galaxy stack) |
| **Reflection symmetry** | Mirror across vertical axis (perpendicular to pair axis) | Mirror across Y-axis (same) |
| **Resolution** | N_side=1024 (stacking), kappa smoothed to 10' FWHM | N_side=2048 (no additional smoothing beyond Planck alm) |
| **Significance** | 1.9 sigma (kappa, between pairs), 3.1 sigma (extended) | TBD |
| **Covariance** | Jackknife (N_sub=250, ~40 deg^2 regions, Hartlap correction) | Jackknife for single stacks; TBD for pairs |
| **Control sample** | Non-physical pairs (40-200 h^-1 Mpc LOS, 13.6M pairs) | Random-point stack subtraction |
| **Cosmology** | h=0.68, Omega_m=0.31, Omega_b=0.049 | Planck18 via astropy |

### Key Differences and Implications

1. **Our r_par is too large.** They use < 5 h^-1 Mpc (with explicit justification in footnote 1). Our current 20 Mpc/h would include many unconnected projections. This is the fourth paper confirming the ~5 Mpc/h consensus.

2. **They rescale pairs; we don't.** Their stacking rescales all pairs to unit separation before stacking, coherently aligning the filament signal regardless of actual separation (6-14 Mpc/h range). Our fixed-grid approach means pairs of different separations within a bin (e.g. 18-22 Mpc/h) produce slightly different filament positions, smearing the signal.

3. **Halo subtraction approaches differ.** Their isotropic profile fitting from the vertical cone is elegant: it uses the data itself to determine the halo contribution, without needing a separate single-galaxy stack. Our control-map approach requires a good single-galaxy stack and assumes the single-galaxy halo profile is representative. However, their method assumes perfect spherical symmetry of the halos, which they verify is robust to ~10% profile width changes.

4. **Their kappa detection is only 1.9 sigma.** This is sobering for our project. They used the full CMASS sample (~1M pairs) and still got marginal kappa significance. Our advantage might be:
   - We use the full-resolution kappa map (N_side=2048) without smoothing to 10 arcmin (they smoothed for consistency with the y-map)
   - We probe different separations (5, 10, 20 Mpc/h) independently rather than combining 6-14 Mpc/h
   - But we have the same fundamental CMB lensing noise limitation

5. **The non-physical pair control sample is a great idea.** Pairs with 40-200 h^-1 Mpc LOS separation have the same projected geometry but no filament. This directly tests for projection effects and uncorrelated systematics. We could construct an analogous control sample.

6. **Baryon budget context.** Their combined tSZ + lensing analysis places the filament gas at rho_b ~ 5.5 * rho_bar at T ~ 2.7 x 10^6 K, consistent with WHIM predictions. Our kappa-only analysis cannot constrain temperature, but can independently confirm the matter overdensity.

7. **Extended filament signal.** They find that the filament extends *beyond* the galaxy pair (significance increases from 2.9 to 3.8 sigma for tSZ, 1.9 to 3.1 for kappa). This is physically expected -- filaments don't terminate at the halo virial radius. Our 100 Mpc/h grid box is large enough to capture this.

### References from this paper relevant to us

- Clampitt et al. (2016), Epps & Hudson (2017) -- pair selection predecessors (already in our list)
- He et al. (2018) -- CMB lensing x filaments via density ridges (reviewed above)
- Tanimura et al. (2019) -- independent tSZ filament detection with LRGs at z < 0.4
- Colberg et al. (2005) -- simulations predicting filaments connect pairs up to 20 h^-1 Mpc
- Jenkins et al. (1998) -- pairwise velocity dispersion (cited for the RSD footnote)
- Hartlap et al. (2007) -- correction factor for inverse covariance matrix bias in jackknife

---

## Tanimura et al. (2020) -- "Density and temperature of cosmic-web filaments on scales of tens of megaparsecs"

**Citation:** Tanimura, Aghanim, Bonjean, Malavasi & Douspis (2020), A&A 637, A41. [arXiv:1911.09706](https://arxiv.org/abs/1911.09706)
**Local PDF:** `notes/references/1911.09706v2.pdf`

### Overview

This paper measures the density and temperature of **large-scale cosmic-web filaments** (lengths 30-100 Mpc) by stacking Planck tSZ y-maps and CMB lensing kappa maps at filament positions identified by the **DisPerSE** algorithm applied to SDSS DR12 galaxy data. Unlike the pair-stacking approach (de Graaff, Kondo, Clampitt), this paper uses algorithmically-detected filament spines as stacking centres.

Key results:
- **tSZ detection: 4.4 sigma** (filament gas pressure)
- **CMB lensing detection: 8.1 sigma** (filament total matter)
- Joint MCMC fit: central overdensity delta_c = 19.0 +/- 4.3, core radius r_c = 1.5 +/- 0.3 Mpc, temperature T_e = (1.4 +/- 0.4) x 10^6 K
- Baryon fraction in filaments: 0.080 +/- 0.025 of total Omega_b

### Data

**Filament catalogue:**
- Galaxy sample: SDSS DR12 LOWZ + CMASS combined (galaxies at 0.2 < z < 0.6)
- Filament finder: **DisPerSE** (Discrete Persistent Structures Extractor, Sousbie 2011)
  - Identifies filaments as field lines of the density gradient connecting critical points (maxima, saddle points)
  - Applied to projected 2D galaxy positions in redshift slices
  - Persistence threshold: 3 sigma (removes noise-dominated features)
- Filament length cuts: 30 < L < 100 Mpc (excludes short filaments that may be bound structures, and very long ones that may be spurious chains)
- **24,544 filaments** after all cuts
- Each filament represented by its spine (series of connected segments)

**Planck tSZ map:**
- MILCA Compton y-parameter map, N_side = 2048
- 40% Galactic mask applied
- Point source and SZ cluster masking (Planck SZ catalogue clusters masked to 3 x R_500)

**Planck CMB lensing map:**
- Planck 2018 lensing convergence (kappa), N_side = 2048
- Lensing mask applied

**Cosmology:** Planck 2018, H_0 = 67.74 km/s/Mpc, Omega_m = 0.3075, Omega_b = 0.0486

### Methodology

#### 1. DisPerSE Filament Finding

DisPerSE works fundamentally differently from pair selection:
- Computes the density field from galaxy positions (Delaunay tessellation)
- Identifies critical points: maxima (nodes/clusters), saddle points, minima (voids)
- Traces field lines of the density gradient connecting saddle points to maxima
- These field lines = filament spines
- Persistence threshold removes features below a significance level (3 sigma here)

The algorithm is applied in 2D angular projections within redshift slices (Delta_z = 0.05, 8 slices from z=0.2 to z=0.6). This is coarser than He+2018's Delta_z=0.005 slices but uses the combined LOWZ+CMASS sample for higher galaxy density.

#### 2. Critical Point and Cluster Masking

**Crucial step to isolate filament signal from cluster/node contamination:**
- DisPerSE critical points (nodes) are masked to a radius of **5 Mpc** on each side of the filament
- Planck SZ catalogue clusters are masked to **3 x R_500**
- This removes the dense endpoints of filaments where clusters/groups reside
- Without this masking, the stacked signal would be dominated by bound structures at filament endpoints

#### 3. Stacking Procedure

For each filament segment:
1. Define the filament orientation from the segment direction
2. Extract a rectangular strip centred on the segment, oriented perpendicular to the filament spine
3. Stack strips from all segments, weighted equally
4. The result is a 1D radial profile of y or kappa as a function of distance from the filament spine

The stacking is done in **physical coordinates** (Mpc), converting angular separations using the filament redshift.

#### 4. Filament Model

Isothermal cylindrical beta-model:
- Electron density: n_e(r) = n_0 / (1 + (r/r_c)^2)^(beta/2), with beta = 2/3
- This is the cylindrical analogue of the isothermal beta-model used for galaxy clusters
- Projected along the line of sight for comparison with observed profiles
- Free parameters: central overdensity delta_c (= n_0/n_bar), core radius r_c, electron temperature T_e

For the tSZ signal: y(r) proportional to integral of n_e(r) * T_e along LOS
For the lensing signal: kappa(r) proportional to integral of rho(r) along LOS, assuming universal baryon fraction

#### 5. MCMC Fitting

- Joint fit to tSZ and lensing profiles simultaneously
- 3 free parameters: delta_c, r_c, T_e
- Covariance matrices from **10,000 bootstrap resamplings** of the filament catalogue
- Separate fits to tSZ-only, kappa-only, and joint tSZ+kappa

### Results

#### Detection Significance

| Signal | S/N | Mean value |
|--------|-----|------------|
| tSZ (y) | 4.4 sigma | y_bar = (0.6 +/- 0.2) x 10^-8 |
| CMB lensing (kappa) | 8.1 sigma | kappa_bar = (1.31 +/- 0.22) x 10^-3 |

The 8.1 sigma kappa detection is **far stronger** than de Graaff+2019's 1.9 sigma, despite both using Planck CMB lensing. The main reason: DisPerSE identifies actual filament locations (where the overdensity truly exists), while pair-stacking assumes a filament exists between every pair -- many pairs may not actually be connected by a filament.

#### Fitted Filament Parameters

| Parameter | tSZ only | kappa only | Joint fit |
|-----------|----------|------------|-----------|
| delta_c | -- | 14.4 +/- 3.0 | 19.0 +/- 4.3 |
| r_c [Mpc] | -- | 2.3 +/- 0.4 | 1.5 +/- 0.3 |
| T_e [10^6 K] | -- | -- | 1.4 +/- 0.4 |

- Central matter overdensity ~19x the mean cosmic density
- Core radius ~1.5 Mpc (filament "thickness")
- Temperature ~1.4 million K, consistent with WHIM predictions (10^5 - 10^7 K)

#### Baryon Budget

- Gas mass per unit filament length: M_gas/L = (4.2 +/- 1.6) x 10^13 M_sun / Mpc
- Baryon fraction in their filament sample: f_b = 0.080 +/- 0.025 of Omega_b
- Combined with de Graaff+2019 (shorter filaments, pair-based): total filament baryon fraction approaches ~30% of Omega_b
- Still leaves ~30% of baryons unaccounted for (possibly in lower-density regions not captured by either method)

### Null Tests and Systematics

| Test | Result |
|------|--------|
| Random filament orientations | Consistent with null (no signal when filament directions are randomised) |
| Rotated maps (90, 180 deg) | Consistent with null |
| CIB contamination | Negligible (stacking 857 GHz Planck map gives signal consistent with zero after cluster masking) |
| Redshift dependence | Signal consistent across z-bins (no strong evolution) |
| Filament length dependence | Longer filaments give slightly stronger signal per unit length |

### Comparison with Other Studies (from the paper)

| Study | Method | Filaments | tSZ S/N | kappa S/N |
|-------|--------|-----------|---------|-----------|
| de Graaff+2019 | Pair stacking (CMASS) | 1M pairs, 6-14 Mpc/h | 2.9 sigma | 1.9 sigma |
| Tanimura+2019 (earlier) | Pair stacking (LRG, z<0.4) | 262k pairs | 5.3 sigma | -- |
| **This paper** | DisPerSE filaments (LOWZ+CMASS) | 24,544 filaments, 30-100 Mpc | **4.4 sigma** | **8.1 sigma** |

The dramatically higher kappa S/N compared to de Graaff is attributed to:
1. Filaments identified by DisPerSE are genuine overdensities (pair-stacking includes many "empty" pairs)
2. Longer filaments (30-100 Mpc) provide more LOS integration path for kappa
3. Cluster masking (3xR_500) is more aggressive, reducing contamination from bound structures

### Contrast with Our Approach

| Aspect | Tanimura+2020 | Our project |
|--------|--------------|-------------|
| **Filament identification** | DisPerSE algorithm (topological, traces density gradient field lines) | Galaxy pair selection (geometric, assumes filament connects each pair) |
| **What is a "filament"** | Algorithmically-detected spine between critical points, 30-100 Mpc long | Matter bridge between two galaxies separated by 5-20 Mpc/h transverse |
| **Filament scale** | Large-scale (30-100 Mpc = ~44-147 Mpc/h) | Small-scale (5-20 Mpc/h transverse separation) |
| **Number of objects** | 24,544 filaments | ~10^5-10^6 pairs |
| **Galaxy sample** | SDSS DR12 LOWZ + CMASS combined (z=0.2-0.6) | BOSS CMASS only (z~0.4-0.7) |
| **CMB lensing map** | Planck 2018 kappa, N_side=2048 | Planck kappa from alm, N_side=2048 |
| **Stacking geometry** | Perpendicular strips along filament spine | 2D grid oriented along pair axis |
| **Output** | 1D radial profile kappa(r_perp) from filament spine | 2D kappa map in pair coordinate system |
| **Halo/cluster removal** | Mask DisPerSE critical points (5 Mpc) + Planck SZ clusters (3xR_500) | Control map subtraction (shifted single-galaxy stack) |
| **Significance (kappa)** | 8.1 sigma | TBD |
| **Filament model** | Isothermal cylindrical beta-profile (3 params: delta_c, r_c, T_e) | No model fitting yet; 2D visual comparison |
| **Covariance** | Bootstrap (10,000 resamplings) | Jackknife for single stacks; TBD for pairs |
| **Cosmology** | Planck 2018 | Planck 2018 (same, via astropy) |

### Key Differences and Implications

1. **Fundamentally different filament scales.** Tanimura+2020 studies filaments 30-100 Mpc long (mega-scale cosmic web), while our pair-stacking probes matter bridges at 5-20 Mpc/h separations (closer to the inter-halo scale). These probe different regimes of the cosmic web: their filaments are the large-scale "highways" connecting clusters, while our filaments are shorter bridges between galaxy-mass halos. The two approaches are complementary.

2. **DisPerSE vs. pair selection.** DisPerSE identifies where filaments actually are by tracing the density field topology. Pair selection assumes every close pair is connected by a filament -- which is statistically true on average but individually uncertain. This is likely why Tanimura achieves 8.1 sigma kappa detection with only 24,544 filaments, while de Graaff gets 1.9 sigma with ~1M pairs: **signal purity matters more than sample size**.

3. **Cluster/node masking strategy.** Tanimura masks filament endpoints aggressively (5 Mpc from critical points + 3xR_500 for SZ clusters). Our control-map subtraction aims to remove the average halo contribution but doesn't explicitly identify and mask galaxy clusters. We could improve our analysis by additionally masking known clusters near our galaxy pairs.

4. **1D radial profile vs. 2D map.** They compress to a 1D profile kappa(r_perp), which maximises S/N for a cylindrically symmetric signal. Our 2D maps retain more spatial information (potential asymmetry along vs. across the pair axis) but at the cost of lower per-pixel S/N.

5. **Model fitting.** Their isothermal cylindrical beta-model provides physically meaningful parameters (overdensity, core radius, temperature). We have not yet attempted model fitting. Fitting a similar cylindrical model to our pair-stacked residual maps would be a natural next step to extract filament properties.

6. **Baryon budget context.** They find 8% of cosmic baryons in their long filaments. de Graaff finds 11% in shorter pair-connected filaments. Combined, filaments may account for ~20-30% of baryons, still leaving a significant "missing" fraction in even lower-density regions (sheets, voids, diffuse IGM).

7. **Temperature measurement.** The joint tSZ+kappa fit constrains T_e = 1.4 x 10^6 K, consistent with WHIM predictions. Our kappa-only analysis cannot constrain temperature, but if we eventually include Planck y-map stacking, we could attempt a similar joint constraint for our shorter inter-halo filaments.

8. **Bootstrap vs. jackknife for covariance.** They use 10,000 bootstrap resamplings (sampling filaments with replacement), while de Graaff used 250 spatial jackknife regions. Both are valid; bootstrap is more computationally expensive but avoids the assumption that sub-regions are independent. For our pair stacks, either approach would work.

### References from this paper relevant to us

- Sousbie (2011) -- DisPerSE algorithm description
- Bonjean et al. (2018) -- earlier Planck y-map filament stacking work
- de Graaff et al. (2019) -- pair-stacking predecessor (reviewed above)
- Planck 2018 lensing paper -- same CMB lensing data we use
- Colberg et al. (2005) -- N-body predictions for filament properties between pairs

---
