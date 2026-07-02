# Simulation Analysis Plan

## Main Goal

Use BigMDPL simulation data to build controlled kappa stacks around halos and halo pairs, then use those stacks to interpret the BOSS/Planck results. The immediate first deliverables are:

1. Stacked kappa map around single halos in a narrow mass bin near `log10(Mvir / h^-1 Msun) = 13`.
2. Stacked kappa map around halo pairs in the same mass bin.
3. Comparison of fixed-separation stacking vs normalized pair-separation stacking.

## Motivation

The updated observational plots still show features that are difficult to interpret, including low-amplitude rings in the random-pair kappa maps and ambiguous filament residuals. Simulation-based analyses should help in two ways:

1. Provide controlled tests that guide improvements to the observational analysis.
2. Produce a theoretical model to compare directly against the observational results.

## Phase 1: Set Up Simulation Inputs

Download or obtain the two BigMDPL files at `z = 0.547`:

- Downsampled particles: `dm_particles_snap_030.dat.bz2`
- Rockstar halo catalog: `hlist_0.64640.list`

Because the halo catalog is very large, first create a compact derived halo file with only:

```text
id, upid, Mvir, x, y, z, vx, vy, vz
```

Then filter to host halos only:

```text
upid = -1
```

Create the initial science sample:

```text
|log10(Mvir / h^-1 Msun) - 13| < 0.05
```

Keep the mass cut configurable so nearby mass bins can be tested later.

## Phase 2: Build the Simulated Kappa Map

Construct a projected 2D kappa map from the downsampled dark matter particles by projecting along the simulation `z` direction.

Use:

```text
Sigma_ij = mp_eff [Np(x_i, y_j) / l^2 - Np_tot / L^2]
kappa_ij = Sigma_ij / Sigma_c
```

with:

```text
L = 2500 h^-1 Mpc
mp = 2.359e10 h^-1 Msun
mp_eff = 200 mp
Sigma_c = 8.88e14 h Msun Mpc^-2
```

Start with a manageable grid before going full resolution:

```text
l = 0.5 h^-1 Mpc   -> 5000 x 5000 map
l = 0.25 h^-1 Mpc  -> 10000 x 10000 map
l = 0.1 h^-1 Mpc   -> 25000 x 25000 final map
```

Store the kappa map as `float32` and use `np.memmap` for the large map products. The final `25000 x 25000` map is already about 2.5 GB in `float32`; avoid `float64` unless a validation test shows it is needed.

If multiprocessing is used for map construction or stacking, do not rely on implicit inherited in-memory arrays unless running on Linux with `fork`. For portability across macOS and spawn-based cluster jobs, workers should open the memmap by path.

Subtracting the global mean surface density is essential. Without this, the stacks will have an arbitrary positive background.

## Phase 3: Single-Halo Stack

Reuse the logic of the current observational single stack, but replace:

```text
galaxy catalog -> halo catalog
Planck kappa map -> simulated kappa map
sky projection -> periodic x-y map sampling
```

For each halo, cut out a square kappa patch centered at `(xh, yh)` using periodic wrapping. Match the current observational physical map size and binning, likely controlled by:

- `BOX_SIZE_HMPC`
- `GRID_SIZE`

All non-grid-centered samples from the simulated kappa map should use periodic bilinear interpolation. Avoid nearest-neighbor sampling because it can create pixel-locking artifacts in radial profiles and rotated pair cutouts.

Suggested outputs:

```text
kappa_single_sim_mass13.csv
kappa_single_sim_mass13.png
radial_profile_single_sim_mass13.csv
radial_profile_single_direct_annuli_mass13.csv
```

Validation checks:

- The stacked map should be circularly symmetric before any symmetrization.
- The radial profile should be positive near the center and approach zero at large radius.
- Changing grid resolution should not significantly change the profile beyond the pixel-scale region.
- Directly compute an independent single-halo radial profile from particle counts in periodic transverse annuli around halos, using the same mean-density subtraction and normalization. This direct-annular profile should agree with the map-sampled profile outside pixel-scale differences.

## Phase 4: Halo Pair Selection

Select pairs from the same host-halo mass bin.

Apply redshift-space distortion only for pair selection:

```text
z_p <- z_p + (1 + z) vz / H(z)
```

Use the BigMDPL cosmology for this step:

```text
Omega_m = 0.307
h = 0.677
H_RSD(z) = 100 sqrt(Omega_m (1 + z)^3 + 1 - Omega_m) km s^-1 (h^-1 Mpc)^-1
```

Do not accidentally inherit a default Planck18 cosmology here. Use the `H0 = 100` convention above, not `100 h`, so `(1 + z) vz / H_RSD(z)` is directly in the simulation's `h^-1 Mpc` coordinate units before periodic wrapping back into `[0, L)`.

Then select pairs using cuts that match the observational analysis as closely as possible:

```text
projected separation r_perp: e.g. 5, 10, 20 h^-1 Mpc
line-of-sight separation: match observational pair cut
```

For each pair, store:

```text
x1, y1, z1_redshiftspace
x2, y2, z2_redshiftspace
M1, M2
r_perp
r_parallel
pair_center_x, pair_center_y
cos_theta, sin_theta
```

Use on-the-fly periodic treatment. For the pair axis, shift one halo coordinate by `L` if needed before computing the pair center and angle. For sampling kappa around the pair center, wrap `x - xc` and `y - yc` into `[-L/2, L/2]`.

## Phase 5: Pair Kappa Stack

Implement the simulation analog of `analysis/boss/scripts/stack_pairs.py`.

For each pair:

1. Define pair center `(xc, yc)`.
2. Define the X axis along the pair.
3. Create the same `101 x 101` pair grid used in the observational stack.
4. Inverse-rotate grid points into simulation `(x, y)`.
5. Sample the simulated kappa map with periodic wrapping and bilinear interpolation.
6. Average all pairs.
7. Apply the same reflection symmetrization as the observational pipeline.

Suggested outputs:

```text
kappa_pairs_sim_mass13_rperp5.csv
kappa_pairs_sim_mass13_rperp10.csv
kappa_pairs_sim_mass13_rperp20.csv
```

Validation checks:

- Two halo peaks should appear at roughly `X = +/- r_perp / 2`.
- The stack should be symmetric under `Y -> -Y`.
- If pair order is randomly swapped, the map should not change after reflection symmetrization.
- If halo positions are randomized, the pair signal should mostly vanish.

## Phase 6: Test Normalized Pair Separation

This directly addresses the suggestion to stretch pairs to a common normalized separation before stacking.

Run two versions.

### A. Fixed Physical Stacking

Use current observational-style coordinates:

```text
X, Y in h^-1 Mpc
```

Pairs with different separations are stacked without rescaling.

### B. Normalized Stacking

Before stacking, rescale each pair so the two halos land at common coordinates:

```text
X_norm = X / r_perp
Y_norm = Y / r_perp
```

or stretch the kappa sampling grid so each pair maps to endpoints:

```text
halo 1 at X = -1/2
halo 2 at X = +1/2
```

Compare:

- Halo peak sharpness.
- Filament or inter-pair bridge contrast.
- Residual ring-like features.
- Radial profiles around the pair center.
- Transverse profile at `X = 0`.

This is likely the most useful simulation experiment for deciding whether normalized stacking should be used on the observational data.

## Phase 7: Understand the Filament Signal

After pair stacks are working, split the pair map into components:

```text
raw pair stack
minus stacked single-halo model at both halo positions
equals residual / filament-like map
```

For each pair stack, compute:

```text
kappa_pair_raw(X, Y)
kappa_two_halo_model(X, Y) = kappa_single around halo 1 + kappa_single around halo 2
kappa_residual = kappa_pair_raw - kappa_two_halo_model
```

This will test whether the apparent bridge is:

- A genuine correlated matter filament.
- Overlap of two halo outskirts.
- An artifact of pair-selection geometry.
- Caused by smoothing.
- Caused by stacking pairs with a broad separation distribution.

## Phase 8: Investigate Ring Features

Use simulation to isolate whether rings can arise without Planck noise, masks, or random catalogs.

Run a controlled sequence:

1. Stack around true halos.
2. Stack around random points in the periodic box.
3. Stack around true halo pairs.
4. Stack around random pairs with matched separation distribution.
5. Stack around shuffled halo pairs.
6. Repeat with and without smoothing.
7. Repeat with fixed vs normalized separation.

If rings appear in simulation random-pair stacks, likely causes include interpolation, map pixelization, pair-coordinate rescaling, annular averaging, or symmetrization.

If rings only appear in observational randoms, the cause is more likely survey mask, Planck mask, sky projection, or random-catalog geometry.

## Phase 9: Smoothing Tests

Do smoothing in simulation after constructing the high-resolution kappa map.

Convert the observational smoothing scale into comoving transverse size at `z = 0.547`:

```text
R_smooth = D_l * theta_smooth
```

With `D_l = 1424.50 h^-1 Mpc`:

```text
8 arcmin ~= 3.31 h^-1 Mpc
4 arcmin ~= 1.66 h^-1 Mpc
2 arcmin ~= 0.83 h^-1 Mpc
```

Run stacks at:

```text
no smoothing
2 arcmin equivalent
4 arcmin equivalent
8 arcmin equivalent
```

This will show how much the bridge and ring-like features are washed out by the current Planck smoothing.

## Phase 10: Observation Comparison

Only after the simulation behavior is understood, compare to BOSS/Planck.

Use the same map size, binning, pair cuts, mass scale, and smoothing. The comparison should include:

```text
single-halo radial profile: observation vs simulation
pair stack: observation vs simulation
two-halo-subtracted residual: observation vs simulation
normalized vs unnormalized pair stack
```

At that point, decide whether to revisit:

- Planck PR4/NPIPE.
- Smaller smoothing.
- Improved pair weighting.
- HOD or abundance-matched galaxy-halo modeling.

ACT DR6 can wait or likely be dropped for now because the BOSS overlap is too small to help much at this stage.

## Immediate Coding Milestones

1. Add `analysis/sim/prepare_halos.py`.
   - Reads Rockstar catalog.
   - Filters `upid = -1`.
   - Writes compact mass-selected halo file.

2. Add `analysis/sim/make_kappa_map.py`.
   - Reads downsampled particles.
   - Projects to 2D.
   - Subtracts mean surface density.
   - Writes kappa map.

3. Add `analysis/sim/stack_single_sim.py`.
   - Produces single-halo kappa stack.

4. Add `analysis/sim/find_pairs_sim.py`.
   - Applies redshift-space `z` shift.
   - Performs periodic pair selection.
   - Writes pair catalog.

5. Add `analysis/sim/stack_pairs_sim.py`.
   - Produces pair kappa stacks.
   - Include options:

```text
--normalize-separation
--smooth none|2arcmin|4arcmin|8arcmin
--rperp-bin
--mass-bin
```

## First Target Figure

The first useful result should be a four-panel figure:

```text
single-halo stack
pair stack, fixed physical separation
pair stack, normalized separation
pair residual after subtracting two single-halo profiles
```

This figure should give a concrete basis for deciding how to update the observational analysis.
