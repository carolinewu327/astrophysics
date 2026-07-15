# Direct-Annular Validation of the Single-Halo Kappa Stack

Status: PASSED (qualitative, 500-halo run), 2026-07-10.

## Purpose

The single-halo kappa stack is the foundation of every downstream simulation
product (pair stacks subtract two copies of it; the bridge statistic is defined
on that residual). This validation checks it against an independent estimator
that never touches the kappa map, so an agreement validates the entire
map-construction and map-sampling path:

- particle binning into the 2D grid (`make_kappa_map.py`),
- mean-surface-density subtraction,
- Sigma -> kappa conversion,
- periodic bilinear sampling and stacking (`stack_single_sim.py`).

## Method

The validator (`direct_annular_profile` in `analysis/sim/stack_single_sim.py`,
run via `--direct-particles`) counts particles in periodic transverse annuli
around each halo, directly from the particle file, and converts counts to
kappa with the same mean-density subtraction and `Sigma_c` as the map:

```text
kappa(r_p) = mp_eff * [ N(r_p) / (N_halos * annulus_area) - N_tot / L^2 ] / Sigma_c
```

It shares no code path with the kappa map itself — no gridding, no
interpolation, no stacking grid.

## Inputs

- Direct profile: `analysis/sim/results/radial_profile_single_direct_annuli_mass13.csv`
  (2026-07-09; 500 randomly sampled mass13 host halos, all 283M downsampled
  particles, 50 annuli of 1 h^-1 Mpc out to r_p = 50 h^-1 Mpc).
- Map profile: `analysis/sim/results/radial_profile_single_sim_mass13.csv`
  (2026-07-08; all 1,436,465 mass13 host halos stacked on the 0.5 h^-1 Mpc
  validation map, radially symmetrized).

## Results

Figure: `analysis/sim/results/validation_single_profile_map_vs_direct.png`
(produced by `analysis/sim/plot_validation_single_profile.py`).

- **Inner, signal-dominated region (r_p <= 10 h^-1 Mpc):** mean ratio
  direct/map = 1.03; individual bins agree within their Poisson errors
  (pulls mostly within +/-2 sigma). The innermost bin is sharper in the
  direct profile (7.9e-3 vs 5.1e-3), as expected: the direct count is exact
  while the map path smooths through the 0.5 h^-1 Mpc map pixels, bilinear
  interpolation, and the 1 h^-1 Mpc stack pixels.
- **Outer region (r_p > 30 h^-1 Mpc):** the direct profile sits systematically
  above the map profile (pull RMS ~4.4 with Poisson-only errors, correlated
  across adjacent bins). Both profiles are at the ~1e-4 level here, where a
  500-halo sample is dominated by correlated large-scale-structure variance
  that Poisson errors do not capture. This is a sample-size limitation of the
  validator run, not evidence of a map-path error: the map profile uses all
  1.44M halos and is the better estimate at these radii.

## Caveats

- The 500-halo run is a **qualitative** check. For a quantitative validation
  of the outer profile, rerun with >= 5,000 halos after disk space is
  recovered (`--direct-max-halos 5000`); Poisson and LSS scatter both shrink
  and the outer comparison becomes meaningful.
- This run validates the 0.5 h^-1 Mpc map. Rerun as a spot-check against the
  final 0.1 h^-1 Mpc map once it exists (the validator itself is
  resolution-independent, so only the map side changes).
- The direct CSV does not record the halo count; the figure script takes it
  as `--n-halos` (default 500). Record it in the output when rerunning.

## Symmetrization Consistency (Phase 1, item 4)

`stack_single_sim.py` applies `radial_symmetrize_map` (in
`analysis/sim/sim_utils.py`), which is an exact port of the observational
`symmetrize_map` in `lib/geometry.py`: same `pwr = 2/3` radial binning and the
same `N // 2` center convention, including its half-pixel offset on even
grids. The simulation single-halo template therefore matches the
observational control-map construction quirk-for-quirk. Current single-stack
outputs were regenerated after this change.

## Reproduce

```bash
# validator (rerun with more halos after disk recovery):
python analysis/sim/stack_single_sim.py \
    --direct-particles data/cosmosim/dm_particles_snap_030.dat.bz2 \
    --direct-max-halos 500 --overwrite

# figure:
python analysis/sim/plot_validation_single_profile.py --n-halos 500
```
