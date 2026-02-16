# Outstanding Questions

Tracking open methodology questions that need discussion with adviser or further literature review.

---

## 1. Choice of r_par (line-of-sight separation) for pair selection

**Status:** Open

**Question:** Should `r_par` scale with the transverse separation being probed (e.g. r_par=5 for 5 Mpc/h pairs), or remain fixed at a large value (~5-20 Mpc/h) for all separations?

### Context

In filament pair-stacking studies, galaxy pairs are selected within a cylindrical shell defined by:
- **r_perp** (transverse/perpendicular separation): the projected on-sky distance between the two galaxies. This is the "separation" we vary (5, 10, 20 Mpc/h).
- **r_par** (parallel/line-of-sight separation): the radial comoving distance difference between the two galaxies. This absorbs redshift-space distortions (RSD).

Redshift-space distortions arise because galaxy peculiar velocities (~300-600 km/s for cluster-scale environments) add a Doppler shift on top of the Hubble flow, making galaxies appear displaced along the line of sight by several Mpc/h. This is the "Fingers of God" effect (for virial motions within clusters) and the Kaiser squashing effect (for coherent infall). As a result, two galaxies that are truly close in 3D space can appear separated by ~5-20 Mpc/h in the radial direction purely from velocity offsets.

### What the literature uses

| Study | r_perp range | r_par (LOS) criterion | Notes |
|---|---|---|---|
| Clampitt et al. 2016 | 6-14 h^-1 Mpc | Delta_z < 0.004 (~12 Mpc/h at z~0.3) | SDSS LRG pairs, weak lensing |
| Epps & Hudson 2017 | 6-10 h^-1 Mpc | Delta_z < 0.002 (~5 h^-1 Mpc comoving) | BOSS LOWZ+CMASS, weak lensing |
| de Graaff et al. 2019 | 6-14 h^-1 Mpc | < 5 h^-1 Mpc (comoving) | ~1M CMASS pairs, tSZ stacking |

**Key observation:** These foundational studies use r_par in the range 5-6 Mpc/h, which is much smaller than our current default of 20 Mpc/h. They do NOT scale r_par with r_perp -- it stays fixed as a "physically connected" criterion.

### What the papers say about RSD and r_par uncertainty

All three papers explicitly acknowledge that the line-of-sight separation is uncertain due to peculiar velocities, but they still use a relatively tight r_par cut:

**Clampitt et al. 2016** (r_los < 6 Mpc/h):
> "Note that this line-of-sight separation assumes the LRG velocity is only due to Hubble flow; in other words, the redshift difference can arise from the difference of line-of-sight peculiar velocities (Delta_v = 1200 km/s for r_los = 6 Mpc/h) even if the two LRGs are in the same distance. This is the so-called redshift space distortion (RSD)."
>
> Systematics section (p10): "Redshift space distortions: the line of sight separation of the LRGs is uncertain owing to their relative peculiar velocity. We have attempted to account for it by adding RSD to the simulations."

They explicitly note that their 6 Mpc/h cut corresponds to Delta_v = 1200 km/s of peculiar velocity scatter. They accept the RSD contamination and model it in simulations rather than widening r_par.

**Epps & Hudson 2017** (Delta_z < 0.002, ~5 h^-1 Mpc):
> "parameterized the line-of-sight separation of the two LRGs by a Gaussian distribution with width sigma_LRG = 8 h^-1 Mpc. This separation in redshift space includes both the peculiar velocities of each LRG in the pair and the Hubble flow. The peculiar velocities are difficult to model since they include contributions from relative infall motions as well as 'thermal' motions of the LRGs themselves within their host halos."

They note that the *true* LOS distribution has sigma ~ 8 h^-1 Mpc (wider than their 5 Mpc/h cut), meaning they intentionally use a tight cut knowing it rejects some real pairs. They also emphasize that spectroscopic redshifts are essential:
> "the uncertainty associated with photometric redshifts will scatter true physical pairs away from each other... sigma_z_phot = 0.05 corresponds to ~150 h^-1 Mpc scatter. This is much larger than the physical line-of-sight separation of order 10 h^-1 Mpc."

**de Graaff et al. 2019** (r_par < 5 h^-1 Mpc):
> "Our criterion on line-of-sight separation implicitly assumes that the total redshift differences are purely cosmological; but in reality, peculiar velocities can have an effect -- as can be seen from e.g. Fig. 11 of Jenkins et al. (1998). For filaments of our length that are near to the plane of the sky, the pairwise dispersion in radial velocity is close to 500 km/s, and this has a convolving effect. Thus, some of our pairs will have a true radial separation slightly larger than our limit of 5 h^-1 Mpc, and some pairs with a smaller radial separation will be rejected because the apparent separation in redshift space is above our limit. But our selection still picks pairs that are nearly transverse to the line of sight."

This is the most explicit statement: they acknowledge 500 km/s pairwise velocity dispersion causes misclassification in both directions, but argue the selection still isolates nearly-transverse pairs, which is what matters for filament detection.

### Summary of RSD handling

The literature consensus is:
1. **Use a tight r_par cut** (5-6 Mpc/h), not a wide one
2. **Accept that RSD causes some misclassification** (real pairs rejected, false pairs included)
3. **The key goal is selecting pairs that are nearly transverse to the line of sight** -- a tight r_par ensures this
4. **A wider r_par (like 20 Mpc/h) would include too many chance projections**, diluting the filament signal
5. **Model the RSD effect in the analysis** rather than trying to avoid it by widening the cut

Our current default of r_par = 20 Mpc/h is ~4x larger than what these studies use. This should be discussed with the adviser.

### Trade-offs

- **Smaller r_par (e.g. 5 Mpc/h):** Consistent with literature. Selects pairs more likely to be physically connected and nearly transverse. Reduces projection contamination. Rejects some real pairs due to peculiar velocities (acceptable).
- **Larger r_par (e.g. 20 Mpc/h):** Much more inclusive than literature standard. Captures pairs with large RSD offsets, but also includes many unconnected projections. Dilutes the filament signal.
- **Scaling r_par with r_perp:** Not done in the literature. All three studies use a fixed r_par regardless of r_perp range.

### Decision needed

1. Should we reduce r_par from 20 to ~5 Mpc/h to match the literature?
2. Should we run sensitivity tests with r_par = 5, 10, 20 to see the effect?
3. Should r_par remain the same for all transverse separation bins?

### References

- Clampitt, Miyatake, Jain & Takada (2016), "Detection of stacked filament lensing between SDSS luminous red galaxies", MNRAS 457, 2391. [arXiv:1402.3302](https://arxiv.org/abs/1402.3302) -- local: `notes/references/1402.3302v3.pdf`
- Epps & Hudson (2017), "The weak lensing masses of filaments between luminous red galaxies", MNRAS 468, 2605. [arXiv:1702.08485](https://arxiv.org/abs/1702.08485) -- local: `notes/references/1702.08485v1.pdf`
- de Graaff, Cai, Heymans & Peacock (2019), "Probing the missing baryons with the Sunyaev-Zel'dovich effect from filaments", A&A 624, A48. [arXiv:1709.10378](https://arxiv.org/abs/1709.10378) -- local: `notes/references/1709.10378v3.pdf`
- Jenkins et al. (1998) -- cited by de Graaff for pairwise velocity dispersion (Fig. 11)

---

## 2. Sign convention for galaxy-random subtraction

**Status:** Open

**Question:** When computing the filament excess signal, should we subtract as (pair_stacked - control_pair) or (control_pair - pair_stacked)?

### Context

The filament signal is isolated by subtracting the "control" (what you'd expect from two independent galaxy halos) from the pair-stacked signal. The sign convention determines whether a positive filament signal means excess convergence (overdensity) between the pair members, which is the physical expectation.

### Decision needed

Confirm with adviser which sign convention is used and ensure `plot_results.py` implements it consistently.

---

## 3. Jackknife error estimation for pair stacking

**Status:** Open

**Question:** Should we add jackknife error estimation to `stack_pairs.py` to match `stack_single.py`?

### Context

Currently `stack_single.py` produces a jackknife error map (dropping one HEALPix region at a time), but `stack_pairs.py` produces only the mean stacked kappa map with no error estimate. This means we cannot compute a proper signal-to-noise map for the filament signal (`filament = pair_stacked - control`), since we'd be missing the pair-stack uncertainty.

### What the literature does

- **de Graaff+2019:** Jackknife with 250 spatial sub-regions (~40 deg² each, ~4000 pairs each). Full stacking + halo fitting repeated per realisation. Hartlap correction applied to inverse covariance.
- **Kondo+2019:** 108 mock realisations from full-sky ray-tracing simulations (more expensive but captures sample variance).
- **Tanimura+2020:** 10,000 bootstrap resamplings of the filament catalogue.

### Trade-offs

- Jackknife for pairs is expensive: each realisation re-stacks ~100k+ pairs against the full Planck map. With ~100 regions, that's ~100x the computation of a single stack.
- Bootstrap (resample pairs with replacement) is an alternative that doesn't require spatial region assignment.
- Without any error on the pair stack, S/N estimates are incomplete.

### Decision needed

1. Is this worth the computation cost now, or should we first confirm we see a signal before investing in error estimation?
2. If yes, use jackknife (spatial regions, matching `stack_single.py`) or bootstrap (simpler, no random catalog needed)?

---

## 4. Inverse-variance weighting when combining North/South regions

**Status:** Open

**Question:** Should we use inverse-variance weighting instead of simple averaging when combining North and South kappa maps?

### Context

`plot_results.py` currently combines per-region maps with a simple `np.mean()`. However, BOSS CMASS North has ~2x more galaxies than South (~560k vs ~290k), so the maps have different noise levels. A simple average gives equal weight to both, which is suboptimal.

The proper approach is inverse-variance weighting:
```
combined = (map_N / σ_N² + map_S / σ_S²) / (1/σ_N² + 1/σ_S²)
```

We already have per-pixel jackknife error maps for single-galaxy stacks (`error_single_{dataset}_{region}.csv`), so the information is available for Maps (1) and (3). For pair-stacked maps (2) and (4), we'd first need error estimates (see question #3).

### Decision needed

1. Implement inverse-variance weighting in `plot_results.py` once error maps are available for all components?
2. Is simple averaging acceptable for a first-look analysis?
