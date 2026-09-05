# Paper: Filaments between Dark Matter Halos from CMB Lensing

Manuscript sources, final figures, and the record of which script makes
which figure. Prose is edited in Overleaf (project
`Filaments_from_CMB_Lensing`); this directory holds the snapshots that are
committed, the figures the scripts produce, and the bibliography.

## Layout

```
paper/
  README.md                 this file
  styau/
    main.tex                aastex701 preamble; \input's the numbered sections
    0_abstract_draft.tex
    1_intro_draft.tex
    2_data_draft.tex
    3_methods_results_draft.tex   Sections 3 (Methods), 4 (Expected signal), 5 (Results)
    4_summary_draft.tex
    5_acknowledgments_draft.tex
    refs.bib
    Figure/                 final PDFs, referenced via \graphicspath{{Figure/}}
```

Compile locally with `latexmk -pdf main.tex` from `paper/styau/`. On
Overleaf, upload the section files and `Figure/` and set `main.tex` as the
main document.

## Conventions inside the .tex files

- `% CHECK` -- a number recorded in the repo (notes, docstrings, figures) but
  not reproducible from a checked-in result file. Verify before submission.
- `% TODO` -- a quantity the pipeline produces but whose value is not in the
  repo. Fill from the run outputs.
- `% UNVERIFIED` -- describes code that is not in the repo.

## Scoring the mock with the fixed physical bands

**Problem.** The mock pair maps are already stacked at fixed separation
(`jackknife_pair_stack.py` calls `stack_pairs(...,
normalize_separation=False)`), so the stacking is fine. The mismatch is in
the scoring: that script passes `axis / rperp_center` to
`summarize_sim_sensitivity.map_stats`, so its central and off-centre Y bands
scale with separation (`|Y| <= 0.15 r_perp` and `0.45-0.85 r_perp`). The
BOSS estimator uses fixed physical bands, `|Y| <= 1.5` and
`1.5 <= |Y| <= 10.5 h^-1 Mpc`, with only the bridge X window tied to the
separation, `|X| <= 0.35 r_perp`. The mock numbers quoted so far (6, 3.4,
1.5 x 1e-4) are therefore not the same statistic as the BOSS numbers, and the
apparent factor-of-four decline with separation is partly the band definition
changing, not the physics.

**Decision.** Keep the BOSS estimator and change the mock scoring. This is
required for a like-for-like data--mock comparison: the transverse aperture
and smoothing scale are defined in physical units in the measurement, while
the bridge length follows the pair separation. Separation-scaled Y bands
remain useful as a sensitivity test, but answer a different question. Fixed
Y bands do not by themselves remove halo leakage; the bridge-window edge is
`0.15 r_perp` from each halo, so leakage can still vary between bins.

**Fix.** In `analysis/sim/jackknife_pair_stack.py`, add
`bridge_excess(map, axis, rperp_center)` from `lib/geometry.py` for both the raw
and template-subtracted maps, including every leave-one-out map. Retain the
existing `map_stats(..., norm_axis, ...)` outputs as explicitly legacy,
scaled-band diagnostics, but make production readers require the fixed-band
result. Keep `summarize_sim_sensitivity.py` unchanged. Correct the
matched-single path to
`kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv`.

Provenance rules for the output JSON, so old and new files cannot be
confused:

- Write the fixed-band result under **new keys**
  (`residual_bridge_excess_fixedband_kappa` and
  `raw_bridge_excess_fixedband_kappa`),
  not by redefining `residual_bridge_excess_kappa`. Old JSONs keep their old
  meaning.
- Add a `band_definition` block (`central_half_y_hmpc`, `off_lo_hmpc`,
  `off_hi_hmpc`, `bridge_half_x_frac`), the configured r_perp selection bounds
  (`rperp_min_hmpc`, `rperp_max_hmpc`), and the observed extrema
  (`rperp_observed_min_hmpc`, `rperp_observed_max_hmpc`). Add required
  `--rperp-min` and `--rperp-max` arguments to the stacker, validate every pair
  against them, and reject a catalogue whose observed extent is inconsistent
  with the declared bin. Longer term, copy the configured bounds from metadata
  written by `find_pairs_sim.py` rather than relying on manually repeated
  command-line values.
- Make every reader that compares mock to BOSS (`sim_diff` in
  `plot_observed_figures.py`, `plot_obs_vs_sim_filament.py`) **assert** that
  the band definition equals the `lib/geometry.py` constants and that the
  r_perp bin equals the BOSS bin, alongside the existing `rpar_max` check.
  Nothing asserted the r_perp bin before, which is how the 19-21 stack got
  into a committed figure.

**Also in scope: two scripts the plan above did not name.** The scaled
convention is hard-coded in `analysis/sim/compare_observed_sim.py:109-113` and
`analysis/sim/check_boss_consistency.py:110-114`, both tracked. Left alone
they would become a *fourth* band definition, and the first is an
observed-versus-mock comparison whose filename advertises nothing. They are
not migrated -- doing so would silently change numbers that have been quoted
-- but they are now labelled: a warning block at the top of each module, a
`band_convention` column in the CSV they write, and a `_scaledband` suffix on
the output filenames (`boss_vs_sim_bridge_stats_scaledband_*.csv`,
`boss_consistency_stats_scaledband.{csv,png}`). Nothing read those files
programmatically, so the rename is safe.

**Order of operations. DONE.** The 20 h^-1 Mpc mock bin was fixed
(19-21 -> 18-22) *before* any number went into the text, so the fill was a
single pass with final values and no footnote. Original reasoning:
fix the bin *before* filling any numbers into the text. Doing it
afterwards means the table, abstract, summary, Section 4 sentence, the
quadrupole figure, and the blocks file all change twice. While there, note
that the manifest row for `sim_quadrupole_diff_maps.pdf` uses
`stack_rperp5_rpar10_matched.csv` (r_par <= 10) whereas `tab:bridge` and
`fig3` use the matched r_par <= 5 stack at 5 h^-1 Mpc; the Section 4 map and
its numbers should come from the same sample.

**Recompute.** For 5 and 10 h^-1 Mpc, central values require no pair
re-stacking: score the existing matched maps and two-halo controls with the
shared estimator. The current 20 h^-1 Mpc value can be scored only as a
provisional regression check; its final value requires the corrected 18--22
h^-1 Mpc pair catalogue and a new stack. Correct jackknife errors require the
per-block accumulators `stack_rperp{5,10,20}_blocks.npz`, which are external and
may not be present in a fresh checkout. If an accumulator is missing, rerun
only `jackknife_pair_stack.py` for that bin with `--blocks-output`; do not
regenerate the HOD catalogue or convergence map, and do not regenerate the 5
or 10 h^-1 Mpc pair catalogues. Budget about one hour for the code and checks,
plus roughly five minutes of stacking per missing bin on the machine used for
the original runs.

**Acceptance checks.**

1. The full map reconstructed from block sums matches the saved stack (the
   existing files agree to 3e-8 of peak).
2. ~~`bridge_excess` equals `band_profile` averaged over the bridge
   window.~~ **Dropped: it cannot fail.** `geometry.bridge_excess` *is*
   `band_profile(arr, axis)[window].mean()` -- one line, deliberately, so the
   plotted curve and the quoted number cannot drift. Asserting it tests a
   function against its own body. The real cross-check is check 3, which
   compares two genuinely different code paths (the JSON route builds its
   template with `sim_utils.make_two_halo_template` and jackknifes over the
   block sums; `fig3` uses `geometry.two_halo_template` and rebuilds the
   leave-one-out maps itself).
3. Regression targets for the currently saved stacks, **measured**:
   `6.32 +/- 0.36`, `4.91 +/- 0.21`, `2.18 +/- 0.12` (x 1e-4).

   The central values are as predicted. The errors at 10 and 20 are *not* the
   `+/- 0.23` and `+/- 0.17` quoted earlier: those were the mean of the
   per-bin jackknife errors across the window, which is not the error on the
   window mean. The correct quantity jackknifes the windowed mean itself. The
   two coincide at 5 h^-1 Mpc (0.357 against 0.356) because the window is
   narrow and its bins are almost perfectly correlated, and diverge as the
   window widens -- by 6 per cent at 10 and 31 per cent at 20. Use the
   error-of-the-mean; the JSON now reports it.

   The 20 h^-1 Mpc value is **final**: the bin fix is done. The corrected
   18-22 catalogue (`pairs_hodnmatch50_rperp20wide_rsd50.csv`, 7,884,478
   pairs, 3.6 min) and its stack (1,856,389 pairs after r_par <= 10, 18 min)
   replaced the 19-21 products in place; the superseded ones are archived
   under `hod_pairs/archive_bin19to21/`. The value moved 1.94 -> 2.18, a 12
   per cent shift and larger than either error -- which is why the bins had
   to match rather than be caveated.
4. Comparison statistics, to replace the "N bins beyond 2 sigma" count
   (1 h^-1 Mpc bins are correlated over the 3.3 h^-1 Mpc beam, so exceedances
   cluster and the count is not a 5 per cent Gaussian test): add
   `chi2_A1 = (BOSS - mock)^T C^-1 (BOSS - mock)` on the 22-element compressed
   vector to `fit_band_covariance.py`, alongside its existing BOSS-versus-zero
   `chi2_A0` and best-fit-amplitude statistics. Also report the fitted
   amplitude A of the mock profile to BOSS in the bridge window with its
   error. The covariance is the BOSS jackknife covariance; state that the much
   smaller mock statistical covariance is neglected.

**Propagate the regenerated values to** (one pass, after the bin fix):

- `stack_rperp*_matched.json`, `band_profiles_obs_vs_sim.pdf`,
  `sim_quadrupole_diff_maps.pdf`, `sim_band_profiles_all_seps.pdf`, and
  `obs_vs_sim_quadrupole_maps.pdf`.
- `3_methods_results_draft.tex`: Mock column of `tab:bridge`; the Section 4
  sentence quoting "approximately 6, 3.4, 1.5"; the trend wording (the
  fixed-band sequence is a factor ~3 overall with a small 5 -> 10 drop, not a
  factor of four).
- `0_abstract_draft.tex` and `4_summary_draft.tex`: the sentences quoting the
  mock prediction and its decline.
- **BOSS numbers in the same pass.** The production values in
  `filament_jackknife_BOSS_North_South.csv` are 5.38 +/- 9.08,
  11.91 +/- 5.10, 1.24 +/- 2.22 x 1e-4 (0.6, 2.3, 0.6 sigma). The drafts,
  the abstract, and the S/N column of `tab:bridge` still carry the older
  0.7 / 2.2 / -1.1 sigma from the August lit-review note; the 20 h^-1 Mpc
  sign has flipped. That CSV is the authoritative source.
- Caption of `fig:band_profiles_obs_sim` and the text: the mock band is
  statistical only (3-5 per cent of the BOSS error); it excludes the
  uncertainty in the HOD and the single mock realisation.

**Optional, recommended.** The 6.32 at 5 h^-1 Mpc should be interpreted as
bridge plus any *residual* halo leakage, since the window edge is 0.75 h^-1
Mpc from the halo centre, well inside the beam. The already-subtracted
superposed-single template cannot measure leakage caused by its own mismatch
to paired haloes. Quantifying that term needs a clearly defined no-filament
null that preserves the paired-halo mass, environment, selection, and
smoothing--for example, a paired-halo-only simulation with the inter-halo
matter removed. `4_summary_draft.tex` has a TODO for this.

## After submission

Deliberately deferred. Each is understood and none affects a number in the
paper; they are structural debt that is cheaper to pay once the deadline is
past. Recorded here so the reasoning is not lost.

1. **Put the transverse bin in the filename.**
   `stack_rperp20_matched.json` does not say which r_perp bin produced it,
   which is why the 18-22 run overwrote the 19-21 one in place and needed a
   hand-made archive (`hod_pairs/archive_bin19to21/`, see its README).
   `stack_rperp20_bin18to22_matched.json` makes a new bin a new file and
   clobbering structurally impossible. Touches `SIM_STACK`/`SIM_BLOCKS` in
   `plot_observed_figures.py`, `PAIR_STACK` in `deconvolve_pair_profile.py`,
   the `stack_rperp{key}_matched` pattern in `plot_obs_vs_sim_filament.py`,
   and the figure manifest above. Mechanical, but it touches every reader at
   once -- and those readers are producing the figures being shipped.
   Mitigated meanwhile: file *contents* are self-describing (`rperp_min_hmpc`
   etc. in every JSON and `.npz`) and `assert_comparable()` checks the bin
   before plotting, so the dangerous failure mode -- a wrong comparison -- is
   already closed. What remains is the inconvenience.

2. **`chi2_A1` and the fitted mock amplitude** (acceptance check 4 above).
   Not implemented. Until it is, `Obs.-Mock` in `tab:bridge` is a
   per-separation z-score, which ignores the bin-to-bin covariance.

3. **`plot_theory_figures.py` has no `assert_comparable`.** The Section 4
   figures bypass the band/bin guard that the Section 5 figures now enforce.

4. **Section 4 and Section 5 use different samples at 5 h^-1 Mpc.**
   `deconvolve_pair_profile.PAIR_STACK["5"]` is the r_par <= 10 stack;
   `tab:bridge` and `fig3` use r_par <= 5. The Section 4 map and its numbers
   should come from one sample.

5. **Read the r_perp bounds from `find_pairs_sim.py` metadata** rather than
   from repeated command-line values (noted as "longer term" above).

6. **Migrate or retire the two legacy scaled-band scripts.**
   `compare_observed_sim.py` and `check_boss_consistency.py` are labelled, not
   migrated. Labelling stops a wrong quote; it does not remove the fourth
   convention from the tree.

7. **Deposit the external inputs with a DOI.** The ~45 MB BOSS accumulators
   and the `*_blocks.npz` are not in git, so no figure marked **external**
   above is reproducible from a fresh clone.

8. **`numbers.tex`.** Generated macros written by the pipeline and `\input` by
   `main.tex`, so text numbers cannot drift from the results -- the standing
   fix for the `% CHECK` / `% TODO` problem described at the end of this file.

9. **`.gitignore` covers only `data/`.** `analysis/*/results/`, `*.npz`,
   `*.float32`, `*.bin`, `.DS_Store` and `.idea/` are all one `git add -A`
   away from the index; the tree holds ~27 GB of such files.

10. **~20 analysis scripts remain untracked**, including
    `boss_band_core_tail.py`, `band_chi2_significance.py`, `band_ratio_Q.py`
    and `quadrupole_moments.py`. Pure code, no data question.

## Figure manifest

Every figure in the paper should have a row here: the file, the script that
writes it, its inputs, and the exact command. Run commands from the repo
root in the WSL `astro` environment. Regenerate figures before each
Overleaf upload.

| Figure file (`styau/Figure/`) | Section | Script | Inputs | Command |
|---|---|---|---|---|
| `footprint.pdf` | 2, Fig. 1 left | `analysis/boss/scripts/plot_footprint.py` | `data/BOSS/galaxy_DR12v5_CMASS_{North,South}.fits.gz`, `data/planck/COM_Lensing_4096_R3.00/mask.fits` | `python analysis/boss/scripts/plot_footprint.py --output-dir paper/styau/Figure` |
| `zdist.pdf` | 2, Fig. 1 right | same run as above | same | same |
| `footprint_stats.txt` | 2 (numbers in text) | same run as above | same | same; not a figure, the counts and redshift statistics quoted in the Data section |
| `sim_pair_vs_control_maps_sep10.pdf` | 4 | `analysis/sim/plot_theory_figures.py` | `analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv`, `analysis/sim/results/hod_pairs/stack_rperp10_matched.csv` | `PYTHONPATH=lib:analysis/sim python analysis/sim/plot_theory_figures.py --paper --only 1 --output-dir paper/styau/Figure` |
| `sim_quadrupole_diff_maps.pdf` | 4 | same script, `--only 2` | same single, plus `stack_rperp{5_rpar10,10,20}_matched.csv` | `PYTHONPATH=lib:analysis/sim python analysis/sim/plot_theory_figures.py --paper --only 2 --output-dir paper/styau/Figure` |
| `sim_band_profiles_all_seps.pdf` | 4 | same script, `--only 3` | same as above | `PYTHONPATH=lib:analysis/sim python analysis/sim/plot_theory_figures.py --paper --only 3 --output-dir paper/styau/Figure` |
| `boss_single_and_pair_maps.pdf` | 5.1 | `analysis/boss/scripts/plot_observed_figures.py` | **external** `analysis/boss/results/jk/acc_single_{galaxy_scw,random_scw_frac100}_BOSS_North_South.npz` and `acc_pairs_galaxy_{5,10,20}_*.npz` | `PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python analysis/boss/scripts/plot_observed_figures.py --paper --only 1 --output-dir paper/styau/Figure` |
| `obs_vs_sim_quadrupole_maps.pdf` | 5.1 | same script, `--only 2` | the same accumulators, plus `analysis/sim/results/hod_pairs/stack_rperp{5,10,20}_matched.{csv,json}` and the mock single | `PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python analysis/boss/scripts/plot_observed_figures.py --paper --only 2 --output-dir paper/styau/Figure` |
| `band_profiles_obs_vs_sim.pdf` | 5.2 | `analysis/boss/scripts/plot_observed_figures.py`, `fig3` | **external** BOSS accumulators as above, plus `analysis/sim/results/hod_pairs/stack_rperp{5,10,20}_matched.{csv,json}` and the mock single; **external** `stack_rperp{5,10,20}_blocks.npz` (~0.9 MB each, kept out of git with the other accumulators) | `PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python analysis/boss/scripts/plot_observed_figures.py --paper --only 3 --output-dir paper/styau/Figure` |
| `filament_snr_vs_smoothing.pdf` | 5.4 | `analysis/boss/scripts/plot_smoothing_scan.py` | **external** `acc_pairs_galaxy_{5,10,20}{,_fwhm6,_fwhm4}_*.npz` and the matching single accumulators; smoothings are discovered from whatever is on disk | `PYTHONPATH=lib python analysis/boss/scripts/plot_smoothing_scan.py --paper --output-dir paper/styau/Figure` |

Paper figures should carry no explanatory text inside the image: no
`fig.suptitle` notes and no in-panel annotation paragraphs. Panel titles are
labels only, for example `(a) Mock pair stack`, and the explanation goes in
the caption. Diagnostic versions of the same plots can keep their notes;
give the plotting script a `--paper` flag that switches them off.

Inputs marked **external** are not in the repository. The BOSS jackknife
accumulators are ~45 MB each and several are needed per figure, which is past
what belongs in git; they live only on the analysis machine today. Those
figures are therefore reproducible from a checked-in script and a checked-in
command, but not yet from a checked-in input. Depositing the accumulators
with a DOI and citing it here is the intended fix.

The plotting scripts write `.png` and `.pdf` side by side. Only the `.pdf`
belongs in `styau/Figure/`; the `.png` is for reading in a terminal or a
browser and is left in `output/plots/`. Run the diagnostic version (no
`--paper`) into `output/plots/` and the paper version into
`styau/Figure/`, then delete the stray `.png` the second run leaves behind.

## Producing `band_profiles_obs_vs_sim.pdf` (Section 5.2)

**What it shows.** One column per separation (5, 10, 20 h^-1 Mpc). Upper
panel: the observed band profile B(X) = central-band mean minus off-centre
band mean of the BOSS filament map, as points with jackknife error bars,
overlaid on the mock's B(X) as a line with a shaded box-jackknife band.
Dashed verticals at X = +/- d/2, shaded bridge window |X| <= 0.35 d, zero
line. Lower strip: (BOSS - mock) / sigma_BOSS. Units 10^-4 on the y axis.
No suptitle, no in-panel text; the caption in
`3_methods_results_draft.tex` (label `fig:band_profiles_obs_sim`) carries the
explanation. This is panel (d) of `plot_obs_vs_sim_filament.py` laid out
across the three separations with a mock error band and a residual strip
added.

**Why it is needed.** It is the only figure that compares data with the
LCDM prediction as a function of position along the pair axis, so it
separates a bridge centred at X = 0 from halo leakage near X = +/- d/2, and
it shows the separation dependence that Section 1 promises.

**Step 1 -- mock per-block stacks.** Mostly already done. The
`stack_rperp*_matched.json` files hold only the bridge statistic's jackknife
error, but the per-block maps exist on disk for two of the three separations
and match their stacks exactly:

| separation | blocks file | n_pairs / r_par vs the matched JSON |
|---|---|---|
| 5 | `stack_rperp5_blocks.npz` | regenerated here; the pre-existing `stack_rperp5_rpar10_blocks.npz` is the r_par <= 10 cut (487,057 pairs), not the r_par <= 5 the matched stack uses (319,772) |
| 10 | `stack_rperp10_blocks.npz` | 655,957 / r_par <= 10 -- identical |
| 20 | `stack_rperp20_blocks.npz` | 929,157 / r_par <= 10 -- identical |

Both existing files reconstruct their saved stack to 3e-8 of peak
(`sums.sum(0) / counts.sum()`, then `reflect_symmetrize_map`), so they are the
same run and can be used directly. Note the naming: `<stack stem>_blocks.npz`,
following what `jackknife_pair_stack.py` already writes -- not
`blocks_rperp*_matched.npz`.

Only the 5 h^-1 Mpc case needs a re-run (about 5 minutes: 25 blocks at ~8 s
each, plus reading the 346 MB pair catalogue). Send `--output` and
`--stack-output` somewhere scratch unless you mean to overwrite the committed
stack; only `--blocks-output` is new:

```
PYTHONPATH=lib:analysis/sim python analysis/sim/jackknife_pair_stack.py \
    --pairs   analysis/sim/results/hod_pairs/pairs_hodnmatch50_rperp5_rsd50.csv \
    --kappa-map analysis/sim/results/kappa_map_l0p1_s8arcmin.float32 \
    --single  analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv \
    --rperp-center 5 --rperp-min 4 --rperp-max 6 \
    --rpar-max 5 --rpar-space redshift \
    --output       /tmp/stack_rperp5_matched.json \
    --stack-output /tmp/stack_rperp5_matched.csv \
    --blocks-output analysis/sim/results/hod_pairs/stack_rperp5_blocks.npz
```

**The r_perp = 20 mock bin does not match BOSS.** Measured from the pair
catalogues themselves (`r_perp` column, min/max):

| separation | mock bin | BOSS bin | |
|---|---|---|---|
| 5 | 4-6 | 4-6 | match |
| 10 | 9-11 | 9-11 | match |
| 20 | **19-21** | **18-22** | **mismatch** |

This is not only a risk for the new figure: `stack_rperp20_matched.csv` is
built on that same 19-21 catalogue, is already committed, and already feeds
the committed `obs_vs_sim_quadrupole_maps.pdf`. `sim_diff()` asserts that
`rpar_max` matches BOSS but nothing asserts the r_perp bin, which is how it
got through. Fixing it means regenerating
`pairs_hodnmatch50_rperp20_rsd50.csv` with `find_pairs_sim.py --rperp-min 18
--rperp-max 22`, re-stacking, and re-issuing both that stack and the
quadrupole figure. Until then the 20 h^-1 Mpc column compares a narrower mock
bin against the wider BOSS one, and the caption says so.

**Step 2 -- add `fig3` to `plot_observed_figures.py`.** Follow the pattern of
`fig1`/`fig2` and the `--only` / `--paper` switches already there. Per
separation key:

1. BOSS side: `t = build_terms(...)` as in `fig2`. Profile
   `boss_prof = band_profile(t["filament"], axis)`; error from the
   leave-one-out maps exactly as `plot_obs_vs_sim_filament.py` lines
   127-131: `loo = [band_profile(m.reshape(shape), axis) for m in
   t["loo"]["filament"]]`, then `jackknife_error(loo)`. Do not derive the
   error from per-pixel errors; the 8 arcmin beam correlates neighbouring
   pixels (README section above, and Section 3.6 of the paper).
2. Mock side: `sim_fil, axis = sim_diff(key)` (already in the script -- it
   returns the map *and* its axis) gives the mock filament map; `sim_prof = band_profile(sim_fil, axis)`. For the
   band: load `stack_rperp{key}_blocks.npz`, form each leave-one-out
   map as `(sum_total - sum_b) / (n_total - n_b)`, subtract the *same*
   `two_halo_template` control used for the full mock stack (held fixed, as
   in `jackknife_pair_stack.py`), take `band_profile` of each, and apply
   `jackknife_error` with the (K-1)/K factor. K is the number of non-empty
   blocks recorded in the JSON as `n_blocks`.
3. Residual strip: `(boss_prof - sim_prof) / boss_err`.
4. Layout: `plt.subplots(2, 3, sharex="col", height_ratios=[3, 1])`,
   x range +/- 2.5 d per column (the `zoom` used elsewhere), y axis in
   units of 10^-4, `figure*` width (about 7 in wide). BOSS as points with
   error bars, mock as a line with `fill_between` band, colours as in the
   other figures (BLUE for BOSS, ORANGE for mock).
5. Save as `band_profiles_obs_vs_sim.{pdf,png}` via the script's `save()`.

**Step 3 -- consistency checks before using the figure.**

- The mean of `boss_prof` over the bridge window must equal the observed
  bridge excess in `analysis/boss/results/filament_jackknife_BOSS_North_South.csv`
  from `combine_filament_jackknife.py` (columns `filament_bridge_excess`,
  `filament_bridge_err`, `significance`), and its windowed jackknife error must
  match the table error to within a few per cent (the 2D box mean and the
  1D window differ slightly in pixel weighting; a warning for this already
  exists in `plot_separation_summary.py`).
- The windowed mean of `sim_prof` is the fixed-band mock prediction. Copy
  it, with its jackknife error from step 2, into the Mock column of
  `tab:bridge` in `3_methods_results_draft.tex`, replacing the `% CHECK`
  values (6, 3.4, 1.5), which came from the separation-scaled bands of
  `summarize_sim_sensitivity.py` and are not the same statistic. Also
  update the two sentences quoting them in `0_abstract_draft.tex` and
  `4_summary_draft.tex`.
- The residual strip should be within about +/- 2 everywhere if the paper's
  "consistent with the prediction" statement holds; note any bin outside
  that in the text.

**Step 4 -- done.** Implemented as `fig3`; the manifest row above carries
the command. The three `stack_rperp*_blocks.npz` are about 0.9 MB each and
are **not** committed -- they are per-block accumulators and belong with the
other accumulators in whatever deposit those go to, rather than being the one
exception in git. Regenerating them is cheap: `jackknife_pair_stack.py` with
`--blocks-output`, about 5 minutes per separation (25 blocks at ~8 s each plus
reading the pair catalogue). `fig3` raises a `FileNotFoundError` naming the
missing file and pointing back to step 1 if they are absent.

Measured when the figure was first produced, for reference:

- Step 3's first check is exact, not approximate. The bridge-window mean of
  `boss_prof` reproduces `filament_jackknife_BOSS_North_South.csv` to the
  printed digits at all three separations: 5.378 +/- 9.081, 11.908 +/- 5.101,
  1.236 +/- 2.221 (x 1e-4).
- The mock's box-jackknife error is 3-5 per cent of the BOSS error, so the
  shaded band is drawn but is far thinner than the observed error bars. That
  is the real ratio, not a plotting failure.
- Fixed-band mock bridge predictions, for the Mock column of `tab:bridge`:
  **6.32 +/- 0.36**, **4.91 +/- 0.21**, **2.18 +/- 0.12** (x 1e-4) at 5, 10,
  20 (errors corrected, and the 20 value is now the 18-22 bin -- see
  acceptance check 3).  Applied to `tab:bridge`, the abstract and the
  summary. The `% CHECK` values now in the table (6, 3.4, 1.5) came from
  `summarize_sim_sensitivity.py`'s separation-scaled bands, a different
  statistic; 10 h^-1 Mpc differs most (4.91 against 3.4).
- Residual bins beyond 2 sigma: 0 of 33 at 5, **6 of 45 at 10**, 2 of 65 at
  20. These are plotting diagnostics only: adjacent 1 h^-1 Mpc bins are
  correlated by the 3.3 h^-1 Mpc beam, so the count has no binomial Gaussian
  interpretation. Use the covariance-based comparison specified above for
  inference.

## Numbers quoted in the text

`Figure/footprint_stats.txt` is the source for the galaxy counts, median and
weighted-mean redshift, Planck mask sky fraction, and unmasked fractions in
Section 2. The bridge values, jackknife errors, template amplitudes,
line-of-sight sweep numbers, and mock predictions in Sections 4 to 6 are not
yet backed by committed result files; see the `% CHECK` and `% TODO`
markers. The intended fix is a generated `numbers.tex` of macros written by
the pipeline and `\input` by `main.tex`, so that text numbers cannot drift
from the results.
