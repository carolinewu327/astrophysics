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
| `band_profiles_obs_vs_sim.pdf` | 5.2 | planned as `fig3` in `analysis/boss/scripts/plot_observed_figures.py`; see "Producing band_profiles_obs_vs_sim.pdf" below | the BOSS accumulators above, plus `analysis/sim/results/hod_pairs/stack_rperp{5,10,20}_matched.{csv,json}`, `blocks_rperp{5,10,20}_matched.npz` (new), and the mock single | `PYTHONPATH=lib:analysis/boss/scripts:analysis/sim python analysis/boss/scripts/plot_observed_figures.py --paper --only 3 --output-dir paper/styau/Figure` |
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

**Step 1 -- mock per-block stacks (run once per separation, on the machine
with the mock maps).** The existing `stack_rperp*_matched.json` holds only
the bridge statistic's jackknife error, not per-block maps, so the mock
profile currently has no error band. Re-run the matched stacks with
`--blocks-output`:

```
PYTHONPATH=lib:analysis/sim python analysis/sim/jackknife_pair_stack.py \
    --pairs   analysis/sim/results/hod_pairs/pairs_hodnmatch50_rperp10_rsd50.csv \
    --kappa-map analysis/sim/results/kappa_map_l0p1_s8arcmin.float32 \
    --single  analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv \
    --rperp-center 10 --rpar-max 10 --rpar-space redshift \
    --output       analysis/sim/results/hod_pairs/stack_rperp10_matched.json \
    --stack-output analysis/sim/results/hod_pairs/stack_rperp10_matched.csv \
    --blocks-output analysis/sim/results/hod_pairs/blocks_rperp10_matched.npz
```

Repeat for 5 (`--rpar-max 5`) and 20 (`--rpar-max 10`). Check the input
file names against what is on disk; the pair-catalogue and map names above
are the ones the scripts default to. Two things to confirm in the JSON
before going on: the 20 h^-1 Mpc mock bin must be 18-22 to match the BOSS
catalogue (`run_sim_sensitivity.py` used 19-21; if the pair file was made
with 19-21, regenerate it with `find_pairs_sim.py --rperp-min 18
--rperp-max 22`), and `rpar_max_hmpc` must read 5, 10, 10.

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
2. Mock side: `sim_fil = sim_diff(key)` (already in the script) gives the
   mock filament map; `sim_prof = band_profile(sim_fil, axis)`. For the
   band: load `blocks_rperp{key}_matched.npz`, form each leave-one-out
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
  bridge excess in `analysis/boss/results/filament_bridge_*.csv` from
  `combine_filament_jackknife.py`, and its windowed jackknife error must
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

**Step 4 -- run and file.** Diagnostic version (with notes) to
`output/plots/`, paper version with `--paper` to `paper/styau/Figure/`,
delete the stray `.png` there, then update the manifest row above from
"planned" to the real command and commit script, PDF, and the three
`blocks_*.npz` files if they are small (about 2 MB each per the
`jackknife_pair_stack.py` docstring).

## Numbers quoted in the text

`Figure/footprint_stats.txt` is the source for the galaxy counts, median and
weighted-mean redshift, Planck mask sky fraction, and unmasked fractions in
Section 2. The bridge values, jackknife errors, template amplitudes,
line-of-sight sweep numbers, and mock predictions in Sections 4 to 6 are not
yet backed by committed result files; see the `% CHECK` and `% TODO`
markers. The intended fix is a generated `numbers.tex` of macros written by
the pipeline and `\input` by `main.tex`, so that text numbers cannot drift
from the results.
