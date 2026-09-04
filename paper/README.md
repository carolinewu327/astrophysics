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
| `band_profiles_obs_vs_sim.pdf` | 5.2 | not yet committed | BOSS band profile with jackknife errors, mock band profile with box-jackknife band | -- |
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

## Numbers quoted in the text

`Figure/footprint_stats.txt` is the source for the galaxy counts, median and
weighted-mean redshift, Planck mask sky fraction, and unmasked fractions in
Section 2. The bridge values, jackknife errors, template amplitudes,
line-of-sight sweep numbers, and mock predictions in Sections 4 to 6 are not
yet backed by committed result files; see the `% CHECK` and `% TODO`
markers. The intended fix is a generated `numbers.tex` of macros written by
the pipeline and `\input` by `main.tex`, so that text numbers cannot drift
from the results.
