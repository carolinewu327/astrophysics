# Paper, wide-bin version (`paper/sts/`)

The manuscript redone with wide r_perp bins, for journal publication. The
submitted version is `paper/styau/`; it is frozen and documented in
`paper/README.md`. Nothing here writes to `paper/styau/`.

## What differs from the submitted version

| | Submitted (`styau/`) | This version (`sts/`) |
|---|---|---|
| r_perp bins | 4-6, 9-11, 18-22 | 4-6, 6-15, 15-25 |
| r_par cut | 5, 10, 10 | 5 everywhere |
| Pair weight | w1 w2 | w1 w2 / Sigma_crit^2 |
| Control | single stack shifted to +/- R/2 | same, averaged over each pair's actual separation |
| Randoms | one run per bin | one run in 1 h^-1 Mpc sub-bins over 3-25, pooled per bin |

How the measurement was run, and how to redo any part of it:
`analysis/boss/results/widebin_rpar5/README.md`.

## Layout

```
paper/sts/
  README.md     this file
  main.tex      the whole manuscript in one file (started from the styau sections)
  ref.bib
  Figure/       figure PDFs, same file names as styau/Figure/
```

Rule: anything that appears in the paper has a script, or at least a written
command, that produces it. Figures are in the manifest below; numbers are in
the two tables after it.

## Figure manifest

Run from the repository root with
`PYTHONPATH=lib:analysis/boss/scripts:analysis/sim` and the `astro`
environment (`/opt/miniconda3/envs/astro/bin/python`).

Figures 2-7 come from one command:

```bash
python analysis/boss/scripts/widebin_paper_figures.py --paper --overwrite
```

| File (`sts/Figure/`) | Figure | Shows | Status |
|---|---|---|---|
| `footprint.pdf`, `zdist.pdf`, `footprint_stats.txt` | 1 | Footprint and n(z); do not depend on the binning | copied from `styau/Figure/` (command in `paper/README.md`) |
| `sim_pair_vs_control_maps_sep10.pdf` | 2 | Mock pair stack and its separation-averaged control, 6-15 | done |
| `sim_quadrupole_diff_maps.pdf` | 3 | Mock residual maps, three bins | done |
| `sim_band_profiles_all_seps.pdf` | 4 | Mock central band, off-centre band and difference | done |
| `boss_single_and_pair_maps.pdf` | 5 | BOSS single stack and three pair stacks | done |
| `obs_vs_sim_quadrupole_maps.pdf` | 6 | Observed and mock residual maps | done |
| `band_profiles_obs_vs_sim.pdf` | 7 | B(X), BOSS against the mock, with residual strip | done |
| `filament_snr_vs_smoothing.pdf` | 8 | Bridge S/N at 8, 6 and 4 arcmin | done (separate command, below) |

Figure 8 has its own script. It needs the 6 and 4 arcmin galaxy-pair stacks,
which are made once (see "Smoothing scan" in the wide-bin README):

```bash
python analysis/boss/scripts/widebin_smoothing_scan.py --overwrite
```

It also writes the numbers to
`analysis/boss/results/widebin_rpar5/smoothing_scan_wide_BOSS.csv`. The random
stacks are the 8 arcmin ones at every smoothing, as in the submitted figure;
the caption or text should say so.

File names are kept from the submitted version so the `\includegraphics`
lines carry over unchanged. `sim_pair_vs_control_maps_sep10.pdf` now shows the
6-15 bin despite the `sep10` in its name.

Figures are paper style only: panel titles are bare labels and the
explanation belongs in the caption. The script writes PDF only here.

Inputs the script reads are listed in the wide-bin README above. The galaxy
pair accumulators (45 MB per bin) are not in git, so the figures can be
redrawn only on the analysis machine.

## Numbers for the text and Table 1

All in `analysis/boss/results/widebin_rpar5/`:

| Quantity | File |
|---|---|
| Bridge excess, error, S/N, pair counts, mean r_perp | `filament_jackknife_wide_avgctrl_BOSS_North_South.csv` |
| Mock prediction, observed minus mock, chi^2, fitted amplitudes, noise floor | `comparison/bridge_study_archived_vs_wide.csv` (rows with `version = wide`) |
| Pooled bridge excess | `comparison/bridge_study_pooled.csv` |
| Bridge excess as the window moves away from the galaxies | `validation/halo_cut/bridge_window_scan.csv` |

## Tests quoted in the text

One row per test. Commands run from the repository root with the setup above.
Only the small result tables are in git; inputs, intermediate files and how to
regenerate each test are in `analysis/boss/results/widebin_rpar5/README.md`.

| Test | Where in `main.tex` | Command | Result file | Status |
|---|---|---|---|---|
| Non-physical pair null test (is the control unbiased?) | "Expected signal from the simulation", the sentence "The nonphysical mock pairs show no bridge..." and its `% TODO` | `python analysis/sim/widebin_null_test.py --bin g6_15` | `analysis/sim/results/hod_pairs/widebin_rpar5/nulltest/nonphysical_pair_test_g6_15.csv` | 6-15 done; 4-6 and 15-25 not run |
| Smoothing scan | "Robustness tests", Figure 8 | `python analysis/boss/scripts/widebin_smoothing_scan.py` | `analysis/boss/results/widebin_rpar5/smoothing_scan_wide_BOSS.csv` | done |
| Halo-cut scan (bridge as the window moves away from the galaxies) | summary, halo-leakage paragraph | none yet; run by hand | `analysis/boss/results/widebin_rpar5/validation/halo_cut/bridge_window_scan.csv` | numbers exist, script missing |
| Line-of-sight sweep | "Expected signal", `% TODO` | not re-run for the wide bins | - | open |
| Sensitivity grid | Methods, mock tests | not re-run for the wide bins | - | open |

### Non-physical pair null test, 6-15 bin (2026-10-10)

Mock pairs in the 6-15 bin, split by real-space line-of-sight separation, with
the separation-averaged control. Bridge excess in units of 1e-4, errors from
the 25-block jackknife. Uses a seeded random 25 per cent of the mock galaxies
(one sixteenth of the pairs); the method and the reason are in the script's
docstring.

| Line-of-sight band (h^-1 Mpc) | Pairs | Before control | After control |
|---|---:|---:|---:|
| 0-20 | 294,403 | 21.35 | 4.09 +/- 0.34 (12.2 sigma) |
| 20-40 | 225,799 | 18.43 | 1.91 +/- 0.37 (5.2 sigma) |
| 40-60 | 215,670 | 16.34 | -0.16 +/- 0.31 (-0.5 sigma) |
| 60-90 | 318,052 | 16.30 | -0.20 +/- 0.29 (-0.7 sigma) |
| 90-120 | 317,547 | 16.15 | -0.33 +/- 0.25 (-1.3 sigma) |

Reading: beyond 40 h^-1 Mpc, where pairs are not physically associated, the
control leaves a residual consistent with zero. The three null bands are all
slightly negative, as was the 40-60 band of the old control at 9-11
(-0.33 +/- 0.28); they share one simulated map, so they are not independent.
Any control bias is below about 0.3e-4, against a BOSS error of 3.9e-4.

## Still to do

1. ~~Figure 8.~~ Done 2026-10-09.
2. Copy the `.tex` sections and the Figure 1 files from `styau/`, then update
   methods, captions, Table 1, abstract and summary for the wide bins.
3. Mock tests: the non-physical pair null test is done for 6-15 (above);
   decide whether to run 4-6 and 15-25. The line-of-sight sweep and the
   sensitivity grid are still at the old bins: re-run or say so in the text.
4. Section 5.1's "2e-3 observed against 0.9e-3 mock at the galaxies": the mock
   value could not be reproduced, even for the old bins.
5. Report the halo-cut result (6-15 falls from 2.4 to 0.9 sigma at the pair
   centre) and give it a script; it was run by hand.
