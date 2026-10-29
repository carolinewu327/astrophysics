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
| `band_fit_profiles_obs_vs_sim.pdf` | new (number set in `main.tex`) | 22-element central and off-centre band vectors, BOSS against the mock, with the full-covariance amplitude fit | done (separate command, below) |

Figure 8 has its own script. It needs the 6 and 4 arcmin galaxy-pair stacks,
which are made once (see "Smoothing scan" in the wide-bin README):

```bash
python analysis/boss/scripts/widebin_smoothing_scan.py --overwrite
```

It also writes the numbers to
`analysis/boss/results/widebin_rpar5/smoothing_scan_wide_BOSS.csv`. The random
stacks are the 8 arcmin ones at every smoothing, as in the submitted figure;
the caption or text should say so.

The band-fit figure (added 2026-10-29 at Zheng's request) has its own script.
It redoes the fit with the same code as `widebin_vs_archived.py` and stops if
A, sigma_A or either chi^2 differs from
`widebin_rpar5/comparison/bridge_study_archived_vs_wide.csv`:

```bash
python analysis/boss/scripts/widebin_band_fit_figure.py --paper --overwrite
```

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
| Mock prediction, observed minus mock, chi^2 before and after the amplitude fit, fitted amplitudes, the free constant with its error, noise floor | `comparison/bridge_study_archived_vs_wide.csv` (rows with `version = wide`) |
| Pooled bridge excess: inverse-variance mean with its error neglecting and including the covariance between bins, the bin-to-bin correlations, and the minimum-variance combination | `comparison/bridge_study_pooled.csv` |
| Residual at the mean galaxy position and at the pair midpoint, with jackknife errors, data and mock; pair-stack convergence at the midpoint | `comparison/paper_figures/residual_at_galaxies.csv` (`widebin_paper_figures.py --residuals-only`) |
| Jackknife regions: 287, with 10.5 per cent RMS scatter in occupancy | `analysis/boss/results/jk/regions_BOSS_North_South_nside10.npz` (18 KB; read `region_of_pix` and `footprint_counts`) |
| Independent-pixel error factor, 2.24 / 2.25 / 2.62 for single-galaxy band means at 10-20, 15-25, 20-40 | `analysis/boss/results/profile_bands_corrected_scw_rand_scw_frac100_BOSS_North_South_joint_linear.csv`, columns `kappa_err` over `kappa_err_diag_only` (from `combine_jackknife.py`, command in `CLAUDE.md`) |
| Bridge excess as the window moves away from the galaxies | `validation/halo_cut/halo_cut_scan.csv` |

## Tests quoted in the text

One row per test. Commands run from the repository root with the setup above.
Only the small result tables are in git; inputs, intermediate files and how to
regenerate each test are in `analysis/boss/results/widebin_rpar5/README.md`.

| Test | Where in `main.tex` | Command | Result file | Status |
|---|---|---|---|---|
| Non-physical pair null test (is the control unbiased?) | "Expected signal from the simulation", the paragraph "To test whether the control itself generates a ridge..."; appendix, the control item | `python analysis/sim/widebin_null_test.py --bin g6_15` (also `g4_6`, `g15_25`) | `analysis/sim/results/hod_pairs/widebin_rpar5/nulltest/nonphysical_pair_test_<bin>.csv` | done for all three bins |
| Smoothing scan | "Robustness tests", Figure 8 | `python analysis/boss/scripts/widebin_smoothing_scan.py` | `analysis/boss/results/widebin_rpar5/smoothing_scan_wide_BOSS.csv` | done |
| Halo-cut scan (bridge as the window moves away from the galaxies) | "Robustness tests", the paragraph "Bridge window"; one sentence each in the Summary and in the appendix item on beam leakage | `python analysis/boss/scripts/widebin_halo_cut_scan.py` | `analysis/boss/results/widebin_rpar5/validation/halo_cut/halo_cut_scan.csv` | done, both binnings |
| Line-of-sight cut, part A: completeness and interloper fraction of the selection | "Expected signal from the simulation", the paragraph "We also use the mock ... to assess the line-of-sight selection"; appendix, the line-of-sight item; Introduction, fourth design point | `python analysis/sim/widebin_los_purity.py` (after the null test, which builds its input) | `analysis/sim/results/hod_pairs/widebin_rpar5/los_purity_wide.csv` | done |
| Line-of-sight cut, part B: bridge signal versus the cut | "Expected signal", one sentence labelled as the narrower-bin analysis, with a `% TODO (wide bins)`; a second `% TODO (wide bins)` in the appendix | `python analysis/sim/widebin_los_signal.py` (after the null test, which builds its input; about 40 minutes) | `analysis/sim/results/hod_pairs/widebin_rpar5/los_signal/los_signal_wide.csv` | done 2026-10-11; in the text only as two commented-out `% OPTIONAL` paragraphs (Expected signal and appendix) for Zheng to decide |
| Sensitivity grid | not in the text: the paragraph was in a block removed from the submitted version | not re-run for the wide bins | - | open only if the paragraph is restored |

### Non-physical pair null test (2026-10-10)

Mock pairs in each bin, split by real-space line-of-sight separation, with the
separation-averaged control. Bridge excess after the control in units of 1e-4,
errors from the 25-block jackknife. Uses a seeded random 25 per cent of the
mock galaxies (one sixteenth of the pairs); the method and the reason are in
the script's docstring.

| Line-of-sight band (h^-1 Mpc) | 4-6 | 6-15 | 15-25 |
|---|---:|---:|---:|
| 0-20 (associated) | 6.19 +/- 0.78 | 4.09 +/- 0.34 | 1.55 +/- 0.22 |
| 20-40 | 3.36 +/- 1.43 | 1.91 +/- 0.37 | 1.17 +/- 0.17 |
| 40-60 | 0.89 +/- 1.26 | -0.16 +/- 0.31 | 0.16 +/- 0.12 |
| 60-90 | -0.96 +/- 0.77 | -0.20 +/- 0.29 | 0.16 +/- 0.13 |
| 90-120 | -1.62 +/- 0.84 | -0.33 +/- 0.25 | -0.10 +/- 0.15 |
| Largest null-band deviation | -1.9 sigma | -1.3 sigma | +1.4 sigma |
| BOSS bridge error, for scale | 9.0 | 3.9 | 2.3 |

Pair counts and the bridge before the control are in the three CSVs.

Reading: beyond 40 h^-1 Mpc, where pairs are not physically associated, the
control leaves a residual consistent with zero in every bin; no null band
reaches 2 sigma. The sign is not the same across bins, so there is no
evidence of a common bias. In each bin the residual does fall slightly with
line-of-sight distance; the bands share one simulated map, so that is not
three independent hints. Any control bias is at most about a fifth of the
BOSS error (4-6) and under a tenth at 6-15 and 15-25. The 20-40 band still
carries real signal in all three bins, so "no bridge" holds beyond 40, not
beyond 20.

### Halo-cut scan (2026-10-10)

Bridge excess as the window |X| <= w shrinks from the default (0.35 R) to the
pair centre alone. Units 1e-4; errors are the full jackknife error of each
windowed mean. "Gap" is the distance from the window edge to the nearest
galaxy in the bin. The default rows reproduce the main result tables exactly.

| Bin | Window, gap | BOSS | S/N | Mock |
|---|---|---:|---:|---:|
| paper 9-11 | default 3.5, gap 1.0 | 11.91 +/- 5.10 | 2.3 | 4.91 |
| | 2, gap 2.5 | 9.65 +/- 5.72 | 1.7 | 5.14 |
| | 1, gap 3.5 | 8.26 +/- 6.66 | 1.2 | 5.21 |
| | 0, gap 4.5 | 6.22 +/- 7.47 | 0.8 | 4.92 |
| wide 6-15 | default 3.5, gap -0.5 | 9.23 +/- 3.92 | 2.4 | 4.91 |
| | 2, gap 1.0 | 7.55 +/- 4.32 | 1.7 | 4.86 |
| | 1, gap 2.0 | 5.69 +/- 4.82 | 1.2 | 4.87 |
| | 0, gap 3.0 | 4.62 +/- 5.25 | 0.9 | 4.77 |
| paper 18-22 | default 7, gap 2.0 | 1.24 +/- 2.22 | 0.6 | 2.18 |
| | 0, gap 9.0 | 2.36 +/- 4.59 | 0.5 | 2.21 |
| wide 15-25 | default 7, gap 0.5 | 1.63 +/- 2.30 | 0.7 | 2.18 |
| | 0, gap 7.5 | -1.64 +/- 4.31 | -0.4 | 2.10 |

The 4-6 bin (both binnings), the intermediate windows and the 6-10 / 10-15
split are in the CSV.

Reading: at R = 10 the BOSS excess falls by about half as the window closes
to the pair centre, in the paper's bin and the wide bin alike, while the mock
stays flat at about 4.9. The part of the 2.3-2.4 sigma that lies near the
galaxies is therefore not something the mock's filament produces; what is
left at the centre matches the mock in size but is under 1 sigma. The paper's
9-11 bin shows this even though its window ends 1 h^-1 Mpc short of the
nearest galaxy, which is inside the 3.3 h^-1 Mpc beam. The R = 20 bins are
consistent with zero and with the mock at every window.

### Line-of-sight selection: completeness and interlopers (2026-10-10)

Mock, same seeded quarter of the galaxies as the null test. A pair counts as
truly associated if its true line-of-sight separation is at most 20 h^-1 Mpc.
Completeness: share of the associated pairs that pass the redshift-space cut.
Interlopers: share of the selected pairs that are not associated.

| Cut on redshift-space r_par | 4-6 | 6-15 | 15-25 |
|---|---:|---:|---:|
| <= 5 (used): completeness | 48.1% | 35.8% | 28.9% |
| <= 5 (used): interlopers | 0.1% | 0.1% | 0.1% |
| <= 10: completeness / interlopers | 73.0% / 0.5% | 63.4% / 0.6% | 56.0% / 0.6% |
| <= 20: completeness / interlopers | 97.7% / 8.2% | 96.5% / 10.0% | 94.8% / 10.6% |
| <= 30: completeness / interlopers | 100% / 25.9% | 100% / 31.1% | 100% / 34.1% |
| Selected at <= 5: median / 90th / 99th percentile of the true separation | 3.3 / 8.5 / 14.9 | 3.6 / 9.0 / 15.4 | 3.6 / 8.9 / 15.2 |

Checks: pair counts at <= 5 times 16 match the full-sample mock stacks to 1 per
cent; no selected pair has a true separation beyond 60, well inside the
catalogue's 120.

Reading: under the 20 h^-1 Mpc line-of-sight criterion, the r_par <= 5
selection is pure but incomplete. Almost none of the selected pairs have a
true line-of-sight separation above 20 h^-1 Mpc, and it keeps between a half
and under a third of the pairs that lie within 20 h^-1 Mpc in real space. This
shows the selected pairs are close along the line of sight; it does not show
that they share a filament. It is a
conservative cut that follows earlier work; it is not what the 90 per cent
completeness / 10 per cent interloper rule in
`notes/phase3_los_preregistration.md` would pick (that rule lands near 20), so
the text must not say the cut was chosen by that rule. The completeness also
depends on the 20 h^-1 Mpc definition: at 15-25 the associated sample is close
to what unclustered pairs would give (5/20 = 25 per cent), so part of it is
chance alignment within the definition, not lost filaments.

### Line-of-sight cut: mock bridge signal versus the cut (2026-10-11)

Mock, same seeded quarter of the galaxies, separation-averaged control. S is
the bridge excess in units of 1e-4. The ratio to the cut of 5 has a paired
jackknife error (nested cuts, same 25 blocks).

| Bin | Cut | Pairs (ratio to cut 5) | S | S / S(5) | Upper bound on S/N gain |
|---|---|---:|---:|---:|---:|
| 4-6 | 5 | 20,184 | 4.84 +/- 1.39 | 1 | 1 |
| | 10 | 30,775 (1.52) | 6.61 +/- 0.94 | 1.37 +/- 0.26 | 1.69 |
| | 20 | 44,631 (2.21) | 6.96 +/- 0.73 | 1.44 +/- 0.36 | 2.14 |
| 6-15 | 5 | 105,618 | 5.14 +/- 0.52 | 1 | 1 |
| | 10 | 187,880 (1.78) | 5.12 +/- 0.42 | 1.00 +/- 0.06 | 1.33 |
| | 20 | 315,530 (2.99) | 4.35 +/- 0.31 | 0.85 +/- 0.07 | 1.47 |
| 15-25 | 5 | 149,296 | 1.83 +/- 0.29 | 1 | 1 |
| | 10 | 290,375 (1.94) | 1.67 +/- 0.26 | 0.91 +/- 0.09 | 1.27 |
| | 20 | 546,516 (3.66) | 1.69 +/- 0.21 | 0.92 +/- 0.10 | 1.76 |

The last column is S/S(5) times the square root of the pair ratio: the S/N
gain if the data error fell as 1/sqrt(N). It is an upper bound, because pairs
added by deepening the cut lie on the same sky as the existing ones and share
their Planck noise; the mock has no such noise and cannot measure the true gain.

Checks: the cut-5 values agree with the full-sample mock stacks (5.97, 4.91,
2.18) within the quarter-sample errors; for every stack the bridge rebuilt from
the saved blocks matches the stacker's own value.

Reading: widening the cut from 5 to 10 does not dilute the mock signal at 6-15
or 15-25 and raises it at 4-6, while the pair count grows by 1.5 to 1.9. At 20
the 6-15 signal falls by about 15 per cent. So r_par <= 5 is safe but not
optimal: a cut of 10 could gain up to about 30 per cent in S/N at 6-15, and
realistically less. Testing it on the data needs new galaxy-pair stacks and a
new random-pair run at r_par <= 10.

## Still to do

1. ~~Figure 8.~~ Done 2026-10-09.
2. ~~Copy the `.tex` sections and the Figure 1 files from `styau/`, then
   update the text for the wide bins.~~ Done 2026-10-10. `main.tex` is a copy
   of `styau/submitted.tex` with the dead blocks stripped; every section
   (abstract, Introduction, Methods, expected signal, Results, robustness
   tests, Summary, appendix) has been updated for the wide bins. It has not
   been compiled since.
3. Mock tests: the non-physical pair null test and the line-of-sight purity
   table are done for all three bins and are in the text. The
   signal-versus-cut sweep is run for the wide bins (above); its two paragraphs
   are commented out in the text as `% OPTIONAL`, since showing them invites
   the question of why the cut is 5.
   Open question for Zheng: is a data run at
   r_par <= 10 worth a second GCP random run? The fitted-control test is
   still at the narrower bins, with a `% TODO (wide bins)` in the text.
   The sensitivity grid is not in the text.
4. ~~Residual at the galaxy positions.~~ Resolved 2026-10-11.
   `residual_at_galaxies.csv` now has jackknife errors and a defined position
   (half the mean pair separation, interpolated on the pair axis); regenerate
   with `python analysis/boss/scripts/widebin_paper_figures.py --residuals-only --overwrite`.
   Wide bins, observed (mock), units 1e-3: 4-6 -0.34 +/- 0.78 (0.15);
   6-15 0.68 +/- 0.57 (0.12); 15-25 0.55 +/- 0.51 (0.07). All within 1.2
   sigma of zero, so the sentence in Section 5.1 stands. The earlier mismatch
   with the `kappa_filament_*` maps was a convention difference, not an error:
   those maps are not reflection-symmetrized, so +X and -X differ (0.64 and
   0.87 at 6-15); the table uses the symmetrized map, which is their mean.
5. ~~Halo-cut script.~~ Done 2026-10-10 (above). The "Bridge window"
   paragraph in `main.tex` quotes `halo_cut_scan.csv`, including the mock
   in the same windows; its `% CHECK` and `% TODO` are removed. Still open:
   the significance of the fall from the full window to the pair centre is
   not quantified (the nested windows are correlated).
