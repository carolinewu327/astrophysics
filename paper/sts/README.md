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
| Bridge excess as the window moves away from the galaxies | `validation/halo_cut/halo_cut_scan.csv` |

## Tests quoted in the text

One row per test. Commands run from the repository root with the setup above.
Only the small result tables are in git; inputs, intermediate files and how to
regenerate each test are in `analysis/boss/results/widebin_rpar5/README.md`.

| Test | Where in `main.tex` | Command | Result file | Status |
|---|---|---|---|---|
| Non-physical pair null test (is the control unbiased?) | "Expected signal from the simulation", the sentence "The nonphysical mock pairs show no bridge..." and its `% TODO` | `python analysis/sim/widebin_null_test.py --bin g6_15` (also `g4_6`, `g15_25`) | `analysis/sim/results/hod_pairs/widebin_rpar5/nulltest/nonphysical_pair_test_<bin>.csv` | done for all three bins |
| Smoothing scan | "Robustness tests", Figure 8 | `python analysis/boss/scripts/widebin_smoothing_scan.py` | `analysis/boss/results/widebin_rpar5/smoothing_scan_wide_BOSS.csv` | done |
| Halo-cut scan (bridge as the window moves away from the galaxies) | summary, halo-leakage paragraph | `python analysis/boss/scripts/widebin_halo_cut_scan.py` | `analysis/boss/results/widebin_rpar5/validation/halo_cut/halo_cut_scan.csv` | done, both binnings |
| Line-of-sight sweep | "Expected signal", `% TODO` | not re-run for the wide bins | - | open |
| Sensitivity grid | Methods, mock tests | not re-run for the wide bins | - | open |

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

## Still to do

1. ~~Figure 8.~~ Done 2026-10-09.
2. Copy the `.tex` sections and the Figure 1 files from `styau/`, then update
   methods, captions, Table 1, abstract and summary for the wide bins.
3. Mock tests: the non-physical pair null test is done for all three bins
   (above). The line-of-sight sweep and the
   sensitivity grid are still at the old bins: re-run or say so in the text.
4. Section 5.1's "2e-3 observed against 0.9e-3 mock at the galaxies": the mock
   value could not be reproduced, even for the old bins.
5. ~~Halo-cut script.~~ Done 2026-10-10 (above). `main.tex` still has a
   `% CHECK` saying the scan has no script and a `% TODO` for the mock
   counterpart; both are now answered by `halo_cut_scan.csv`.
