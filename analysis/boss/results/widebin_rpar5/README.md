# Wide-bin pair study: how it was run and how to redo it

The submitted paper used narrow r_perp bins (4-6, 9-11, 18-22). This study
repeats it with wide bins to see whether the signal-to-noise improves.

| | Paper | Wide-bin study |
|---|---|---|
| r_perp bins | 4-6, 9-11, 18-22 | 4-6, 6-15, 15-25 (splits 6-10, 10-15) |
| r_par cut | 5, 10, 10 | 5 everywhere |
| Pair weight | w1 w2 | w1 w2 / Sigma_crit^2 |
| Control | single stack shifted to +/- R/2 | same, averaged over each pair's actual separation |
| Randoms | one run per bin | one run in 1 h^-1 Mpc sub-bins over 3-25, pooled per bin |

Everything lives in directories named `widebin_rpar5/`. No step writes to the
archived results, and every script stops rather than replace an existing file
unless you pass `--overwrite`.

## Setup

Run every command from the repository root.

```bash
alias python=/opt/miniconda3/envs/astro/bin/python
export PYTHONPATH=lib:analysis/boss/scripts:analysis/sim
W=analysis/boss/results/widebin_rpar5              # BOSS products (this directory)
GW=data/paircatalogs/BOSS/widebin_rpar5            # galaxy pair catalogs
HW=analysis/sim/results/hod_pairs/widebin_rpar5    # mock products
TESS=analysis/boss/results/jk/regions_BOSS_North_South_nside10.npz
RLABEL=random_w3_25_rpar5_scw_frac100
```

Tests: `python -m unittest discover -s tests`

## Part 1: the measurement

Steps 1-4 are independent of each other. Step 5 needs all of them.

### Step 1. Galaxy pair catalogs (local)

One catalog per hemisphere covering r_perp 3-25, r_par <= 5. Every bin is cut
from these two files.

```bash
for REG in North South; do
  python analysis/boss/scripts/find_pairs.py \
    --dataset BOSS --region $REG --catalog-type galaxy \
    --rpar 5 --rperp-min 3 --rperp-max 25 --n-processes 6 --output-dir $GW
done
```

Output: `$GW/galaxy_pairs_BOSS_{North,South}_5.0_3.0_25.0hmpc.csv` (3,121,653 pairs in total).

### Step 2. Galaxy pair stacks, one per bin (local, a few minutes each)

```bash
for G in 4_6 6_15 15_25 6_10 10_15; do
  LO=${G%_*}; HI=${G#*_}
  HALF=$([ "$HI" = 25 ] || echo --rperp-half-open)   # [lo, hi) except the last bin
  python analysis/boss/scripts/stack_pairs_jk.py \
    --pair-catalogs $GW/galaxy_pairs_BOSS_North_5.0_3.0_25.0hmpc.csv,$GW/galaxy_pairs_BOSS_South_5.0_3.0_25.0hmpc.csv \
    --label galaxy_g${G}_rpar5_scw --rperp-min $LO --rperp-max $HI $HALF \
    --sigma-crit-weight --fine-bin-width 0.5 --n-processes 6 \
    --output-dir $W --jk-regions-path $TESS \
    2>&1 | tee $W/logs/galaxy_g${G}_rpar5_scw.log
done
```

Output per bin: `$W/jk/acc_pairs_galaxy_g{lo}_{hi}_rpar5_scw_BOSS_North_South.npz`
(per-jackknife-region sums, about 45 MB) and the stacked map
`$W/kappa_pairs_galaxy_g{lo}_{hi}_rpar5_scw_BOSS_North_South_joint.csv`.

To try a new bin, add it to the list (for example `3_7 7_15`, which were also run).

### Step 3. Random pairs (GCP, about 40 hours; do not rerun unless needed)

One run per hemisphere, full random catalog, 22 sub-bins of 1 h^-1 Mpc over
r_perp 3-25. Any bin with integer edges inside 3-25 can be built from it, so a
new bin does not need a new random run. Full procedure and checks:
`notes/gcp_migration_plan.md`.

```bash
python analysis/boss/scripts/find_and_stack_pairs.py \
  --dataset BOSS --region South --catalog-type random --fraction 1.0 \
  --rpar 5 --rperp-bin-edges "$(seq -s, 3 25)" --sigma-crit-weight \
  --shard 0/1 --label $RLABEL --output-dir $W \
  --checkpoint-path $W/checkpoints/${RLABEL}_BOSS_South.npz \
  --checkpoint-interval 10 --n-processes 40 --chunk-size 5000
```

Repeat with `--region North`. Runtime on a 44-vCPU machine with 40 workers:
South 8.8 h, North 31.8 h.

Output: `$W/kappa_pairs_${RLABEL}_BOSS_{North,South}.npz` (3.5 MB each, the
product everything downstream reads), plus a `.csv` quick-look map and a
`.meta.json` with the full provenance.

### Step 4. Mock (local, about 75 minutes for all bins)

```bash
python analysis/sim/find_pairs_sim.py \
  --halos analysis/sim/results/galaxies_hod_nmatch_frac50.csv \
  --output $HW/pairs_hodnmatch50_rperp3_25_rsd5.csv \
  --rperp-min 3 --rperp-max 25 --rpar-max 5 --rpar-space redshift

# bin : bridge centre
for GC in 4_6:5 6_15:10 15_25:20 6_10:8 10_15:12.5; do
  G=${GC%:*}; C=${GC#*:}; LO=${G%_*}; HI=${G#*_}
  python analysis/sim/jackknife_pair_stack.py \
    --pairs $HW/pairs_hodnmatch50_rperp3_25_rsd5.csv \
    --kappa-map analysis/sim/results/kappa_map_l0p1_s8arcmin.float32 \
    --single analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv \
    --rperp-min $LO --rperp-max $HI $([ "$HI" = 25 ] || echo --rperp-half-open) \
    --rperp-center $C --rpar-max 5 --rpar-space redshift --control separation_averaged \
    --output $HW/stack_g${G}_rpar5_avgctrl.json \
    --stack-output $HW/stack_g${G}_rpar5_avgctrl_map.csv \
    --blocks-output $HW/stack_g${G}_rpar5_avgctrl_blocks.npz
done
```

Output per bin: `.json` (bridge excess and error), `_map.csv` (stacked map),
`_blocks.npz` (25 jackknife blocks). Summary of the run: `$HW/D3_summary.txt`.

### Step 5. Combine: filament maps and the bridge excess (local, under a minute)

```bash
python analysis/boss/scripts/combine_filament_jackknife.py \
  --groupings g4_6,g6_15,g15_25,g6_10,g10_15 \
  --random-subbins-label $RLABEL --control separation_averaged \
  --single-tag _scw --single-random-tag _scw_frac100 \
  --results-dir analysis/boss/results --output-dir $W
```

It also needs the two archived single-galaxy accumulators,
`analysis/boss/results/jk/acc_single_{galaxy_scw,random_scw_frac100}_BOSS_North_South.npz`.
These were not re-made for this study.

Output:
- `$W/filament_jackknife_wide_avgctrl_BOSS_North_South.csv`: **the main result**,
  one row per bin with pair count, mean r_perp, bridge excess, jackknife error and S/N.
- `$W/kappa_filament_*` and `$W/error_filament_*`: the residual map and its
  per-pixel error for each bin.

The script refuses to replace an existing summary. To run other bins beside the
main ones, give the run a tag, as was done for 3-7 and 7-15:
`--groupings g3_7,g7_15 --tag alt37` writes `filament_jackknife_wide_alt37_avgctrl_*.csv`.

The bridge window is `|X| <= 0.35 x centre`. Centres are set in
`BRIDGE_CENTERS` in `combine_filament_jackknife.py` (5, 10, 20, and 8 and 12.5
for the splits); `--bridge-centers` overrides them. The mock in step 4 must use
the same centres.

## Part 2: tables and charts in the write-up

Three scripts in `analysis/boss/scripts/` write to `$W/comparison/`. Each
recomputes the paper's bins and the wide bins through one code path, so the
"paper" column is a check: it must reproduce the manuscript.

```bash
python analysis/boss/scripts/widebin_vs_archived.py              # tables
python analysis/boss/scripts/widebin_control_profiles.py         # pair stack vs control
python analysis/boss/scripts/widebin_paper_figures.py --wide-only  # Figs 2-7
```

The same script makes the manuscript's figures. With `--paper` it writes
Figs 2-7 for the wide bins, in paper style, to `paper/sts/Figure/`
(see `paper/sts/README.md`):

```bash
python analysis/boss/scripts/widebin_paper_figures.py --paper --overwrite
```

| Output (in `$W/comparison/`) | In the write-up | Script |
|---|---|---|
| `bridge_study_archived_vs_wide.csv` | Table 1 paper vs wide; profile tests (chi^2, fitted amplitudes) | `widebin_vs_archived.py` |
| `bridge_study_pooled.csv` | Pooled bridge excess | same run |
| `band_profiles_archived_vs_wide.pdf` | (first version of the band-profile chart; replaced by Fig 7) | same run |
| `control_profiles_archived_vs_wide.pdf` | Pair stack vs control | `widebin_control_profiles.py` |
| `paper_figures/fig2_mock_pair_vs_control_wide.pdf` | Fig 2 | `widebin_paper_figures.py --wide-only` |
| `paper_figures/fig3_mock_residual_maps_wide.pdf` | Fig 3 | same run |
| `paper_figures/fig4_mock_band_profiles_wide.pdf` | Fig 4 | same run |
| `paper_figures/fig5_boss_single_and_pair_maps_wide.pdf` | Fig 5 | same run |
| `paper_figures/fig6_obs_vs_mock_residual_maps_wide.pdf` | Fig 6 | same run |
| `paper_figures/fig7_band_profiles_obs_vs_mock.pdf` | Fig 7 (both binnings) | same run |
| `paper_figures/residual_at_galaxies.csv` | Residual at the galaxy positions | same run |
| `paper_figures/fig{2..6}_*.pdf` (no `_wide`) | Not in the write-up: paper bins and wide bins in one image | `widebin_paper_figures.py` (no flag) |

Each script writes a `.png` beside every `.pdf`. Add `--overwrite` to redraw a
figure that already exists.

What these scripts read (none of it is re-stacked):

| Input | From | In git |
|---|---|---|
| `$W/jk/acc_pairs_galaxy_g*_rpar5_scw_*.npz` | step 2 | no (45 MB each) |
| `$W/kappa_pairs_${RLABEL}_BOSS_{North,South}.npz` | step 3 | yes |
| `$HW/stack_g*_rpar5_avgctrl{.json,_map.csv}` | step 4 | yes |
| `$HW/stack_g*_rpar5_avgctrl_blocks.npz` | step 4 | yes (the five main bins) |
| `analysis/boss/results/jk/acc_single_*_scw*_BOSS_North_South.npz` | archived | no |
| `analysis/boss/results/jk/acc_pairs_galaxy_{5,10,20_rpar10}_*.npz` | archived (paper bins) | no |
| `analysis/sim/results/hod_pairs/stack_rperp{5,10,20}_matched.*`, `stack_rperp{5,10,20}_blocks.npz` | archived (paper bins) | csv/json yes, blocks no |
| `analysis/sim/results/kappa_single_sim_hodnmatch_8arcmin_centered_g101.csv` | archived | yes |

## Smoothing scan (manuscript Figure 8)

Bridge S/N at 8, 6 and 4 arcmin. Only the galaxy-pair stacks are re-made at
each smoothing (about 13 minutes for all six). The single-galaxy stacks at 6
and 4 arcmin already exist in `analysis/boss/results/jk/`. The random stacks
are the 8 arcmin ones at every smoothing, as in the submitted figure.

```bash
for F in 6 4; do for G in 4_6 6_15 15_25; do
  LO=${G%_*}; HI=${G#*_}
  HALF=$([ "$HI" = 25 ] || echo --rperp-half-open)
  python analysis/boss/scripts/stack_pairs_jk.py \
    --pair-catalogs $GW/galaxy_pairs_BOSS_North_5.0_3.0_25.0hmpc.csv,$GW/galaxy_pairs_BOSS_South_5.0_3.0_25.0hmpc.csv \
    --label galaxy_g${G}_rpar5_scw_fwhm${F} --rperp-min $LO --rperp-max $HI $HALF \
    --sigma-crit-weight --fine-bin-width 0.5 --n-processes 6 --fwhm-arcmin $F \
    --output-dir $W --jk-regions-path $TESS \
    2>&1 | tee $W/logs/galaxy_g${G}_rpar5_scw_fwhm${F}.log
done; done

python analysis/boss/scripts/widebin_smoothing_scan.py
```

Output: `paper/sts/Figure/filament_snr_vs_smoothing.pdf` and
`$W/smoothing_scan_wide_BOSS.csv`. The 8 arcmin rows reproduce the main result
table exactly, which is the check that the script matches step 5.

## If you change something, what to rerun

| Change | Rerun |
|---|---|
| A figure's look (labels, colours, layout) | the one script in Part 2, with `--overwrite` |
| The bridge window or centres | step 4 (mock, new centres), step 5 with a new `--tag`, then Part 2 |
| A new bin with integer edges inside 3-25 | step 2 and step 4 for that bin, step 5 with a new `--tag`; add the bin to `ROWS` in `widebin_paper_figures.py` and `widebin_control_profiles.py` and to the bin list in `widebin_vs_archived.py` if it should appear in the charts |
| A bin outside 3-25, a different r_par cut, or different weights | everything, including the 40-hour random run |
| The control (for example back to single-separation) | step 5 with `--control nominal` and a new `--tag`; step 4 with `--control nominal` |
| The paper's own bins | `combine_filament_jackknife.py --separations 5,10,20` still reproduces the archived numbers exactly |

## Not yet scripted

`validation/halo_cut/` holds the halo-cut scan: the bridge excess recomputed
with the window shrunk away from the galaxies (6-15 falls from 2.4 to 0.9
sigma). It was run by hand and has no script in the repository. A script that
covers both the wide bins and the paper's bins is still to be written.
