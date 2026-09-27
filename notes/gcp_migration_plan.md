# GCP random-pair stacking runbook

Last verified: 2026-09-26 (wide-bin commands rehearsed locally at 1%; the
August commands below are the verified record of that run)

Project: `astrophysics-cwu` (`1025419119877`)

VM: `astro-pairs`, `us-central1-a`

This is the operational runbook for full BOSS random-pair stacks with
`analysis/boss/scripts/find_and_stack_pairs.py`. That script finds pairs and
stacks them directly; it must be used for full random catalogs because writing
billions of pairs to an intermediate CSV is not practical.

The most recent production run measured the `r_perp = 5 Mpc/h` bin
(`4 <= r_perp <= 6 Mpc/h`) with `r_par <= 10 Mpc/h` in BOSS North and South.
Both jobs completed, were bundled on GCP, downloaded, extracted into the local
repository, and checked.

**Next run (2026-09): the wide-bin random run.** One pass per region over
`r_perp = 3-25 Mpc/h` in 1 Mpc/h sub-bins, `r_par <= 5 Mpc/h`, with 1/Sigma_crit^2
pair weights, optionally sharded across VMs. Its steps are in the sections
marked **Wide-bin**: pilots (3), production and shard merge (4), completion
checks (6), bundle (7). Design and test record: `notes/wide_bin_implementation_plan.md`.

---

## 1. Current retained GCP setup

The VM is stopped, not deleted. Starting it again preserves the repository,
Conda environment, catalogs, Planck files, checkpoints, and results on its
persistent boot disk.

Verified configuration from the August 2026 run:

| Item | Configuration |
|---|---|
| Project | `astrophysics-cwu` |
| Zone | `us-central1-a` |
| VM | `astro-pairs` |
| Machine | `c3-standard-22` |
| CPU / memory | 22 vCPU / 88 GiB RAM |
| OS | Debian 13 cloud image; kernel reported `6.12.105+deb13-cloud-amd64` |
| Boot disk | 200 GB persistent disk |
| Production workers | 16 |
| Chunk size | 5,000 catalog objects |
| Python environment | Miniforge environment `astro` |
| Repository | `~/astrophysics` |

At the August launch the regional C3 CPU quota was 24, so a machine larger than
22 vCPUs was not available. The 22-vCPU VM was deliberately run with 16
workers to leave memory and system headroom. **The quota was raised to 128 C3
vCPUs in us-central1 on 2026-09-26**, enough for about five 22-vCPU VMs (e.g.
one per shard) or larger machines.

### Restart the retained VM

In Google Cloud Console:

1. Select project `astrophysics-cwu`.
2. Open **Compute Engine -> VM instances**.
3. Select `astro-pairs` and click **Start/Resume**.
4. When its status is **Running**, click **SSH**.

Do not create a new VM unless `astro-pairs` or its disk has been deleted. An
ephemeral external IP may change after a stop/start; Browser SSH is unaffected.

### Verify the retained environment

```bash
cd ~/astrophysics
git status --short
git log -1 --oneline

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate astro
python -c "import healpy, astropy, numpy, pandas, scipy; print('Environment OK')"

ls -lh \
  data/BOSS/random0_DR12v5_CMASS_North.fits.gz \
  data/BOSS/random0_DR12v5_CMASS_South.fits.gz \
  data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits \
  data/planck/COM_Lensing_4096_R3.00/mask.fits
```

Before pulling new code, commit and push relevant local changes. Then on GCP:

```bash
cd ~/astrophysics
git pull --ff-only
git rev-parse HEAD
python -m py_compile analysis/boss/scripts/find_and_stack_pairs.py \
  analysis/boss/scripts/merge_pair_shards.py lib/catalog.py
```

Record the commit hash in the run log or experiment notes.

---

## 2. Rebuilding the VM if necessary

Only use this section if the retained VM/disk no longer exists.

### VM

The tested size is `c3-standard-22` with a 200 GB persistent disk in
`us-central1-a`. The exact image family offered by the console can change; the
successful VM ran Debian 13. A current Debian or Ubuntu x86-64 image is fine.

The equivalent command is approximately:

```bash
gcloud compute instances create astro-pairs \
  --project=astrophysics-cwu \
  --zone=us-central1-a \
  --machine-type=c3-standard-22 \
  --image-family=debian-13 \
  --image-project=debian-cloud \
  --boot-disk-size=200GB \
  --boot-disk-type=pd-balanced
```

Confirm the currently available image family before using the command. In the
console, select a 200 GB balanced persistent disk unless the next job has a
demonstrated need for SSD throughput.

### Packages and Miniforge

```bash
sudo apt-get update
sudo apt-get install -y build-essential git wget curl tmux htop

wget -O /tmp/Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash /tmp/Miniforge3.sh -b -p "$HOME/miniforge3"
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda init bash

conda create -n astro -c conda-forge \
  python=3.12 numpy matplotlib astropy healpy pandas scipy tqdm -y
conda activate astro
python -c "import healpy, astropy, numpy, pandas, scipy; print('Environment OK')"
```

### Repository and directories

```bash
cd ~
git clone git@github.com:carolinewu327/astrophysics.git
cd ~/astrophysics
mkdir -p \
  data/BOSS \
  data/planck/COM_Lensing_4096_R3.00/MV \
  analysis/boss/results/checkpoints \
  analysis/boss/results/logs \
  analysis/boss/results/validation
```

If GitHub SSH authentication is not configured, clone over HTTPS instead.

### BOSS catalogs

```bash
cd ~/astrophysics/data/BOSS
wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_North.fits.gz
wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_South.fits.gz
```

The production logs should report approximately 30,113,342 North objects and
10.74 million South objects after preprocessing. A materially different count
means the catalog or selection should be checked before continuing.

### Planck files

Two files are required at these exact paths:

```text
data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits
data/planck/COM_Lensing_4096_R3.00/mask.fits
```

They can be uploaded independently with Browser SSH's **Upload file** button,
then moved into place, or copied from the local computer:

```bash
gcloud compute scp \
  /local/path/dat_klm.fits \
  astro-pairs:~/astrophysics/data/planck/COM_Lensing_4096_R3.00/MV/ \
  --project=astrophysics-cwu --zone=us-central1-a

gcloud compute scp \
  /local/path/mask.fits \
  astro-pairs:~/astrophysics/data/planck/COM_Lensing_4096_R3.00/ \
  --project=astrophysics-cwu --zone=us-central1-a
```

---

## 3. Mandatory preflight pilot

Always run the South 0.1% pilot after a code, environment, map, catalog, or VM
change. It exercises catalog loading, Planck ALM conversion and 8-arcmin
smoothing, pair finding, multiprocessing, stacking, checkpointing, and output
serialization.

```bash
cd ~/astrophysics
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate astro

mkdir -p \
  analysis/boss/results/validation \
  analysis/boss/results/checkpoints \
  analysis/boss/results/logs

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
  --dataset BOSS \
  --region South \
  --catalog-type random \
  --fraction 0.001 \
  --seed 12345 \
  --rpar 10 \
  --rperp-min 4 \
  --rperp-max 6 \
  --label random_5_rpar10_pilot_frac001 \
  --output-dir analysis/boss/results/validation \
  --checkpoint-path analysis/boss/results/checkpoints/random_5_rpar10_pilot_frac001_BOSS_South.npz \
  --checkpoint-interval 1 \
  --n-processes 4 \
  --chunk-size 500 \
  --overwrite
```

Verified pilot result from 2026-08-29:

- 10,740 catalog objects
- 22 chunks
- 92 pairs
- zero skipped pairs
- about 29 seconds total, including Planck loading
- final CSV and metadata JSON written successfully

The small number of pilot pairs makes its map noisy. This pilot tests the
pipeline, not the scientific signal.

### Wide-bin smoke test (R2): production flags at 0.1%

Same flags as production except fraction, seed and output location. Checks the
sub-bin routing, the 1/Σ_crit² weight and the new `.npz` output on the VM.

```bash
W=analysis/boss/results/widebin_rpar5
mkdir -p $W/validation $W/checkpoints $W/logs
EDGES=$(seq -s, 3 25)
LABEL=random_w3_25_rpar5_scw_pilot_frac001
STEM=${LABEL}_BOSS_South

PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
  --dataset BOSS --region South --catalog-type random \
  --fraction 0.001 --seed 12345 \
  --rpar 5 --rperp-bin-edges "$EDGES" --sigma-crit-weight \
  --label "$LABEL" --output-dir $W/validation \
  --checkpoint-path "$W/checkpoints/${STEM}.npz" \
  --checkpoint-interval 1 --n-processes 16 --chunk-size 5000 \
  2>&1 | tee "$W/logs/${STEM}.log"
```

Pass criteria: the log shows `1/Sigma_crit^2 pair weight: ON`, `r_perp sub-bins: 22`,
and `Pairs per r_perp sub-bin: [3,4): …` at the end; counts grow roughly with r⊥ (annulus area); the
`.npz` exists next to the CSV. No `--overwrite`: a rerun under the same name
stops instead of replacing it.

### Wide-bin performance pilot (R3)

As above with `--fraction 0.05 --label random_w3_25_rpar5_scw_pilot_frac05`
and `--checkpoint-interval 10`. Record the stage times from the log and peak
memory (`peak_worker_rss_mb`, `parent_rss_mb` in the metadata). Pairs scale as
fraction², so the full run has ~400× the pilot's pairs; the search stage does
not scale that way. Use the pilot to decide how many shards/VMs to use.

Rehearsed locally on 2026-09-26 at 1% (South unsharded, North as two shards
with 4 processes, then merged): identical to serial runs to 1e-15, all log
lines present, section-6 check and section-7 bundle work as written. On macOS
`seq -s,` leaves a trailing comma (`...,25,`); the parser ignores it.

### Optional checkpoint/resume smoke test

Interrupt a pilot only after at least one checkpoint message appears, then use
the same science parameters, label, chunk size, and checkpoint path with
`--resume-checkpoint`. Do not use `--overwrite` when resuming.

```bash
PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
  --dataset BOSS \
  --region South \
  --catalog-type random \
  --fraction 0.001 \
  --rpar 10 \
  --rperp-min 4 \
  --rperp-max 6 \
  --label random_5_rpar10_pilot_frac001 \
  --output-dir analysis/boss/results/validation \
  --checkpoint-path analysis/boss/results/checkpoints/random_5_rpar10_pilot_frac001_BOSS_South.npz \
  --resume-checkpoint \
  --n-processes 4 \
  --chunk-size 500
```

For a subsample, omit `--seed` on resume: the script reloads the saved seed
from the checkpoint.

---

## 4. Production launch

### Start and use tmux

```bash
tmux new -s rperp5_south
```

Inside tmux:

```bash
cd ~/astrophysics
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate astro

mkdir -p analysis/boss/results/checkpoints analysis/boss/results/logs
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
set -o pipefail
```

Detach with `Ctrl-b`, then `d`. Reattach later with:

```bash
tmux attach -t rperp5_south
```

Closing Browser SSH does not stop a tmux job. Stopping the VM does.

### Verified South command: `r_perp=4-6`, `r_par<=10`

```bash
PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
  --dataset BOSS \
  --region South \
  --catalog-type random \
  --fraction 1.0 \
  --rpar 10 \
  --rperp-min 4 \
  --rperp-max 6 \
  --label random_5_rpar10_frac100 \
  --output-dir analysis/boss/results \
  --checkpoint-path analysis/boss/results/checkpoints/random_5_rpar10_frac100_BOSS_South.npz \
  --checkpoint-interval 10 \
  --n-processes 16 \
  --chunk-size 5000 \
  --overwrite \
  2>&1 | tee analysis/boss/results/logs/random_5_rpar10_frac100_BOSS_South.log
```

### Verified North command

Use a separate session, for example `tmux new -s rperp5_north`, and the same
environment setup as above.

```bash
PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
  --dataset BOSS \
  --region North \
  --catalog-type random \
  --fraction 1.0 \
  --rpar 10 \
  --rperp-min 4 \
  --rperp-max 6 \
  --label random_5_rpar10_frac100 \
  --output-dir analysis/boss/results \
  --checkpoint-path analysis/boss/results/checkpoints/random_5_rpar10_frac100_BOSS_North.npz \
  --checkpoint-interval 10 \
  --n-processes 16 \
  --chunk-size 5000 \
  --overwrite \
  2>&1 | tee analysis/boss/results/logs/random_5_rpar10_frac100_BOSS_North.log
```

Do not launch North and South simultaneously on the 22-vCPU VM. They would
oversubscribe CPU and memory. Complete one, validate and bundle it, then run the
other.

### Wide-bin production (R4): r⊥ 3–25 in 1 h⁻¹Mpc sub-bins, r∥ ≤ 5, 1/Σ_crit²

```bash
W=analysis/boss/results/widebin_rpar5
mkdir -p $W/checkpoints $W/logs
EDGES=$(seq -s, 3 25)
RLABEL=random_w3_25_rpar5_scw_frac100
REGION=South            # then North
SHARD=0/1               # e.g. 0/2 and 1/2 on two VMs for North
SUFFIX=$([ "$SHARD" = 0/1 ] || echo "_shard${SHARD%/*}of${SHARD#*/}")
STEM=${RLABEL}_BOSS_${REGION}${SUFFIX}
LOGS=$W/logs

{ git rev-parse HEAD; date -u; hostname; nproc; free -h; df -h .; } > $LOGS/${STEM}_env.txt
vmstat -t 300 > $LOGS/${STEM}.vmstat & VMSTAT_PID=$!

PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
  --dataset BOSS --region $REGION --catalog-type random --fraction 1.0 \
  --rpar 5 --rperp-bin-edges "$EDGES" --sigma-crit-weight \
  --shard $SHARD --label "$RLABEL" --output-dir $W \
  --checkpoint-path "$W/checkpoints/${STEM}.npz" \
  --checkpoint-interval 10 --n-processes 16 --chunk-size 5000 \
  2>&1 | tee $LOGS/${STEM}.log

kill $VMSTAT_PID
```

Notes:
- Every shard of one region must use identical flags apart from `--shard`
  (same edges, chunk size, fraction, label). The merge refuses anything else.
- No `--overwrite` anywhere. New outputs live only under `widebin_rpar5/`.
- `--n-processes 16` is for `c3-standard-22`; scale with the machine.

### Merge shards (only when a region was sharded)

On the VM (or locally after downloading all shard bundles):

```bash
PYTHONPATH=lib python analysis/boss/scripts/merge_pair_shards.py \
  --inputs $W/kappa_pairs_${RLABEL}_BOSS_North_shard0of2.npz,$W/kappa_pairs_${RLABEL}_BOSS_North_shard1of2.npz \
  --output $W/kappa_pairs_${RLABEL}_BOSS_North.npz
```

It checks that the shards are exactly 0..N−1, share one configuration, and
together cover every chunk once; it writes the unsharded `.npz`, `.csv` and
`.meta.json` (with each shard's SHA-256 and provenance).

### Resume a production job

Use the identical dataset, region, catalog fraction, separation cuts, label,
output directory, checkpoint path, and chunk size. Replace `--overwrite` with
`--resume-checkpoint`. The number of workers may be changed if necessary.

If a final output CSV already exists, the script intentionally refuses to
resume. Confirm whether that CSV is complete before moving it aside.

Wide-bin runs: resume keeps each earlier invocation's code version in the
metadata (`invocations`, `mixed_code_versions`), and may rewrite a partial
`.npz`/`.meta.json` but never a finished CSV. A sharded run resumes each shard
with its own `--shard` and checkpoint path.

---

## 5. Observed production performance

### Latest run on `c3-standard-22`, 16 workers

| Region | Objects | Chunks | Selected pairs | Skipped geometry | Runtime |
|---|---:|---:|---:|---:|---:|
| South | 10.74 M | 2,149 | 84,077,920 | 0 | 4.14 h |
| North | 30,113,342 | 6,023 | 252,484,638 | 798,798 | 36.00 h |

The North skipped count was recorded by the completed run and did not prevent
all chunks or the final map from completing. It represents pairs rejected by
the map-coordinate geometry guard, not incomplete chunks.

North was about 8.7 times slower than South even though it produced only about
3 times as many pairs. Do not estimate North solely by multiplying the South
runtime by catalog size or pair count. For scheduling on the retained VM, use
the measured region-specific runtimes.

### Other completed full-catalog benchmarks

These were produced with 24-28 workers in an earlier environment and are
useful for scale, but are not directly comparable to the latest 16-worker run.

| Cut | Region | Pairs | Workers | Runtime |
|---|---|---:|---:|---:|
| `r_perp=4-6`, `r_par<=5` | South | 42.1 M | 24 | 1.37 h |
| `r_perp=4-6`, `r_par<=5` | North | 126.3 M | 28 | 9.92 h |
| `r_perp=9-11`, `r_par<=10` | South | 166.5 M | 28 | 3.86 h |
| `r_perp=9-11`, `r_par<=10` | North | 502.4 M | 28 | 23.20 h |
| `r_perp=18-22`, `r_par<=10` | South | 651.9 M | 24 | 8.97 h |
| `r_perp=18-22`, `r_par<=10` | North | 1.987 B | 28 | 36.33 h |

### Planning a larger run

For a first estimate, pair counts scale approximately with

```text
r_par_max * (rperp_max^2 - rperp_min^2)
```

relative to the latest baseline
`10 * (6^2 - 4^2) = 200`. This is only a sizing estimate; measured wall time is
not linear enough to justify a precise forecast.

Examples at `r_par<=10`:

| New perpendicular bin | Pair-count factor vs `4-6` | Conservative scheduling note |
|---|---:|---|
| `3-7` | 2x | Roughly 8-12 h South and 72-96 h North on the retained VM |
| `8-12` | 4x | Use the measured `9-11` benchmark as the better reference |
| `18-22` | 8x | A prior full run measured 9.0 h South / 36.3 h North on 24-28 workers |

Doubling `r_par_max` from 10 to 20 approximately doubles the candidate and
valid pair counts again.

The time window for `3-7` is deliberately broad. It applies only to the
retained 22-vCPU VM and should be replaced by a measured performance-pilot
estimate before launch.

Before any materially larger run:

1. Run the 0.1% South smoke test with the exact proposed cuts.
2. Run a 5-10% South performance pilot with the same 16 workers and chunk size
   as production.
3. Record objects, chunks, pair count, peak memory, and wall time.
4. Extrapolate South and North separately.
5. Check free disk space with `df -h` and memory with `free -h`/`htop`.
6. Confirm C3 quota and expected cost before changing machine type.

For a random catalog fraction `f`, the number of possible pairs scales roughly
as `f^2`. Thus a 10% catalog pilot has about 1% of the full-run pairs, not 10%.
Subtract fixed startup time (catalog and Planck loading) before scaling the
pair-processing portion. Treat the result as a range: the latest North and
South jobs demonstrate that throughput also changes with footprint and load.

The existing 22-vCPU/88-GiB VM can run a larger pair selection; it will mainly
take longer. With the 128-vCPU quota (2026-09-26) the options are more VMs (one
shard each, `--shard K/N`, then `merge_pair_shards.py`) or a larger machine. On
`c3-standard-44`, begin with 32-36 workers rather than all 44, and repeat the
pilot: multiprocessing speedup is not guaranteed to
be linear, and each worker has its own candidate-pair arrays and 101x101
accumulators.

---

## 6. Completion checks

A job is complete only when the log contains all of the following:

- `Saved stacked map to ...csv`
- `Saved metadata to ...meta.json`
- `Completed in ...`

Then inspect the metadata:

```bash
python - <<'PY'
import json
from pathlib import Path

for region in ("South", "North"):
    path = Path(
        f"analysis/boss/results/"
        f"kappa_pairs_random_5_rpar10_frac100_BOSS_{region}.meta.json"
    )
    if not path.exists():
        continue
    meta = json.loads(path.read_text())
    assert meta["n_chunks_completed"] == meta["n_chunks_total"]
    assert meta["rpar"] == 10.0
    assert meta["rperp_min"] == 4.0
    assert meta["rperp_max"] == 6.0
    print(
        region,
        f"chunks={meta['n_chunks_completed']}/{meta['n_chunks_total']}",
        f"pairs={meta['total_pairs']:,}",
        f"hours={meta['runtime_seconds']/3600:.2f}",
        f"skipped={meta['total_skipped']:,}",
    )
PY
```

Validate the numerical map using the `astro` environment:

```bash
python - <<'PY'
from pathlib import Path
import numpy as np
import pandas as pd

for region in ("South", "North"):
    path = Path(
        f"analysis/boss/results/"
        f"kappa_pairs_random_5_rpar10_frac100_BOSS_{region}.csv"
    )
    if not path.exists():
        continue
    arr = pd.read_csv(path, index_col=0).to_numpy()
    assert arr.shape == (101, 101), arr.shape
    assert np.isfinite(arr).all()
    print(region, arr.shape, arr.min(), arr.max(), arr.mean())
PY
```

Expected latest-run metadata:

| Region | Completed chunks | Total pairs | Runtime |
|---|---:|---:|---:|
| South | 2,149 / 2,149 | 84,077,920 | 14,887.7 s |
| North | 6,023 / 6,023 | 252,484,638 | 129,587.6 s |

### Wide-bin completion checks

Log must contain `Saved stacked map to …csv`, `Saved metadata to …meta.json`,
`Completed in …`, `Saved raw per-bin sums to …npz`, `1/Sigma_crit^2 pair weight: ON`,
and `Pairs per r_perp sub-bin: …`.

```bash
python - <<'PY'
import json
from pathlib import Path
import numpy as np

W = Path("analysis/boss/results/widebin_rpar5")
for meta_path in sorted(W.glob("kappa_pairs_random_w3_25_rpar5_scw_frac100_BOSS_*.meta.json")):
    meta = json.loads(meta_path.read_text())
    if "merged_from" in meta:                      # merged shards
        assert meta["n_chunks_completed"] == meta["n_chunks_total"]
        print(meta_path.name, "merged from", meta["n_shards"], "shards",
              f"pairs={meta['total_pairs']:,}", "mixed code:", meta["mixed_code_versions"])
    else:
        if meta["shard"] == "0/1":                 # a shard holds only its share
            assert meta["n_chunks_completed"] == meta["n_chunks_catalog"]
        assert meta["rpar"] == 5.0 and meta["sigma_crit_weight"] is True
        assert meta["rperp_bin_edges"] == [float(x) for x in range(3, 26)]
        print(meta_path.name, f"shard={meta['shard']}",
              f"chunks={meta['n_chunks_completed']}/{meta['n_chunks_catalog']}",
              f"pairs={meta['total_pairs']:,}", f"skipped={meta['total_skipped']:,}",
              f"hours={meta['runtime_seconds_all_invocations']/3600:.2f}",
              f"commit={meta['git_commit'][:8]} dirty={meta['git_dirty']}")
    sums = np.load(str(meta_path).replace(".meta.json", ".npz"))
    assert sums["sum_wk"].shape == (22, 101, 101)
    assert np.isfinite(sums["sum_wk"]).all() and (sums["sum_w"].sum(0) > 0).all()
    counts = sums["n_pairs"]
    print("  pairs per bin (should rise ~linearly with r_perp):", counts.tolist())
PY
```

The glob also picks up shard files (`…_shard0of2.meta.json`); their chunk
completeness is checked by `merge_pair_shards.py`, which refuses an unfinished
shard. `git_dirty` should be `false` on the VM.

Tested locally (2026-09-26) on the 1% N+S runs and a merged 2-shard run.

---

## 7. Bundle and download

August runs kept four artifacts per region (wide-bin runs keep seven; see
"Wide-bin bundle" below):

1. final map CSV;
2. metadata JSON;
3. final checkpoint NPZ;
4. complete log.

### Reusable bundle command

Set the same region and label used for the production command:

```bash
RUN_REGION=South
RUN_LABEL=random_5wide_rpar10_frac100
RUN_STEM="${RUN_LABEL}_BOSS_${RUN_REGION}"

tar -czf "$HOME/${RUN_STEM}_bundle.tar.gz" \
  "analysis/boss/results/kappa_pairs_${RUN_STEM}.csv" \
  "analysis/boss/results/kappa_pairs_${RUN_STEM}.meta.json" \
  "analysis/boss/results/checkpoints/${RUN_STEM}.npz" \
  "analysis/boss/results/logs/${RUN_STEM}.log"

tar -tzf "$HOME/${RUN_STEM}_bundle.tar.gz"
sha256sum "$HOME/${RUN_STEM}_bundle.tar.gz"
```

### Bundle South

From `~/astrophysics`:

```bash
tar -czf "$HOME/random_5_rpar10_frac100_BOSS_South_bundle.tar.gz" \
  analysis/boss/results/kappa_pairs_random_5_rpar10_frac100_BOSS_South.csv \
  analysis/boss/results/kappa_pairs_random_5_rpar10_frac100_BOSS_South.meta.json \
  analysis/boss/results/checkpoints/random_5_rpar10_frac100_BOSS_South.npz \
  analysis/boss/results/logs/random_5_rpar10_frac100_BOSS_South.log
```

### Bundle North

```bash
tar -czf "$HOME/random_5_rpar10_frac100_BOSS_North_bundle.tar.gz" \
  analysis/boss/results/kappa_pairs_random_5_rpar10_frac100_BOSS_North.csv \
  analysis/boss/results/kappa_pairs_random_5_rpar10_frac100_BOSS_North.meta.json \
  analysis/boss/results/checkpoints/random_5_rpar10_frac100_BOSS_North.npz \
  analysis/boss/results/logs/random_5_rpar10_frac100_BOSS_North.log
```

Inspect before downloading:

```bash
tar -tzf "$HOME/random_5_rpar10_frac100_BOSS_North_bundle.tar.gz"
sha256sum "$HOME/random_5_rpar10_frac100_BOSS_North_bundle.tar.gz"
```

### Wide-bin bundle: seven files per region or shard

Per `STEM` (region or shard), from `~/astrophysics` with `W` and `STEM` set as in the launch:

```bash
tar -czf "$HOME/${STEM}_bundle.tar.gz" \
  "$W/kappa_pairs_${STEM}.npz" \
  "$W/kappa_pairs_${STEM}.csv" \
  "$W/kappa_pairs_${STEM}.meta.json" \
  "$W/checkpoints/${STEM}.npz" \
  "$W/logs/${STEM}.log" \
  "$W/logs/${STEM}_env.txt" \
  "$W/logs/${STEM}.vmstat"
sha256sum "$HOME/${STEM}_bundle.tar.gz" | tee "$HOME/${STEM}_bundle.tar.gz.sha256"
tar -tzf "$HOME/${STEM}_bundle.tar.gz"
```

`tar` stops if a listed file is missing: `_env.txt` and `.vmstat` come from the
two launch lines around the run (Linux `free`/`vmstat -t`), so if either failed
on the VM, drop that file from the list rather than skipping the bundle.

The `.npz` next to the CSV is the **main product** (raw per-sub-bin sums,
unsymmetrized); the CSV is a quick-look map only. The final checkpoint is now
current (written after the last chunk); older runs' checkpoints were stale.
Keep shard bundles until the merge is verified locally.

### Download option A: Browser SSH

Click **Download file** at the top of Browser SSH and enter the absolute path,
for example:

```text
/home/carolinewu327/random_5_rpar10_frac100_BOSS_North_bundle.tar.gz
```

### Download option B: `gcloud` from the local computer

```bash
gcloud compute scp \
  astro-pairs:~/random_5_rpar10_frac100_BOSS_North_bundle.tar.gz \
  ~/Downloads/ \
  --project=astrophysics-cwu --zone=us-central1-a
```

### Extract locally into the repository

The archive stores repository-relative paths, so extract it from the repository
root:

```bash
cd /Users/carolinewu/PycharmProjects/astrophysics
tar -xzf ~/Downloads/random_5_rpar10_frac100_BOSS_North_bundle.tar.gz
```

List the archive before extraction if the bundle was made differently:

```bash
tar -tzf ~/Downloads/random_5_rpar10_frac100_BOSS_North_bundle.tar.gz
```

The verified August 2026 North and South archives each contained exactly the
four artifacts listed above, and their extracted CSV/metadata files matched the
archive copies byte-for-byte.

---

## 8. Stop without deleting the setup

After the bundle is downloaded and validated locally:

1. Confirm no production process remains: `pgrep -af find_and_stack_pairs`.
2. Exit or kill the finished tmux session if desired.
3. In **Compute Engine -> VM instances**, select `astro-pairs` and click
   **Stop**, not **Delete**.

Stopping the VM ends CPU/RAM billing but retains the persistent disk and the
full environment. A tmux session does not reduce billing; only stopping the VM
does. The stopped VM's tmux sessions will not survive a restart, but its files
will.

---

## 9. Pipeline implementation notes

- The parent loads the Planck kappa map and mask before creating the process
  pool. Linux `fork` shares these read-only arrays by copy-on-write.
- Workers return 101x101 weighted accumulator grids, not pair tables.
- `imap_unordered` balances chunks; checkpoints therefore store the complete
  set of finished chunk IDs, not only a last-chunk number.
- Checkpoints are written atomically and contain accumulators, pair/skipped
  counts, chunk layout, completed chunk IDs, run configuration, output path,
  and subsampling seed.
- Resume validates the saved run configuration. Never combine a checkpoint
  with changed cuts, fraction, label, output path, chunk size, or seed.
- `--fraction < 1` must use a reproducible seed. A generated seed is stored and
  automatically reused on resume.
- `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `MKL_NUM_THREADS` are set to 1
  so every multiprocessing worker does not start its own numerical thread
  pool.
- Wide-bin additions (2026-09):
  - `--rperp-bin-edges` routes each pair into `[e_i, e_i+1)` sub-bins (top one
    closed) in one pass; accumulators are `(n_bins, 101, 101)`. Search cost
    depends on r_par, not on the r_perp range.
  - `--sigma-crit-weight` multiplies each pair weight by 1/Sigma_crit^2 at the
    pair's mean redshift (4096-point table over the catalog z range).
  - `--shard K/N` keeps chunks with `id % N == K`; the chunk layout is the same
    for every shard. Fractional multi-shard runs require `--seed`.
  - A final checkpoint is written after the last chunk (August checkpoints
    are stale: North 6,020 of 6,023 chunks).
  - Existing outputs are refused unless `--overwrite`.
- Full random production must use `find_and_stack_pairs.py`. The older
  `find_pairs.py -> CSV -> stack_pairs.py` route remains appropriate only for
  smaller pair catalogs and validation.

Relevant tracked files:

| File | Role |
|---|---|
| `analysis/boss/scripts/find_and_stack_pairs.py` | combined full-random pair finder and stacker |
| `analysis/boss/scripts/merge_pair_shards.py` | validate and sum shard outputs |
| `analysis/boss/scripts/find_pairs.py` | CSV pair finder for smaller runs |
| `analysis/boss/scripts/stack_pairs.py` | stack an existing pair catalog |
| `lib/catalog.py` | catalog and Planck-map loading |
