# Plan: Migrate Pair-Finding Pipeline to GCP c3-standard-32

## Context

The full BOSS random catalog (~15M+ objects North, ~5M+ South) produces too many pairs (~2B) to materialize as CSV or hold in memory. The current two-step pipeline (`find_pairs.py -> CSV -> stack_pairs.py`) cannot scale for full random pairs. This migration provides:
1. `find_and_stack_pairs.py`, a combined find-and-stack script that never materializes all pairs
2. performance fixes to the existing CSV-writing `find_pairs.py` path for galaxy pairs and smaller validation runs
3. a complete VM setup and run guide for the full-random production run

---

## Part 1: VM Setup

### 1.1 Create GCP VM
```bash
gcloud compute instances create astro-pairs \
    --zone=us-central1-a \
    --machine-type=c3-standard-32 \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --boot-disk-size=200GB \
    --boot-disk-type=pd-ssd
```
200 GB SSD covers data (~5 GB), conda env (~5 GB), and working space.

### 1.2 SSH + system packages
```bash
gcloud compute ssh astro-pairs --zone=us-central1-a
sudo apt-get update && sudo apt-get install -y build-essential git wget curl tmux htop
```
`tmux` is essential — the run will take many hours and must survive SSH disconnects.

### 1.3 Install Miniforge (mamba)
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda init bash && source ~/.bashrc
```

### 1.4 Create conda environment
```bash
mamba create -n astro python=3.12 -y
mamba activate astro
mamba install -c conda-forge numpy matplotlib astropy healpy pandas scipy tqdm -y
```

### 1.5 Verify
```bash
python -c "import healpy, astropy, numpy, pandas, scipy; print('OK')"
```

---

## Part 2: Data Download

### 2.1 Clone repo + create directories
```bash
cd ~ && git clone git@github.com:carolinewu327/astrophysics.git && cd astrophysics
mkdir -p data/BOSS data/planck/COM_Lensing_4096_R3.00/MV data/paircatalogs/BOSS
mkdir -p analysis/boss/results analysis/boss/results/checkpoints output/plots
```

### 2.2 BOSS catalogs from SDSS SAS (~4.4 GB total)
```bash
cd ~/astrophysics/data/BOSS

# Galaxy catalogs (~183 MB)
wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_CMASS_North.fits.gz
wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_CMASS_South.fits.gz

# Random catalogs (~4.2 GB)
wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_North.fits.gz
wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASS_South.fits.gz
```

### 2.3 Planck lensing data from ESA PLA
The Planck Legacy Archive uses a JavaScript UI that makes direct URL construction unreliable. Best approach:
1. Navigate to https://pla.esac.esa.int/#maps on your local browser
2. Search for "COM_Lensing_4096_R3.00" in Cosmology products
3. Download the .tgz tarball (~370 MB)
4. Upload to VM: `gcloud compute scp COM_Lensing_4096_R3.00.tgz astro-pairs:~/astrophysics/data/planck/ --zone=us-central1-a`

On the VM:
```bash
cd ~/astrophysics/data/planck
tar -xzf COM_Lensing_4096_R3.00.tgz
ls COM_Lensing_4096_R3.00/MV/dat_klm.fits COM_Lensing_4096_R3.00/mask.fits
```

Or, if you already have the data locally:
```bash
# From your local machine (WSL):
gcloud compute scp /mnt/c/Users/ECM-Guest/dev/astrophysics/data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits \
    astro-pairs:~/astrophysics/data/planck/COM_Lensing_4096_R3.00/MV/ --zone=us-central1-a
gcloud compute scp /mnt/c/Users/ECM-Guest/dev/astrophysics/data/planck/COM_Lensing_4096_R3.00/mask.fits \
    astro-pairs:~/astrophysics/data/planck/COM_Lensing_4096_R3.00/ --zone=us-central1-a
```

### 2.4 Upload existing results (optional)
If galaxy pair catalogs and single-stack results already exist locally, upload them to avoid re-running:
```bash
gcloud compute scp /mnt/c/Users/ECM-Guest/dev/astrophysics/data/paircatalogs/BOSS/galaxy_pairs_*.csv \
    astro-pairs:~/astrophysics/data/paircatalogs/BOSS/ --zone=us-central1-a
gcloud compute scp /mnt/c/Users/ECM-Guest/dev/astrophysics/analysis/boss/results/*.csv \
    astro-pairs:~/astrophysics/analysis/boss/results/ --zone=us-central1-a
```

### 2.5 Set PYTHONPATH
```bash
echo 'export PYTHONPATH=$HOME/astrophysics/lib:$PYTHONPATH' >> ~/.bashrc
source ~/.bashrc
```

### 2.6 Preflight check
```bash
cd ~/astrophysics
python -m py_compile \
    analysis/boss/scripts/find_pairs.py \
    analysis/boss/scripts/find_and_stack_pairs.py \
    analysis/boss/scripts/stack_pairs.py \
    lib/catalog.py
python - <<'PY'
from catalog import resolve_catalog_path, resolve_planck_paths
print(resolve_catalog_path("data", "BOSS", "North", "galaxy"))
print(resolve_planck_paths("data"))
PY
```

---

## Part 3: Combined Script — `find_and_stack_pairs.py`

**Implemented file:** `analysis/boss/scripts/find_and_stack_pairs.py`

### Architecture

```
Main process:
  1. Load catalog → l, b, D, z, weights
  2. Load Planck kappa map + mask into module globals BEFORE fork
  3. Build sorted_indices/sorted_D and chunk_ranges = [(chunk_id, start, end), ...]
  4. Create reducer state:
       total_sum_wk : float64[101,101]
       total_sum_w  : float64[101,101]
       total_pairs  : int64
       total_skipped: int64
       completed_chunks : set[int]
  5. Pool context = multiprocessing.get_context("fork")
  6. pool.imap_unordered(process_chunk_and_stack, pending_chunks)
  7. For each worker result:
       - add partial grids to reducer
       - add chunk_id to completed_chunks
       - checkpoint atomically every N completed chunks
  8. Compute mean = total_sum_wk / total_sum_w where total_sum_w > 0
  9. Apply reflect_symmetrize_map
 10. Save final CSV + final metadata JSON
```

### Key design decisions implemented

**Workers return grids, not pairs.** Each worker's `process_chunk_and_stack()` replaces the dict-building loop in `process_chunk_optimized()` with inline stacking adapted from `stack_pairs.py`. It returns:
```python
{
    "chunk_id": int,
    "sum_wk": np.ndarray(shape=(101, 101), dtype=np.float64),
    "sum_w": np.ndarray(shape=(101, 101), dtype=np.float64),
    "n_pairs": int,
    "n_skipped": int,
}
```
That is a few hundred KB per completed chunk instead of millions of pair rows.

**Kappa map shared via fork() COW.** The kmap (~400 MB) and mask (~400 MB) are loaded as module-level globals *before* `Pool()` creation. On Linux, `fork()` shares them read-only via copy-on-write. Do NOT pass them through `partial()` or `initargs` (that triggers pickle/copy). Total physical memory: ~6 GB (not 32 × 800 MB).

**Force fork context explicitly.** Multiprocess production mode uses:
```python
ctx = multiprocessing.get_context("fork")
with ctx.Pool(processes=n_processes) as pool:
    ...
```
Do not rely on the platform default. The memory-sharing design is intended for Linux production. Local validation should use `--n-processes 1`.

**Vectorized Dmid.** The script replaces per-pair `cosmo.comoving_distance()` calls with one vectorized call per chunk on `z_mid = (z[i_final] + z[j_final]) / 2`.

**`imap_unordered`** for load balancing across workers. This means checkpoint/resume must track the full set of completed chunk IDs, not a single "last chunk".

**Chunk size is a memory control.** Start with `--chunk-size 5000` for random catalogs. The candidate-pair arrays (`i_list`, `j_list`) are still the main per-worker memory risk, even though the final pair list is never materialized.

### Module globals and helper functions

The script defines these module globals near the top:
```python
KMAP = None
MASK = None
NSIDE = None
X_GRID = None
Y_GRID = None
```

It sets them once in the parent process before `Pool()` creation:
```python
X_vals = np.linspace(-HALF_SIZE, HALF_SIZE, GRID_RES)
Y_vals = np.linspace(-HALF_SIZE, HALF_SIZE, GRID_RES)
X_GRID, Y_GRID = np.meshgrid(X_vals, Y_vals)
KMAP, MASK, NSIDE = load_kappa_map(...)
```

The implementation uses these helper functions:
- `process_chunk_and_stack(chunk_meta)` where `chunk_meta = (chunk_id, start, end)` and data arrays are module globals
- `stack_pairs_inline(i_final, j_final) -> (sum_wk, sum_w, n_skipped)`
- `save_checkpoint_atomic(path, reducer_state, args_dict)`
- `load_checkpoint(path) -> reducer_state`

### Stacking logic per worker (embedded from stack_pairs.py:95-161)

Inside one worker:

1. Run the current two-pass pair search from `find_pairs.py`:
   - rough filter by `r_par`
   - rough angular filter
   - exact `angular_separation()`
   - final `r_perp` mask
2. Vectorize `Dmid_arr = cosmo.comoving_distance(z_mid).value * cosmo.h`
3. Allocate worker-local accumulators:
```python
sum_wk = np.zeros((GRID_RES, GRID_RES), dtype=np.float64)
sum_w = np.zeros((GRID_RES, GRID_RES), dtype=np.float64)
n_skipped = 0
```
4. For each valid pair found in the chunk:
   - convert `l1,b1,l2,b2` to radians
   - enforce consistent longitude ordering
   - compute `(lc, bc)` midpoint and `(cos_theta, sin_theta)`
   - map `X_GRID, Y_GRID` to sky coordinates
   - sample `KMAP` and `MASK`
   - accumulate:
```python
pair_weight = w1 * w2
sum_wk += pair_weight * kappa_vals
sum_w += pair_weight * valid_mask.astype(np.float64)
```
5. Return the result dict shown above.

This is the same geometry as `stack_pairs.py`, just embedded into the worker so that pairs never become a DataFrame.

### Checkpoint strategy

Use one checkpoint file per production run, for example:
```text
analysis/boss/results/checkpoints/kappa_pairs_random_20_frac100_BOSS_North.npz
```

Checkpoint contents:
```python
np.savez(
    checkpoint_tmp_path,
    total_sum_wk=total_sum_wk,
    total_sum_w=total_sum_w,
    total_pairs=np.int64(total_pairs),
    total_skipped=np.int64(total_skipped),
    completed_chunks=np.array(sorted(completed_chunks), dtype=np.int32),
    chunk_ranges=np.array(chunk_ranges, dtype=np.int64),
    args_json=np.array(json.dumps(serialized_args)),
    output_path=np.array(output_path),
    saved_at=np.array(time.time()),
)
os.replace(checkpoint_tmp_path, checkpoint_path)
```

Also write a sidecar metadata JSON next to the final CSV, for example:
```text
analysis/boss/results/kappa_pairs_random_20_frac100_BOSS_North.meta.json
```

Minimum metadata fields:
- `label`
- `dataset`
- `region`
- `catalog_type`
- `fraction`
- `rpar`
- `rperp_min`
- `rperp_max`
- `n_processes`
- `chunk_size`
- `checkpoint_interval`
- `n_chunks_total`
- `n_chunks_completed`
- `total_pairs`
- `total_skipped`
- `grid_res`
- `started_at`
- `finished_at`
- `runtime_seconds`
- `output_csv`
- `checkpoint_path`
- `seed`

Important implemented corrections:
- It does **not** save only `last_chunk`. That would be wrong under `imap_unordered`.
- Resume rebuilds `pending_chunks = [c for c in chunk_ranges if c[0] not in completed_chunks]`.
- Resume aborts if the checkpoint metadata does not match the current CLI args.
- For `--fraction < 1.0`, the script uses a reproducible subsampling seed. If no `--seed` is provided, it auto-generates one before loading the catalog, stores it in checkpoint metadata, and reuses it on resume.
- Checkpoints are written atomically with `os.replace()` so an interrupted save does not corrupt the file.

### CLI interface
Use the same science flags as `find_pairs.py`:
- `--dataset`
- `--region`
- `--catalog-type`
- `--fraction`
- `--rpar`
- `--rperp-min`
- `--rperp-max`
- `--n-processes`
- `--chunk-size`
- `--data-dir`
- `--output-dir`
- `--overwrite`

Additional flags:
- `--label` for output naming
- `--checkpoint-path` optional explicit path
- `--seed` for reproducible subsampling; auto-generated and checkpointed when `--fraction < 1.0`
- `--checkpoint-interval` number of completed chunks between checkpoint writes (default: 10)
- `--resume-checkpoint` flag to resume from checkpoint

Recommended output names for provenance:
- `--label random_20_frac100`
- `--label random_10_frac100`
- `--label random_20_test_frac001`

Output: `analysis/boss/results/kappa_pairs_{label}_{dataset}_{region}.csv`

---

## Part 4: Modifications to `find_pairs.py`

These have been applied to the existing CSV-writing script. This path remains useful for galaxy-pair catalogs and smaller validation runs; full random production should use `find_and_stack_pairs.py`.

**File:** `analysis/boss/scripts/find_pairs.py`

### 4.1 Remove 12-process cap
```python
# Before:
n_processes = min(n_processes, 12)
# After: remove this line entirely
```

### 4.2 Switch `imap` -> `imap_unordered`
```python
# Before:
for i, result in enumerate(pool.imap(process_func, chunk_ranges)):
# After:
for i, result in enumerate(pool.imap_unordered(process_func, chunk_ranges)):
```
Adjust logging to track completed chunk count and actual `chunk_id`, since completion order is no longer sequential.

### 4.3 Vectorize `Dmid` and replace list-of-dicts output
The per-pair Python loop has been replaced with a vectorized `Dmid` calculation and per-chunk DataFrames:
```python
z_mid = (z[i_final] + z[j_final]) / 2.0
Dmid_arr = cosmo.comoving_distance(z_mid).value * cosmo.h
return pd.DataFrame({...}, columns=PAIR_COLUMNS)
```
`build_pair_catalog()` concatenates chunk DataFrames at the end instead of extending one large Python list of dicts.

### 4.4 Keep checkpoint cleanup changes
The current "remove previous checkpoint / remove final checkpoint" behavior is still useful for the CSV-writing version of `find_pairs.py`; keep it unless it interferes with debugging.

---

## Part 5: Run Plan

### 5.1 Local deterministic validation
First validate the combined script on a deterministic galaxy catalog. This avoids random-subset mismatch and verifies that the embedded find+stack path matches the existing two-step geometry up to the known precision difference from using float64 `Z` in the new script.

```bash
cd ~/astrophysics

PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
    --dataset BOSS --region South --catalog-type galaxy \
    --rpar 5 --rperp-min 4 --rperp-max 6 \
    --label galaxy_5_validation \
    --n-processes 1 --chunk-size 5000 \
    --output-dir analysis/boss/results/validation \
    --overwrite
```

Compare against the existing two-step output, or regenerate the reference:
```bash
PYTHONPATH=lib python analysis/boss/scripts/stack_pairs.py \
    --pair-catalog data/paircatalogs/BOSS/galaxy_pairs_BOSS_South_5.0_4.0_6.0hmpc.csv \
    --label galaxy_5_reference --dataset BOSS --region South \
    --output-dir analysis/boss/results/validation \
    --overwrite

python -c "
import pandas as pd, numpy as np
new = pd.read_csv('analysis/boss/results/validation/kappa_pairs_galaxy_5_validation_BOSS_South.csv', index_col=0).values
old = pd.read_csv('analysis/boss/results/validation/kappa_pairs_galaxy_5_reference_BOSS_South.csv', index_col=0).values
diff = np.abs(new - old)
print(f'Max diff: {diff.max():.3e}')
print(f'RMS diff: {np.sqrt(np.mean(diff**2)):.3e}')
"
```

Pass criterion:
- `max diff ~5e-6` is expected when comparing the new float64-`Z` path to the legacy float32-`Z` pair catalog.
- Larger discrepancies or visibly different morphology should be investigated before the full run.

### 5.2 Linux pilot and resume smoke test
Run a small random pilot on the VM to exercise the production path: Linux `fork`, multiprocessing, checkpoint write, and checkpoint resume. Use an explicit seed if you want the run to be exactly reproducible from scratch.

```bash
PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
    --dataset BOSS --region South --catalog-type random \
    --fraction 0.05 --seed 12345 \
    --rpar 5 --rperp-min 4 --rperp-max 6 \
    --label random_5_pilot_frac005 \
    --output-dir analysis/boss/results/validation \
    --checkpoint-path analysis/boss/results/checkpoints/random_5_pilot_frac005_BOSS_South.npz \
    --checkpoint-interval 2 \
    --n-processes 4 --chunk-size 5000 \
    --overwrite
```

To test resume behavior, interrupt after a checkpoint is written, remove any partial final CSV if present, and rerun:

```bash
PYTHONPATH=lib python analysis/boss/scripts/find_and_stack_pairs.py \
    --dataset BOSS --region South --catalog-type random \
    --fraction 0.05 \
    --rpar 5 --rperp-min 4 --rperp-max 6 \
    --label random_5_pilot_frac005 \
    --output-dir analysis/boss/results/validation \
    --checkpoint-path analysis/boss/results/checkpoints/random_5_pilot_frac005_BOSS_South.npz \
    --resume-checkpoint \
    --n-processes 4 --chunk-size 5000
```

### 5.3 Full production run
```bash
tmux new -s astro

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# South full random (smaller, finishes first)
python analysis/boss/scripts/find_and_stack_pairs.py \
    --dataset BOSS --region South --catalog-type random \
    --fraction 1.0 --label random_20_frac100 \
    --output-dir analysis/boss/results \
    --checkpoint-path analysis/boss/results/checkpoints/random_20_frac100_BOSS_South.npz \
    --n-processes 24 --chunk-size 5000 --overwrite

# North full random (larger)
python analysis/boss/scripts/find_and_stack_pairs.py \
    --dataset BOSS --region North --catalog-type random \
    --fraction 1.0 --label random_20_frac100 \
    --output-dir analysis/boss/results \
    --checkpoint-path analysis/boss/results/checkpoints/random_20_frac100_BOSS_North.npz \
    --n-processes 24 --chunk-size 5000 --overwrite
```

Estimated runtime: South ~8-12 hours, North ~30-50 hours.
Use `--n-processes 24` (75% of 32 vCPUs) to leave headroom.

If using spot/preemptible capacity, resume with:
```bash
python analysis/boss/scripts/find_and_stack_pairs.py \
    --dataset BOSS --region North --catalog-type random \
    --fraction 1.0 --label random_20_frac100 \
    --output-dir analysis/boss/results \
    --checkpoint-path analysis/boss/results/checkpoints/random_20_frac100_BOSS_North.npz \
    --resume-checkpoint \
    --n-processes 24 --chunk-size 5000
```

### 5.4 Download results + teardown
```bash
# From local machine:
gcloud compute scp astro-pairs:~/astrophysics/analysis/boss/results/*.csv \
    /mnt/c/Users/ECM-Guest/dev/astrophysics/analysis/boss/results/ --zone=us-central1-a

gcloud compute instances delete astro-pairs --zone=us-central1-a
```

**Cost note:** check current GCP pricing for `c3-standard-32` before launch. The total cost is:
```text
hourly_instance_price * total_runtime_hours
```
For spot/preemptible runs, checkpoint/resume is mandatory.

---

## Files to create/modify

| File | Action |
|------|--------|
| `analysis/boss/scripts/find_and_stack_pairs.py` | **Created** — combined find+stack script |
| `analysis/boss/scripts/find_pairs.py` | **Modified** — remove cap, switch to `imap_unordered`, vectorize `Dmid`, use chunk DataFrames |
| `lib/catalog.py` | **Modified** — add optional deterministic subsampling seed |
| `CLAUDE.md` | **Updated** — document new script and pipeline behavior |
