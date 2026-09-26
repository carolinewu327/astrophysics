#!/usr/bin/env python
"""Merge ``find_and_stack_pairs.py --shard K/N`` outputs into one run.

Each shard processes the chunks with ``id % N == K`` out of the same chunk
layout, so the raw sums are additive: the merged sums equal those of a single
unsharded run (up to floating-point summation order).  Before adding anything
this refuses a set of shards that could not have come from one run:

- every shard's run configuration must match, apart from its shard id and
  output path (cuts, sub-bin edges, weighting, fraction, seed, chunk size,
  catalog, region);
- the shard ids must be exactly ``0..N-1``, each once;
- each shard's completed chunks must be exactly ``K, K+N, K+2N, ...`` below
  the catalog's chunk count, with no repeats -- so an unfinished, duplicated
  or foreign shard cannot be merged by mistake.

Writes the unsharded names: ``<output>.npz`` (raw per-bin sums, same keys as an
unsharded run), ``<output>.csv`` (reflection-symmetrized all-bin map) and
``<output>.meta.json`` (provenance, including each shard file's SHA-256).

Usage
-----
    PYTHONPATH=lib python analysis/boss/scripts/merge_pair_shards.py \\
        --inputs run_BOSS_North_shard0of2.npz,run_BOSS_North_shard1of2.npz \\
        --output analysis/boss/results/widebin_rpar5/kappa_pairs_run_BOSS_North.npz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os

import numpy as np
import pandas as pd

from catalog import setup_logging
from find_and_stack_pairs import finalize_map

logger = logging.getLogger(__name__)

# Config keys that legitimately differ between shards of one run.
SHARD_SPECIFIC_KEYS = ("shard", "output_path")


def sha256_of(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_shard(path: str) -> dict[str, object]:
    with np.load(path, allow_pickle=False) as data:
        if "shard_index" not in data.files:
            raise ValueError(f"{path} has no shard information; was it written by --shard?")
        shard = {key: data[key] for key in data.files}
    shard["path"] = path
    shard["config"] = json.loads(str(shard["config_json"]))
    return shard


def validate_shards(shards: list[dict[str, object]]) -> None:
    """Raise unless the shards partition one run exactly."""
    first = shards[0]
    count = int(first["shard_count"])
    shared_config = {k: v for k, v in first["config"].items() if k not in SHARD_SPECIFIC_KEYS}

    for shard in shards:
        path = shard["path"]
        config = {k: v for k, v in shard["config"].items() if k not in SHARD_SPECIFIC_KEYS}
        if config != shared_config:
            diff = sorted(
                k for k in set(config) | set(shared_config)
                if config.get(k) != shared_config.get(k)
            )
            raise ValueError(f"{path} was run with a different configuration: {diff}")
        if int(shard["shard_count"]) != count:
            raise ValueError(f"{path} is one of {int(shard['shard_count'])} shards, not {count}.")
        for key in ("rperp_edges", "grid_res", "box_size_hmpc", "n_chunks_total", "sigma_crit_weight"):
            if not np.array_equal(shard[key], first[key]):
                raise ValueError(f"{path} differs from {first['path']} in {key}.")

    indices = sorted(int(s["shard_index"]) for s in shards)
    if indices != list(range(count)):
        raise ValueError(f"Need shards 0..{count - 1} exactly once each; got {indices}.")

    # Each shard must hold exactly the chunks it owns, no more, no fewer and no
    # repeats; together the shards then cover the catalog exactly once.
    n_total = int(first["n_chunks_total"])
    for shard in shards:
        k = int(shard["shard_index"])
        ids = np.asarray(shard["completed_chunks"], dtype=np.int64)
        unique_ids, counts = np.unique(ids, return_counts=True)
        if np.any(counts > 1):
            raise ValueError(f"{shard['path']} repeats chunks {unique_ids[counts > 1][:5].tolist()}")
        expected = np.arange(k, n_total, count)
        if not np.array_equal(unique_ids, expected):
            missing = np.setdiff1d(expected, unique_ids)
            extra = np.setdiff1d(unique_ids, expected)
            raise ValueError(
                f"{shard['path']} (shard {k}/{count}) should hold chunks {k}, {k + count}, ... "
                f"of {n_total}; missing {missing[:10].tolist()}, unexpected {extra[:10].tolist()}. "
                "Is the shard unfinished?"
            )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge sharded find_and_stack_pairs.py outputs.")
    parser.add_argument("--inputs", required=True,
                        help="Comma-separated shard .npz files (the raw-sum outputs, not checkpoints).")
    parser.add_argument("--output", required=True,
                        help="Merged .npz path; the .csv and .meta.json are written beside it.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    args.inputs = [p.strip() for p in args.inputs.split(",") if p.strip()]
    if not args.output.endswith(".npz"):
        parser.error("--output must end in .npz")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()

    sums_path = args.output
    csv_path = sums_path[: -len(".npz")] + ".csv"
    meta_path = sums_path[: -len(".npz")] + ".meta.json"
    existing = [p for p in (sums_path, csv_path, meta_path) if os.path.exists(p)]
    if existing and not args.overwrite:
        raise FileExistsError(f"Output exists (use --overwrite): {', '.join(existing)}")
    if os.path.abspath(sums_path) in {os.path.abspath(p) for p in args.inputs}:
        raise ValueError("--output must not be one of the inputs.")

    shards = sorted((load_shard(p) for p in args.inputs), key=lambda s: int(s["shard_index"]))
    validate_shards(shards)
    count = len(shards)
    logger.info("Validated %d shards covering %d chunks.", count, int(shards[0]["n_chunks_total"]))

    sum_wk = np.sum([s["sum_wk"] for s in shards], axis=0)
    sum_w = np.sum([s["sum_w"] for s in shards], axis=0)
    n_pairs = np.sum([s["n_pairs"] for s in shards], axis=0).astype(np.int64)
    total_skipped = int(sum(int(s["total_skipped"]) for s in shards))
    completed = np.sort(np.concatenate([s["completed_chunks"] for s in shards])).astype(np.int32)

    config = {k: v for k, v in shards[0]["config"].items() if k not in SHARD_SPECIFIC_KEYS}
    config["output_path"] = csv_path
    edges = shards[0]["rperp_edges"]
    shard_files = [{"path": s["path"], "sha256": sha256_of(s["path"])} for s in shards]

    tmp = f"{sums_path}.tmp"
    with open(tmp, "wb") as handle:
        np.savez(
            handle,
            rperp_edges=edges,
            sum_wk=sum_wk,
            sum_w=sum_w,
            n_pairs=n_pairs,
            total_skipped=np.int64(total_skipped),
            completed_chunks=completed,
            n_chunks_total=shards[0]["n_chunks_total"],
            shard_index=np.int32(0),
            shard_count=np.int32(1),
            grid_res=shards[0]["grid_res"],
            box_size_hmpc=shards[0]["box_size_hmpc"],
            sigma_crit_weight=shards[0]["sigma_crit_weight"],
            config_json=np.array(json.dumps(config, sort_keys=True)),
            merged_from=np.array(json.dumps(shard_files, sort_keys=True)),
        )
    os.replace(tmp, sums_path)

    final_map = finalize_map({"total_sum_wk": sum_wk, "total_sum_w": sum_w})
    pd.DataFrame(final_map).to_csv(csv_path, index=True)

    # One record per shard, keyed by shard index, so provenance can never be
    # attributed to the wrong shard even when a sidecar is missing.  The code
    # history comes from the .npz itself; the sidecar adds runtime only.
    shard_records = []
    for s, f in zip(shards, shard_files):
        sidecar = s["path"][: -len(".npz")] + ".meta.json"
        meta = None
        if os.path.exists(sidecar):
            with open(sidecar, encoding="utf-8") as handle:
                meta = json.load(handle)
        else:
            logger.warning("Shard %d: sidecar %s is missing.", int(s["shard_index"]), sidecar)
        history = (
            json.loads(str(s["provenance_history_json"]))
            if "provenance_history_json" in s else None
        )
        complete = history is not None and not bool(s.get("prior_unrecorded", False))
        if not complete:
            logger.warning("Shard %d: code provenance is incomplete.", int(s["shard_index"]))
        shard_records.append({
            "shard_index": int(s["shard_index"]),
            "path": f["path"],
            "sha256": f["sha256"],
            "sidecar_path": sidecar,
            "sidecar_found": meta is not None,
            "runtime_seconds": None if meta is None else meta.get("runtime_seconds_all_invocations"),
            "provenance_complete": complete,
            "invocations": history,
        })
    all_invocations = [inv for r in shard_records for inv in (r["invocations"] or [])]
    commits = sorted({str(inv.get("git_commit")) for inv in all_invocations})
    provenance_complete = all(r["provenance_complete"] for r in shard_records)
    metadata = {
        "merged_from": shard_files,
        "n_shards": count,
        "rperp_bin_edges": [float(x) for x in edges],
        "n_pairs_per_bin": [int(x) for x in n_pairs],
        "total_pairs": int(n_pairs.sum()),
        "total_skipped": total_skipped,
        "n_chunks_completed": int(len(completed)),
        "n_chunks_total": int(shards[0]["n_chunks_total"]),
        "sigma_crit_weight": bool(shards[0]["sigma_crit_weight"]),
        "config": config,
        "output_sums": sums_path,
        "output_csv": csv_path,
        "shards": shard_records,
        "provenance_complete": provenance_complete,
        "git_commits_all": commits,
        "mixed_code_versions": len(commits) > 1 or not provenance_complete,
        "any_git_dirty": any(bool(inv.get("git_dirty")) for inv in all_invocations),
    }
    if len(commits) > 1:
        logger.warning("Shards ran different code versions: %s", commits)
    if metadata["any_git_dirty"]:
        logger.warning("At least one shard ran with uncommitted changes to tracked files.")
    with open(meta_path, "w", encoding="ascii") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")

    logger.info("Merged %d pairs (%d skipped) -> %s", int(n_pairs.sum()), total_skipped, sums_path)
    logger.info("Saved stacked map to %s", csv_path)
    logger.info("Saved metadata to %s", meta_path)


if __name__ == "__main__":
    main()
