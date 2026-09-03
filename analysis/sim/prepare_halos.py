#!/usr/bin/env python
"""Prepare a compact BigMDPL host-halo catalog for simulation stacking.

Two selection modes:

- Mass bin (default): ``|log10(Mvir) - log10_mvir| < half_width`` -> CSV,
  as used for the original mass13 science sample.
- Mass floor (``--min-log10-mvir``): ``log10(Mvir) >= floor`` -> intended for
  the one-time compaction of the full hlist before deleting it. Combine with
  a ``.bin`` output for a streamed structured-array binary (float32 fields,
  int64 ids) plus a JSON sidecar, far smaller per row than CSV.

``--halo-catalog`` may be either a raw Rockstar hlist (text) or a previously
written compact ``.bin`` catalog; the latter is the only source once the raw
hlist has been deleted. Selections reaching below a ``.bin`` input's stored
mass floor are rejected rather than silently truncated.

Load a ``.bin`` catalog with::

    info = json.load(open(path + ".json"))
    halos = np.fromfile(path, dtype=np.dtype([tuple(f) for f in info["dtype"]]))
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time

import numpy as np
import pandas as pd

from sim_utils import ensure_parent, setup_logging


logger = logging.getLogger(__name__)

# Units (from the hlist header): Mvir/Mpeak in Msun/h, Rvir/Rs in comoving
# kpc/h, Vmax/Vpeak in physical km/s. Rvir/Rs/Vmax/Mpeak/Vpeak are carried so
# future HOD or abundance-matching work survives deletion of the raw hlist.
ROCKSTAR_USECOLS = [1, 6, 10, 11, 12, 16, 17, 18, 19, 20, 21, 22, 59, 61]
HALO_COLUMNS = [
    "id", "upid", "Mvir", "Rvir", "Rs", "Vmax",
    "x", "y", "z", "vx", "vy", "vz", "Mpeak", "Vpeak",
]

# float32 positions are exact to ~1.5e-4 h^-1 Mpc over the 2500 box; Mvir to
# ~1e-7 dex. Both far below any physical requirement here.
COMPACT_DTYPE = np.dtype(
    [("id", np.int64), ("upid", np.int64)]
    + [(name, np.float32) for name in HALO_COLUMNS[2:]]
)


def read_bin_sidecar(path: str) -> dict:
    with open(f"{path}.json", "r", encoding="utf-8") as f:
        return json.load(f)


def iter_bin_chunks(
    halo_catalog: str,
    chunksize: int,
    max_rows: int | None = None,
):
    info = read_bin_sidecar(halo_catalog)
    dtype = np.dtype([tuple(field) for field in info["dtype"]])
    arr = np.memmap(halo_catalog, dtype=dtype, mode="r")
    n = len(arr) if max_rows is None else min(len(arr), max_rows)
    for start in range(0, n, chunksize):
        block = arr[start : min(start + chunksize, n)]
        yield pd.DataFrame({name: np.asarray(block[name]) for name in dtype.names})


def iter_halo_chunks(
    halo_catalog: str,
    chunksize: int,
    max_rows: int | None = None,
):
    if halo_catalog.endswith(".bin"):
        yield from iter_bin_chunks(halo_catalog, chunksize, max_rows)
        return

    rows_seen = 0
    reader = pd.read_csv(
        halo_catalog,
        sep=r"\s+",
        comment="#",
        header=None,
        usecols=ROCKSTAR_USECOLS,
        names=HALO_COLUMNS,
        chunksize=chunksize,
        engine="c",
    )

    for chunk in reader:
        if max_rows is not None:
            remaining = max_rows - rows_seen
            if remaining <= 0:
                break
            chunk = chunk.iloc[:remaining]
        rows_seen += len(chunk)
        yield chunk
        if max_rows is not None and rows_seen >= max_rows:
            break


def select_hosts(
    chunk: pd.DataFrame,
    log10_mvir: float,
    half_width: float,
    min_log10_mvir: float | None,
) -> pd.DataFrame:
    mvir = chunk["Mvir"].to_numpy(dtype=np.float64, copy=False)
    mask = (chunk["upid"].to_numpy() == -1) & np.isfinite(mvir) & (mvir > 0)
    if min_log10_mvir is not None:
        mask &= np.log10(np.maximum(mvir, 1.0)) >= min_log10_mvir
    else:
        mask &= np.abs(np.log10(np.maximum(mvir, 1.0)) - log10_mvir) < half_width
    return chunk.loc[mask, HALO_COLUMNS]


def chunk_to_compact(selected: pd.DataFrame) -> np.ndarray:
    arr = np.empty(len(selected), dtype=COMPACT_DTYPE)
    arr["id"] = selected["id"].to_numpy(dtype=np.int64)
    arr["upid"] = selected["upid"].to_numpy(dtype=np.int64)
    for name in HALO_COLUMNS[2:]:
        arr[name] = selected[name].to_numpy(dtype=np.float32)
    return arr


def write_sidecar(
    output: str,
    n_rows: int,
    n_seen: int,
    halo_catalog: str,
    selection: dict,
) -> None:
    meta = {
        "path": os.path.abspath(output),
        "n_rows": n_rows,
        "rows_read_from_source": n_seen,
        "source": os.path.abspath(halo_catalog),
        "columns": HALO_COLUMNS,
        "dtype": [[name, np.dtype(COMPACT_DTYPE[name]).str] for name in COMPACT_DTYPE.names],
        "row_bytes": COMPACT_DTYPE.itemsize,
        "selection": selection,
    }
    with open(f"{output}.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)
        f.write("\n")


def prepare_halos(
    halo_catalog: str,
    output: str,
    log10_mvir: float,
    half_width: float,
    chunksize: int,
    max_rows: int | None = None,
    min_log10_mvir: float | None = None,
) -> int:
    ensure_parent(output)
    if os.path.exists(output):
        os.remove(output)

    if halo_catalog.endswith(".bin"):
        source_floor = read_bin_sidecar(halo_catalog)["selection"].get("min_log10_mvir")
        requested_low = min_log10_mvir if min_log10_mvir is not None else log10_mvir - half_width
        if source_floor is not None and requested_low < source_floor - 1e-9:
            raise ValueError(
                f"Requested selection reaches down to log10(Mvir)={requested_low}, but the "
                f"compact input catalog only contains halos above its floor of {source_floor}. "
                "The result would be silently incomplete."
            )

    binary = output.endswith(".bin")
    selection = (
        {"mode": "floor", "min_log10_mvir": min_log10_mvir, "upid": -1}
        if min_log10_mvir is not None
        else {"mode": "bin", "log10_mvir": log10_mvir, "half_width": half_width, "upid": -1}
    )
    logger.info("Selection: %s -> %s (%s)", selection, output, "binary" if binary else "csv")

    n_written = 0
    n_seen = 0
    header = True
    t0 = time.time()
    out_handle = open(output, "wb") if binary else None

    try:
        for chunk_id, chunk in enumerate(iter_halo_chunks(halo_catalog, chunksize, max_rows), start=1):
            n_seen += len(chunk)
            selected = select_hosts(chunk, log10_mvir, half_width, min_log10_mvir)

            if len(selected) > 0:
                if binary:
                    out_handle.write(chunk_to_compact(selected).tobytes())
                else:
                    selected = selected.copy()
                    selected["id"] = selected["id"].astype("int64")
                    selected["upid"] = selected["upid"].astype("int64")
                    selected.to_csv(output, mode="a", index=False, header=header)
                    header = False
                n_written += len(selected)

            if chunk_id == 1 or chunk_id % 10 == 0:
                elapsed = time.time() - t0
                logger.info(
                    "Processed %.3g rows, selected %d halos, %.0f rows/s",
                    n_seen,
                    n_written,
                    n_seen / max(elapsed, 1e-9),
                )
    finally:
        if out_handle is not None:
            out_handle.close()

    if not binary and header:
        pd.DataFrame(columns=HALO_COLUMNS).to_csv(output, index=False)
    if binary:
        write_sidecar(output, n_written, n_seen, halo_catalog, selection)

    elapsed = time.time() - t0
    logger.info(
        "Wrote %d host halos to %s after reading %d rows in %.1f s (%.0f rows/s)",
        n_written,
        output,
        n_seen,
        elapsed,
        n_seen / max(elapsed, 1e-9),
    )
    return n_written


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a compact mass-selected host-halo catalog from a BigMDPL Rockstar hlist.",
    )
    parser.add_argument(
        "--halo-catalog",
        default="analysis/sim/results/hosts_minmass12p5.bin",
        help="Input catalog: raw Rockstar hlist (text) or compact .bin catalog.",
    )
    parser.add_argument(
        "--output",
        default="analysis/sim/results/halos_mass13.csv",
        help="Output compact halo CSV.",
    )
    parser.add_argument(
        "--log10-mvir",
        type=float,
        default=13.0,
        help="Center of log10(Mvir / h^-1 Msun) mass cut.",
    )
    parser.add_argument(
        "--half-width",
        type=float,
        default=0.05,
        help="Half-width of log10(Mvir) selection.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=1_000_000,
        help="Rows per pandas chunk.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional row cap for smoke tests.",
    )
    parser.add_argument(
        "--min-log10-mvir",
        type=float,
        default=None,
        help=(
            "Floor-mode selection: keep all hosts with log10(Mvir) >= this value, "
            "ignoring --log10-mvir/--half-width. Use with a .bin output for the "
            "one-time hlist compaction."
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()
    prepare_halos(
        halo_catalog=args.halo_catalog,
        output=args.output,
        log10_mvir=args.log10_mvir,
        half_width=args.half_width,
        chunksize=args.chunksize,
        max_rows=args.max_rows,
        min_log10_mvir=args.min_log10_mvir,
    )


if __name__ == "__main__":
    main()

