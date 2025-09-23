#!/usr/bin/env python3
"""
99_file_info_summary.py

Summarize files in a directory by prefix without loading big datasets into memory.
- Uses h5py/loompy directly (avoids reading into memory with anndata).
- Attempts to peek at metadata (age, sex) in a memory-safe way.
"""

import os, sys
import h5py, loompy
import pandas as pd
from datetime import datetime

def safe_counts(values, max_unique=20):
    """Return value counts only if #unique is small, otherwise summarize."""
    series = pd.Series(values)
    uniq = series.unique()
    if len(uniq) <= max_unique:
        return series.value_counts().to_dict()
    else:
        return {"unique_values": len(uniq)}

def file_info(path):
    ctime = os.path.getctime(path)
    created = datetime.fromtimestamp(ctime).strftime("%Y-%m-%d %H:%M:%S")

    shape, columns, metadata = None, None, {}

    try:
        if path.endswith(".csv"):
            with open(path, "r") as f:
                header = f.readline().strip().split(",")
            columns = header
            try:
                row_count = sum(1 for _ in open(path)) - 1
                shape = (row_count, len(columns))
            except Exception:
                shape = ("unknown", len(columns))

        elif path.endswith(".h5ad"):
            with h5py.File(path, "r") as f:
                if "obs" in f and "var" in f:
                    n_obs = f["obs"].shape[0]
                    n_var = f["var"].shape[0]
                    shape = (n_obs, n_var)
                    columns = ["obs", "var"]
                    # peek at obs keys
                    for key in ["age", "Age", "sex", "Sex"]:
                        if f"obs/{key}" in f:
                            dset = f[f"obs/{key}"]
                            # only grab first 10k to avoid mem blowup
                            values = dset[:min(10000, dset.shape[0])]
                            metadata[key] = safe_counts(values)

        elif path.endswith(".loom"):
            ds = loompy.connect(path, mode="r")
            shape = ds.shape
            columns = list(ds.ra.keys())
            for key in ["age", "Age", "sex", "Sex"]:
                if key in ds.ca.keys():
                    values = ds.ca[key][:min(10000, ds.shape[1])]
                    metadata[key] = safe_counts(values)
            ds.close()

    except Exception as e:
        columns = [f"[Error reading file: {e}]"]

    return created, shape, columns, metadata

def main(folder, prefix):
    for fname in sorted(os.listdir(folder)):
        if not fname.startswith(prefix):
            continue
        path = os.path.join(folder, fname)
        if not os.path.isfile(path):
            continue
        created, shape, columns, metadata = file_info(path)
        print(f"File: {fname}")
        print(f"  Created: {created}")
        if shape:
            print(f"  Shape: {shape}")
        if columns:
            print(f"  Columns/Features: {columns[:20]}")
        if metadata:
            print(f"  Metadata:")
            for k, v in metadata.items():
                print(f"    {k}: {v}")
        print("-" * 60)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python 99_file_info_summary.py <folder> <prefix>")
        sys.exit(1)
    folder, prefix = sys.argv[1], sys.argv[2]
    main(folder, prefix)
