#!/usr/bin/env python3
"""
25_subset_adults.py

Subset liftover .h5ad files to adults only (age > 1).
Does the age filtering with h5py (lightweight) to avoid loading the full matrix.
"""

import argparse
import h5py
import numpy as np
import scanpy as sc
from pathlib import Path

BASE = Path("/scratch/easmit31/GRN_copy/scenic/h5ad_files")

parser = argparse.ArgumentParser()
parser.add_argument("--celltype", required=True,
                    choices=["gaba", "astros"])
args = parser.parse_args()

if args.celltype == "gaba":
    infile = BASE / "gaba_allLv_liftover.h5ad"
    outfile = BASE / "gaba_adults_allLv.h5ad"
else:
    infile = BASE / "astros_allLv_nombcb_liftover.h5ad"
    outfile = BASE / "astros_adults_allLv.h5ad"

print(f"[subset] Checking ages in {infile}")
with h5py.File(infile, "r") as f:
    ages = f["obs/age"][:]
    try:
        ages = np.array([float(a.decode("utf-8")) for a in ages])
    except AttributeError:
        ages = ages.astype(float)

mask = ages > 1
keep_idx = np.where(mask)[0]
print(f"[subset] Keeping {len(keep_idx)} of {len(ages)} cells (>{1} years)")

# Now reload with scanpy, but slice only adults
adata = sc.read_h5ad(infile)
adata = adata[keep_idx, :].copy()
adata.write(outfile, compression="gzip")
print(f"[subset] Wrote: {outfile}")
