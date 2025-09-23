#!/usr/bin/env python3
"""
26_subset_adults_full.py

Subset *full liftover* .h5ad files to adults only (age > 1).
Writes to a new file with '_adults_full' in the name.
"""

import argparse, h5py, numpy as np, scanpy as sc
from pathlib import Path

BASE = Path("/scratch/easmit31/GRN_copy/scenic/h5ad_files")

parser = argparse.ArgumentParser()
parser.add_argument("--celltype", required=True, choices=["gaba", "astros"])
args = parser.parse_args()

if args.celltype == "gaba":
    infile = BASE / "gaba_allLv_liftover.h5ad"
    outfile = BASE / "gaba_adults_full.h5ad"
elif args.celltype == "astros":
    infile = BASE / "astros_allLv_nombcb_liftover.h5ad"
    outfile = BASE / "astros_adults_full.h5ad"

print(f"[subset] Reading metadata from {infile}")
with h5py.File(infile, "r") as f:
    ages = f["obs/age"][:]
    try:
        ages = np.array([float(a.decode("utf-8")) for a in ages])
    except AttributeError:
        ages = ages.astype(float)

mask = ages > 1
keep_idx = np.where(mask)[0]
print(f"[subset] Adults: {len(keep_idx)} / {len(ages)} cells")

# Reload & subset properly (need RAM)
adata = sc.read_h5ad(infile)
adata = adata[keep_idx, :].copy()

print(f"[subset] Final shape: {adata.shape}")
adata.write(outfile, compression="gzip")
print(f"[subset] Wrote {outfile}")
