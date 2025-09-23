#!/usr/bin/env python3
"""
30_convert_to_loom.py
Convert a filtered / liftover AnnData object into a `.loom` file for use with pySCENIC.
- Reads the subsetted H5AD or Zarr file
- Prepares row (genes) and column (cells) attributes
- Saves a `.loom` in the working directory
"""

import argparse
import os
import numpy as np
import anndata as ad
import loompy as lp
import scipy.sparse as sp

parser = argparse.ArgumentParser(description="Convert AnnData to loom for pySCENIC")
parser.add_argument("--adata_in", required=True,
                    help="Input AnnData (H5AD or Zarr)")
parser.add_argument("--loom_out", required=True,
                    help="Output loom path")
args = parser.parse_args()

print(f"[loom] Reading: {args.adata_in}")
if args.adata_in.endswith(".zarr") or os.path.isdir(args.adata_in):
    adata = ad.read_zarr(args.adata_in)
else:
    adata = ad.read_h5ad(args.adata_in)

print(f"[loom] Loaded AnnData: {adata.shape}")

row_attrs = {
    "Gene": np.array(adata.var.index, dtype=str),
}
col_attrs = {
    "CellID": np.array(adata.obs.index, dtype=str),
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
}

if os.path.exists(args.loom_out):
    os.remove(args.loom_out)

# ---- Safe patch: ensure float32 matrix ----
X = adata.X
if sp.issparse(X):
    X = X.astype(np.float32)
else:
    X = np.array(X, dtype=np.float32)

print(f"[loom] Writing loom: {args.loom_out}  (dtype={X.dtype}, shape={X.shape})")
lp.create(args.loom_out, X.transpose(), row_attrs, col_attrs)
print("[loom] Done.")
