#!/usr/bin/env python3
"""
00_make_anndata.py

Convert a 10X-style folder (matrix.mtx, barcodes.tsv, genes.tsv)
into an AnnData .h5ad file.
"""

import argparse
import os
import pandas as pd
import scanpy as sc
from scipy import io

parser = argparse.ArgumentParser(description="Convert 10X tsvs to AnnData h5ad")
parser.add_argument("input_dir", help="Directory containing matrix.mtx, barcodes.tsv, genes.tsv")
parser.add_argument("--outdir", default=".", help="Where to write the .h5ad (default=current dir)")
args = parser.parse_args()

indir = args.input_dir.rstrip("/")
folder_name = os.path.basename(indir)

print(f"Reading 10X files from: {indir}")
barcodes = pd.read_csv(os.path.join(indir, "barcodes.tsv"), header=None)
print(f"  Found {len(barcodes)} cells")
genes = pd.read_csv(os.path.join(indir, "genes.tsv"), header=None)
print(f"  Found {len(genes)} genes")
mtx = io.mmread(os.path.join(indir, "matrix.mtx")).tocsr()
print(f"  Matrix shape: {mtx.shape} (genes x cells)")

print("Creating AnnData...")
adata = sc.AnnData(X=mtx.T)
adata.obs_names = barcodes.iloc[:, 0].astype(str)
adata.var_names = genes.iloc[:, 0].astype(str)
print(f"  AnnData: {adata.shape} (cells x genes)")

os.makedirs(args.outdir, exist_ok=True)
out_path = os.path.join(args.outdir, f"{folder_name}_raw.h5ad")
adata.write_h5ad(out_path)
print(f"Saved: {out_path}")
