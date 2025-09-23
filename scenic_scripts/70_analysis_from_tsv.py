#!/usr/bin/env python3
"""
70_analysis_from_tsv.py

Purpose
-------
Post-process AUCell scores (TSV) with metadata from h5ad:
- Merge AUCell regulon scores with obs metadata (e.g., louvain).
- Aggregate by group (default: louvain) to group-by-regulon matrix.
- Save per-cell table, group-mean matrix, and z-score matrix.
"""

import argparse
import numpy as np
import pandas as pd
import h5py

# -------------------
# Parse arguments
# -------------------
ap = argparse.ArgumentParser(description="Summarize AUCell scores TSV + h5ad metadata")
ap.add_argument("--scores_tsv", required=True, help="Path to AUCell scores .tsv")
ap.add_argument("--h5ad_in", required=True, help="Input h5ad (for obs metadata)")
ap.add_argument("--groupby", default="louvain", help="obs column to aggregate by")
ap.add_argument("--out_dir", required=True, help="Output directory")
ap.add_argument("--prefix", required=True, help="Prefix for output files")
args = ap.parse_args()

# -------------------
# Read AUCell scores
# -------------------
print(f"[analysis] Reading AUCell scores TSV: {args.scores_tsv}")
scores = pd.read_csv(args.scores_tsv, sep="\t", index_col=0)
print(f"[analysis] Scores shape: {scores.shape}")

# -------------------
# Read metadata from h5ad
# -------------------
print(f"[analysis] Reading metadata from h5ad (obs/{args.groupby}): {args.h5ad_in}")
with h5py.File(args.h5ad_in, "r") as f:
    cell_ids = f["obs"]["_index"][:].astype(str)

    if args.groupby.encode() not in f["obs"]:
        raise KeyError(f"obs/{args.groupby} not found. Available: {list(f['obs'].keys())}")

    grp = f["obs"][args.groupby]
    if isinstance(grp, h5py.Group) and "codes" in grp and "categories" in grp:
        codes = grp["codes"][:]
        cats  = grp["categories"][:].astype(str)
        groups = pd.Categorical.from_codes(codes, cats).astype(str)
    else:
        groups = grp[:].astype(str)

meta = pd.DataFrame({args.groupby: groups}, index=cell_ids)
print(f"[analysis] Meta shape: {meta.shape}")

# -------------------
# Merge scores + metadata
# -------------------
df = scores.join(meta, how="inner")
print(f"[analysis] Merged shape: {df.shape}")

# Save per-cell table
per_cell_path = f"{args.out_dir}/{args.prefix}_per_cell_selected.tsv.gz"
df.to_csv(per_cell_path, sep="\t", index=True, compression="gzip")
print(f"[analysis] Wrote per-cell selected table: {per_cell_path}")

# -------------------
# Aggregate to group means
# -------------------
if args.groupby not in df.columns:
    raise KeyError(f"Column '{args.groupby}' not found in merged dataframe.")

numeric_cols = [c for c in df.columns if c != args.groupby]
mean_df = df.groupby(args.groupby)[numeric_cols].mean()
result = mean_df.transpose()

out_mat = f"{args.out_dir}/{args.prefix}_{args.groupby}_mean_regulons.tsv.gz"
result.to_csv(out_mat, sep="\t", index=True, compression="gzip")
print(f"[analysis] Wrote {args.groupby}-mean matrix: {out_mat} (shape {result.shape})")

# -------------------
# Z-score by row
# -------------------
means = result.mean(axis=1)
stds  = result.std(axis=1).replace(0, np.nan)
zscore = result.sub(means, axis=0).div(stds, axis=0)

out_z = f"{args.out_dir}/{args.prefix}_{args.groupby}_mean_regulons_zscore.tsv.gz"
zscore.to_csv(out_z, sep="\t", index=True, compression="gzip")
print(f"[analysis] Wrote z-score matrix: {out_z}")
