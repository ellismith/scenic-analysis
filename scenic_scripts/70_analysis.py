#!/usr/bin/env python3
"""
70_analysis.py
Purpose
-------
Post-process pySCENIC AUCELL output:
- Read <CT>_pyscenic_output.h5ad
- Select metadata columns (e.g., 'louvain', 'age') and all 'Regulon' columns
- Aggregate by a grouping column (default: louvain) to get a group-by-regulon matrix
- Save:
    * per-cell table (selected cols)
    * group-mean matrix (regulons x groups)
    * simple z-score standardized version
"""

import argparse, re, numpy as np, pandas as pd, scanpy as sc
from config import CT, WORK_DIR, PYSCENIC_OUT

ap = argparse.ArgumentParser(description="Summarize AUCELL output to group-mean tables")
ap.add_argument("--h5ad_in", default=str(PYSCENIC_OUT), help="Input AUCELL .h5ad")
ap.add_argument("--out_dir", default=str(WORK_DIR), help="Output directory")
ap.add_argument("--meta_keep", nargs="+", default=["louvain", "age"], help="Metadata columns to keep")
ap.add_argument("--groupby", default="louvain", help="obs column to aggregate by (default: louvain)")
ap.add_argument("--prefix", default=CT, help="Prefix for output files (e.g. gaba_adults)")
args = ap.parse_args()

print(f"[analysis] Reading AUCELL h5ad: {args.h5ad_in}")
adata = sc.read_h5ad(args.h5ad_in)
df_obs = adata.obs.copy()

# Grab regulon columns
regulon_cols = [c for c in df_obs.columns if re.match(r'(?i)^regulon', c)]

missing_meta = [m for m in args.meta_keep if m not in df_obs.columns]
if missing_meta:
    print(f"[analysis] Warning: missing meta columns: {missing_meta}")

keep_cols = [m for m in args.meta_keep if m in df_obs.columns] + regulon_cols
df_select = df_obs[keep_cols].copy()

# Write per-cell selected table
per_cell_path = f"{args.out_dir}/{args.prefix}_per_cell_selected.tsv.gz"
df_select.to_csv(per_cell_path, sep="\t", index=True, compression="gzip")
print(f"[analysis] Wrote per-cell selected table: {per_cell_path}")

# Require grouping column
if args.groupby not in df_select.columns:
    raise KeyError(f"Column '{args.groupby}' not found; cannot aggregate.")

# Aggregate to group means
numeric_cols = [c for c in regulon_cols if pd.api.types.is_numeric_dtype(df_select[c])]
mean_df = df_select.groupby(args.groupby)[numeric_cols].mean()
result = mean_df.transpose()

out_mat = f"{args.out_dir}/{args.prefix}_{args.groupby}_mean_regulons.tsv.gz"
result.to_csv(out_mat, sep="\t", index=True, compression="gzip")
print(f"[analysis] Wrote {args.groupby}-mean matrix: {out_mat} (shape {result.shape})")

# Z-score by row
means = result.mean(axis=1)
stds  = result.std(axis=1).replace(0, np.nan)
zscore = result.sub(means, axis=0).div(stds, axis=0)

out_z = f"{args.out_dir}/{args.prefix}_{args.groupby}_mean_regulons_zscore.tsv.gz"
zscore.to_csv(out_z, sep="\t", index=True, compression="gzip")
print(f"[analysis] Wrote z-score matrix: {out_z}")
