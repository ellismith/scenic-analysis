#!/usr/bin/env python3
"""
97_compare_variability.py
Compare top variable regulons between two AUCell outputs (e.g., GABA vs Glut).
"""

import argparse
import scanpy as sc
import pandas as pd

parser = argparse.ArgumentParser(description="Compare regulon variability between two h5ad files.")
parser.add_argument("h5ad1", help="Path to first .h5ad file")
parser.add_argument("h5ad2", help="Path to second .h5ad file")
parser.add_argument("--topN", type=int, default=30, help="Number of top variable regulons to compare")
parser.add_argument("--out", default="regulon_variance_comparison.csv", help="Output CSV file")
args = parser.parse_args()

# --- load two AUCell h5ad outputs ---
adata1 = sc.read_h5ad(args.h5ad1)
adata2 = sc.read_h5ad(args.h5ad2)

# --- identify regulon columns ---
reg_cols1 = [c for c in adata1.obs.columns if c.startswith("Regulon")]
reg_cols2 = [c for c in adata2.obs.columns if c.startswith("Regulon")]
common_regs = list(set(reg_cols1) & set(reg_cols2))

# --- compute variance for each regulon in each dataset ---
var1 = adata1.obs[common_regs].var().sort_values(ascending=False)
var2 = adata2.obs[common_regs].var().sort_values(ascending=False)

# --- rank them ---
rank1 = var1.rank(ascending=False, method="dense").astype(int)
rank2 = var2.rank(ascending=False, method="dense").astype(int)

# --- put into DataFrame for comparison ---
df = pd.DataFrame({
    "var1": var1,
    "rank1": rank1,
    "var2": var2,
    "rank2": rank2,
})

# --- select top N regulons ---
top1 = set(var1.head(args.topN).index)
top2 = set(var2.head(args.topN).index)
overlap = top1 & top2

print(f"Top {args.topN} in {args.h5ad1}: {len(top1)}")
print(f"Top {args.topN} in {args.h5ad2}: {len(top2)}")
print(f"Overlap ({len(overlap)}): {sorted(overlap)}")

# --- save table ---
df.to_csv(args.out)
print(f"Wrote comparison table to {args.out}")
