#!/usr/bin/env python3
"""
97_compare_variability.py
Compare top variable regulons between two AUCell outputs using h5py
"""
import argparse
import h5py
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Compare regulon variability between two h5ad files.")
parser.add_argument("h5ad1", help="Path to first .h5ad file")
parser.add_argument("h5ad2", help="Path to second .h5ad file")
parser.add_argument("--topN", type=int, default=30, help="Number of top variable regulons to compare")
parser.add_argument("--out", default="regulon_variance_comparison.csv", help="Output CSV file")
args = parser.parse_args()

# Load regulon columns and compute variance using h5py
def get_regulon_variances(h5ad_path):
    with h5py.File(h5ad_path, 'r') as f:
        reg_cols = [k for k in f['obs'].keys() if k.startswith('AUC_')]
        variances = {}
        for col in reg_cols:
            vals = f['obs'][col][:]
            variances[col] = np.var(vals)
    return pd.Series(variances).sort_values(ascending=False)

print(f"Loading {args.h5ad1}...")
var1 = get_regulon_variances(args.h5ad1)
print(f"Loading {args.h5ad2}...")
var2 = get_regulon_variances(args.h5ad2)

# Find common regulons
common_regs = list(set(var1.index) & set(var2.index))
print(f"Found {len(common_regs)} common regulons")

# Filter to common and rank
var1_common = var1[common_regs].sort_values(ascending=False)
var2_common = var2[common_regs].sort_values(ascending=False)

rank1 = var1_common.rank(ascending=False, method="dense").astype(int)
rank2 = var2_common.rank(ascending=False, method="dense").astype(int)

# Create comparison dataframe
df = pd.DataFrame({
    "var1": var1_common,
    "rank1": rank1,
    "var2": var2_common,
    "rank2": rank2,
})

# Find top N in each
top1 = set(var1_common.head(args.topN).index)
top2 = set(var2_common.head(args.topN).index)
overlap = top1 & top2

print(f"Top {args.topN} in dataset1: {len(top1)}")
print(f"Top {args.topN} in dataset2: {len(top2)}")
print(f"Overlap ({len(overlap)}): {sorted(overlap)}")

df.to_csv(args.out)
print(f"Wrote comparison table to {args.out}")
