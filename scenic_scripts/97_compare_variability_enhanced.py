#!/usr/bin/env python3
"""
97_compare_variability.py
Compare top variable regulons between two AUCell outputs (e.g., GABA vs Glut).
Now supports both cell-level and cluster-level variance analysis.
"""
import argparse
import scanpy as sc
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="Compare regulon variability between two h5ad files.")
parser.add_argument("h5ad1", help="Path to first .h5ad file")
parser.add_argument("h5ad2", help="Path to second .h5ad file")
parser.add_argument("--topN", type=int, default=30, help="Number of top variable regulons to compare")
parser.add_argument("--out", default="regulon_variance_comparison.csv", help="Output CSV file")
parser.add_argument("--cluster_key", default="louvain", help="Cluster column name for cluster-level variance")
parser.add_argument("--both", action="store_true", help="Analyze both cell-level and cluster-level variance")
args = parser.parse_args()

# --- load two AUCell h5ad outputs ---
adata1 = sc.read_h5ad(args.h5ad1)
adata2 = sc.read_h5ad(args.h5ad2)

# --- identify regulon columns ---
reg_cols1 = [c for c in adata1.obs.columns if c.startswith("Regulon")]
reg_cols2 = [c for c in adata2.obs.columns if c.startswith("Regulon")]
common_regs = list(set(reg_cols1) & set(reg_cols2))

print(f"Found {len(common_regs)} common regulons between datasets")

def compute_cell_variance(adata, regulon_cols):
    """Compute variance across all cells for each regulon"""
    return adata.obs[regulon_cols].var().sort_values(ascending=False)

def compute_cluster_variance(adata, regulon_cols, cluster_key):
    """Compute variance across cluster means for each regulon"""
    if cluster_key not in adata.obs.columns:
        print(f"Warning: {cluster_key} not found in dataset, skipping cluster variance")
        return None
    
    # Calculate mean regulon activity per cluster
    cluster_means = adata.obs.groupby(cluster_key)[regulon_cols].mean()
    # Calculate variance across cluster means
    return cluster_means.var().sort_values(ascending=False)

# --- compute cell-level variance (original behavior) ---
cell_var1 = compute_cell_variance(adata1, common_regs)
cell_var2 = compute_cell_variance(adata2, common_regs)

cell_rank1 = cell_var1.rank(ascending=False, method="dense").astype(int)
cell_rank2 = cell_var2.rank(ascending=False, method="dense").astype(int)

# Initialize results DataFrame
df = pd.DataFrame({
    "cell_var1": cell_var1,
    "cell_rank1": cell_rank1,
    "cell_var2": cell_var2,
    "cell_rank2": cell_rank2,
})

# --- compute cluster-level variance if requested or if both flag is set ---
if args.both or args.cluster_key:
    cluster_var1 = compute_cluster_variance(adata1, common_regs, args.cluster_key)
    cluster_var2 = compute_cluster_variance(adata2, common_regs, args.cluster_key)
    
    if cluster_var1 is not None and cluster_var2 is not None:
        cluster_rank1 = cluster_var1.rank(ascending=False, method="dense").astype(int)
        cluster_rank2 = cluster_var2.rank(ascending=False, method="dense").astype(int)
        
        # Add cluster variance columns
        df["cluster_var1"] = cluster_var1
        df["cluster_rank1"] = cluster_rank1
        df["cluster_var2"] = cluster_var2
        df["cluster_rank2"] = cluster_rank2

# --- analyze top N regulons for cell-level variance ---
cell_top1 = set(cell_var1.head(args.topN).index)
cell_top2 = set(cell_var2.head(args.topN).index)
cell_overlap = cell_top1 & cell_top2

print(f"\n=== CELL-LEVEL VARIANCE ===")
print(f"Top {args.topN} in {args.h5ad1}: {len(cell_top1)}")
print(f"Top {args.topN} in {args.h5ad2}: {len(cell_top2)}")
print(f"Overlap ({len(cell_overlap)}): {sorted(list(cell_overlap))[:10]}{'...' if len(cell_overlap) > 10 else ''}")

# --- analyze top N regulons for cluster-level variance if available ---
if 'cluster_var1' in df.columns:
    cluster_top1 = set(cluster_var1.head(args.topN).index)
    cluster_top2 = set(cluster_var2.head(args.topN).index)
    cluster_overlap = cluster_top1 & cluster_top2
    
    print(f"\n=== CLUSTER-LEVEL VARIANCE ===")
    print(f"Top {args.topN} in {args.h5ad1}: {len(cluster_top1)}")
    print(f"Top {args.topN} in {args.h5ad2}: {len(cluster_top2)}")
    print(f"Overlap ({len(cluster_overlap)}): {sorted(list(cluster_overlap))[:10]}{'...' if len(cluster_overlap) > 10 else ''}")
    
    # --- compare cell vs cluster rankings ---
    print(f"\n=== CELL vs CLUSTER COMPARISON ===")
    cell_cluster_overlap1 = cell_top1 & cluster_top1
    cell_cluster_overlap2 = cell_top2 & cluster_top2
    print(f"Dataset 1 - Cell/Cluster overlap: {len(cell_cluster_overlap1)}")
    print(f"Dataset 2 - Cell/Cluster overlap: {len(cell_cluster_overlap2)}")

# --- save table ---
df.to_csv(args.out)
print(f"\nWrote comparison table to {args.out}")
print(f"Columns: {list(df.columns)}")
