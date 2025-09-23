#!/usr/bin/env python3
"""
96_cluster_regulon_heatmap.py

Purpose
-------
Plot AUCell regulon activity heatmaps.
- If --cluster is numeric → heatmap of cells in that cluster.
- If --cluster all → single summary heatmap of all clusters (regulons x clusters).
"""

import argparse
from pathlib import Path
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="AUCell .h5ad file")
parser.add_argument("--cluster", required=True, help="Cluster ID to plot, or 'all'")
parser.add_argument("--celltype", required=True, help="Cell type prefix for filenames")
parser.add_argument("--groupby", default="louvain", help="obs column to group by (default: louvain)")
parser.add_argument("--max_regulons", type=int, default=50, help="Top regulons to show")
parser.add_argument("--figdir", required=True, help="Directory to save figures")
args = parser.parse_args()

input_path = Path(args.input)
FIG_DIR = Path(args.figdir)
FIG_DIR.mkdir(parents=True, exist_ok=True)

print(f"[heatmap] Reading: {input_path}")
adata = sc.read_h5ad(input_path)

if args.groupby not in adata.obs:
    raise ValueError(f"Groupby '{args.groupby}' not found in obs")

# regulon columns
reg_cols = [c for c in adata.obs.columns if c.startswith("AUC_")]

if args.cluster == "all":
    print("[heatmap] Building summary heatmap across all clusters")
    # average activity per cluster
    cluster_means = adata.obs.groupby(args.groupby)[reg_cols].mean()

    # rank regulons by variance across clusters
    variances = cluster_means.var(axis=0)
    top_regs = variances.sort_values(ascending=False).head(args.max_regulons).index

    plt.figure(figsize=(10, 0.25*len(top_regs)))
    sns.heatmap(cluster_means[top_regs].T, cmap="viridis", cbar_kws={"label": "Mean AUC"})
    plt.title(f"{args.celltype}: Top {args.max_regulons} variable regulons across clusters")
    plt.xlabel("Clusters")
    plt.ylabel("Regulons")
    plt.tight_layout()
    out_file = FIG_DIR / f"{args.celltype}_allclusters_regulons.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"[heatmap] Saved {out_file}")

else:
    # numeric cluster → keep old behavior
    cells = adata[adata.obs[args.groupby] == str(args.cluster), :]
    reg_means = cells.obs[reg_cols].mean(axis=0).sort_values(ascending=False)
    ordered_regs = reg_means.index[:args.max_regulons].tolist()

    plt.figure(figsize=(10, 6))
    sns.heatmap(cells.obs[ordered_regs].T, cmap="viridis", cbar_kws={"label": "AUC"})
    plt.title(f"{args.celltype} cluster {args.cluster} regulons")
    plt.xlabel("Cells")
    plt.ylabel("Regulons (ordered by mean AUC)")
    plt.tight_layout()
    out_file = FIG_DIR / f"{args.celltype}_cluster{args.cluster}_regulons.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"[heatmap] Saved {out_file}")
