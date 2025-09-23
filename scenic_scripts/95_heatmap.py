#!/usr/bin/env python3
"""
95_heatmap.py

Purpose
-------
Plot a heatmap of top variable regulons across groups (default: louvain clusters).

Inputs
------
- Z-score matrix (from 70_analysis_from_tsv.py)
Outputs
-------
- Heatmap PNG file
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------
# Parse arguments
# -------------------
ap = argparse.ArgumentParser(description="Plot regulon heatmap from z-score matrix")
ap.add_argument("--zscore", required=True, help="Path to z-score matrix (TSV.gz)")
ap.add_argument("--out_png", required=True, help="Path to output PNG")
ap.add_argument("--max_rows", type=int, default=80, help="Max regulons to plot (by variance)")
args = ap.parse_args()

print(f"[heatmap] Reading z-score matrix: {args.zscore}")
mat = pd.read_csv(args.zscore, sep="\t", index_col=0)

# -------------------
# Reorder clusters numerically if possible
# -------------------
try:
    mat.columns = mat.columns.astype(int)
    mat = mat[sorted(mat.columns)]
    mat.columns = mat.columns.astype(str)  # keep labels as strings for axis ticks
    print(f"[heatmap] Reordered clusters numerically: {list(mat.columns)}")
except ValueError:
    print("[heatmap] Louvain labels not numeric, keeping original order")

# -------------------
# Select top variable regulons
# -------------------
variances = mat.var(axis=1)
top_regs = variances.sort_values(ascending=False).head(args.max_rows).index
sub = mat.loc[top_regs]

print(f"[heatmap] Submatrix shape: {sub.shape} (top {args.max_rows} regulons)")

# -------------------
# Plot heatmap
# -------------------
plt.figure(figsize=(10, max(6, 0.2 * sub.shape[0])))

im = plt.imshow(
    sub.values,
    aspect="auto",
    cmap="RdBu_r",
    vmin=-2, vmax=2,
)

plt.colorbar(im, label="Regulon activity (z-score)")
plt.yticks(range(sub.shape[0]), sub.index, fontsize=6)
plt.xticks(range(sub.shape[1]), sub.columns, rotation=90)
plt.xlabel("Louvain clusters")
plt.ylabel("Regulons")
plt.title(f"Top {args.max_rows} variable regulons across clusters")

plt.tight_layout()
plt.savefig(args.out_png, dpi=150)
plt.close()
print(f"[heatmap] Wrote {args.out_png}")
