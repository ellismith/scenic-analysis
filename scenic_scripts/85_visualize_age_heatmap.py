#!/usr/bin/env python3
"""
85_visualize_age_heatmap.py
---------------------------
Plot heatmap of regulon activity averaged across groups (e.g. age).
"""

import argparse
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to SCENIC output .h5ad")
    parser.add_argument("--groupby", required=True, help="obs column to group by (e.g. age)")
    parser.add_argument("--out", required=True, help="Path to save PNG")
    args = parser.parse_args()

    print(f"[heatmap] Reading: {args.input}")
    adata = sc.read_h5ad(args.input)

    if args.groupby not in adata.obs.columns:
        raise ValueError(f"{args.groupby} not in obs. Available: {list(adata.obs.columns)}")

    reg_cols = [c for c in adata.obs.columns if c.startswith("Regulon(")]
    if not reg_cols:
        raise ValueError("No regulon columns found in obs.")

    df = adata.obs[reg_cols].copy()
    df[args.groupby] = adata.obs[args.groupby]

    mean_by_group = df.groupby(args.groupby).mean().T
    plt.figure(figsize=(10, max(6, 0.2 * mean_by_group.shape[0])))
    sns.heatmap(mean_by_group, cmap="viridis")
    plt.title(f"Average regulon activity by {args.groupby}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"[heatmap] Saved: {args.out}")

if __name__ == "__main__":
    main()
