#!/usr/bin/env python3
"""
85_visualize_expr_umap.py
-------------------------
Quick UMAP from gene expression in a loom file, colored by a metadata column (e.g. age).
"""

import argparse
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to .loom file")
    parser.add_argument("--colorby", required=True, help="obs column to color by (e.g. age)")
    parser.add_argument("--out", required=True, help="Path to save PNG")
    args = parser.parse_args()

    print(f"[expr-umap] Reading: {args.input}")
    adata = sc.read_loom(args.input)

    if args.colorby not in adata.obs.columns:
        raise ValueError(f"{args.colorby} not in obs. Available: {list(adata.obs.columns)}")

    # Normalize and log-transform expression
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Dimensionality reduction
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Plot and save
    sc.pl.umap(adata, color=args.colorby, show=False)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"[expr-umap] Saved: {args.out}")

if __name__ == "__main__":
    main()
