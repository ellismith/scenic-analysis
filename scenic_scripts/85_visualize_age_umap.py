#!/usr/bin/env python3
"""
85_visualize_age_umap.py
------------------------
Visualize AUCell regulon activity with UMAP colored by a metadata column (e.g. age).
"""

import argparse
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to SCENIC output .h5ad")
    parser.add_argument("--colorby", required=True, help="obs column to color by (e.g. age)")
    parser.add_argument("--out", required=True, help="Path to save PNG")
    args = parser.parse_args()

    print(f"[umap] Reading: {args.input}")
    adata = sc.read_h5ad(args.input)

    if args.colorby not in adata.obs.columns:
        raise ValueError(f"{args.colorby} not in obs. Available: {list(adata.obs.columns)}")

    # Dimensionality reduction
    sc.pp.scale(adata, zero_center=True, max_value=3)
    sc.tl.pca(adata, use_highly_variable=False, svd_solver="arpack")
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata)

    # Plot and save
    sc.pl.umap(adata, color=args.colorby, show=False)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"[umap] Saved: {args.out}")

if __name__ == "__main__":
    main()
