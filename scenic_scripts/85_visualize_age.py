#!/usr/bin/env python3
"""
85_visualize_age.py
-------------------
Visualize AUCell regulon activity with cells colored by a metadata column (e.g. age).
- Loads SCENIC output .h5ad
- Runs PCA/UMAP on regulon activity
- Saves a UMAP colored by the chosen metadata column
"""

import argparse
import scanpy as sc
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Visualize AUCell regulon activity by metadata")
    parser.add_argument("--input", required=True, help="Path to SCENIC output .h5ad")
    parser.add_argument("--colorby", required=True, help="obs column to color by (e.g. age)")
    parser.add_argument("--out", required=True, help="Path to save output PNG")
    args = parser.parse_args()

    print(f"[viz] Reading: {args.input}")
    adata = sc.read_h5ad(args.input)

    if args.colorby not in adata.obs.columns:
        raise ValueError(f"Column {args.colorby} not found in .obs. Available: {list(adata.obs.columns)}")

    # Regulon columns usually in .obs starting with "Regulon("
    reg_cols = [c for c in adata.obs.columns if c.startswith("Regulon(")]
    if not reg_cols:
        raise ValueError("No regulon columns found in .obs. Check your SCENIC output.")

    # Dimensionality reduction
    sc.pp.scale(adata, zero_center=True, max_value=3)
    sc.tl.pca(adata, use_highly_variable=False, svd_solver='arpack')
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata)

    # Plot
    print(f"[viz] Plotting UMAP colored by: {args.colorby}")
    sc.pl.umap(adata, color=args.colorby, save=None, show=False)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    sc.pl.umap(adata, color=args.colorby, show=False)
    import matplotlib.pyplot as plt
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"[viz] Saved to: {args.out}")

if __name__ == "__main__":
    main()
