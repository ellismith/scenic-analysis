#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import umap

def main():
    ap = argparse.ArgumentParser(description="Run UMAP on AUCell matrix.")
    ap.add_argument("--input", required=True, help="AUCell CSV (cells x regulons)")
    ap.add_argument("--out-prefix", required=True, dest="out_prefix", help="prefix for outputs (png, csv)")
    ap.add_argument("--n-regs", type=int, default=100, dest="n_regs", help="top regulons by variance")
    ap.add_argument("--n-pcs", type=int, default=50, dest="n_pcs", help="number of PCs to use before UMAP")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    # Load AUCell (cells x regulons), set first column as index
    df = pd.read_csv(args.input)
    cell_col = df.columns[0]
    df = df.set_index(cell_col)

    # pick top regulons by variance
    vars_ = df.var(axis=0)
    top_regs = vars_.sort_values(ascending=False).head(args.n_regs).index
    X = df[top_regs].to_numpy(dtype=float)

    # PCA
    pca = PCA(n_components=args.n_pcs, random_state=args.seed)
    X_pca = pca.fit_transform(X)

    # UMAP
    reducer = umap.UMAP(random_state=args.seed)
    X_umap = reducer.fit_transform(X_pca)

    # save embedding
    coords = pd.DataFrame(X_umap, index=df.index, columns=["UMAP1","UMAP2"])
    coords.to_csv(f"{args.out_prefix}_umap.csv")
    print("Saved embedding:", coords.shape)

    # plot (random 5000 points if huge for a quick preview)
    plot_df = coords
    if len(plot_df) > 5000:
        plot_df = plot_df.sample(5000, random_state=args.seed)
    plt.figure(figsize=(7,6))
    plt.scatter(plot_df["UMAP1"], plot_df["UMAP2"], s=2, alpha=0.6)
    plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
    plt.title(f"AUCell UMAP ({args.n_regs} regs, {args.n_pcs} PCs)")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}_umap.png", dpi=200)
    print("Saved PNG")

if __name__ == "__main__":
    main()
