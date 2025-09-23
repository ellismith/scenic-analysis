#!/usr/bin/env python3
import argparse, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def find_cell_col(df):
    for cand in ("cell","Cell","barcode","Unnamed: 0","index"):
        if cand in df.columns:
            return cand
    # if no obvious column, assume first col is cell
    return df.columns[0]

def main():
    ap = argparse.ArgumentParser(description="Plot UMAP CSV, colored by metadata column if provided.")
    ap.add_argument("--coords", required=True, help="CSV with columns: cell, UMAP1, UMAP2 (index col name can vary)")
    ap.add_argument("--meta", default="cell_clusters.csv", help="metadata CSV (default: cell_clusters.csv)")
    ap.add_argument("--meta-col", default="louvain", help="metadata column to color by (default: louvain)")
    ap.add_argument("--out", required=True, help="output PNG")
    ap.add_argument("--sample", type=int, default=50000, help="max points to plot (subsample for speed)")
    ap.add_argument("--seed", type=int, default=0, help="random seed")
    args = ap.parse_args()

    coords = pd.read_csv(args.coords)
    cell_col = find_cell_col(coords)
    if "UMAP1" not in coords.columns or "UMAP2" not in coords.columns:
        raise SystemExit("coords CSV must have UMAP1 and UMAP2 columns.")
    if len(coords) > args.sample > 0:
        coords = coords.sample(args.sample, random_state=args.seed)

    if os.path.exists(args.meta):
        meta = pd.read_csv(args.meta, index_col=0)
        # join
        df = coords.rename(columns={cell_col:"cell"}).set_index("cell").join(meta, how="left")
        if args.meta_col not in df.columns:
            raise SystemExit(f"Column '{args.meta_col}' not found in {args.meta}.")
        groups = df[args.meta_col].astype(str)
        uniq = sorted(groups.dropna().unique())
        cmap = plt.cm.get_cmap("tab20", max(1, len(uniq)))
        color_map = {g: cmap(i % cmap.N) for i, g in enumerate(uniq)}
        colors = [color_map.get(g, (0.7,0.7,0.7,0.4)) for g in groups]
        plt.figure(figsize=(7.5,6.5))
        plt.scatter(df["UMAP1"], df["UMAP2"], s=2, c=colors, linewidths=0, alpha=0.9)
        plt.xlabel("UMAP1"); plt.ylabel("UMAP2"); plt.title(f"UMAP colored by {args.meta_col}")
        # tiny legend
        if len(uniq) <= 20:
            from matplotlib.patches import Patch
            handles = [Patch(facecolor=color_map[g], edgecolor='none', label=str(g)) for g in uniq]
            plt.legend(handles=handles, title=args.meta_col, loc="upper right", fontsize=7, frameon=False, markerscale=6)
        plt.tight_layout(); plt.savefig(args.out, dpi=200)
        print("Wrote", args.out)
    else:
        plt.figure(figsize=(7.5,6.5))
        plt.scatter(coords["UMAP1"], coords["UMAP2"], s=2, linewidths=0, alpha=0.7)
        plt.xlabel("UMAP1"); plt.ylabel("UMAP2"); plt.title("UMAP")
        plt.tight_layout(); plt.savefig(args.out, dpi=200)
        print("Wrote", args.out)

if __name__ == "__main__":
    main()
