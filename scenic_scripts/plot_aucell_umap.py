#!/usr/bin/env python3
import argparse, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def try_import_umap():
    try:
        import umap
        return umap
    except Exception:
        sys.stderr.write("ERROR: umap-learn is required. Try: pip install umap-learn\n")
        raise

def main():
    ap = argparse.ArgumentParser(description="UMAP on AUCell matrix (cells x regulons + metadata).")
    ap.add_argument("--input", required=True, help="AUCell scores TSV (tab-delimited, with metadata columns)")
    ap.add_argument("--out-prefix", required=True, help="Prefix for outputs (CSV/PNG)")
    ap.add_argument("--sample-cells", type=int, default=20000, help="Number of cells to sample")
    ap.add_argument("--pca", type=int, default=50, help="PCA components before UMAP")
    ap.add_argument("--neighbors", type=int, default=15, help="UMAP n_neighbors")
    ap.add_argument("--min-dist", type=float, default=0.1, help="UMAP min_dist")
    ap.add_argument("--metric", default="euclidean", help="UMAP metric")
    ap.add_argument("--seed", type=int, default=0, help="Random seed")
    ap.add_argument("--color-regulon", default="", help="Regulon name to color by")
    ap.add_argument("--color-meta", default="", help="Metadata column to color by (e.g. age, sex, louvain)")
    args = ap.parse_args()

    umap = try_import_umap()

    # --- Load TSV (tab-delimited!)
    df = pd.read_csv(args.input, sep="\t")
    df.columns = df.columns.str.strip()

    # --- Identify columns
    cell_col = df.columns[0]
    reg_cols = [c for c in df.columns if c.startswith("AUC_")]
    if not reg_cols:
        raise ValueError("No AUC_ regulon columns found in input TSV")

    # --- Build feature matrix
    X = df[reg_cols].astype("float32")
    cells = df[cell_col].astype(str).values

    # --- Sample cells
    rng = np.random.default_rng(args.seed)
    if args.sample_cells > 0 and args.sample_cells < len(df):
        idx = rng.choice(len(df), size=args.sample_cells, replace=False)
        X = X.iloc[idx].to_numpy(copy=False)
        cells = cells[idx]
        df_color = df.iloc[idx]
    else:
        X = X.to_numpy(copy=False)
        df_color = df

    # --- Optional PCA
    if args.pca and args.pca > 0:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=min(args.pca, X.shape[1]), random_state=args.seed)
        emb_input = pca.fit_transform(X)
    else:
        emb_input = X

    # --- UMAP
    reducer = umap.UMAP(
        n_neighbors=args.neighbors,
        min_dist=args.min_dist,
        metric=args.metric,
        random_state=args.seed,
        n_components=2,
        verbose=False,
    )
    emb = reducer.fit_transform(emb_input).astype("float32")

    # --- Save embedding
    out_csv = f"{args.out_prefix}_umap.csv"
    pd.DataFrame({"cell": cells, "UMAP1": emb[:,0], "UMAP2": emb[:,1]}).to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    # --- Base scatter
    out_png = f"{args.out_prefix}_umap.png"
    plt.figure(figsize=(8,7))
    plt.scatter(emb[:,0], emb[:,1], s=2, linewidths=0, alpha=0.8)
    plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
    plt.title(f"UMAP of AUCell (n={len(cells)} cells)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print(f"Wrote {out_png}")

    # --- Color by regulon
    if args.color_regulon:
        colname = args.color_regulon
        if colname in df_color.columns:
            vals = df_color[colname].astype("float32").values
            plt.figure(figsize=(8,7))
            sc = plt.scatter(emb[:,0], emb[:,1], c=vals, s=2, linewidths=0, alpha=0.9)
            cbar = plt.colorbar(sc); cbar.set_label(f"AUC: {colname}")
            plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
            plt.title(f"UMAP colored by {colname}")
            plt.tight_layout()
            out_col = f"{args.out_prefix}_umap_{colname}.png"
            plt.savefig(out_col, dpi=200)
            print(f"Wrote {out_col}")
        else:
            print(f"WARNING: regulon '{colname}' not found.", file=sys.stderr)

    # --- Color by metadata
    if args.color_meta:
        colname = args.color_meta
        if colname in df_color.columns:
            vals = df_color[colname].values
            plt.figure(figsize=(8,7))
            if np.issubdtype(vals.dtype, np.number):
                sc = plt.scatter(emb[:,0], emb[:,1], c=vals, s=2, linewidths=0, alpha=0.9, cmap="viridis")
                cbar = plt.colorbar(sc); cbar.set_label(colname)
            else:
                cats = pd.Categorical(vals)
                colors = pd.factorize(cats)[0]
                plt.scatter(emb[:,0], emb[:,1], c=colors, s=2, linewidths=0, alpha=0.9, cmap="tab10")
                plt.legend(handles=[
                    plt.Line2D([0],[0], marker='o', color='w',
                               markerfacecolor=plt.cm.tab10(i/10),
                               markersize=6, label=cat)
                    for i, cat in enumerate(cats.categories)
                ], bbox_to_anchor=(1.05,1), loc='upper left')
            plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
            plt.title(f"UMAP colored by {colname}")
            plt.tight_layout()
            out_col = f"{args.out_prefix}_umap_{colname}.png"
            plt.savefig(out_col, dpi=200)
            print(f"Wrote {out_col}")
        else:
            print(f"WARNING: metadata '{colname}' not found.", file=sys.stderr)

if __name__ == "__main__":
    main()
