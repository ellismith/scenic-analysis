#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser(description="Plot a compact AUCell heatmap (subset of cells × regulons).")
    ap.add_argument("--input", required=True, help="Path to auc_mtx.csv")
    ap.add_argument("--out",   required=True, help="Output image path, e.g. auc_heatmap.png")
    ap.add_argument("--n-cells", type=int, default=200, help="Number of cells to sample for the heatmap")
    ap.add_argument("--n-regs",  type=int, default=30,  help="Number of regulons (columns) to plot")
    ap.add_argument("--regulons", default="", help="Comma-separated regulon names to plot (overrides --n-regs)")
    ap.add_argument("--seed", type=int, default=0, help="Random seed for sampling cells")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    if df.shape[1] < 2:
        raise SystemExit("auc_mtx.csv looks empty or malformed.")
    cell_col = df.columns[0]
    df = df.set_index(cell_col)

    # Choose regulons to display
    if args.regulons.strip():
        regs = [r.strip() for r in args.regulons.split(",") if r.strip() in df.columns]
        if not regs:
            raise SystemExit("None of the requested regulons found in the file.")
    else:
        # most variable regulons to make the heatmap informative
        variances = df.var(axis=0)
        regs = list(variances.sort_values(ascending=False).head(args.n_regs).index)

    # Sample cells
    rng = np.random.default_rng(args.seed)
    if args.n_cells < len(df):
        idx = rng.choice(df.index.values, size=args.n_cells, replace=False)
        sub = df.loc[idx, regs]
    else:
        sub = df.loc[:, regs]

    # Normalize each regulon column to [0,1] for contrast (viz only)
    X = sub.values.astype(float)
    col_min = X.min(axis=0, keepdims=True)
    col_max = X.max(axis=0, keepdims=True)
    denom = np.where((col_max - col_min) == 0, 1.0, (col_max - col_min))
    Xn = (X - col_min) / denom

    # Plot (pure matplotlib)
    plt.figure(figsize=(min(22, 2 + 0.3*len(regs)), min(16, 2 + 0.12*len(sub))))
    im = plt.imshow(Xn, aspect="auto", interpolation="nearest")
    cbar = plt.colorbar(im, fraction=0.02, pad=0.02)
    cbar.set_label("AUCell (scaled 0–1 per regulon)")
    plt.xticks(ticks=np.arange(len(regs)), labels=regs, rotation=60, ha="right", fontsize=8)
    if len(sub) <= 60:
        plt.yticks(ticks=np.arange(len(sub)), labels=sub.index, fontsize=7)
    else:
        plt.yticks([])
    plt.title(f"AUCell activity — {sub.shape[0]} cells × {len(regs)} regulons")
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
