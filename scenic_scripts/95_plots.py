#!/usr/bin/env python3
"""
95_plots.py

Purpose
-------
Create summary plots from SCENIC outputs:
1) Heatmap of z-scored louvain-mean regulon activity (from step 70)
2) Partial dependence plots for top-N regulons by GAM R² (from step 90)

Inputs (defaults from config.py)
--------------------------------
- <CT>_pyscenic_output.h5ad               (AUCELL output; for PD plots)
- <CT>_louvain_mean_regulons_zscore.tsv.gz (z-scored matrix; for heatmap)
- <CT>_gam_r2.tsv.gz                       (GAM results; for picking top regulons)

Outputs
-------
- <CT>_heatmap_regulons.png
- <CT>_pdep_<sanitized_regulon>.png  (one per top regulon)

Usage
-----
python 95_plots.py
# or:
python 95_plots.py --top_n 6
"""

import os, re, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from pygam import LinearGAM, s
from config import CT, WORK_DIR, PYSCENIC_OUT

def sanitize(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', name)

def plot_heatmap(zscore_path: str, out_png: str, max_rows: int = 80):
    mat = pd.read_csv(zscore_path, sep="\t", index_col=0)
    # pick top rows by variance to keep figure readable
    variances = mat.var(axis=1)
    top = variances.sort_values(ascending=False).head(max_rows).index
    sub = mat.loc[top]
    plt.figure(figsize=(min(16, 2 + 0.25*sub.shape[1]), min(20, 2 + 0.25*sub.shape[0])))
    im = plt.imshow(sub.values, aspect="auto", interpolation="nearest")
    plt.colorbar(im, fraction=0.03, pad=0.04, label="z-score")
    plt.yticks(range(len(sub.index)), sub.index, fontsize=6)
    plt.xticks(range(len(sub.columns)), sub.columns, rotation=90, fontsize=8)
    plt.title(f"{CT}: Regulon activity (z-scored) by louvain")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_partial_dependence(adata_path: str, gam_table: str, out_dir: str, top_n: int = 5, meta_col: str = "age"):
    adata = sc.read_h5ad(adata_path)
    obs = adata.obs.copy()
    if meta_col not in obs.columns:
        raise KeyError(f"obs['{meta_col}'] not found in {adata_path}")
    y = pd.to_numeric(obs[meta_col], errors="coerce")
    valid = y.notna()
    y = y[valid]

    gam_df = pd.read_csv(gam_table, sep="\t")
    gam_df = gam_df.sort_values("r2", ascending=False)
    regs = [r for r in gam_df["regulon"].head(top_n).tolist() if r in obs.columns]
    if not regs:
        print("[plots] No overlapping regulons between GAM table and obs columns.")
        return

    for i, reg in enumerate(regs, 1):
        x = pd.to_numeric(obs.loc[valid, reg], errors="coerce")
        mask = x.notna()
        X = x[mask].values.reshape(-1, 1)
        yy = y[mask].values
        if len(yy) < 10:
            print(f"[plots] Skipping {reg}: too few points ({len(yy)})")
            continue

        gam = LinearGAM(s(0)).fit(X, yy)
        XX = np.linspace(np.nanmin(X), np.nanmax(X), 200).reshape(-1, 1)
        pdp = gam.partial_dependence(XX)
        lo, hi = gam.prediction_intervals(XX, width=0.95)

        plt.figure(figsize=(6,4))
        plt.plot(XX[:,0], pdp, label="Partial effect")
        plt.fill_between(XX[:,0], lo, hi, alpha=0.3, label="95% CI")
        plt.xlabel(reg)
        plt.ylabel(f"Effect on predicted {meta_col}")
        plt.title(f"{CT}: Partial dependence — {reg}")
        plt.legend()
        plt.grid(True, alpha=0.2)
        plt.tight_layout()
        out_png = os.path.join(out_dir, f"{CT}_pdep_{sanitize(reg)}.png")
        plt.savefig(out_png, dpi=300)
        plt.close()
        print(f"[plots] Wrote PD plot: {out_png}")

def main():
    ap = argparse.ArgumentParser(description="Make heatmaps and partial dependence plots")
    ap.add_argument("--zscore_mat", default=str(WORK_DIR / f"{CT}_louvain_mean_regulons_zscore.tsv.gz"))
    ap.add_argument("--gam_table",  default=str(WORK_DIR / f"{CT}_gam_r2.tsv.gz"))
    ap.add_argument("--h5ad_in",    default=str(PYSCENIC_OUT))
    ap.add_argument("--out_dir",    default=str(WORK_DIR))
    ap.add_argument("--top_n",      type=int, default=5)
    args = ap.parse_args()

    # Heatmap
    heatmap_png = os.path.join(args.out_dir, f"{CT}_heatmap_regulons.png")
    if os.path.exists(args.zscore_mat):
        print(f"[plots] Heatmap from: {args.zscore_mat}")
        plot_heatmap(args.zscore_mat, heatmap_png)
        print(f"[plots] Wrote heatmap: {heatmap_png}")
    else:
        print(f"[plots] Skipping heatmap (not found): {args.zscore_mat}")

    # Partial dependence plots
    if os.path.exists(args.gam_table):
        print(f"[plots] PD plots from GAM table: {args.gam_table} (top_n={args.top_n})")
        plot_partial_dependence(args.h5ad_in, args.gam_table, args.out_dir, top_n=args.top_n, meta_col="age")
    else:
        print(f"[plots] Skipping PD plots (GAM table not found): {args.gam_table}")

if __name__ == "__main__":
    main()
