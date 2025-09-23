#!/usr/bin/env python3
"""
85_visualize_aucell.py
Memory-safe summary plots of AUCell results.
- Barplots of most variable and most active regulons
- Computes stats one column at a time to avoid 'Killed'
"""

import argparse
from pathlib import Path
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="AUCell .h5ad file")
parser.add_argument("--figdir", required=False, help="Output figures directory (default: figures_all)")
args = parser.parse_args()

# derive prefix from input filename
input_path = Path(args.input)
stem = input_path.stem.replace("_pyscenic_output_merged", "").replace("_pyscenic_output", "")
prefix = stem

# figure directory
FIG_DIR = Path(args.figdir) if args.figdir else (input_path.parent.parent / "figures_all")
FIG_DIR.mkdir(parents=True, exist_ok=True)

print(f"[viz] Reading (backed mode): {input_path}")
adata = sc.read_h5ad(input_path, backed="r")

# regulon columns
reg_cols = [c for c in adata.obs.columns if c.startswith("AUC_")]

means = []
variances = []

print(f"[viz] Calculating mean/var for {len(reg_cols)} regulons (streaming)...")
for c in reg_cols:
    col = adata.obs[c].to_numpy()  # load 1 regulon column into memory
    means.append(np.mean(col))
    variances.append(np.var(col))

means = np.array(means)
variances = np.array(variances)

# --- Top variable ---
top_var_idx = np.argsort(variances)[::-1][:30]
top_var_regs = [reg_cols[i] for i in top_var_idx]

plt.figure(figsize=(8, 6))
sns.barplot(x=variances[top_var_idx], y=top_var_regs, orient="h")
plt.title("Top 30 Most Variable Regulons")
plt.tight_layout()
plt.savefig(FIG_DIR / f"{prefix}_summary_most_variable.png", dpi=300)
plt.close()

# --- Top active ---
top_mean_idx = np.argsort(means)[::-1][:30]
top_mean_regs = [reg_cols[i] for i in top_mean_idx]

plt.figure(figsize=(8, 6))
sns.barplot(x=means[top_mean_idx], y=top_mean_regs, orient="h")
plt.title("Top 30 Most Active Regulons (mean AUC)")
plt.tight_layout()
plt.savefig(FIG_DIR / f"{prefix}_summary_most_active.png", dpi=300)
plt.close()

print(f"[viz] Saved figures to {FIG_DIR}")
