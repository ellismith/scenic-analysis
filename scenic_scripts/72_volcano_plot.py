#!/scratch/easmit31/conda_envs/pyscenic/bin/python
"""
72_volcano_plot.py

Make a volcano plot from regulon linear model results.

Usage:
    ./72_volcano_plot.py LM_TSV_PATH TERM FIGURES_DIR
Example:
    ./72_volcano_plot.py /path/to/results.tsv age figures_all
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

if len(sys.argv) != 4:
    sys.exit("Usage: ./72_volcano_plot.py LM_TSV_PATH TERM FIGURES_DIR")

lm_file = sys.argv[1]
term = sys.argv[2]
figures_dir = sys.argv[3]

# Output directory
outdir = Path(figures_dir)
outdir.mkdir(parents=True, exist_ok=True)
outfile = outdir / f"volcano_{term}.png"

# Load results
df = pd.read_csv(lm_file, sep="\t")

# Subset by term
sub = df[df["term"] == term].copy()
if sub.empty:
    sys.exit(f"No rows found for term {term}")

sub["-log10q"] = -np.log10(sub["qval"] + 1e-300)

# Plot
plt.figure(figsize=(7,5))
plt.scatter(sub["coef"], sub["-log10q"], s=10, alpha=0.6)

# Threshold line
plt.axhline(-np.log10(0.05), color="red", linestyle="--", label="FDR=0.05")

# Label top 15 hits
top_hits = sub.sort_values("qval").head(15)
for _, row in top_hits.iterrows():
    plt.text(row["coef"], row["-log10q"], row["regulon"],
             fontsize=7, ha="right" if row["coef"] < 0 else "left")

plt.xlabel(f"Coefficient ({term} effect)")
plt.ylabel("-log10(FDR q-value)")
plt.title(f"Volcano plot: {term} effect on regulon activity")
plt.legend()
plt.tight_layout()
plt.savefig(outfile, dpi=150)
print(f"[done] Saved volcano plot → {outfile}")
