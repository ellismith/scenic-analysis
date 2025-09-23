#!/scratch/easmit31/conda_envs/pyscenic/bin/python
"""
71_plot_regulon_vs_age.py

Plot regulon activity vs age for a single regulon.

Usage:
    ./71_plot_regulon_vs_age.py H5AD_PATH REGULON_NAME FIGURES_DIR
Example:
    ./71_plot_regulon_vs_age.py /path/to/data.h5ad TFAP2C figures_all
"""

import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

if len(sys.argv) != 4:
    sys.exit("Usage: ./71_plot_regulon_vs_age.py H5AD_PATH REGULON_NAME FIGURES_DIR")

h5ad_file = sys.argv[1]
regulon = sys.argv[2]
figures_dir = sys.argv[3]

# Output directory
outdir = Path(figures_dir)
outdir.mkdir(parents=True, exist_ok=True)
outfile = outdir / f"{regulon}_age_scatter.png"

print(f"Loading {h5ad_file}...")
adata = sc.read_h5ad(h5ad_file, backed='r')

# Check for age column and regulon column
if 'age' not in adata.obs.columns:
    sys.exit("Error: 'age' column not found in obs")

regulon_col = f"AUC_{regulon}_regulon"
if regulon_col not in adata.obs.columns:
    regulon_col = f"AUC_{regulon}"
    if regulon_col not in adata.obs.columns:
        sys.exit(f"Error: {regulon} not found in AUC columns")

print(f"Plotting {regulon_col} vs age...")

# Convert to memory for plotting
adata_mem = adata.to_memory()

# Get data
age = pd.to_numeric(adata_mem.obs['age'], errors='coerce')
activity = adata_mem.obs[regulon_col]

# Plot
df = pd.DataFrame({"age": age, "activity": activity})
plt.figure(figsize=(6,4))
sns.regplot(data=df, x="age", y="activity", scatter_kws={"s":10, "alpha":0.4})
plt.title(f"{regulon} regulon activity vs. Age")
plt.tight_layout()
plt.savefig(outfile, dpi=150)
print(f"[done] Saved plot → {outfile}")
