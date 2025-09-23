#!/scratch/easmit31/conda_envs/pyscenic/bin/python
"""
73_age_sex_scatter.py

Make a 4-way scatter plot: age vs sex coefficients.
- Regulon labels cleaned (drop "AUC_" prefix and "_regulon" suffix).

Usage:
    ./73_age_sex_scatter.py LM_TSV_PATH FIGURES_DIR
Example:
    ./73_age_sex_scatter.py /path/to/results.tsv figures_all
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

if len(sys.argv) != 3:
    sys.exit("Usage: ./73_age_sex_scatter.py LM_TSV_PATH FIGURES_DIR")

lm_file = sys.argv[1]
figures_dir = sys.argv[2]

# Derive prefix from filename (strip suffix like _lm.tsv)
lm_path = Path(lm_file)
prefix = lm_path.stem.replace("_lm", "")

# Output directory
outdir = Path(figures_dir)
outdir.mkdir(parents=True, exist_ok=True)
outfile = outdir / f"{prefix}_age_sex_scatter.png"

# Load results
df = pd.read_csv(lm_file, sep="\t")

# Clean regulon names
df["regulon_clean"] = df["regulon"].str.replace("^AUC_", "", regex=True).str.replace("_regulon$", "", regex=True)

# Pivot to get both age and sex coefs per regulon
age_df = df[df["term"] == "age"][["regulon_clean", "coef", "qval"]].rename(columns={"coef":"age_coef","qval":"age_q"})
sex_df = df[df["term"] == "sex[T.M]"][["regulon_clean", "coef", "qval"]].rename(columns={"coef":"sex_coef","qval":"sex_q"})
merged = pd.merge(age_df, sex_df, on="regulon_clean", how="inner")

# Plot
plt.figure(figsize=(7,7))
plt.scatter(merged["age_coef"], merged["sex_coef"], s=12, alpha=0.6)

# Axes lines
plt.axhline(0, color="black", linestyle="--")
plt.axvline(0, color="black", linestyle="--")

# Labels
plt.xlabel("Age coefficient")
plt.ylabel("Sex coefficient (M vs F)")
plt.title("Regulon age vs sex effects")

# Label top hits
merged["combined_q"] = merged[["age_q","sex_q"]].min(axis=1)
top_hits = merged.sort_values("combined_q").head(15)
for _, row in top_hits.iterrows():
    plt.text(row["age_coef"], row["sex_coef"], row["regulon_clean"],
             fontsize=7, ha="right" if row["age_coef"] < 0 else "left")

plt.tight_layout()
plt.savefig(outfile, dpi=150)
print(f"[done] Saved scatter plot → {outfile}")
