#!/scratch/easmit31/conda_envs/pyscenic/bin/python
"""
70_regulon_linear_models.py

Run linear models of regulon activity ~ age + sex (from AUCell results).
- Each regulon’s AUCell score is modeled as a function of age and sex.
- Uses statsmodels OLS regression.
- Outputs coefficients, p-values, R², and FDR-adjusted q-values.

Usage:
    ./70_regulon_linear_models.py input.h5ad output.tsv
"""

import sys
import scanpy as sc
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

def main(in_h5ad, out_tsv):
    # 1. Load the h5ad file (contains AUCell results)
    adata = sc.read_h5ad(in_h5ad)

    # 2. Extract AUCell scores
    if "AUCell" in adata.obsm:
        print("[lm] Using AUCell scores from adata.obsm['AUCell']")
        regulons = pd.DataFrame(
            adata.obsm["AUCell"],
            index=adata.obs_names
        )
    else:
        auc_cols = [c for c in adata.obs.columns if c.startswith("AUC_")]
        if not auc_cols:
            raise ValueError("No AUCell scores found in obsm['AUCell'] or obs['AUC_*']")
        print(f"[lm] Using AUCell scores from {len(auc_cols)} obs columns")
        regulons = adata.obs[auc_cols].copy()
        regulons.index = adata.obs_names

    # 3. Get metadata (must have 'age' and 'sex')
    meta = adata.obs.copy()
    if not {"age", "sex"} <= set(meta.columns):
        raise ValueError("adata.obs must contain 'age' and 'sex' columns")

    results = []

    # 4. Fit linear model for each regulon
    for regulon in regulons.columns:
        df = pd.DataFrame({
            "activity": regulons[regulon],
            "age": meta["age"],
            "sex": meta["sex"],
        })
        df["sex"] = df["sex"].astype("category")

        try:
            model = smf.ols("activity ~ age + sex", data=df).fit()
            for term in model.params.index:
                results.append({
                    "regulon": regulon,
                    "term": term,
                    "coef": model.params[term],
                    "pval": model.pvalues[term],
                    "r2": model.rsquared
                })
        except Exception as e:
            print(f"[warn] {regulon} failed: {e}", file=sys.stderr)

    results_df = pd.DataFrame(results)

    # 5. Multiple testing correction (FDR, Benjamini–Hochberg)
    reject, qvals, _, _ = multipletests(results_df["pval"], method="fdr_bh")
    results_df["qval"] = qvals

    # 6. Save results
    results_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[done] Saved results to {out_tsv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: ./70_regulon_linear_models.py input.h5ad output.tsv")
    main(sys.argv[1], sys.argv[2])
