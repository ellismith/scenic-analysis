#!/usr/bin/env python3
"""
90_gam_models.py

Purpose
-------
Fit univariate GAMs (pyGAM) for each regulon vs. a numeric metadata column (default: 'age').
Compute R² values for each model as a measure of fit.

Inputs
------
- PYSCENIC_OUT (AnnData .h5ad) from config.py
  Expects:
    - obs[meta_col] numeric (default: 'age')
    - regulon columns in obs (Regulon* from AUCELL)

Outputs
-------
- <WORK_DIR>/<CT>_gam_r2.tsv.gz  (columns: regulon, r2, n)

Usage
-----
python 90_gam_models.py
# or override:
python 90_gam_models.py --meta_col age
"""
import re, argparse
import numpy as np
import pandas as pd
import scanpy as sc
from pygam import LinearGAM, s
from sklearn.metrics import r2_score
from config import CT, WORK_DIR, PYSCENIC_OUT

def main():
    ap = argparse.ArgumentParser(description="Fit GAMs per regulon vs metadata")
    ap.add_argument("--h5ad_in", default=str(PYSCENIC_OUT), help="Input AUCELL .h5ad")
    ap.add_argument("--meta_col", default="age", help="Metadata column in obs to regress against")
    ap.add_argument("--out", default=str(WORK_DIR / f"{CT}_gam_r2.tsv.gz"),
                    help="Output TSV (gz)")
    args = ap.parse_args()

    print(f"[gam] Reading: {args.h5ad_in}")
    adata = sc.read_h5ad(args.h5ad_in)
    if args.meta_col not in adata.obs.columns:
        raise KeyError(f"obs['{args.meta_col}'] not found in AUCELL output.")

    obs = adata.obs.copy()
    y = pd.to_numeric(obs[args.meta_col], errors="coerce")
    valid = y.notna()

    regulon_cols = [c for c in obs.columns if re.match(r'(?i)^regulon', c)]
    if not regulon_cols:
        raise RuntimeError("No 'Regulon*' columns found in adata.obs")

    results = []
    for reg in regulon_cols:
        x = pd.to_numeric(obs.loc[valid, reg], errors="coerce")
        mask = x.notna()
        if mask.sum() < 10:
            results.append({"regulon": reg, "r2": np.nan, "n": int(mask.sum())})
            continue

        X = x[mask].values.reshape(-1, 1)
        y_reg = y[mask].values

        try:
            gam = LinearGAM(s(0)).fit(X, y_reg)
            y_pred = gam.predict(X)
            r2 = r2_score(y_reg, y_pred)
            results.append({"regulon": reg, "r2": float(r2), "n": int(len(y_reg))})
        except Exception as e:
            print(f"[gam] Error with {reg}: {e}")
            results.append({"regulon": reg, "r2": np.nan, "n": int(len(y_reg))})

    df = pd.DataFrame(results).sort_values("r2", ascending=False)
    df.to_csv(args.out, sep="\t", index=False, compression="gzip")
    print(f"[gam] Wrote: {args.out} (rows={len(df)})")

if __name__ == "__main__":
    main()
