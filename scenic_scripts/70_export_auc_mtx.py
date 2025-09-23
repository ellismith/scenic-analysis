#!/scratch/easmit31/conda_envs/pyscenic/bin/python
"""
70_export_auc_mtx.py

Export AUCell regulon activity + metadata from .h5ad into a .tsv.

Usage:
    ./70_export_auc_mtx.py --input INPUT_H5AD --out OUTPUT_TSV
"""

import argparse
import scanpy as sc
import pandas as pd
from pathlib import Path

def main(input_h5ad, output_tsv):
    print(f"[export] Reading: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)

    # --- regulon columns ---
    auc_cols = [c for c in adata.obs.columns if c.startswith("AUC_")]
    if not auc_cols:
        raise ValueError("No AUC_* regulon columns found in adata.obs")
    print(f"[export] Found {len(auc_cols)} regulon columns")

    df_auc = adata.obs[auc_cols].copy()

    # --- metadata columns ---
    meta_cols = [c for c in ["louvain", "age", "sex"] if c in adata.obs.columns]
    if meta_cols:
        print(f"[export] Adding metadata columns: {meta_cols}")
        df_meta = adata.obs[meta_cols].copy()
    else:
        print("[export] No metadata columns found")
        df_meta = pd.DataFrame(index=adata.obs_names)

    # --- combine ---
    df_combined = pd.concat([df_meta, df_auc], axis=1)
    df_combined.insert(0, "cell", adata.obs_names)

    # --- write ---
    out_path = Path(output_tsv)
    df_combined.to_csv(out_path, sep="\t", index=False)
    print(f"[export] Wrote {out_path} with shape {df_combined.shape}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export AUCell matrix + metadata to TSV")
    parser.add_argument("--input", required=True, help="Input .h5ad with AUCell scores in obs")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()
    main(args.input, args.out)
