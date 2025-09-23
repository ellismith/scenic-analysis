#!/usr/bin/env python3
"""
55_fix_varnames.py

Purpose
-------
Fix gene identifiers in AnnData .h5ad files so that AUCell can match regulons.
Specifically, replace Ensembl IDs (e.g. ENSMMUG00000000001) with human_gene_name
from .var if present.

Usage
-----
python 55_fix_varnames.py --input <input.h5ad> --output <output_fixed.h5ad>
"""

import argparse
import scanpy as sc

def main():
    ap = argparse.ArgumentParser(description="Fix var gene identifiers for AUCell compatibility")
    ap.add_argument("--input", required=True, help="Input .h5ad (with Ensembl IDs in .var)")
    ap.add_argument("--output", required=True, help="Output .h5ad with fixed gene names")
    args = ap.parse_args()

    print(f"[fix-varnames] Reading: {args.input}", flush=True)
    adata = sc.read_h5ad(args.input)

    if "human_gene_name" not in adata.var:
        raise KeyError("Expected var['human_gene_name'] column not found.")

    # Replace var_names
    adata.var_names = adata.var["human_gene_name"].astype(str)
    adata.var_names_make_unique()

    print(f"[fix-varnames] Writing: {args.output}", flush=True)
    adata.write(args.output)
    print("[fix-varnames] Done.", flush=True)

if __name__ == "__main__":
    main()
