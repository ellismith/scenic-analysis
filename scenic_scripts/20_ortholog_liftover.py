#!/usr/bin/env python3
"""
20_ortholog_liftover.py

Purpose
-------
Map macaque genes to one-to-one human orthologs and collapse counts to
human gene symbols, producing a liftover AnnData (.h5ad) for downstream SCENIC.

Inputs
------
- SUBSET_H5AD  (from config.py): cells filtered by Louvain clusters
- MART_CSV     (from config.py): Ensembl mapping table with columns:
    ['Gene stable ID', 'Gene name', 'Human gene name', 'Human homology type', ...]
    Only rows with Human homology type == 'ortholog_one2one' are used.

Outputs
-------
- LIFTOVER_H5AD: AnnData with cells x human genes.
  adata.var has 'n_mm_geneid' = number of macaque Ensembl IDs collapsed per human gene.

CLI overrides
-------------
You can override config paths/params:
    python 20_ortholog_liftover.py \
        --h5ad_in /path/to/subset.h5ad \
        --mart_csv /path/to/ensembl113_mmul10_macaque_human.csv \
        --h5ad_out /path/to/liftover.h5ad \
        --min_cells 100 \
        --downsample 100000

Notes
-----
- Requires: anndata, numpy, pandas, scanpy (only for read, not required strictly),
  but we use anndata directly to avoid heavy deps.
- Memory usage depends on matrix size; downsample if needed.
"""

import argparse
import numpy as np
import pandas as pd
import anndata as ad

# Pull defaults from your shared config
from config import SUBSET_H5AD, LIFTOVER_H5AD, MART_CSV, NCELLS_FILTER

def main():
    p = argparse.ArgumentParser(description="Macaque→Human ortholog liftover and gene collapse")
    p.add_argument("--h5ad_in",  default=str(SUBSET_H5AD), help="Input subset .h5ad (cells x macaque genes)")
    p.add_argument("--mart_csv", default=str(MART_CSV),    help="Ensembl macaque-human mapping CSV")
    p.add_argument("--h5ad_out", default=str(LIFTOVER_H5AD), help="Output liftover .h5ad (cells x human genes)")
    p.add_argument("--min_cells", type=int, default=NCELLS_FILTER, help="Keep genes expressed in at least this many cells")
    p.add_argument("--downsample", type=int, default=100000,
                   help="If >0 and n_cells > this, randomly downsample cells to this number")
    args = p.parse_args()

    print("=== Liftover configuration ===")
    print("Input subset h5ad:", args.h5ad_in)
    print("Mapping CSV       :", args.mart_csv)
    print("Output h5ad       :", args.h5ad_out)
    print("Min cells per gene:", args.min_cells)
    print("Downsample cells  :", args.downsample)

    # --- Load subset data
    adata = ad.read_h5ad(args.h5ad_in)
    print(f"Loaded subset: {adata.shape} (cells x genes)")

    # --- Load mapping and restrict to 1:1 orthologs present in data
    mart = pd.read_csv(args.mart_csv)
    before = len(mart)
    mart = mart.loc[mart["Human homology type"] == "ortholog_one2one"].dropna(subset=["Gene name","Human gene name"])
    print(f"Mapping rows: {before} -> {len(mart)} after 1:1 + dropna")

    # The input AnnData must have Ensembl IDs in adata.var['ensembl_gene_id']
    if "ensembl_gene_id" not in adata.var.columns:
        raise KeyError("Expected 'ensembl_gene_id' in adata.var. Please ensure your input has this column.")

    # Optional prefilter to only mapped Ensembl IDs
    mapped_ids = set(mart["Gene stable ID"])
    print("Intersecting genes with mapping table...")
    adata.var["nCell"] = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    keep_vars = (adata.var["nCell"] >= args.min_cells) & (adata.var["ensembl_gene_id"].isin(mapped_ids))
    kept = int(keep_vars.sum())
    print(f"Keeping {kept} genes after min_cells & mapping filter.")
    adata = adata[:, keep_vars].copy()

    # --- Optional downsample cells
    if args.downsample and adata.n_obs > args.downsample:
        print(f"Downsampling cells: {adata.n_obs} -> {args.downsample}")
        rng = np.random.default_rng(42)
        keep_idx = rng.choice(adata.obs_names, size=args.downsample, replace=False)
        adata = adata[keep_idx].copy()

    # --- Build macaque Ensembl ID -> human gene name map
    ens2human = dict(zip(mart["Gene stable ID"], mart["Human gene name"]))
    adata.var["human_gene_name"] = adata.var["ensembl_gene_id"].map(ens2human)

    # Safety check
    if adata.var["human_gene_name"].isna().any():
        missing = int(adata.var["human_gene_name"].isna().sum())
        print(f"Warning: {missing} genes in subset lacked a human mapping; they will be dropped in collapse.")

    # --- Collapse counts to human gene names
    # Strategy: X (cells x macaque-ens) dot one-hot(human_gene_name) => (cells x human genes)
    print("Collapsing counts to human gene names...")
    one_hot = pd.get_dummies(adata.var["human_gene_name"]).values.astype(int)  # (mac_ens x unique human genes)
    X = adata.X
    # result is (cells x human genes)
    pb = X.dot(one_hot)

    # Build var with n_mm_geneid (how many macaque IDs collapsed into each human gene)
    counts_per_human = adata.var.groupby("human_gene_name").size().rename("n_mm_geneid")
    pb_var = counts_per_human.reset_index().set_index("human_gene_name")

    pb_adata = ad.AnnData(pb, obs=adata.obs.copy(), var=pb_var.copy())
    print(f"Liftover output: {pb_adata.shape} (cells x human genes)")

    # --- Write output
    pb_adata.write_h5ad(args.h5ad_out)
    print(f"Saved liftover AnnData: {args.h5ad_out}")

if __name__ == "__main__":
    main()
