#!/scratch/easmit31/conda_envs/pyscenic_final/bin/python

import argparse
import datetime
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd   # for parsing ages safely

print(f"[{datetime.datetime.now()}] Starting 10_select_cells")

parser = argparse.ArgumentParser(description="Subset AnnData by metadata fields")
parser.add_argument("--h5ad_in",  required=True,
                    help="Input AnnData file (.h5ad)")
parser.add_argument("--h5ad_out", required=True,
                    help="Output subset AnnData file (.h5ad)")
parser.add_argument("--clusters", nargs="+", default=None,
                    help="Louvain clusters to keep (strings)")
parser.add_argument("--regions", nargs="+", default=None,
                    help="Regions to keep (strings)")
parser.add_argument("--min_individuals", type=int, default=None,
                    help="Minimum number of unique individuals required per cell type")
parser.add_argument("--age_min", type=float, default=None,
                    help="Minimum age to keep (numeric, e.g. 1.0)")

args = parser.parse_args()

print(f"[select] Reading: {args.h5ad_in}")
adata = sc.read_h5ad(args.h5ad_in)

mask = np.ones(adata.n_obs, dtype=bool)

# Cluster filter
if args.clusters is not None:
    if "louvain" not in adata.obs.columns:
        raise KeyError("obs['louvain'] not found in input AnnData.")
    clusters = [str(c) for c in args.clusters]
    mask &= adata.obs["louvain"].astype(str).isin(clusters)
    print(f"[select] Keeping clusters: {clusters}")

# Region filter
if args.regions is not None:
    if "region" not in adata.obs.columns:
        raise KeyError("obs['region'] not found in input AnnData.")
    mask &= adata.obs["region"].astype(str).isin(args.regions)
    print(f"[select] Keeping regions: {args.regions}")

# Age filter
if args.age_min is not None:
    if "age" not in adata.obs.columns:
        raise KeyError("obs['age'] not found in input AnnData.")
    ages = pd.to_numeric(adata.obs["age"], errors="coerce")
    mask &= ages >= args.age_min
    print(f"[select] Keeping cells with age >= {args.age_min}")

adata_subset = adata[mask].copy()

# Debugging: report before/after
print(f"[select] Original: {adata.shape} → Subset: {adata_subset.shape}")
if "age" in adata_subset.obs.columns:
    ages = pd.to_numeric(adata_subset.obs["age"], errors="coerce")
    print(f"[select] Subset age range: min={ages.min()}, max={ages.max()}")

# Minimum individuals filter
if args.min_individuals is not None:
    if "individual" not in adata_subset.obs.columns:
        raise KeyError("obs['individual'] not found in input AnnData.")
    counts = adata_subset.obs.groupby("cell_type")["individual"].nunique()
    print("[select] Unique individuals per cell_type after filters:")
    print(counts.to_string())

    valid_types = counts[counts >= args.min_individuals].index
    adata_subset = adata_subset[adata_subset.obs["cell_type"].isin(valid_types)].copy()

    if adata_subset.n_obs == 0:
        print(f"[WARNING] No cells passed the min_individuals filter (≥{args.min_individuals}).")
        print("[WARNING] Subset not written. Exiting without creating output file.")
        sys.exit(0)

    print(f"[select] Keeping cell types with ≥{args.min_individuals} individuals: {list(valid_types)}")

adata_subset.write_h5ad(args.h5ad_out)
print(f"[select] Wrote subset to: {args.h5ad_out}")
