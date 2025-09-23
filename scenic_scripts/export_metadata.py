#!/usr/bin/env python3
import scanpy as sc
import sys

# Read your h5ad file
adata = sc.read_h5ad(sys.argv[1])

# Export obs (metadata) to CSV
adata.obs.to_csv(sys.argv[2])
print(f"Exported metadata to {sys.argv[2]}")
print(f"Available columns: {list(adata.obs.columns)}")
