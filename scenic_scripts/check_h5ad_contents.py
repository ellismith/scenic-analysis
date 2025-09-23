#!/usr/bin/env python3
import scanpy as sc
import sys

adata = sc.read_h5ad(sys.argv[1])
print(f"Shape: {adata.shape}")
print(f"obs keys: {list(adata.obs.columns)}")
print(f"obsm keys: {list(adata.obsm.keys())}")
print(f"var keys: {list(adata.var.columns)}")

# Check for AUCell specifically
if "AUCell" in adata.obsm:
    print(f"AUCell shape: {adata.obsm['AUCell'].shape}")
    print(f"AUCell type: {type(adata.obsm['AUCell'])}")
else:
    print("No 'AUCell' found in obsm")
    
# Look for regulon-related data
regulon_obs = [col for col in adata.obs.columns if 'regulon' in col.lower()]
regulon_var = [col for col in adata.var.columns if 'regulon' in col.lower()]
print(f"Regulon columns in obs: {regulon_obs}")
print(f"Regulon columns in var: {regulon_var}")
