#!/bin/bash
#SBATCH --job-name=test_sparse
#SBATCH --partition=public
#SBATCH --time=00:10:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/test_sparse_%j.out
#SBATCH --error=logs/test_sparse_%j.err

module load mamba/latest
PYTHON="/scratch/easmit31/conda_envs/pyscenic/bin/python"

cd /scratch/easmit31/GRN_copy/scenic/scenic_scripts || exit 1

# Step 1: make a toy h5ad + mapping
$PYTHON <<'PY'
import scanpy as sc
import pandas as pd
import scipy.sparse as sp

# Load a slice of your real data (small)
adata = sc.read_h5ad("/scratch/easmit31/GRN_copy/scenic/h5ad_files/gaba_adults_allLv_ready.h5ad")
adata = adata[:1000, :100].copy()  # just 1000 cells × 100 genes
adata.write("test_small.h5ad")

# Make fake mapping for those 100 genes
map_df = pd.DataFrame({
    "Gene stable ID": adata.var_names,
    "Human gene name": ["HG"+str(i) for i in range(len(adata.var_names))],
    "Human homology type": ["ortholog_one2one"]*len(adata.var_names),
})
map_df.to_csv("test_map.csv", index=False)
print("Wrote test_small.h5ad and test_map.csv")
PY

# Step 2: run the patched sparse liftover
$PYTHON 20_ortholog_liftover_sparse.py test_small.h5ad test_map.csv test_out.h5ad

# Step 3: sanity-check the output
$PYTHON <<'PY'
import scanpy as sc
adata = sc.read_h5ad("test_out.h5ad")
print("Output shape:", adata.shape)
print("First few genes:", list(adata.var.index)[:5])
PY
