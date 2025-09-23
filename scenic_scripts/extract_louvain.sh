#!/bin/bash
# run with: sbatch --mem=32G --cpus-per-task=4 extract_louvain.sh

PYBIN=/scratch/easmit31/conda_envs/pyscenic/bin/python

$PYBIN <<'PY'
import scanpy as sc

# load h5ad
adata = sc.read_h5ad("Res1_glutamatergic-neurons_update.h5ad")

# list metadata keys
print("Available obs columns:", list(adata.obs.columns))

# if 'louvain' exists, export
if "louvain" in adata.obs.columns:
    adata.obs[["louvain"]].to_csv("cell_clusters.csv")
    print("Wrote cell_clusters.csv with shape", adata.obs[["louvain"]].shape)
else:
    print("No 'louvain' column found. Available:", list(adata.obs.columns))
PY
