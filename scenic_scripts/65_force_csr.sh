#!/bin/bash
#SBATCH --job-name=force_csr
#SBATCH --output=logs/65_force_csr_%j.out
#SBATCH --error=logs/65_force_csr_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=8
#SBATCH --partition=highmem

set -euo pipefail
source ~/.bashrc
conda activate /scratch/easmit31/conda_envs/pyscenic

PY=/scratch/easmit31/conda_envs/pyscenic/bin/python
IN=/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv_liftover_sparse.h5ad
OUT=/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv_liftover_ready.h5ad

$PY - <<'PYCODE'
import scanpy as sc
from scipy.sparse import csr_matrix

infile = "'''$IN'''"
outfile = "'''$OUT'''"

print("[force_csr] Reading full AnnData")
adata = sc.read_h5ad(infile, backed=None)
print("[force_csr] Shape:", adata.shape)

adata.X = csr_matrix(adata.X.astype("float32"))
print("[force_csr] Converted X to CSR float32")

adata.write(outfile)
print("[force_csr] Wrote", outfile)
PYCODE
