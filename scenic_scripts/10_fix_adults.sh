#!/bin/bash
#SBATCH --job-name=fix_adults
#SBATCH --partition=highmem
#SBATCH --time=04:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/10_fix_adults_%j.out
#SBATCH --error=logs/10_fix_adults_%j.err

set -euo pipefail

PREFIX=$1
IN_FILE="/scratch/easmit31/GRN_copy/scenic/h5ad_files/${PREFIX}_adults_allLv.h5ad"
OUT_FILE="/scratch/easmit31/GRN_copy/scenic/h5ad_files/${PREFIX}_adults_allLv_fixed.h5ad"

echo "[fix] Input:  $IN_FILE"
echo "[fix] Output: $OUT_FILE"

/scratch/easmit31/conda_envs/pyscenic/bin/python - <<'PY'
import sys, scipy.sparse as sp, anndata as ad
prefix = sys.argv[1]
in_file  = f"/scratch/easmit31/GRN_copy/scenic/h5ad_files/{prefix}_adults_allLv.h5ad"
out_file = f"/scratch/easmit31/GRN_copy/scenic/h5ad_files/{prefix}_adults_allLv_fixed.h5ad"

adata = ad.read_h5ad(in_file)   # full load into memory
print("[fix] Loaded:", adata.shape, "Matrix type:", type(adata.X))

# ensure .X is CSR
if not sp.isspmatrix_csr(adata.X):
    adata.X = sp.csr_matrix(adata.X)

print("[fix] Converted matrix type:", type(adata.X))
adata.write(out_file)
print("[fix] Wrote:", out_file)

# reopen to confirm
adata2 = ad.read_h5ad(out_file, backed="r")
print("[fix] Re-loaded:", adata2.shape, "Matrix type:", type(adata2.X))
adata2.file.close()
PY $PREFIX
