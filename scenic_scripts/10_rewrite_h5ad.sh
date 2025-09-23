#!/bin/bash
#SBATCH --job-name=rewrite_h5ad
#SBATCH --output=logs/10_rewrite_h5ad_%j.out
#SBATCH --error=logs/10_rewrite_h5ad_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=highmem
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4

set -euo pipefail

INPUT=$1   # e.g. gaba_adults_allLv_fixed.h5ad
OUTPUT=$2  # e.g. gaba_adults_allLv_ready.h5ad

echo "[rewrite_h5ad] Input:  $INPUT"
echo "[rewrite_h5ad] Output: $OUTPUT"

source /packages/apps/mamba/2.0.8/etc/profile.d/conda.sh
conda activate /scratch/easmit31/conda_envs/pyscenic

python3 - <<PY
import scanpy as sc
print("[rewrite_h5ad] Reading:", "$INPUT")
adata = sc.read_h5ad("$INPUT")  # full load into memory
if hasattr(adata.X, "toarray"):
    adata.X = adata.X.toarray()
print("[rewrite_h5ad] Writing:", "$OUTPUT")
adata.write("$OUTPUT")
print("[rewrite_h5ad] Done.")
PY
