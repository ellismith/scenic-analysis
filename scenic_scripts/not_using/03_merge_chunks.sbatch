#!/bin/bash
#SBATCH --job-name=aucell_merge
#SBATCH --output=logs/03_merge_chunks_%j.out
#SBATCH --error=logs/03_merge_chunks_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail

PYBIN="/scratch/easmit31/conda_envs/pyscenic/bin/python"
TAG="${1:-manual}"

echo "[merge] starting (tag=$TAG) on $(hostname) at $(date)"
"$PYBIN" 03_merge_chunks.py --pattern "chunk_*_auc.csv" --out "aucell_full_auc_mtx.csv"
echo "[merge] finished at $(date)"
