#!/bin/bash
#SBATCH --job-name=chunk_h5ad
#SBATCH --output=logs/chunk_%j.out
#SBATCH --error=logs/chunk_%j.err
#SBATCH --time=3:59:00
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc

set -euo pipefail

PY=/scratch/easmit31/conda_envs/pyscenic/bin/python
PROJECT_DIR=/scratch/easmit31/GRN_copy/scenic/scenic_scripts

cd "$PROJECT_DIR"

echo "Starting chunking job..."
$PY chunk_h5ad.py "$@"
echo "Chunking complete!"
