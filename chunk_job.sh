#!/bin/bash
#SBATCH --job-name=chunk_h5ad
#SBATCH --output=logs/chunk_%j.out
#SBATCH --error=logs/chunk_%j.err
#SBATCH --time=3:59:00
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc

module load gcc-11.2.0-gcc-8.5.0
set -euo pipefail

PY=/scratch/easmit31/conda_envs/pyscenic_final/bin/python
PROJECT_DIR=/scratch/easmit31/GRN/scenic/scenic_scripts

cd "$PROJECT_DIR"

echo "Starting chunking job..."
$PY chunk_h5ad.py "$@"
echo "Chunking complete!"
