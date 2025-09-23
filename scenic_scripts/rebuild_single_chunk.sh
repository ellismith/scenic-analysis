#!/bin/bash
#SBATCH --job-name=rebuild_chunk
#SBATCH --output=logs/rebuild_%j.out
#SBATCH --error=logs/rebuild_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=htc

set -euo pipefail

PY=/scratch/easmit31/conda_envs/pyscenic/bin/python
PROJECT_DIR=/scratch/easmit31/GRN_copy/scenic/scenic_scripts

cd "$PROJECT_DIR"

echo "Rebuilding chunk: $1"
$PY rebuild_all_chunks.py "$1"
