#!/bin/bash
#SBATCH --job-name=liftover
#SBATCH --partition=highmem
#SBATCH --time=10:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/20_ortholog_liftover_%j.out

module load gcc-11.2.0-gcc-8.5.0
module load mamba/latest
PROJECT_DIR="/scratch/easmit31/GRN/scenic/scenic_scripts"
PYTHON="/scratch/easmit31/conda_envs/pyscenic_final/bin/python"
cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"
# Usage: sbatch 20_ortholog_liftover.sh <input.h5ad> <mapping.csv> <output.h5ad>
"$PYTHON" 20_ortholog_liftover_sparse.py "$@"
