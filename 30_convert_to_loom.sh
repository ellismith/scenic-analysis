#!/bin/bash
#SBATCH --job-name=loom
#SBATCH --time=03:59:00
#SBATCH --mem=256GB
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/30_convert_to_loom_%j.out
#SBATCH --partition=highmem

module load gcc-11.2.0-gcc-8.5.0
module load mamba/latest
PY=/scratch/easmit31/conda_envs/pyscenic_final/bin/python
PROJECT_DIR="/scratch/easmit31/GRN/scenic/scenic_scripts"
cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
mkdir -p logs
# Forward all sbatch args into Python
"$PY" 30_convert_to_loom.py "$@"
