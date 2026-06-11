#!/bin/bash
#SBATCH --job-name=aucell_final
#SBATCH --output=logs/60_aucell_final_%j.out
#SBATCH --error=logs/60_aucell_final_%j.err
#SBATCH --time=3:59:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
module load gcc-11.2.0-gcc-8.5.0
set -euo pipefail
PY=/scratch/easmit31/conda_envs/pyscenic_final/bin/python
PROJECT_DIR=/scratch/easmit31/GRN/scenic/scenic_scripts
cd "$PROJECT_DIR"
$PY 60_aucell_final.py "$@"
