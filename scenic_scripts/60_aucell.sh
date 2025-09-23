#!/bin/bash
#SBATCH --job-name=aucell
#SBATCH --output=logs/60_aucell_%j.out
#SBATCH --error=logs/60_aucell_%j.err
#SBATCH --time=3:59:00
#SBATCH --mem=216GB
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc

set -euo pipefail

PY=/scratch/easmit31/conda_envs/pyscenic/bin/python
PROJECT_DIR=/scratch/easmit31/GRN_copy/scenic/scenic_scripts
cd "$PROJECT_DIR"

$PY 60_aucell.py \
  "$@"
