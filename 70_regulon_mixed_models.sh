#!/bin/bash
#SBATCH --job-name=mixed_models
#SBATCH --output=logs/70_mixed_models_%j.out
#SBATCH --error=logs/70_mixed_models_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

set -euo pipefail
set -x

module load gcc-11.2.0-gcc-8.5.0

PY=/scratch/easmit31/conda_envs/pyscenic_final/bin/python
IN_H5AD=$1
OUT_TSV=$2

$PY /scratch/easmit31/GRN/scenic/scenic_scripts/70_regulon_mixed_models.py "$IN_H5AD" "$OUT_TSV"
