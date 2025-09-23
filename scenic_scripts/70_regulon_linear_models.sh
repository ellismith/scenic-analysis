#!/bin/bash
#SBATCH --job-name=reg_lm
#SBATCH --output=logs/70_regulon_linear_models_%j.out
#SBATCH --error=logs/70_regulon_linear_models_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail
set -x

# Point to your Python in pyscenic env
PY=/scratch/easmit31/conda_envs/pyscenic/bin/python

# Inputs
IN_H5AD=$1
OUT_TSV=$2

# Run script
$PY 70_regulon_linear_models.py "$IN_H5AD" "$OUT_TSV"
