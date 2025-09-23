#!/bin/bash
#SBATCH --job-name=export_auc
#SBATCH --output=logs/70_export_auc_%j.out
#SBATCH --error=logs/70_export_auc_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail
set -x

PY="/scratch/easmit31/conda_envs/pyscenic/bin/python"
cd /scratch/easmit31/GRN_copy/scenic/scenic_scripts

# Arguments
INPUT_H5AD=$1
OUTPUT_TSV=$2

$PY 70_export_auc_mtx.py --input "$INPUT_H5AD" --out "$OUTPUT_TSV"
