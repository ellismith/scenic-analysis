#!/bin/bash
#SBATCH --job-name=cluster_models
#SBATCH --output=logs/74_cluster_models_%j.out
#SBATCH --error=logs/74_cluster_models_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

set -euo pipefail
set -x

module load gcc-11.2.0-gcc-8.5.0

PY=/scratch/easmit31/conda_envs/pyscenic_final/bin/python
IN_H5AD=$1
SOURCE_H5AD=$2
OUT_TSV=$3
$PY /scratch/easmit31/GRN/scenic/scenic_scripts/74_cluster_mixed_models.py "$IN_H5AD" "$SOURCE_H5AD" "$OUT_TSV"
