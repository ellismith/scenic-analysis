#!/bin/bash
#SBATCH --job-name=fix_h5ad
#SBATCH --output=logs/65_fix_h5ad_%j.out
#SBATCH --error=logs/65_fix_h5ad_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=8
#SBATCH --partition=highmem

set -euo pipefail
source ~/.bashrc
conda activate /scratch/easmit31/conda_envs/pyscenic

PY=/scratch/easmit31/conda_envs/pyscenic/bin/python

IN=/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv_liftover_sparse.h5ad
OUT=/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv_liftover_ready.h5ad

$PY 65_fix_h5ad.py $IN $OUT
