#!/bin/bash
#SBATCH --job-name=analysis70
#SBATCH --output=logs/70_analysis_%j.out
#SBATCH --error=logs/70_analysis_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
set -euo pipefail
set -x

module load mamba/latest
source activate /scratch/easmit31/conda_envs/pyscenic

H5AD_IN="$1"
PREFIX="$2"

/scratch/easmit31/conda_envs/pyscenic/bin/python \
  70_analysis.py \
  --h5ad_in "$H5AD_IN" \
  --prefix "$PREFIX" \
  --groupby louvain \
  --out_dir /scratch/easmit31/GRN_copy/scenic
