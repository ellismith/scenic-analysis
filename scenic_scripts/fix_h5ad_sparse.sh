#!/bin/bash
#SBATCH --job-name=fix_h5ad
#SBATCH --output=logs/fix_h5ad_%j.out
#SBATCH --error=logs/fix_h5ad_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch fix_h5ad_sparse.sh <in.h5ad> <out.h5ad>"
  exit 1
fi

PY=/scratch/easmit31/conda_envs/pyscenic/bin/python
INFILE=$1
OUTFILE=$2

$PY fix_h5ad_sparse.py "$INFILE" "$OUTFILE"
