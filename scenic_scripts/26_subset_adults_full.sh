#!/bin/bash
#SBATCH --job-name=subset_adults_full
#SBATCH --output=logs/26_subset_adults_full_%j.out
#SBATCH --error=logs/26_subset_adults_full_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=8
#SBATCH --partition=highmem

set -euo pipefail
set -x

source /packages/apps/mamba/2.0.8/etc/profile.d/conda.sh
conda activate /scratch/easmit31/conda_envs/pyscenic

PYTHON=/scratch/easmit31/conda_envs/pyscenic/bin/python

$PYTHON 26_subset_adults_full.py "$@"
