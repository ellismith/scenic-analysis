#!/bin/bash
#SBATCH --job-name=fix_env
#SBATCH --output=logs/fix_env_%j.out
#SBATCH --error=logs/fix_env_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -euo pipefail

mamba install -y -p /scratch/easmit31/conda_envs/pyscenic \
  "numpy=1.24.*" \
  "scikit-learn=1.3.*" \
  "pandas<2.1" \
  "dask=2023.3.2" \
  "distributed=2023.3.2" \
  "pyarrow=14.*"
