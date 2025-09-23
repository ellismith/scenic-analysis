#!/bin/bash
#SBATCH --job-name=fix_varnames
#SBATCH --output=logs/55_fix_varnames_%j.out
#SBATCH --error=logs/55_fix_varnames_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -eo pipefail

# Load conda
source ~/.bashrc
conda activate /scratch/easmit31/conda_envs/pyscenic

INPUT="$1"
OUTPUT="$2"

echo "[fix-varnames.sh] Input:  $INPUT"
echo "[fix-varnames.sh] Output: $OUTPUT"

python 55_fix_varnames.py --input "$INPUT" --output "$OUTPUT"
