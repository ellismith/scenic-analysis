#!/bin/bash
#SBATCH --job-name=merge_aucell
#SBATCH --output=logs/merge_aucell_%j.out
#SBATCH --error=logs/merge_aucell_%j.err
#SBATCH --time=03:59:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc

set -euo pipefail

module load gcc-11.2.0-gcc-8.5.0

PY=/scratch/easmit31/conda_envs/pyscenic_final/bin/python
PROJECT_DIR=/scratch/easmit31/GRN_copy/scenic/scenic_scripts

cd "$PROJECT_DIR"

if [ "$#" -lt 2 ]; then
  echo "Usage: sbatch merge_aucell_chunks.sh <chunks_dir> <output_file>"
  exit 1
fi

CHUNKS_DIR=$1
OUTPUT_FILE=$2

echo "Merging AUCell chunks from: $CHUNKS_DIR"
echo "Output will be written to: $OUTPUT_FILE"

$PY merge_aucell_chunks.py "$CHUNKS_DIR" "$OUTPUT_FILE"
