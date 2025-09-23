#!/bin/bash
#SBATCH --job-name=heatmap
#SBATCH --output=logs/95_heatmap_%j.out
#SBATCH --error=logs/95_heatmap_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
set -euo pipefail
set -x

# Load conda
source ~/.bashrc
conda activate /scratch/easmit31/conda_envs/pyscenic

INPUT="$1"
OUT_PNG="$2"
GROUPBY="louvain"
MAX_REGULONS=50

echo "[heatmap.sh] Input file:    $INPUT"
echo "[heatmap.sh] Output PNG:   $OUT_PNG"
echo "[heatmap.sh] Grouping by:  $GROUPBY"
echo "[heatmap.sh] Max regulons: $MAX_REGULONS"

python3 95_heatmap.py \
  --zscore "$INPUT" \
  --out_png "$OUT_PNG" \
  --max_rows "$MAX_REGULONS"
