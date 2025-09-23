#!/bin/bash
#SBATCH --job-name=submit_clusters
#SBATCH --output=logs/submit_clusters_%j.out
#SBATCH --error=logs/submit_clusters_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
# Usage: sbatch 96_submit_clusters.sh <prefix> <suffix>
# Example: sbatch 96_submit_clusters.sh gaba _adults

if [ -z "$1" ] || [ -z "$2" ]; then
  echo "❌ Error: Must provide prefix and suffix (e.g. gaba _adults)"
  exit 1
fi

CELLTYPE="$1"
SUFFIX="$2"
INPUT="/scratch/easmit31/GRN_copy/scenic/h5ad_files/${CELLTYPE}_allLv_pyscenic_output.h5ad"
GROUPBY="louvain"
MAX_REGULONS=40
SAMPLE=1000

# detect clusters
NCLUSTERS=$(python3 - <<PY
import scanpy as sc
adata = sc.read_h5ad("${INPUT}", backed="r")
print(len(adata.obs["${GROUPBY}"].astype(str).unique()))
PY
)

echo "[submit-clusters] Prefix: ${CELLTYPE}, Suffix: ${SUFFIX}"
echo "[submit-clusters] Input: ${INPUT}, Clusters: ${NCLUSTERS}"

for cl in $(seq 0 $((NCLUSTERS-1))); do
  sbatch --job-name="${CELLTYPE}_cl${cl}" \
         --output=logs/${CELLTYPE}${SUFFIX}_cluster${cl}_heatmap_%j.out \
         --error=logs/${CELLTYPE}${SUFFIX}_cluster${cl}_heatmap_%j.err \
         --mem=16G --time=01:00:00 --cpus-per-task=2 \
         --wrap "
  source /packages/apps/mamba/2.0.8/etc/profile.d/conda.sh
  conda activate /scratch/easmit31/conda_envs/pyscenic
  python3 96_cluster_regulon_heatmap.py \
    --input ${INPUT} \
    --celltype ${CELLTYPE} \
    --groupby ${GROUPBY} \
    --cluster ${cl} \
    --max_regulons ${MAX_REGULONS} \
    --sample ${SAMPLE} \
    --suffix ${SUFFIX}
  "
done
