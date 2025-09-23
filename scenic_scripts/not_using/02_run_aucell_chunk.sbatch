#!/bin/bash
#SBATCH --job-name=aucell_chunk
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G

set -euo pipefail
export PYTHONNOUSERSITE=1

# >>> activate pyscenic environment explicitly <<<
export PATH="/scratch/easmit31/conda_envs/pyscenic/bin:$PATH"
PYBIN="/scratch/easmit31/conda_envs/pyscenic/bin/python"
PYS="/scratch/easmit31/conda_envs/pyscenic/bin/pyscenic"

# env passed from launcher
H5AD="${H5AD:?missing H5AD}"
REGULONS="${REGULONS:?missing REGULONS}"
CHUNK_SIZE="${CHUNK_SIZE:?missing CHUNK_SIZE}"
N_OBS="${N_OBS:?missing N_OBS}"

# compute start/end by array index
TASK=${SLURM_ARRAY_TASK_ID}
S=$(( TASK * CHUNK_SIZE ))
E=$(( S + CHUNK_SIZE ))
if (( S >= N_OBS )); then
  echo "[task $TASK] S=$S >= N_OBS=$N_OBS, nothing to do."
  exit 0
fi
if (( E > N_OBS )); then E=$N_OBS; fi

CID=$(printf "%04d" ${TASK})
LOOM="tmp_chunk_${CID}.loom"
AUC="chunk_${CID}_auc.csv"

echo "[task $CID] host=$(hostname) time=$(date)"
echo "[task $CID] slicing rows [$S:$E) from $H5AD  -> loom: $LOOM"

"$PYBIN" - <<PY
import scanpy as sc, loompy
h5ad = r"${H5AD}"
s = int(${S})
e = int(${E})
loom = r"${LOOM}"

# slice and write loom
ad = sc.read_h5ad(h5ad, backed="r")
sub = ad[s:e].to_memory()
sub.write_loom(loom)

# reopen with loompy to rename attributes
with loompy.connect(loom, 'r+') as ds:
    if 'var_names' in ds.ra:
        ds.ra['Gene'] = ds.ra['var_names']
        del ds.ra['var_names']
    if 'obs_names' in ds.ca:
        ds.ca['CellID'] = ds.ca['obs_names']
        del ds.ca['obs_names']
    print("[PY] wrote", loom, "with", sub.n_obs, "cells and fixed attrs")
PY

echo "[task $CID] AUCell scoring -> $AUC"
"$PYS" aucell "${LOOM}" "${REGULONS}" -o "${AUC}" --num_workers "${SLURM_CPUS_PER_TASK}"

echo "[task $CID] removing temp loom"
rm -f "${LOOM}"

echo "[task $CID] done $(date)"
