#!/usr/bin/env bash
set -euo pipefail

# ====== EDIT THESE LINES ======
H5AD="/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/Res1_glutamatergic-neurons_update.h5ad"
REGULONS="regulons_for_aucell.csv"    # or .gmt
CHUNK_SIZE=100000                      # ↓ from 200k → 100k to reduce RAM per job
# =================================

PYBIN="${PYBIN:-/scratch/easmit31/conda_envs/pyscenic/bin/python}"

mkdir -p logs

# get total cells via h5py (no Scanpy)
N_OBS=$("$PYBIN" -c 'import h5py,sys
with h5py.File(sys.argv[1],"r") as f:
    try:
        n = f["obs"]["_index"].shape[0]
    except Exception:
        try: n = f["X"].shape[0]
        except Exception: raise RuntimeError("Could not determine n_obs from H5AD")
    print(n)' "$H5AD")

echo "[prep] n_obs=$N_OBS"
LAST=$(( (N_OBS + CHUNK_SIZE - 1) / CHUNK_SIZE - 1 ))
echo "[prep] chunk_size=$CHUNK_SIZE  -> array indices 0..$LAST"

# submit array; pass params via environment; logs go to logs/
jid=$(sbatch \
  --export=ALL,H5AD="$H5AD",REGULONS="$REGULONS",CHUNK_SIZE="$CHUNK_SIZE",N_OBS="$N_OBS" \
  --array=0-"$LAST" \
  --output=logs/aucell_chunk_%A_%a.out \
  --error=logs/aucell_chunk_%A_%a.err \
  02_run_aucell_chunk.sbatch | awk '{print $4}')

echo "[submit] job array id: $jid"
echo
echo "[watch] squeue -j $jid"
echo "[tail ] tail -f logs/aucell_chunk_${jid}_0000.err"
echo
echo "[merge] when array COMPLETES:"
echo "        sbatch --output=logs/03_merge_chunks_%j.out --error=logs/03_merge_chunks_%j.err 03_merge_chunks.sbatch $jid"
