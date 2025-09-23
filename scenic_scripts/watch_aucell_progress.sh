#!/usr/bin/env bash
set -euo pipefail

# Usage: ./watch_aucell_progress.sh <loom> <auc_csv> [interval_seconds]
LOOM="${1:-aucell_full.loom}"
OUTCSV="${2:-aucell_full_auc_mtx.csv}"
INTERVAL="${3:-30}"

PYBIN="${PYBIN:-python}"

# total cells = number of columns in the loom's /matrix
tot=$("$PYBIN" -c 'import h5py,sys; 
with h5py.File(sys.argv[1],"r") as f: 
    print(f["/matrix"].shape[1])' "$LOOM")

prev_done=0
prev_ts=$(date +%s)

while true; do
  now_ts=$(date +%s)
  # CSV grows by one header + one line per cell
  if [[ -f "$OUTCSV" ]]; then
    lines=$(wc -l < "$OUTCSV")
  else
    lines=0
  fi
  done=$(( lines>0 ? lines-1 : 0 ))

  pct=$("$PYBIN" -c 'import sys; tot=int(sys.argv[1]); done=int(sys.argv[2]); 
print(f"{(done/tot*100):.2f}" if tot>0 else "0.00")' "$tot" "$done")

  dt=$(( now_ts - prev_ts ))
  dcell=$(( done - prev_done ))
  rate=$("$PYBIN" -c 'import sys; dt=float(sys.argv[1]); d=float(sys.argv[2]);
print(f"{(d/dt):.2f}" if dt>0 else "0.00")' "$dt" "$dcell")

  rem=$(( tot - done ))
  eta=$("$PYBIN" - <<'PY' "$rate" "$rem"
import sys, math
rate=float(sys.argv[1]); rem=int(sys.argv[2])
if rate<=0:
    print("NA")
else:
    sec = rem / rate
    h = int(sec // 3600); m = int((sec % 3600) // 60)
    print(f"{h}h{m:02d}m")
PY
)

  ts=$(date "+%Y-%m-%d %H:%M:%S")
  echo "[$ts] $done/$tot cells  (${pct}%)  rate=${rate} cells/s  ETA=${eta}  file=${OUTCSV}"

  prev_done=$done
  prev_ts=$now_ts
  sleep "$INTERVAL"
done
