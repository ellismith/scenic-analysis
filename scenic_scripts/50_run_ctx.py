#!/usr/bin/env python3
"""
50_run_ctx.py

Run pySCENIC ctx to prune the adjacency matrix with cisTarget databases.
This script looks for all *.rankings.feather files (ignores *.scores.feather).
"""

import argparse
import glob
import sys
import subprocess
import datetime

print(f"[{datetime.datetime.now()}] Starting 50_run_ctx")

parser = argparse.ArgumentParser(description="Run pySCENIC ctx with ranking DBs")
parser.add_argument("--adj", required=True, help="Adjacencies CSV from grn step")
parser.add_argument("--expression_mtx_fname", required=True, help="Filtered loom file (.loom)")
parser.add_argument("--db_glob", required=True, help="Glob to cisTarget DB files")
parser.add_argument("--annotations_fname", required=True, help="Motif-to-TF annotations table (.tbl)")
parser.add_argument("--output", required=True, help="Output .csv file with regulons")
parser.add_argument("--workers", type=int, default=20, help="Number of workers")
args = parser.parse_args()

# Expand DB files and filter only rankings
db_files = sorted(glob.glob(args.db_glob))
rank_db_files = [p for p in db_files if p.endswith(".rankings.feather")]

if not rank_db_files:
    print(f"[ctx] ERROR: No rankings DBs found for pattern:\n  {args.db_glob}", file=sys.stderr)
    sys.exit(1)

print(f"[ctx] Using {len(rank_db_files)} ranking DBs")
for p in rank_db_files[:5]:
    print(f"  - {p}")
if len(rank_db_files) > 5:
    print("  ...")

cmd = [
    "pyscenic", "ctx", args.adj,
    *rank_db_files,
    "--annotations_fname", args.annotations_fname,
    "--expression_mtx_fname", args.expression_mtx_fname,
    "--output", args.output,
    "--mask_dropouts",
    "--num_workers", str(args.workers),
]

print("[ctx] Running command:")
print("  " + " ".join(cmd))

res = subprocess.run(cmd)
if res.returncode != 0:
    print(f"[ctx] pySCENIC ctx failed with exit code {res.returncode}", file=sys.stderr)
    sys.exit(res.returncode)
else:
    print("[ctx] pySCENIC ctx finished successfully.")
