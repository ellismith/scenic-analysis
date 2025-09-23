#!/usr/bin/env python3
"""
60_aucell.py — run pySCENIC AUCELL.
Flexible: can pass in paths via CLI, or fall back to config.py defaults.
"""

import os, shlex, subprocess, sys, shutil, argparse
from config import CT, WORK_DIR, PYSCENIC_OUT

# -------------------
# CLI args
# -------------------
parser = argparse.ArgumentParser(description="Run pySCENIC AUCELL")
parser.add_argument("--expr", default=str(WORK_DIR / f"{CT}_liftover_adata.h5ad"),
                    help="Expression matrix (H5AD) with human orthologs")
parser.add_argument("--reg", default=str(WORK_DIR / f"{CT}_reg.csv"),
                    help="Regulons CSV from ctx step")
parser.add_argument("--out", default=str(PYSCENIC_OUT),
                    help="Output h5ad with AUCELL scores")
parser.add_argument("--workers", type=int, default=20,
                    help="Number of worker threads")
args = parser.parse_args()

# -------------------
# pyscenic binary
# -------------------
ENV_PYSCENIC = "/scratch/easmit31/conda_envs/pyscenic/bin/pyscenic"
exe = ENV_PYSCENIC if os.path.exists(ENV_PYSCENIC) else shutil.which("pyscenic")
if not exe:
    print("[aucell] ERROR: pyscenic not found on PATH.", file=sys.stderr)
    sys.exit(127)

print("[aucell] Using pyscenic:", exe)

cmd = [
    exe, "aucell",
    args.expr, args.reg,
    "--output", args.out,
    "--num_workers", str(args.workers),
]

print("[aucell] Running command:\n ", " ".join(shlex.quote(c) for c in cmd))

# -------------------
# Run
# -------------------
ret = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr).returncode
if ret != 0:
    print(f"[aucell] pySCENIC AUCELL failed with exit code {ret}", file=sys.stderr)
    sys.exit(ret)

print(f"[aucell] Done. Output: {args.out}")
