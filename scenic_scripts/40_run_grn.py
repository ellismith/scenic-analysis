#!/usr/bin/env python3
"""
40_run_grn.py
Run pySCENIC GRN step (GRNBoost2) on a loom file.
Now flexible: paths passed via CLI arguments instead of hard-coded defaults.
"""

import argparse, os, sys, time, subprocess, importlib

parser = argparse.ArgumentParser(description="Run pySCENIC GRN (env-forced)")
parser.add_argument("--loom_in", required=True,
                    help="Input loom file (from liftover/convert step)")
parser.add_argument("--tfs", required=True,
                    help="TF list (txt) for GRN inference")
parser.add_argument("--out", default="adj.csv",
                    help="Output adjacency CSV")
parser.add_argument("--workers", type=int, default=20,
                    help="Number of worker threads")
args = parser.parse_args()

os.makedirs("logs", exist_ok=True)
stamp   = time.strftime("%Y%m%d_%H%M%S")
logpath = f"logs/grn_{stamp}.log"

# -------------------
# Preflight checks
# -------------------
pre = []
try:
    import loompy
    ds = loompy.connect(args.loom_in, "r")
    pre.append(f"LOOM shape rows×cols = {ds.shape}")
    pre.append(f"row attrs: {list(ds.ra.keys())[:4]}")
    pre.append(f"col attrs: {list(ds.ca.keys())[:4]}")
    ds.close()
except Exception as e:
    pre.append(f"ERROR opening loom: {e}")

def ver(modname, import_name=None):
    try:
        m = importlib.import_module(import_name or modname)
        return getattr(m, "__version__", "unknown")
    except Exception as e:
        return f"import error: {e}"

pre.append("VERSIONS: pyscenic="    + ver("pyscenic"))
pre.append("VERSIONS: dask="        + ver("dask"))
pre.append("VERSIONS: distributed=" + ver("distributed"))
pre.append("VERSIONS: dask_expr="   + ver("dask._expr", "dask._expr"))

# -------------------
# Build command
# -------------------
time_bin = "/usr/bin/time" if os.path.exists("/usr/bin/time") else None
cmd_core = [
    sys.executable, "-m", "pyscenic.cli.pyscenic", "grn",
    args.loom_in, args.tfs,
    "-o", args.out,
    "--num_workers", str(args.workers),
    "--sparse",
    "--gene_attribute", "Gene",
    "--cell_id_attribute", "CellID",
    "--method", "grnboost2",
]
cmd = ([time_bin, "-v"] + cmd_core) if time_bin else cmd_core

# -------------------
# Environment: pin threads, ignore user-site
# -------------------
env = os.environ.copy()
env.update({
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "PYTHONNOUSERSITE": "1",
})

# -------------------
# Run
# -------------------
with open(logpath, "w") as LOG:
    LOG.write("=== GRN preflight ===\n")
    for line in pre:
        LOG.write(line + "\n")
    LOG.write("\n=== Command ===\n" + " ".join(cmd) + "\n\n")
    LOG.flush()
    print(f"[GRN] logging to {logpath}")
    rc = subprocess.call(cmd, env=env, stdout=LOG, stderr=subprocess.STDOUT)

print(f"[GRN] exit code {rc}  (log: {logpath})")
sys.exit(rc)
