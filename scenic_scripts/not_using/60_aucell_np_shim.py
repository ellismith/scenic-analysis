#!/usr/bin/env python3
import sys, time, traceback
# --- inputs ---
EXPR = "/scratch/easmit31/GRN_copy/glutamatergic-neurons_L2-3_TEST/glutamatergic-neurons_L2-3_TEST_filtered_scenic.loom"
REG  = "/scratch/easmit31/GRN_copy/glutamatergic-neurons_L2-3_TEST/glutamatergic-neurons_L2-3_TEST_reg.csv"
OUT  = "/scratch/easmit31/GRN_copy/glutamatergic-neurons_L2-3_TEST/glutamatergic-neurons_L2-3_TEST_pyscenic_output.h5ad"

print("[AUCELL] start")
# Re-add numpy deprecated aliases if sitecustomize removed them
import numpy as np
for name, val in (("object", object), ("bool", bool), ("int", int), ("float", float)):
    if not hasattr(np, name):
        setattr(np, name, val)

from pyscenic.cli.pyscenic import main as pyscenic_main
argv = ["pyscenic", "aucell", EXPR, REG, "--output", OUT, "--num_workers", "6"]
print("[AUCELL] args:", " ".join(argv))
t0 = time.time()
try:
    sys.argv = argv
    pyscenic_main()
    print(f"[AUCELL] done in {time.time()-t0:.1f}s → {OUT}")
except SystemExit as e:
    code = e.code if isinstance(e.code, int) else 1
    print(f"[AUCELL] EXIT code {code}", file=sys.stderr); raise
except Exception:
    print("[AUCELL] ERROR\n" + traceback.format_exc(), file=sys.stderr); sys.exit(1)
