#!/usr/bin/env python3
import argparse, os, math
import numpy as np
import scanpy as sc

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--chunk-size", type=int, default=200000)
ap.add_argument("--outdir", default="chunks")
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)
adata = sc.read_h5ad(args.h5ad, backed="r")
n = adata.n_obs
k = args.chunk_size
num_chunks = math.ceil(n / k)
print(f"[info] cells={n}, chunk_size={k}, chunks={num_chunks}")

# Just read the index (obs_names)
names = np.array(adata.obs_names)

for i in range(num_chunks):
    s, e = i*k, min((i+1)*k, n)
    out = os.path.join(args.outdir, f"cell_ids_chunk_{i:04d}.txt")
    with open(out, "w") as f:
        for cell in names[s:e]:
            f.write(f"{cell}\n")
print("[done] wrote:", num_chunks, "files in", args.outdir)
