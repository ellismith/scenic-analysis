#!/usr/bin/env python3
"""
20_filter_adults.py

Safely filter an AnnData .h5ad file to keep only adult cells (age >= 1.0).
Does not load whole matrix into memory; copies obs/var and X row-by-row.
"""

import h5py
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--in_file", required=True, help="Input .h5ad")
parser.add_argument("--out_file", required=True, help="Output .h5ad (adults only)")
parser.add_argument("--age_key", default="age", help="obs column with age values")
args = parser.parse_args()

with h5py.File(args.in_file, "r") as fin:
    obs = pd.DataFrame({k: fin["obs"][k][:] for k in fin["obs"].keys()})
    obs.index = fin["obs"]["_index"][:].astype(str)
    obs.index.name = "_index"

    # build mask
    obs[args.age_key] = obs[args.age_key].astype(float)
    mask = obs[args.age_key] >= 1.0
    adult_ids = obs.index[mask]

    print(f"Total cells: {obs.shape[0]}, Adult cells: {mask.sum()}")

    # Now copy structure into new file
    with h5py.File(args.out_file, "w") as fout:
        # copy var as-is
        fin.copy("var", fout)

        # subset obs
        g_obs = fout.create_group("obs")
        for k in fin["obs"].keys():
            dset = fin["obs"][k]
            g_obs.create_dataset(k, data=dset[mask.to_numpy()])

        # copy X in chunks
        X = fin["X"]
        n_cells, n_genes = X.shape
        gX = fout.create_dataset("X", shape=(mask.sum(), n_genes), dtype=X.dtype)

        row = 0
        for i, keep in enumerate(mask.to_numpy()):
            if keep:
                gX[row, :] = X[i, :]
                row += 1

print(f"Wrote filtered h5ad: {args.out_file}")
