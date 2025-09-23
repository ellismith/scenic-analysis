#!/usr/bin/env python3
import argparse, glob, pandas as pd, sys

ap = argparse.ArgumentParser()
ap.add_argument("--pattern", default="chunk_*_auc.csv")
ap.add_argument("--out", default="aucell_full_auc_mtx.csv")
args = ap.parse_args()

files = sorted(glob.glob(args.pattern))
if not files:
    sys.exit("[err] no chunk files found")

print("[merge] files:", len(files))
df0 = pd.read_csv(files[0])
cols = df0.columns
parts = [df0]

for f in files[1:]:
    df = pd.read_csv(f)
    df = df.reindex(columns=cols)  # align columns/order
    parts.append(df)

out = pd.concat(parts, axis=0, ignore_index=False)
out.to_csv(args.out, index=False)
print("[merge] wrote", args.out, "shape=", out.shape)
