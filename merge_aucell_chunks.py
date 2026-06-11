#!/usr/bin/env python3
import sys, glob, os, anndata as ad

if len(sys.argv) != 3:
    print("Usage: python merge_aucell_chunks.py <chunks_dir> <output_file>")
    sys.exit(1)

chunks_dir = sys.argv[1]
out_file = sys.argv[2]

# Build glob pattern inside the script
pattern = os.path.join(chunks_dir, "*_aucell.h5ad")
files = sorted(glob.glob(pattern))

if not files:
    print(f"ERROR: No chunk files found in: {chunks_dir}")
    sys.exit(1)

print(f"Found {len(files)} chunks to merge.")
print("First few:", files[:3])

# Load first chunk as base
print("Loading base:", files[0])
adata_merged = ad.read_h5ad(files[0])

# Iteratively append others
for f in files[1:]:
    print("Appending:", f)
    adata_next = ad.read_h5ad(f)
    adata_merged = ad.concat([adata_merged, adata_next], axis=0, join="outer", merge="same")

print("Final shape:", adata_merged.shape)
adata_merged.write(out_file, compression="gzip")
print("✅ Wrote merged file:", out_file)
