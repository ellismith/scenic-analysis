#!/usr/bin/env python3
import os
import glob
import subprocess
import sys

def rebuild_chunk(chunk_file, output_dir):
    """Rebuild one chunk into CSR format and save"""
    base_name = os.path.basename(chunk_file).replace('.h5ad', '')
    rebuilt_file = os.path.join(output_dir, f"{base_name}_rebuilt.h5ad")

    rebuild_script = f'''
/scratch/easmit31/conda_envs/pyscenic/bin/python <<'PY'
import scanpy as sc
import scipy.sparse as sp

print("👉 Loading {chunk_file}")
adata = sc.read_h5ad("{chunk_file}")

print(f"Original X type: {{type(adata.X)}}")
print(f"Shape: {{adata.shape}}")

# Ensure X is CSR
if not isinstance(adata.X, sp.csr_matrix):
    print("Converting to CSR...")
    adata.X = sp.csr_matrix(adata.X)
else:
    print("Already CSR")

print(f"Final X type: {{type(adata.X)}}")
adata.write("{rebuilt_file}", compression="gzip")
print("✅ Saved rebuilt: {rebuilt_file}")
PY
'''
    return rebuild_script, rebuilt_file

def main():
    if len(sys.argv) != 2:
        print("Usage: python rebuild_all_chunks.py <chunks_dir_or_file>")
        sys.exit(1)

    target = sys.argv[1]

    if os.path.isdir(target):
        chunks = sorted(glob.glob(os.path.join(target, "*.h5ad")))
        rebuilt_dir = target.replace("_chunks", "_rebuilt")
    else:
        chunks = [target]
        rebuilt_dir = os.path.dirname(target).replace("_chunks", "_rebuilt")

    os.makedirs(rebuilt_dir, exist_ok=True)
    print(f"Rebuilding {len(chunks)} chunk(s)")
    print(f"Output → {rebuilt_dir}")

    for i, chunk_file in enumerate(chunks):
        print(f"\n=== Rebuilding {i+1}/{len(chunks)}: {os.path.basename(chunk_file)} ===")
        rebuild_script, rebuilt_file = rebuild_chunk(chunk_file, rebuilt_dir)
        result = subprocess.run(rebuild_script, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"❌ ERROR rebuilding {chunk_file}")
            print(result.stderr)
        else:
            print(result.stdout)
            print(f"✅ Rebuilt: {rebuilt_file}")

if __name__ == "__main__":
    main()
