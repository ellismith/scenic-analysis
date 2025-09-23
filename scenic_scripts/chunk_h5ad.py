import scanpy as sc
import sys
import numpy as np
import os

input_file = sys.argv[1]
chunk_size = int(sys.argv[2])  # cells per chunk
output_prefix = sys.argv[3]

print(f"Loading {input_file}...")
adata = sc.read_h5ad(input_file)
print(f"Total cells: {adata.n_obs}")

n_chunks = int(np.ceil(adata.n_obs / chunk_size))
print(f"Creating {n_chunks} chunks of ~{chunk_size} cells each")

# Create output directory based on prefix
output_dir = f"/scratch/easmit31/GRN_copy/scenic/h5ad_files/{output_prefix}_chunks"
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory: {output_dir}")

for i in range(n_chunks):
    start_idx = i * chunk_size
    end_idx = min((i + 1) * chunk_size, adata.n_obs)
    
    chunk = adata[start_idx:end_idx].copy()
    output_file = os.path.join(output_dir, f"{output_prefix}_{i:03d}.h5ad")
    
    print(f"Writing chunk {i+1}/{n_chunks}: {chunk.n_obs} cells -> {output_file}")
    chunk.write(output_file)

print("Done!")
