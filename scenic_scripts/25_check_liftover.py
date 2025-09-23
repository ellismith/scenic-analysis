import os, h5py

files = [
    "/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv_liftover_sparse.h5ad",
    "/scratch/easmit31/GRN_copy/scenic/h5ad_files/gaba_adults_allLv_liftover_sparse.h5ad",
]

for path in files:
    print("\n=== Checking:", path)
    if not os.path.exists(path):
        print("⚠️ File not found (maybe job still running)")
        continue

    with h5py.File(path, "r") as f:
        print("Keys:", list(f.keys()))

        if isinstance(f["X"], h5py.Dataset):
            print("Shape of X:", f["X"].shape)
        else:
            # Sparse: shape is an attribute
            print("Sparse X shape:", tuple(f["X"].attrs["shape"]))

        genes = [g.decode() for g in f["var"]["_index"][:5]]
        cells = [c.decode() for c in f["obs"]["_index"][:5]]
        print("First 5 genes:", genes)
        print("First 5 cells:", cells)

print("\n✅ Done")
