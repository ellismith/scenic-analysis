#!/usr/bin/env python3
"""
65_fix_h5ad.py
--------------
Force AnnData .h5ad to have CSR matrix with correct attrs for pySCENIC AUCell.
"""

import h5py

in_file  = "/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv_liftover_sparse.h5ad"
out_file = "/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv_liftover_ready.h5ad"

print("[fix_h5ad] Opening:", in_file)
with h5py.File(in_file, "r") as fin, h5py.File(out_file, "w") as fout:
    # Copy everything
    fin.copy("/", fout)

    # Make sure /X has shape + attrs
    X = fout["X"]
    if "shape" not in X:
        shape = [fin["X/indptr"].shape[0]-1, int(fin["X/indices"][:].max())+1]
        X.create_dataset("shape", data=shape, dtype="int64")
        print("[fix_h5ad] Added missing shape:", shape)

    # Add attrs
    X.attrs["encoding-type"] = np.string_("csr_matrix")
    X.attrs["encoding-version"] = np.string_("0.1.0")
    print("[fix_h5ad] Added attrs:", dict(X.attrs))

print("[fix_h5ad] Wrote:", out_file)
