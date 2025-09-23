#!/usr/bin/env python3
"""
96_all_clusters_regulon_heatmaps.py

Purpose
-------
Loop through all clusters in <groupby> and call 96_cluster_regulon_heatmap.py
for each one. Uses backed mode to avoid loading the full h5ad into memory.
"""

import argparse
import scanpy as sc
import subprocess
import sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="AUCell h5ad file")
    ap.add_argument("--celltype", required=True, help="Cell type prefix (e.g. gaba, astros)")
    ap.add_argument("--groupby", default="louvain", help="obs column to group by (default: louvain)")
    ap.add_argument("--max_regulons", type=int, default=50, help="Top regulons to show")
    ap.add_argument("--figdir", required=True, help="Directory to save figures")
    args = ap.parse_args()

    # Ensure we call the correct python
    project_dir = Path(__file__).resolve().parent
    python_bin = "/scratch/easmit31/conda_envs/pyscenic/bin/python"
    cluster_script = project_dir / "96_cluster_regulon_heatmap.py"

    print(f"[all-clusters] Reading in backed mode: {args.input}")
    adata = sc.read_h5ad(args.input, backed="r")

    if args.groupby not in adata.obs:
        raise ValueError(f"Groupby '{args.groupby}' not found in obs")

    clusters = sorted(adata.obs[args.groupby].astype(str).unique())
    print(f"[all-clusters] Found {len(clusters)} clusters: {clusters}")

    # Loop over clusters, call cluster heatmap script with correct python
    for cl in clusters:
        print(f"[all-clusters] Running cluster {cl}")
        cmd = [
            python_bin, str(cluster_script),
            "--input", args.input,
            "--celltype", args.celltype,
            "--groupby", args.groupby,
            "--cluster", str(cl),
            "--max_regulons", str(args.max_regulons),
            "--figdir", args.figdir,
        ]
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
