#!/usr/bin/env python3
"""
Convert a cleaned pySCENIC ctx output (reg_clean.csv) into:
1. A GMT file (for reference/compatibility)
2. An AUCell-ready CSV file (TF,target pairs)

Collapses multiple motif rows per TF into one regulon.
"""

import sys
import pandas as pd
import re
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 3:
    sys.exit("Usage: python make_regulons_gmt.py <reg_clean.csv> <out_prefix>")

in_csv, out_prefix = sys.argv[1], sys.argv[2]
df = pd.read_csv(in_csv)
targets_col = df.columns[-2]

# collect all targets per TF
tf_targets = defaultdict(set)
for _, row in df.iterrows():
    tf = row["TF"]
    raw_targets = str(row[targets_col])
    genes = re.findall(r"\('([^']+)'", raw_targets)
    tf_targets[tf].update(genes)

# write GMT
gmt_path = Path(out_prefix + "_regulons.gmt")
with open(gmt_path, "w") as f:
    for tf, genes in tf_targets.items():
        if genes:
            line = [f"{tf}_regulon", "na"] + sorted(genes)
            f.write("\t".join(line) + "\n")
print(f"Wrote {gmt_path} with {len(tf_targets)} regulons")

# write CSV
csv_path = Path(out_prefix + "_regulons.csv")
with open(csv_path, "w") as f:
    f.write("regulon,target\n")
    for tf, genes in tf_targets.items():
        for g in sorted(genes):
            f.write(f"{tf},{g}\n")
print(f"Wrote {csv_path} with {len(tf_targets)} regulons")
