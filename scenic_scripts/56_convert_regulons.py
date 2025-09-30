#!/usr/bin/env python3
import pandas as pd
import sys
import ast

if len(sys.argv) != 3:
    sys.exit("Usage: ./56_convert_regulons.py input_regulons.csv output_simple.csv")

input_csv, output_csv = sys.argv[1:]

print(f"Reading: {input_csv}")
df = pd.read_csv(input_csv, skiprows=3, header=None)

regulon_target_pairs = []

for idx, row in df.iterrows():
    tf_name = row[0]
    target_str = row[8]
    
    try:
        # It's a string, so parse it
        targets = ast.literal_eval(target_str)
        for gene, weight in targets:
            regulon_target_pairs.append({"regulon": tf_name, "target": gene})
    except Exception as e:
        print(f"Row {idx} ({tf_name}): {e}")
        continue

result_df = pd.DataFrame(regulon_target_pairs)
result_df.to_csv(output_csv, index=False)
print(f"✓ Wrote {len(result_df)} regulon-target pairs")
print(f"✓ Unique regulons: {result_df['regulon'].nunique()}")
