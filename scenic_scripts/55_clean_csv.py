#!/usr/bin/env python3
"""
55_clean_csv.py - Clean pySCENIC CTX output files for AUCell input
This recreates the EXACT format that pySCENIC CTX originally outputs.
"""

import pandas as pd
import sys
import os
import argparse

def clean_scenic_csv(cell_type, input_dir="/scratch/easmit31/GRN_copy/scenic/h5ad_files"):
    """Recreate the exact pySCENIC CTX output format"""
    
    input_file = os.path.join(input_dir, f"{cell_type}_allLv_reg.csv")
    output_file = os.path.join(input_dir, f"{cell_type}_allLv_reg_clean.csv")
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    print(f"[clean_csv] Processing: {input_file}")
    
    # Read the data (skip the malformed headers, start from actual data)
    df = pd.read_csv(input_file, skiprows=3, header=None)
    
    print(f"[clean_csv] Read {df.shape[0]} regulons with {df.shape[1]} columns")
    
    # Recreate the EXACT header format that pySCENIC expects
    # This is what the original CTX output should look like
    header_lines = [
        ",".join([""] * 2 + ["Enrichment"] * 8),  # First header row
        ",".join([""] * 2 + ["AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax"]),  # Second header row  
        ",".join(["TF", "MotifID"] + [""] * 8)  # Third header row
    ]
    
    # Write the file with the exact format
    with open(output_file, 'w') as f:
        # Write headers
        for header in header_lines:
            f.write(header + '\n')
        
        # Write data
        for _, row in df.iterrows():
            f.write(','.join([f'"{str(val)}"' if ',' in str(val) else str(val) for val in row]) + '\n')
    
    print(f"[clean_csv] SUCCESS: Created properly formatted file: {output_file}")
    
    # Verify by checking first few lines
    with open(output_file, 'r') as f:
        lines = [f.readline().strip() for _ in range(7)]
    
    print("[clean_csv] Output file first 7 lines:")
    for i, line in enumerate(lines):
        print(f"  {i}: {line[:80]}...")

def main():
    parser = argparse.ArgumentParser(description="Clean pySCENIC CTX output for AUCell")
    parser.add_argument("cell_type", help="Cell type prefix (e.g., 'astros', 'gaba')")
    parser.add_argument("--input_dir", default="/scratch/easmit31/GRN_copy/scenic/h5ad_files",
                       help="Directory containing input files")
    
    args = parser.parse_args()
    clean_scenic_csv(args.cell_type, args.input_dir)

if __name__ == "__main__":
    main()
