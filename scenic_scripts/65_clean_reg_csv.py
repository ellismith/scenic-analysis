#!/usr/bin/env python3
"""
65_clean_reg_csv.py

Purpose
-------
Clean pySCENIC regulon CSV outputs by skipping metadata lines
and saving a proper CSV with only the real table.

Usage
-----
python 65_clean_reg_csv.py --infile <reg.csv> --outfile <reg_clean.csv>
"""

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Clean pySCENIC regulon CSV.")
    parser.add_argument("--infile", required=True, help="Path to input reg.csv")
    parser.add_argument("--outfile", required=True, help="Path to output reg_clean.csv")
    args = parser.parse_args()

    # assume the first valid header is the first line that starts with "TF,"
    header_line = None
    with open(args.infile) as f:
        for i, line in enumerate(f):
            if line.strip().startswith("TF,"):
                header_line = i
                break

    if header_line is None:
        raise RuntimeError(f"Could not find a valid header in {args.infile}")

    df = pd.read_csv(args.infile, skiprows=header_line)
    df.to_csv(args.outfile, index=False)
    print(f"[clean] Wrote cleaned file to {args.outfile} with shape {df.shape}")

if __name__ == "__main__":
    main()
