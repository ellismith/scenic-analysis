#!/usr/bin/env python3
import argparse, pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Summarize AUCell matrix by regulon (continuous, no binarization).")
    ap.add_argument("--input", required=True, help="Path to auc_mtx.csv")
    ap.add_argument("--out",   required=True, help="Output CSV, e.g. regulon_summary.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    cell_col = df.columns[0]
    mat = df.drop(columns=[cell_col])

    summary = pd.DataFrame({
        "mean":   mat.mean(axis=0),
        "std":    mat.std(axis=0),
        "var":    mat.var(axis=0),
        "min":    mat.min(axis=0),
        "q25":    mat.quantile(0.25, axis=0),
        "median": mat.median(axis=0),
        "q75":    mat.quantile(0.75, axis=0),
        "max":    mat.max(axis=0),
    }).sort_values("var", ascending=False)

    print(f"Cells: {df.shape[0]} | Regulons: {mat.shape[1]}")
    print(f"Writing summary to {args.out}")
    summary.to_csv(args.out, index=True)

if __name__ == "__main__":
    main()
