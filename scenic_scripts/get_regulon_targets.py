#!/usr/bin/env python3
import argparse, sys, csv

def load_gmt(path):
    # GMT: set_name \t description \t gene1 \t gene2 ...
    d = {}
    with open(path, newline='') as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                name = parts[0]
                genes = [g for g in parts[2:] if g]
                d[name] = genes
    return d

def load_tf_target_csv(path):
    # CSV: TF,target
    d = {}
    with open(path, newline='') as f:
        r = csv.DictReader(f)
        if not set(["TF", "target"]).issubset(r.fieldnames or []):
            sys.exit(f"[error] {path} must have columns TF,target")
        for row in r:
            tf = str(row["TF"]).strip()
            tg = str(row["target"]).strip()
            if tf and tg:
                d.setdefault(tf, []).append(tg)
    # de-duplicate while preserving order
    for k, v in d.items():
        seen = set(); out=[]
        for g in v:
            if g not in seen:
                out.append(g); seen.add(g)
        d[k] = out
    return d

def main():
    ap = argparse.ArgumentParser(description="Lookup regulon targets by name.")
    ap.add_argument("--regulons", required=True,
                    help="Path to regulons file (.gmt from AUCell or TF,target CSV).")
    ap.add_argument("--name", help="Regulon/TF name to fetch (must match the set name/TF).")
    ap.add_argument("--list", action="store_true", help="List all regulon names and exit.")
    ap.add_argument("--out", default="", help="Optional output file to write targets (one per line).")
    args = ap.parse_args()

    if args.regulons.lower().endswith(".gmt"):
        mapping = load_gmt(args.regulons)
    else:
        mapping = load_tf_target_csv(args.regulons)

    if args.list:
        for k in sorted(mapping.keys()):
            print(k)
        return

    if not args.name:
        sys.exit("[error] pass --name REGULON (or use --list to see available names)")

    name = args.name
    if name not in mapping:
        # try a lenient match (case-insensitive exact)
        alt = {k.lower(): k for k in mapping}
        if name.lower() in alt:
            name = alt[name.lower()]
        else:
            sys.exit(f"[error] regulon '{args.name}' not found. Try --list to see options.")

    genes = mapping[name]
    if not genes:
        print(f"[warn] regulon '{name}' has no targets in {args.regulons}", file=sys.stderr)

    if args.out:
        with open(args.out, "w") as f:
            for g in genes:
                f.write(g + "\n")
        print(f"Wrote {len(genes)} targets for '{name}' to {args.out}")
    else:
        for g in genes:
            print(g)

if __name__ == "__main__":
    main()
