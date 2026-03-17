#!/usr/bin/env python3
import argparse
import pathlib


def read_sscore(path):
    with open(path) as f:
        header = f.readline().strip().split()
        idx = {k: i for i, k in enumerate(header)}
        needed = ["IID", "SCORE1_SUM"]
        for k in needed:
            if k not in idx:
                raise ValueError(f"{path} missing {k}")

        rows = []
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            rows.append((parts[idx["IID"]], float(parts[idx["SCORE1_SUM"]])))
    return rows


def percentile_rank(values, x):
    n = len(values)
    lt = sum(v < x for v in values)
    eq = sum(v == x for v in values)
    return 100.0 * (lt + 0.5 * eq) / n


def main():
    ap = argparse.ArgumentParser(description="Summarize PRS percentile for a target sample from PLINK .sscore files")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--target-iid", required=True, help="IID for the true sample (not synthetic)")
    args = ap.parse_args()

    results_dir = pathlib.Path(args.results_dir)
    files = sorted(results_dir.glob("*.sscore"))
    if not files:
        raise SystemExit("No .sscore files found")

    print("PGS_ID\tN\tTARGET_SCORE\tPERCENTILE")
    for f in files:
        rows = read_sscore(f)
        vals = [v for _, v in rows]
        lookup = dict(rows)
        if args.target_iid not in lookup:
            continue
        x = lookup[args.target_iid]
        p = percentile_rank(vals, x)
        pgs_id = f.stem
        print(f"{pgs_id}\t{len(vals)}\t{x:.6f}\t{p:.2f}")


if __name__ == "__main__":
    main()
