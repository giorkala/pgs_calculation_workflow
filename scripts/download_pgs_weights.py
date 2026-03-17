#!/usr/bin/env python3
import argparse
import gzip
import pathlib
import re
import sys
import urllib.request

REST_BASE = "https://www.pgscatalog.org/rest/score/"


def read_ids(ids_file: pathlib.Path):
    ids = []
    with ids_file.open() as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                line = line.split("\t", 1)[0]
            ids.append(line)
    return ids


def fetch_json(url: str):
    with urllib.request.urlopen(url) as resp:
        return resp.read().decode("utf-8")


def extract_ftp_url(json_text: str):
    m = re.search(r'"ftp_scoring_file"\s*:\s*"([^"]+)"', json_text)
    if m:
        return m.group(1)
    m = re.search(
        r'"ftp_harmonized_scoring_files"\s*:\s*\{[^}]*"GRCh37"\s*:\s*\{[^}]*"positions"\s*:\s*"([^"]+)"',
        json_text,
        flags=re.S,
    )
    if m:
        return m.group(1)
    return None


def open_text_maybe_gz(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def convert_to_plink_scorefile(in_file: pathlib.Path, out_file: pathlib.Path):
    # Output columns: ID ALLELE WEIGHT
    with open_text_maybe_gz(in_file) as fin, out_file.open("w") as fout:
        header = None
        idx = {}
        fout.write("ID\tALLELE\tWEIGHT\n")
        for line in fin:
            if line.startswith("#"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                idx = {k: i for i, k in enumerate(header)}
                if not any(k in idx for k in ("rsID", "chr_name")):
                    raise ValueError(f"Cannot find variant ID columns in {in_file}")
                for k in ("effect_allele", "effect_weight"):
                    if k not in idx:
                        raise ValueError(f"Missing column {k} in {in_file}")
                continue

            parts = line.rstrip("\n").split("\t")
            var_id = None
            if "rsID" in idx and idx["rsID"] < len(parts) and parts[idx["rsID"]]:
                var_id = parts[idx["rsID"]]
            elif all(k in idx for k in ("chr_name", "chr_position", "effect_allele", "other_allele")):
                chrom = parts[idx["chr_name"]]
                pos = parts[idx["chr_position"]]
                ea = parts[idx["effect_allele"]]
                oa = parts[idx["other_allele"]]
                var_id = f"{chrom}:{pos}:{oa}:{ea}"

            if not var_id:
                continue

            allele = parts[idx["effect_allele"]]
            weight = parts[idx["effect_weight"]]
            if allele and weight:
                fout.write(f"{var_id}\t{allele}\t{weight}\n")


def main():
    p = argparse.ArgumentParser(description="Download PGS scoring files and emit PLINK2 scorefiles.")
    p.add_argument("--ids", required=True, help="Text/TSV file with one PGS ID per line in first column.")
    p.add_argument("--outdir", default="scores", help="Output directory.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    args = p.parse_args()

    outdir = pathlib.Path(args.outdir)
    rawdir = outdir / "raw"
    plinkdir = outdir / "plink2"
    rawdir.mkdir(parents=True, exist_ok=True)
    plinkdir.mkdir(parents=True, exist_ok=True)

    pgs_ids = read_ids(pathlib.Path(args.ids))
    if not pgs_ids:
        print("No PGS IDs found in --ids file", file=sys.stderr)
        sys.exit(1)

    for pgs_id in pgs_ids:
        rest_url = f"{REST_BASE}{pgs_id}"
        print(f"[INFO] Querying {rest_url}")
        try:
            js = fetch_json(rest_url)
        except Exception as e:
            print(f"[WARN] REST query failed for {pgs_id}: {e}", file=sys.stderr)
            continue

        ftp_url = extract_ftp_url(js)
        if not ftp_url:
            print(f"[WARN] Could not find scoring file URL for {pgs_id}", file=sys.stderr)
            continue

        raw_name = pathlib.Path(ftp_url).name
        raw_path = rawdir / raw_name
        if raw_path.exists() and not args.overwrite:
            print(f"[INFO] Exists, skip download: {raw_path}")
        else:
            print(f"[INFO] Downloading {ftp_url}")
            try:
                urllib.request.urlretrieve(ftp_url, raw_path)
            except Exception as e:
                print(f"[WARN] Download failed for {pgs_id}: {e}", file=sys.stderr)
                continue

        out_score = plinkdir / f"{pgs_id}.score.tsv"
        try:
            convert_to_plink_scorefile(raw_path, out_score)
            print(f"[OK] Wrote {out_score}")
        except Exception as e:
            print(f"[WARN] Convert failed for {pgs_id}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
