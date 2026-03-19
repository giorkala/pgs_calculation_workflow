#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import math
import pathlib
import re
import shutil
import subprocess
import sys
import urllib.request
from collections import defaultdict

REST_BASE = "https://www.pgscatalog.org/rest/score/"
PROJECT_ROOT = pathlib.Path(".").resolve()
MERGED_PFILE = PROJECT_ROOT / "merge_1000g/eur_plus_imputed"
MERGED_PVAR = PROJECT_ROOT / "merge_1000g/eur_plus_imputed.pvar.zst"
RESULTS_DIR = PROJECT_ROOT / "results/plink2_batch"
SCORES_DIR = PROJECT_ROOT / "scores/batch"
RAW_DIR = SCORES_DIR / "raw"
PLINK_DIR = SCORES_DIR / "plink2"
FAKE_RE = re.compile(r"^FAKE_[0-9]+_FAKE_[0-9]+$")


def read_panel(path: pathlib.Path):
    out = []
    with path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            out.append((row["PGS_ID"], row["TRAIT"]))
    return out


def fetch_text(url: str):
    with urllib.request.urlopen(url) as resp:
        return resp.read().decode("utf-8")


def extract_ftp_url(js: str):
    data = json.loads(js)
    if data.get("ftp_harmonized_scoring_files"):
        g37 = data["ftp_harmonized_scoring_files"].get("GRCh37")
        if g37 and g37.get("positions"):
            return g37["positions"]
    return data.get("ftp_scoring_file")


def download(url: str, dest: pathlib.Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)


def open_text(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_header_and_rows(path: pathlib.Path):
    with open_text(path) as f:
        data_lines = (line for line in f if not line.startswith('#'))
        reader = csv.reader(data_lines, delimiter='	', quotechar='"')
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"No header found in {path}")
        yield header
        for parts in reader:
            if not parts:
                continue
            yield parts


def choose_cols(header):
    idx = {k: i for i, k in enumerate(header)}
    col = {}
    for key, options in {
        "chrom": ["hm_chr", "chr_name"],
        "pos": ["hm_pos", "chr_position"],
        "ea": ["effect_allele"],
        "oa": ["other_allele", "hm_inferOtherAllele"],
        "weight": ["effect_weight"],
        "rsid": ["rsID"],
    }.items():
        for opt in options:
            if opt in idx:
                col[key] = idx[opt]
                break
    return col


def build_direct_scorefile(raw_path: pathlib.Path, out_path: pathlib.Path):
    rows = parse_header_and_rows(raw_path)
    header = next(rows)
    col = choose_cols(header)
    required = ["chrom", "pos", "ea", "oa", "weight"]
    if not all(k in col for k in required):
        return False, 0, "needs_remap"
    n = 0
    with out_path.open("w") as out:
        out.write("ID\tALLELE\tWEIGHT\n")
        for parts in rows:
            chrom = parts[col["chrom"]]
            pos = parts[col["pos"]]
            ea = parts[col["ea"]]
            oa = parts[col["oa"]]
            wt = parts[col["weight"]]
            if not all([chrom, pos, ea, oa, wt]):
                continue
            if oa == "?":
                continue
            out.write(f"{chrom}:{pos}:{oa}:{ea}\t{ea}\t{wt}\n")
            n += 1
    return True, n, "direct"


def build_remapped_scorefile(raw_path: pathlib.Path, out_path: pathlib.Path):
    zstdcat = shutil.which("zstdcat")
    if not zstdcat:
        raise RuntimeError("zstdcat not found")
    rows = parse_header_and_rows(raw_path)
    header = next(rows)
    col = choose_cols(header)
    required = ["chrom", "pos", "ea", "weight"]
    if not all(k in col for k in required):
        return False, 0, "unmappable"
    by_pos = defaultdict(list)
    total = 0
    for parts in rows:
        chrom = parts[col["chrom"]]
        pos = parts[col["pos"]]
        ea = parts[col["ea"]]
        wt = parts[col["weight"]]
        rsid = parts[col["rsid"]] if "rsid" in col else ""
        if not all([chrom, pos, ea, wt]):
            continue
        by_pos[(chrom, pos)].append((ea, wt, rsid))
        total += 1
    proc = subprocess.Popen([zstdcat, str(MERGED_PVAR)], stdout=subprocess.PIPE, text=True)
    idx = None
    hdr = None
    matched = []
    for line in proc.stdout:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").lstrip("#").split("\t")
            idx = {c: i for i, c in enumerate(hdr)}
            continue
        parts = line.rstrip("\n").split("\t")
        if idx is None or len(parts) < len(hdr):
            continue
        key = (parts[idx["CHROM"]], parts[idx["POS"]])
        if key not in by_pos:
            continue
        vid = parts[idx["ID"]]
        ref = parts[idx["REF"]]
        alt = parts[idx["ALT"]]
        if "," in alt:
            continue
        for ea, wt, rsid in by_pos[key]:
            if ea == ref or ea == alt:
                matched.append((vid, ea, wt))
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"zstdcat failed with code {ret}")
    with out_path.open("w") as out:
        out.write("ID\tALLELE\tWEIGHT\n")
        for vid, ea, wt in matched:
            out.write(f"{vid}\t{ea}\t{wt}\n")
    return True, len(matched), f"remapped_from_{total}"


def run(cmd):
    res = subprocess.run(cmd, text=True, capture_output=True)
    if res.returncode != 0:
        raise RuntimeError(res.stderr or res.stdout)
    return res


def plink_score(score_path: pathlib.Path, out_prefix: pathlib.Path):
    cmd = [
        "bash", "-lc",
        "module load cellgen/plink/2.00 >/dev/null 2>&1 && "
        f"plink2 --pfile {MERGED_PFILE} vzs --score {score_path} 1 2 3 header list-variants --out {out_prefix}"
    ]
    return run(cmd)


def parse_log(log_path: pathlib.Path):
    txt = log_path.read_text()
    proc = None
    miss = 0
    m = re.search(r"--score: ([0-9]+) variants processed", txt)
    if m:
        proc = int(m.group(1))
    m = re.search(r"Warning: --score: ([0-9]+) entries .* skipped due to missing variant IDs", txt)
    if m:
        miss = int(m.group(1))
    return proc, miss


def summarize_sscore(path: pathlib.Path, target_iid: str):
    with path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
    score_col = "SCORE1_AVG" if "SCORE1_AVG" in rows[0] else "SCORE1_SUM"
    refs = [float(r[score_col]) for r in rows if r["IID"] != target_iid and not FAKE_RE.match(r["IID"])]
    target = next(float(r[score_col]) for r in rows if r["IID"] == target_iid)
    fake_vals = [float(r[score_col]) for r in rows if FAKE_RE.match(r["IID"])]
    le = sum(1 for x in refs if x <= target)
    mean = sum(refs) / len(refs)
    sd = math.sqrt(sum(x * x for x in refs) / len(refs) - mean * mean) if refs else 0.0
    z = (target - mean) / sd if sd > 0 else 0.0
    return {
        "target_score": target,
        "percentile": 100.0 * le / len(refs) if refs else 0.0,
        "ref_n": len(refs),
        "z": z,
        "fake_min": min(fake_vals) if fake_vals else "",
        "fake_max": max(fake_vals) if fake_vals else "",
    }


def parse_original_count(prep_mode, prepared):
    m = re.match(r"remapped_from_([0-9]+)", str(prep_mode))
    if m:
        return int(m.group(1))
    return int(prepared or 0)


def technical_confidence(prep_mode, prepared, processed):
    if not processed or processed == 0 or not prepared:
        return "failed"
    original = parse_original_count(prep_mode, prepared)
    recovery = prepared / original if original else 0.0
    processed_frac = processed / prepared if prepared else 0.0

    # Sparse scores are inherently less robust to missingness/noise.
    if prepared < 100:
        return "low"
    if prepared < 1000:
        return "medium"

    if processed_frac < 0.95:
        return "low"
    if recovery >= 0.9 and prepared >= 100000:
        return "high"
    if recovery >= 0.6 and prepared >= 10000:
        return "medium"
    return "low"


def overall_confidence(prep_mode, prepared, processed):
    tech = technical_confidence(prep_mode, prepared, processed)
    if tech == "failed":
        return "failed"
    # We are not fetching score-level validation metrics (e.g. external R2/AUC),
    # so cap overall interpretive confidence at medium.
    if tech == "high":
        return "medium"
    return tech


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default=str(PROJECT_ROOT / "config/pgs_panel.tsv"))
    ap.add_argument("--target-iid", required=True, help="IID of the real target sample in the merged pfile.")
    args = ap.parse_args()

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PLINK_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for pgs_id, trait in read_panel(pathlib.Path(args.panel)):
        rest_url = f"{REST_BASE}{pgs_id}"
        try:
            js = fetch_text(rest_url)
            ftp_url = extract_ftp_url(js)
            if not ftp_url:
                raise RuntimeError("no scoring file URL")
            raw_path = RAW_DIR / pathlib.Path(ftp_url).name
            if not raw_path.exists():
                download(ftp_url, raw_path)
            score_path = PLINK_DIR / f"{pgs_id}.score.tsv"
            ok, prepared, prep_mode = build_remapped_scorefile(raw_path, score_path)
            if not ok:
                ok, prepared, prep_mode = build_direct_scorefile(raw_path, score_path)
            out_prefix = RESULTS_DIR / pgs_id
            plink_score(score_path, out_prefix)
            processed, missing_ids = parse_log(pathlib.Path(f"{out_prefix}.log"))
            stats = summarize_sscore(pathlib.Path(f"{out_prefix}.sscore"), args.target_iid)
            summary_rows.append({
                "PGS_ID": pgs_id,
                "TRAIT": trait,
                "PREP_MODE": prep_mode,
                "ORIGINAL_SCORE_VARIANTS": parse_original_count(prep_mode, prepared),
                "PROCESSED_VARIANTS": processed,
                "RECOVERY_FRACTION": round(processed / parse_original_count(prep_mode, prepared), 4) if processed is not None and parse_original_count(prep_mode, prepared) else "",
                "SKIPPED_MISSING_ID": missing_ids,
                "TARGET_PERCENTILE": round(stats["percentile"], 2),
                "TARGET_Z": round(stats["z"], 3),
                "TARGET_SCORE": stats["target_score"],
                "FAKE_MIN": stats["fake_min"],
                "FAKE_MAX": stats["fake_max"],
                "TECHNICAL_CONFIDENCE": technical_confidence(prep_mode, prepared, processed),
                "OVERALL_CONFIDENCE": overall_confidence(prep_mode, prepared, processed),
            })
            print(f"[OK] {pgs_id} {trait}: {processed}/{prepared} scored via {prep_mode}")
        except Exception as e:
            summary_rows.append({
                "PGS_ID": pgs_id,
                "TRAIT": trait,
                "PREP_MODE": "failed",
                "ORIGINAL_SCORE_VARIANTS": "",
                "PROCESSED_VARIANTS": "",
                "RECOVERY_FRACTION": "",
                "SKIPPED_MISSING_ID": "",
                "TARGET_PERCENTILE": "",
                "TARGET_Z": "",
                "TARGET_SCORE": "",
                "FAKE_MIN": "",
                "FAKE_MAX": "",
                "TECHNICAL_CONFIDENCE": "failed",
                "OVERALL_CONFIDENCE": "failed",
                "ERROR": str(e).replace("\n", " "),
            })
            print(f"[WARN] {pgs_id} {trait}: {e}", file=sys.stderr)

    out = RESULTS_DIR / "summary.tsv"
    fields = [
        "PGS_ID", "TRAIT", "PREP_MODE", "ORIGINAL_SCORE_VARIANTS", "PROCESSED_VARIANTS",
        "RECOVERY_FRACTION", "SKIPPED_MISSING_ID", "TARGET_PERCENTILE", "TARGET_Z",
        "TARGET_SCORE", "FAKE_MIN", "FAKE_MAX", "TECHNICAL_CONFIDENCE", "OVERALL_CONFIDENCE", "ERROR"
    ]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"[OK] wrote {out}")

if __name__ == "__main__":
    main()
