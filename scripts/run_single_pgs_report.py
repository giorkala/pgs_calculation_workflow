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

import matplotlib.pyplot as plt

REST_BASE = "https://www.pgscatalog.org/rest/score/"
DEFAULT_PROJECT = pathlib.Path("/nfs/users/nfs_g/gk18/prs_application")
DEFAULT_TARGET_IID = "genome1402_genome1402"
FAKE_RE = re.compile(r"^FAKE_[0-9]+_FAKE_[0-9]+$")


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


def extract_trait_reported(js: str):
    data = json.loads(js)
    return data.get("trait_reported") or data.get("name") or data.get("pgs_id")


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
        reader = csv.reader(data_lines, delimiter='\t', quotechar='"')
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
        "ea": ["effect_allele", "ALLELE"],
        "oa": ["other_allele", "hm_inferOtherAllele"],
        "weight": ["effect_weight", "WEIGHT"],
        "id": ["ID"],
    }.items():
        for opt in options:
            if opt in idx:
                col[key] = idx[opt]
                break
    return col


def count_original_variants(raw_path: pathlib.Path):
    with open_text(raw_path) as fh:
        for line in fh:
            if line.startswith('#variants_number='):
                return int(line.strip().split('=', 1)[1])
            if not line.startswith('#'):
                break
    return None


def is_plink_scorefile(path: pathlib.Path):
    rows = parse_header_and_rows(path)
    header = next(rows)
    return header[:3] == ["ID", "ALLELE", "WEIGHT"]


def copy_plink_scorefile(src: pathlib.Path, dest: pathlib.Path):
    with open_text(src) as fin, dest.open('w') as fout:
        for line in fin:
            fout.write(line)
    with dest.open() as f:
        prepared = max(0, sum(1 for _ in f) - 1)
    return prepared, "direct_plink"


def build_remapped_scorefile(raw_path: pathlib.Path, out_path: pathlib.Path, merged_pvar: pathlib.Path):
    zstdcat = shutil.which("zstdcat")
    if not zstdcat:
        raise RuntimeError("zstdcat not found")
    rows = parse_header_and_rows(raw_path)
    header = next(rows)
    col = choose_cols(header)

    if "id" in col and "ea" in col and "weight" in col and "chrom" not in col and "pos" not in col:
        raise RuntimeError("Local score file has PLINK-style IDs but not genomic positions; use it directly only if IDs already match the dataset.")

    required = ["chrom", "pos", "ea", "weight"]
    if not all(k in col for k in required):
        raise RuntimeError(f"Cannot remap score file {raw_path}; missing genomic position/effect columns.")

    by_pos = {}
    total = 0
    for parts in rows:
        chrom = parts[col["chrom"]]
        pos = parts[col["pos"]]
        ea = parts[col["ea"]]
        wt = parts[col["weight"]]
        if not all([chrom, pos, ea, wt]):
            continue
        by_pos.setdefault((chrom, pos), []).append((ea, wt))
        total += 1

    proc = subprocess.Popen([zstdcat, str(merged_pvar)], stdout=subprocess.PIPE, text=True)
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
        for ea, wt in by_pos[key]:
            if ea == ref or ea == alt:
                matched.append((vid, ea, wt))
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"zstdcat failed with code {ret}")

    with out_path.open('w') as out:
        out.write("ID\tALLELE\tWEIGHT\n")
        for vid, ea, wt in matched:
            out.write(f"{vid}\t{ea}\t{wt}\n")
    return len(matched), f"remapped_from_{total}"


def run(cmd):
    res = subprocess.run(cmd, text=True, capture_output=True)
    if res.returncode != 0:
        raise RuntimeError(res.stderr or res.stdout)
    return res


def plink_score(pfile_prefix: pathlib.Path, score_path: pathlib.Path, out_prefix: pathlib.Path):
    cmd = [
        "bash", "-lc",
        "module load cellgen/plink/2.00 >/dev/null 2>&1 && "
        f"plink2 --pfile {pfile_prefix} vzs --score {score_path} 1 2 3 header list-variants --out {out_prefix}"
    ]
    run(cmd)


def parse_log(log_path: pathlib.Path):
    txt = log_path.read_text()
    m = re.search(r"--score: ([0-9]+) variants processed", txt)
    processed = int(m.group(1)) if m else None
    m = re.search(r"Warning: --score: ([0-9]+) entries .* skipped due to missing variant IDs", txt)
    skipped = int(m.group(1)) if m else 0
    return processed, skipped


def summarize_sscore(path: pathlib.Path, target_iid: str):
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    score_col = 'SCORE1_AVG' if 'SCORE1_AVG' in rows[0] else 'SCORE1_SUM'
    refs = [float(r[score_col]) for r in rows if r['IID'] != target_iid and not FAKE_RE.match(r['IID'])]
    target = next(float(r[score_col]) for r in rows if r['IID'] == target_iid)
    fakes = [float(r[score_col]) for r in rows if FAKE_RE.match(r['IID'])]
    mean = sum(refs) / len(refs)
    sd = math.sqrt(sum(x * x for x in refs) / len(refs) - mean * mean) if refs else 0.0
    z = (target - mean) / sd if sd > 0 else 0.0
    percentile = 100.0 * sum(1 for x in refs if x <= target) / len(refs) if refs else 0.0
    if percentile <= 25:
        quartile = 'Q1'
    elif percentile >= 75:
        quartile = 'Q4'
    elif percentile <= 50:
        quartile = 'Q2'
    else:
        quartile = 'Q3'
    return {
        'refs': refs,
        'target': target,
        'fake_min': min(fakes) if fakes else None,
        'fake_max': max(fakes) if fakes else None,
        'ref_n': len(refs),
        'percentile': percentile,
        'quartile': quartile,
        'outlier_q1_q4': quartile in {'Q1', 'Q4'},
        'z': z,
    }


def technical_confidence(original, processed):
    if not processed or processed == 0:
        return 'failed'
    recovery = processed / original if original else 0.0
    if processed < 100:
        return 'low'
    if processed < 1000:
        return 'medium'
    if recovery >= 0.9 and processed >= 100000:
        return 'high'
    if recovery >= 0.6 and processed >= 10000:
        return 'medium'
    return 'low'


def overall_confidence(original, processed):
    tech = technical_confidence(original, processed)
    if tech == 'failed':
        return 'failed'
    return 'medium' if tech == 'high' else tech


def make_plot(refs, target, fake_min, fake_max, title, png_path=None, pdf_path=None):
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.hist(refs, bins=30, color='#b8c9d9', edgecolor='#4a6572', linewidth=0.6)
    if fake_min is not None and fake_max is not None:
        ax.axvspan(fake_min, fake_max, color='#e74c3c', alpha=0.3)
    ax.axvline(target, color='#c0392b', linewidth=2.0)
    ax.set_title(title)
    ax.set_xlabel('PRS')
    ax.set_ylabel('EUR count')
    ax.grid(alpha=0.2, linewidth=0.5)
    if png_path or pdf_path:
        if png_path:
            fig.savefig(png_path, dpi=200, bbox_inches='tight')
        if pdf_path:
            fig.savefig(pdf_path, bbox_inches='tight')
    else:
        print("No output path for plot provided; skipping save.")
    plt.close(fig)


def write_report_txt(path: pathlib.Path, report: dict):
    lines = [
        f"PGS_ID: {report.get('pgs_id','')}",
        f"Trait: {report['trait_name']}",
        f"PFILE: {report['pfile_prefix']}",
        f"Original score variants: {report.get('original_score_variants','')}",
        f"Processed variants: {report.get('processed_variants','')}",
        f"Recovery fraction: {report.get('recovery_fraction','')}",
        f"Target IID: {report['target_iid']}",
        f"Percentile: {report['target_percentile']:.2f}",
        f"Quartile: {report['target_quartile']}",
        f"Outlier Q1/Q4: {'yes' if report['outlier_q1_q4'] else 'no'}",
        f"Target z: {report['target_z']:.3f}",
        f"Target score: {report['target_score']}",
        f"Technical confidence: {report['technical_confidence']}",
        f"Overall confidence: {report['overall_confidence']}",
        f"Score prep mode: {report['prep_mode']}",
        f"Skipped missing IDs: {report['skipped_missing_id']}",
        f"Histogram PNG: {report['histogram_png']}",
        f"Histogram PDF: {report['histogram_pdf']}",
    ]
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description='Single-PGS scorer/report generator.')
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--pgs-id', help='PGS Catalog ID to fetch and score, e.g. PGS002204')
    src.add_argument('--score-file', help='Existing local allele-weight file. Can be raw PGS Catalog file or PLINK score file.')
    ap.add_argument('--trait-name', help='Trait/display name. Required with --score-file unless it can be inferred.')
    ap.add_argument('--name', help='Short run name for output prefix. Defaults to PGS ID or score-file stem.')
    ap.add_argument('--pfile', default=str(DEFAULT_PROJECT / 'merge_1000g/eur_plus_imputed'), help='PLINK2 pfile prefix')
    ap.add_argument('--merged-pvar', default=str(DEFAULT_PROJECT / 'merge_1000g/eur_plus_imputed.pvar.zst'), help='Merged pvar.zst path for remapping')
    ap.add_argument('--target-iid', default=DEFAULT_TARGET_IID)
    ap.add_argument('--outdir', default=str(DEFAULT_PROJECT / 'results/single_pgs'))
    args = ap.parse_args()

    pfile_prefix = pathlib.Path(args.pfile)
    merged_pvar = pathlib.Path(args.merged_pvar)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pgs_id = args.pgs_id
    trait_name = args.trait_name
    original_variants = None

    if pgs_id:
        print('Fetching scoring file from the PGS catalog...')
        js = fetch_text(f"{REST_BASE}{pgs_id}")
        ftp_url = extract_ftp_url(js)
        if not ftp_url:
            raise RuntimeError(f'Could not find scoring file URL for {pgs_id}')
        trait_name = trait_name or extract_trait_reported(js) or pgs_id
        raw_path = outdir / pathlib.Path(ftp_url).name
        if not raw_path.exists():
            download(ftp_url, raw_path)
        original_variants = count_original_variants(raw_path)
        base_name = args.name or pgs_id
        print('Done with fetching.')
    else:
        raw_path = pathlib.Path(args.score_file)
        if not raw_path.exists():
            raise FileNotFoundError(raw_path)
        trait_name = trait_name or raw_path.stem
        base_name = args.name or raw_path.stem
        original_variants = count_original_variants(raw_path)

    score_path = outdir / f'{base_name}.score.tsv'
    if is_plink_scorefile(raw_path):
        prepared, prep_mode = copy_plink_scorefile(raw_path, score_path)
        if original_variants is None:
            original_variants = prepared
    else:
        prepared, prep_mode = build_remapped_scorefile(raw_path, score_path, merged_pvar)
        if original_variants is None:
            m = re.match(r'remapped_from_([0-9]+)', prep_mode)
            original_variants = int(m.group(1)) if m else prepared

    print('Calling PLINK for scoring...')
    out_prefix = outdir / base_name
    plink_score(pfile_prefix, score_path, out_prefix)
    print('Done with PLINK scoring.')
    processed, skipped = parse_log(out_prefix.with_suffix('.log'))
    stats = summarize_sscore(out_prefix.with_suffix('.sscore'), args.target_iid)

    recovery = processed / original_variants if original_variants and processed else None
    report = {
        'pgs_id': pgs_id or '',
        'trait_name': trait_name,
        'pfile_prefix': str(pfile_prefix),
        'original_score_variants': original_variants,
        'processed_variants': processed,
        'recovery_fraction': round(recovery, 4) if recovery is not None else '',
        'target_iid': args.target_iid,
        'target_percentile': stats['percentile'],
        'target_quartile': stats['quartile'],
        'outlier_q1_q4': stats['outlier_q1_q4'],
        'target_z': stats['z'],
        'target_score': stats['target'],
        'technical_confidence': technical_confidence(original_variants, processed),
        'overall_confidence': overall_confidence(original_variants, processed),
        'prep_mode': prep_mode,
        'skipped_missing_id': skipped,
    }

    png_path = outdir / f'{base_name}.hist.png'
    pdf_path = outdir / f'{base_name}.hist.pdf'
    title = f"{trait_name}\n{base_name} | pct {stats['percentile']:.2f} | z {stats['z']:.3f}"
    make_plot(stats['refs'], stats['target'], stats['fake_min'], stats['fake_max'], title, png_path, pdf_path=None)
    report['histogram_png'] = str(png_path)
    report['histogram_pdf'] = str(pdf_path)

    json_path = outdir / f'{base_name}.report.json'
    txt_path = outdir / f'{base_name}.report.txt'
    # json_path.write_text(json.dumps(report, indent=2) + "\n")
    write_report_txt(txt_path, report)

    print(f"score_tsv\t{score_path}")
    print(f"sscore\t{out_prefix.with_suffix('.sscore')}")
    print(f"report_txt\t{txt_path}")
    # print(f"report_json\t{json_path}")
    # print(f"hist_png\t{png_path}")
    # print(f"hist_pdf\t{pdf_path}")
    print(f"Target Percentile: {report['target_percentile']}")


if __name__ == '__main__':
    main()
