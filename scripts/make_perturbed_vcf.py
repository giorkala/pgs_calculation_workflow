#!/usr/bin/env python3
import argparse
import gzip
import pathlib
import random
import shutil
import subprocess
import sys

import pysam


def open_maybe_gzip(path: str, mode: str):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def mutate_gt(gt: str, max_allele_index: int, mutation_rate: float, rng: random.Random) -> str:
    if gt in {'.', './.', '.|.'}:
        return gt

    sep = '/'
    if '|' in gt:
        sep = '|'

    alleles = gt.replace('|', '/').split('/')
    mutated = []

    for allele in alleles:
        if allele == '.':
            mutated.append(allele)
            continue

        if rng.random() < mutation_rate:
            try:
                current = int(allele)
            except ValueError:
                mutated.append(allele)
                continue

            if max_allele_index <= 0:
                mutated.append(allele)
                continue

            choices = [str(i) for i in range(max_allele_index + 1) if i != current]
            if choices:
                mutated.append(rng.choice(choices))
            else:
                mutated.append(allele)
        else:
            mutated.append(allele)

    return sep.join(mutated)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Create perturbed genotype copies and optionally prepare Michigan-imputation VCFs.'
    )
    parser.add_argument('-i', '--input', required=True, help='Input VCF/VCF.GZ')
    parser.add_argument('-o', '--output', help='Output VCF/VCF.GZ for perturbation mode')
    parser.add_argument('-n', '--copies', type=int, default=10, help='Number of synthetic copies (default: 10)')
    parser.add_argument('-r', '--mutation-rate', type=float, default=0.01, help='Per-allele mutation probability (default: 0.01)')
    parser.add_argument('--seed', type=int, default=1402, help='Random seed (default: 1402)')
    parser.add_argument('--sample-index', type=int, default=0, help='0-based sample column index to clone (default: 0)')
    parser.add_argument('--prefix', default='SYN', help='Prefix for synthetic sample IDs (default: SYN)')
    parser.add_argument(
        '--michigan-prep',
        action='store_true',
        help='Prepare Michigan-ready bgzip/indexed files: one all-chrom VCF plus chr1..chr22 VCFs',
    )
    parser.add_argument(
        '--from-23andme-raw',
        action='store_true',
        help='Convert 23andMe raw txt/txt.gz to a clean VCF using reference FASTA (drops hom-ref/unsupported sites)',
    )
    parser.add_argument('--ref-fasta', default=None, help='Reference FASTA path (required with --from-23andme-raw)')
    parser.add_argument(
        '--panel-pvar',
        default=None,
        help='Optional PLINK2 .pvar/.pvar.zst with REF/ALT (e.g., 1000G b37) to retain hom-ref calls',
    )
    parser.add_argument('--sample-id', default='SAMPLE', help='Sample ID used in VCF when converting 23andMe raw')
    parser.add_argument('--michigan-outdir', default='michigan_input', help='Output directory for Michigan files')
    parser.add_argument(
        '--michigan-prefix',
        default=None,
        help='Prefix for Michigan output files (default: input basename without .vcf/.vcf.gz)',
    )
    parser.add_argument('--keep-plain', action='store_true', help='Keep intermediate plain chr*.vcf files')
    return parser.parse_args()


def bgzip_from_any(in_path: pathlib.Path, out_bgz: pathlib.Path):
    with open_maybe_gzip(str(in_path), 'rb') as fin, pysam.BGZFile(str(out_bgz), 'wb') as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)


def tabix_index_vcf(vcf_gz: pathlib.Path):
    pysam.tabix_index(str(vcf_gz), preset='vcf', force=True)


def infer_base_prefix(input_path: str) -> str:
    name = pathlib.Path(input_path).name
    if name.endswith('.vcf.gz'):
        return name[:-7]
    if name.endswith('.vcf'):
        return name[:-4]
    return pathlib.Path(input_path).stem


def run_michigan_prep(args):
    outdir = pathlib.Path(args.michigan_outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    base_prefix = args.michigan_prefix or infer_base_prefix(args.input)

    # 1) Full multi-chromosome file, bgzip + tbi index.
    allchr_bgz = outdir / f'{base_prefix}.allchr.vcf.gz'
    bgzip_from_any(pathlib.Path(args.input), allchr_bgz)
    tabix_index_vcf(allchr_bgz)
    sys.stderr.write(f'[OK] Wrote {allchr_bgz} and {allchr_bgz}.tbi\n')

    # 2) Autosome-specific files (chr1..chr22), then bgzip + tbi index each.
    plain_paths = {}
    writers = {}
    header_lines = []
    records_seen = 0

    for i in range(1, 23):
        plain = outdir / f'{base_prefix}.chr{i}.vcf'
        plain_paths[i] = plain
        writers[i] = plain.open('w')

    with open_maybe_gzip(args.input, 'rt') as fin:
        for line in fin:
            if line.startswith('##'):
                header_lines.append(line)
                continue

            if line.startswith('#CHROM'):
                header_lines.append(line)
                for i in range(1, 23):
                    for h in header_lines:
                        writers[i].write(h)
                continue

            records_seen += 1
            chrom = line.split('\t', 1)[0]
            if chrom.startswith('chr'):
                chrom = chrom[3:]

            if chrom.isdigit():
                c = int(chrom)
                if 1 <= c <= 22:
                    writers[c].write(line)

    for fh in writers.values():
        fh.close()

    if records_seen == 0:
        raise ValueError('No variant records found in input VCF.')

    for i in range(1, 23):
        plain = plain_paths[i]
        bgz = outdir / f'{base_prefix}.chr{i}.vcf.gz'
        bgzip_from_any(plain, bgz)
        tabix_index_vcf(bgz)
        if not args.keep_plain:
            plain.unlink(missing_ok=True)
        sys.stderr.write(f'[OK] Wrote {bgz} and {bgz}.tbi\n')


def run_perturbation(args):
    if not args.output:
        raise ValueError('--output is required unless using --michigan-prep')

    rng = random.Random(args.seed)

    with open_maybe_gzip(args.input, 'rt') as fin, open_maybe_gzip(args.output, 'wt') as fout:
        source_sample_name = None

        for line in fin:
            if line.startswith('##'):
                fout.write(line)
                continue

            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                if len(header) < 10:
                    raise ValueError('Input VCF has no sample columns.')

                sample_names = header[9:]
                if args.sample_index < 0 or args.sample_index >= len(sample_names):
                    raise ValueError(f'sample-index {args.sample_index} out of range for {len(sample_names)} samples.')

                source_sample_name = sample_names[args.sample_index]
                synthetic_names = [f"{args.prefix}_{i+1}" for i in range(args.copies)]

                new_header = header[:9] + [source_sample_name] + synthetic_names
                fout.write('\t'.join(new_header) + '\n')
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 10:
                fout.write(line)
                continue

            ref = fields[3]
            alts = fields[4]
            alt_list = [] if alts == '.' else alts.split(',')
            max_allele_index = len(alt_list)

            fmt = fields[8].split(':')
            try:
                gt_idx = fmt.index('GT')
            except ValueError:
                raise ValueError('FORMAT column does not contain GT for at least one variant.')

            sample_values = fields[9:]
            src_sample = sample_values[args.sample_index]
            src_parts = src_sample.split(':')
            src_gt = src_parts[gt_idx] if gt_idx < len(src_parts) else '.'

            synth_gts = [mutate_gt(src_gt, max_allele_index, args.mutation_rate, rng) for _ in range(args.copies)]

            out_fields = fields[:8] + ['GT', src_gt] + synth_gts
            fout.write('\t'.join(out_fields) + '\n')

        if source_sample_name is None:
            raise ValueError('No VCF header (#CHROM) found in input.')

    sys.stderr.write(
        f"Wrote {args.output} with source sample '{source_sample_name}' and {args.copies} perturbed copies at mutation rate {args.mutation_rate}.\n"
    )


def normalize_chrom_for_fasta(chrom: str) -> str:
    return chrom.strip()


def build_fasta_chrom_map(fasta: pysam.FastaFile):
    refs = list(fasta.references)
    refset = set(refs)
    mapping = {}
    for i in range(1, 23):
        raw = str(i)
        if raw in refset:
            mapping[raw] = raw
        elif f'chr{raw}' in refset:
            mapping[raw] = f'chr{raw}'
        elif f'Chr{raw}' in refset:
            mapping[raw] = f'Chr{raw}'
    return mapping


def is_acgt_base(a: str) -> bool:
    return len(a) == 1 and a in 'ACGT'


def normalize_panel_chrom(chrom: str) -> str:
    c = chrom.strip()
    if c.startswith('chr'):
        c = c[3:]
    if c.startswith('Chr'):
        c = c[3:]
    return c


def load_23andme_snp_keys(raw_path: str):
    # Build a lookup set from raw autosomal SNP-like rows to limit pvar scanning memory.
    keys = set()
    with open_maybe_gzip(raw_path, 'rt') as fin:
        for line in fin:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            _rsid, chrom, pos_s, gt = parts[0], parts[1], parts[2], parts[3]
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            if len(gt) != 2 or any(a not in 'ACGT' for a in gt):
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            keys.add((chrom, pos))
    return keys


def load_panel_candidates(panel_pvar: str, needed_keys):
    # Map (chrom,pos) -> list[(id,ref,alt)] for biallelic SNPs only.
    by_pos = {}
    cmd = ['zstdcat', panel_pvar] if panel_pvar.endswith('.zst') else ['cat', panel_pvar]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    if proc.stdout is None:
        raise RuntimeError('Failed to open panel pvar stream')

    for line in proc.stdout:
        if not line or line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 5:
            continue
        chrom = normalize_panel_chrom(parts[0])
        if chrom not in {str(i) for i in range(1, 23)}:
            continue
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        key = (chrom, pos)
        if key not in needed_keys:
            continue
        var_id = parts[2]
        ref = parts[3].upper()
        alt = parts[4].upper()
        if ',' in alt:
            continue
        if not is_acgt_base(ref) or not is_acgt_base(alt):
            continue
        by_pos.setdefault(key, []).append((var_id, ref, alt))

    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f'Error streaming panel pvar (exit {rc})')
    return by_pos


def choose_candidate(rsid: str, gt: str, cands):
    a1, a2 = gt[0], gt[1]
    if not cands:
        return None

    # Priority 1: exact rsID match and genotype-compatible.
    for cid, ref, alt in cands:
        if cid == rsid and a1 in {ref, alt} and a2 in {ref, alt}:
            return cid, ref, alt
    # Priority 2: genotype-compatible candidate.
    for cid, ref, alt in cands:
        if a1 in {ref, alt} and a2 in {ref, alt}:
            return cid, ref, alt
    # Priority 3: rsID match (even if genotype not compatible, usually strand/array issue).
    for cid, ref, alt in cands:
        if cid == rsid:
            return cid, ref, alt
    return None


def run_from_23andme_raw(args):
    if not args.output:
        raise ValueError('--output is required with --from-23andme-raw')
    if not args.ref_fasta:
        raise ValueError('--ref-fasta is required with --from-23andme-raw')

    fasta = pysam.FastaFile(args.ref_fasta)
    out_path = pathlib.Path(args.output)
    chrom_map = build_fasta_chrom_map(fasta)
    panel_by_pos = None
    if args.panel_pvar:
        sys.stderr.write('[INFO] Collecting target loci from 23andMe raw file...\n')
        needed = load_23andme_snp_keys(args.input)
        sys.stderr.write(f'[INFO] Raw autosomal SNP-like loci: {len(needed)}\n')
        sys.stderr.write(f'[INFO] Loading panel variants from {args.panel_pvar}...\n')
        panel_by_pos = load_panel_candidates(args.panel_pvar, needed)
        sys.stderr.write(f'[INFO] Panel loci matched: {len(panel_by_pos)}\n')

    kept = 0
    skipped_non_autosome = 0
    skipped_bad_gt = 0
    skipped_ref_missing = 0
    skipped_hom_ref = 0
    skipped_multialt = 0
    skipped_panel_miss = 0
    skipped_panel_mismatch = 0

    with open_maybe_gzip(args.input, 'rt') as fin, pysam.BGZFile(str(out_path), 'w') as fout:
        # Header
        fout.write(b'##fileformat=VCFv4.2\n')
        fout.write(b'##source=make_perturbed_vcf.py:from-23andme-raw\n')
        for i in range(1, 23):
            try:
                ln = fasta.get_reference_length(str(i))
            except Exception:
                ln = 0
            if ln > 0:
                fout.write(f'##contig=<ID={i},length={ln}>\n'.encode())
        fout.write(b'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fout.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample_id}\n'.encode())

        for line in fin:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue

            rsid, chrom, pos_s, gt = parts[0], parts[1], parts[2], parts[3]
            # Michigan input here should be autosomal only.
            if chrom not in {str(i) for i in range(1, 23)}:
                skipped_non_autosome += 1
                continue

            if len(gt) != 2 or any(a not in 'ACGT' for a in gt):
                skipped_bad_gt += 1
                continue

            try:
                pos = int(pos_s)
            except ValueError:
                skipped_bad_gt += 1
                continue

            fasta_chrom = chrom_map.get(chrom, normalize_chrom_for_fasta(chrom))
            try:
                ref = fasta.fetch(fasta_chrom, pos - 1, pos).upper()
            except Exception:
                skipped_ref_missing += 1
                continue

            if len(ref) != 1 or ref not in 'ACGT':
                skipped_ref_missing += 1
                continue

            a1, a2 = gt[0], gt[1]
            if panel_by_pos is not None:
                cands = panel_by_pos.get((chrom, pos), [])
                chosen = choose_candidate(rsid, gt, cands)
                if chosen is None:
                    skipped_panel_miss += 1
                    continue
                out_id, pref, alt = chosen
                allele_to_gt = {pref: '0', alt: '1'}
                gt_num = f"{allele_to_gt.get(a1, '.')}/{allele_to_gt.get(a2, '.')}"
                if '.' in gt_num:
                    skipped_panel_mismatch += 1
                    continue
                out_ref = pref
            else:
                non_ref = sorted({a for a in (a1, a2) if a != ref})
                if len(non_ref) == 0:
                    # Cannot infer ALT for hom-ref from raw 23andMe alone; skip.
                    skipped_hom_ref += 1
                    continue
                if len(non_ref) > 1:
                    # Would be multi-allelic in a single sample (rare for SNPs); skip for clean biallelic VCF.
                    skipped_multialt += 1
                    continue

                alt = non_ref[0]
                allele_to_gt = {ref: '0', alt: '1'}
                gt_num = f"{allele_to_gt.get(a1, '.')}/{allele_to_gt.get(a2, '.')}"
                if '.' in gt_num:
                    skipped_bad_gt += 1
                    continue
                out_id = rsid
                out_ref = ref

            rec = f'{chrom}\t{pos}\t{out_id}\t{out_ref}\t{alt}\t.\tPASS\t.\tGT\t{gt_num}\n'
            fout.write(rec.encode())
            kept += 1

    tabix_index_vcf(out_path)
    sys.stderr.write(
        f'[OK] Wrote {out_path} and {out_path}.tbi\n'
        f'[STATS] kept={kept} skipped_non_autosome={skipped_non_autosome} skipped_bad_gt={skipped_bad_gt} '
        f'skipped_ref_missing={skipped_ref_missing} skipped_hom_ref={skipped_hom_ref} skipped_multialt={skipped_multialt} '
        f'skipped_panel_miss={skipped_panel_miss} skipped_panel_mismatch={skipped_panel_mismatch}\n'
    )


def main():
    args = parse_args()
    if args.from_23andme_raw:
        run_from_23andme_raw(args)
    elif args.michigan_prep:
        run_michigan_prep(args)
    else:
        run_perturbation(args)


if __name__ == '__main__':
    main()
