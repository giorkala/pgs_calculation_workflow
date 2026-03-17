#!/usr/bin/env python3
import csv
import math
import pathlib
import re

import matplotlib.pyplot as plt

PROJECT = pathlib.Path('/nfs/users/nfs_g/gk18/prs_application')
SUMMARY = PROJECT / 'results/plink2_batch/summary_broad.tsv'
RESULTS_DIR = PROJECT / 'results/plink2_batch'
OUT_PNG = PROJECT / 'results/plink2_batch/prs_histograms_grid.png'
OUT_PDF = PROJECT / 'results/plink2_batch/prs_histograms_grid.pdf'
TARGET_IID = 'genome1402_genome1402'
FAKE_RE = re.compile(r'^FAKE_[0-9]+_FAKE_[0-9]+$')


def read_summary(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def load_scores(pgs_id):
    path = RESULTS_DIR / f'{pgs_id}.sscore'
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    score_col = 'SCORE1_AVG' if 'SCORE1_AVG' in rows[0] else 'SCORE1_SUM'
    refs = [float(r[score_col]) for r in rows if r['IID'] != TARGET_IID and not FAKE_RE.match(r['IID'])]
    target = next(float(r[score_col]) for r in rows if r['IID'] == TARGET_IID)
    fakes = [float(r[score_col]) for r in rows if FAKE_RE.match(r['IID'])]
    return refs, target, fakes


def main():
    summary = read_summary(SUMMARY)
    n = len(summary)
    ncols = 4
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.4 * nrows), constrained_layout=True)
    axes = axes.flatten()

    for ax, row in zip(axes, summary):
        refs, target, fakes = load_scores(row['PGS_ID'])
        ax.hist(refs, bins=30, color='#b8c9d9', edgecolor='#4a6572', linewidth=0.6)
        if fakes:
            ax.axvspan(min(fakes), max(fakes), color='#e74c3c', alpha=0.3)
        ax.axvline(target, color='#c0392b', linewidth=2.0)
        title_color = '#b03a2e' if row['TARGET_QUARTILE'] == 'Q1' else '#1f618d' if row['TARGET_QUARTILE'] == 'Q4' else '#222222'
        ax.set_title(f"{row['TRAIT']}\n{row['PGS_ID']} | pct {row['TARGET_PERCENTILE']} | z {row['TARGET_Z']}", fontsize=10, color=title_color)
        ax.set_xlabel('PRS')
        ax.set_ylabel('EUR count')
        ax.grid(alpha=0.2, linewidth=0.5)

    for ax in axes[n:]:
        ax.axis('off')

    fig.suptitle('PRS Distributions in 1000G EUR with genome1402 highlighted', fontsize=16)
    fig.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
    fig.savefig(OUT_PDF, bbox_inches='tight')
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == '__main__':
    main()
