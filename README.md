# PRS Application Workflow

Simple pipeline for calculating polygenic risk scores (PGS; aka PRS) for an individual's genotypes, using variant weights from the PGS Catalog. To provide context, it first merges the user's genotypes with those from the 1000 Genomes Project (1kGP) and calculates PGS for all samples.

It covers:
- converting raw 23andMe text into a valid VCF using 1kGP Phase 3 alleles
- creating synthetic copies to satisfy Michigan Imputation Server minimum sample count (also serving as probes for sensitivity analysis), saving to chromosome-based VCFs
- merging imputed data with 1kGP EUR reference samples
- fetching or supplying PGS weights
- scoring PRS against the merged EUR reference cohort
- generating summary tables and visualizations

The current working sample is `user123` and the merged reference cohort contains:
- `503` 1kGP EUR samples
- `1` real target sample
- `10` synthetic perturbed copies

## Project Layout

Important files and directories:
- `user123.ttam.txt.gz`: raw 23andMe text file
- `merge_1kgp/`: merged reference/imputed cohort files
- `scores/`: raw and converted score files
- `results/`: PRS outputs, summaries, figures
- `scripts/make_perturbed_vcf.py`: conversion and perturbation utility
- `scripts/merge_imputed_with_1kgp.sh`: merge helper
- `scripts/run_pgs_batch.py`: batch PGS fetch/score/QC
- `scripts/run_single_pgs_report.py`: single-trait wrapper
- `scripts/plot_prs_histograms.py`: histogram grid plotting

## Requirements

Expected environment:
- Linux shell
- Python 3
- `plink2`
- `bcftools`
- `tabix`
- `zstdcat`
- Python package `matplotlib` (only for last step)

## 1. Start From Raw 23andMe Data

The raw 23andMe file looks like:

```text
# rsid  chromosome  position  genotype
rs548049170  1  69869   TT
rs9326622    1  567092  CC
rs116587930  1  727841  AG
```

Important limitation:
- raw 23andMe genotypes do not tell you the alternate allele for hom-ref sites
- therefore a FASTA alone is not enough to build a useful VCF
- we use the 1kGP Phase 3 panel to assign valid `REF/ALT`

Build a clean autosomal VCF anchored to the Phase 3 panel:

```bash
python3 scripts/make_perturbed_vcf.py \
  --from-23andme-raw \
  -i user123.ttam.txt.gz \
  -o user123.raw_autosomal_phase3panel.vcf.gz \
  --sample-id user123 \
  --ref-fasta FIXTHIS/Homo_sapiens.GRCh37.dna.all.fa \
  --panel-pvar FIXTHIS/1kGPenomes/b37/all_phase3.pvar.zst
```

Output:
- `user123.raw_autosomal_phase3panel.vcf.gz`
- `user123.raw_autosomal_phase3panel.vcf.gz.tbi`

## 2. Create Synthetic Copies For Michigan

Michigan requires at least five samples. We satisfy that by creating perturbed copies of the target sample.

Generate `10` synthetic copies with `1%` random allele perturbation:

```bash
python3 scripts/make_perturbed_vcf.py \
  -i user123.raw_autosomal_phase3panel.vcf.gz \
  -o user123.raw_autosomal_phase3panel.plus10.vcf.gz \
  -n 10 \
  -r 0.01 \
  --seed 1402 \
  --sample-index 0 \
  --prefix FAKE
```

Output:
- `user123.raw_autosomal_phase3panel.plus10.vcf.gz`

## 3. Make Michigan-Ready Chromosome Files

Split into chromosome-specific bgzipped VCFs plus an all-chromosome file:

```bash
python3 scripts/make_perturbed_vcf.py \
  -i user123.raw_autosomal_phase3panel.plus10.vcf.gz \
  --michigan-prep \
  --michigan-outdir michigan_input_phase3panel_plus10 \
  --michigan-prefix user123.raw_autosomal_phase3panel.plus10
```

Outputs go to:
- `michigan_input_phase3panel_plus10/`

If you want one compressed bundle for easier handling:

```bash
tar -czf michigan_input_phase3panel_plus10.tar.gz michigan_input_phase3panel_plus10
```

## 4. Upload To Michigan Imputation Server

Use the chromosome VCFs in `michigan_input_phase3panel_plus10/`.

Key assumptions for this project:
- build: `GRCh37`
- upload chromosome-specific `vcf.gz` files
- input contains the real sample and the synthetic copies

After the job completes, download the imputed output files into:
- `imputation/`

We worked with the `chr*.dose.vcf.gz` files.

Notes:
- `chr*.dose.vcf.gz` = full imputed output
- `chr*.empiricalDose.vcf.gz` = typed-only leave-one-out style output
- use the `dose` VCFs for this workflow

## 5. Prepare The 1kGP EUR Reference Cohort

Extract EUR from Phase 3 and normalize IDs.

The final EUR reference dataset used here was (not provided):
- `merge_1kgp/1kgp_eur.snpidfix.{pgen,pvar.zst,psam}`

A critical lesson from this project:
- variant IDs must be normalized consistently
- use quoted `@:#:$r:$a` when setting PLINK IDs
- shell expansion will break unquoted `$r` and `$a`

## 6. Merge Imputed Data With 1kGP EUR

PLINK does not support sample-level merging of different datasets, so we convert the pgen data to VCFs to work with bcftools.

Main helper script:

```bash
bash scripts/merge_imputed_with_1kgp.sh
```

That script does the following:
1. subsets 1kGP EUR to common SNPs
2. imports the concatenated imputed dose VCF
3. exports both datasets as VCF
4. merges samples with `bcftools merge`

Final merged cohort:
- `merge_1kgp/eur_plus_imputed.vcf.gz`
- `merge_1kgp/eur_plus_imputed.vcf.gz.tbi`

Sanity checks from this project:
- samples: `514`
- variants: about `30.16M`

Convert merged VCF to PLINK2:

```bash

plink2 \
  --vcf merge_1kgp/eur_plus_imputed.vcf.gz dosage=DS \
  --double-id \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 200 missing \
  --make-pgen vzs \
  --out merge_1kgp/eur_plus_imputed
```

Outputs:
- `merge_1kgp/eur_plus_imputed.pgen`
- `merge_1kgp/eur_plus_imputed.pvar.zst`
- `merge_1kgp/eur_plus_imputed.psam`

## 7. Score A Single PGS

Use the single-trait wrapper.

Fetch from PGS Catalog by PGS ID:

```bash
python3 scripts/run_single_pgs_report.py \
  --pgs-id PGS002204 \
  --pfile merge_1kgp/eur_plus_imputed \
  --outdir results/single_pgs/PGS002204
```

Use an existing local allele-weight file:

```bash
python3 scripts/run_single_pgs_report.py \
  --score-file /path/to/weights.tsv \
  --trait-name "My trait" \
  --name my_trait \
  --pfile merge_1kgp/eur_plus_imputed \
  --outdir results/single_pgs/my_trait
```

The wrapper:
- fetches or reads the score file
- remaps by `chr`, `pos`, and `effect_allele` onto the merged cohort
- runs `plink2 --score`
- makes a histogram
- writes short text and JSON reports

Outputs per run:
- `<name>.score.tsv`
- `<name>.sscore`
- `<name>.hist.png`
- `<name>.hist.pdf`
- `<name>.report.txt`
- `<name>.report.json`

## 8. Score A Batch Of PGS

Batch scoring script:

```bash
python3 scripts/run_pgs_batch.py --panel pgs_candidate_pool.tsv
```

Important result from this project:
- score preparation must be done by remapping against the merged cohort
- using `other_allele/effect_allele` directly as `REF/ALT` is not reliable
- the remap-first approach rescued multiple scores and removed a systematic mismatch problem

Current broad summary:
- `results/plink2_batch/summary.tsv`

Columns include:
- `ORIGINAL_SCORE_VARIANTS`
- `PROCESSED_VARIANTS`
- `RECOVERY_FRACTION`
- `TARGET_PERCENTILE`
- `TARGET_QUARTILE`
- `OUTLIER_Q1_OR_Q4`
- `TARGET_Z`
- `TECHNICAL_CONFIDENCE`
- `OVERALL_CONFIDENCE`

Current broad panel config:
- `pgs_candidate_pool.tsv`

## 9. Visualization and Interpretation

Make one histogram per trait in a grid:

```bash
python3 scripts/plot_prs_histograms.py --target-iid user123_user123
```

Outputs:
- `results/plink2_batch/prs_histograms_grid.png`
- `results/plink2_batch/prs_histograms_grid.pdf`

Figure features:
- one histogram per trait
- distribution from the `503` EUR reference samples
- red vertical line for `user123_user123`
- translucent red band for the fake-genotype score range
- Q1/Q4 title coloring

### Percentile
Your percentile is relative to the EUR reference cohort used in this project.

Examples:
- `5` percentile: lower than almost all EUR reference samples
- `50` percentile: near the middle
- `90` percentile: higher than most EUR reference samples

### Target z
`TARGET_Z` is:

```text
(your_score - mean_reference_score) / sd_reference_score
```

Interpretation:
- `0`: near the EUR mean
- positive: above the mean
- negative: below the mean
- magnitude shows how unusual the score is relative to the reference

### Confidence
This project uses two confidence notions:
- `TECHNICAL_CONFIDENCE`: overlap and harmonization quality
- `OVERALL_CONFIDENCE`: capped because we are not automatically ingesting published score performance metrics such as external `R2` or AUC

Important:
- PRS is relative ranking, not diagnosis
- direction must be checked against the underlying phenotype definition
- sparse scores and scores with lower recovery fractions should be interpreted more cautiously

## 10. Practical Lessons From This Project

1. Raw 23andMe data must be anchored to a real variant panel, not just a FASTA.
2. Michigan input needs valid `REF/ALT`; hom-ref calls are still informative and should not be dropped.
3. Consistent variant ID construction is critical.
4. For PGS files, remapping by position and effect allele against the final merged cohort is much more reliable than trusting the score file's `other_allele` as `REF`.
5. Sparse scores can still be useful, but they deserve lower confidence.
6. Indels were intentionally dropped in the current merged SNP-only pipeline; recovering them would require a more complex branch.

## 11. Future work

Natural future extensions:
- append single-score runs automatically into a master summary TSV
- add trait metadata and publication links to summaries
- add trait-direction interpretation labels where phenotype coding is unambiguous
- build an indel-inclusive branch of the merge pipeline for older sparse scores

## 12. Acknowledgements
This workflow was developed with help from OpenAI's GPT-5.4 (Codex).
