#!/usr/bin/env bash
set -euo pipefail

# Merge imputed Michigan output with 1000G EUR reference samples.
# Assumes the following already exist:
#   - merge_1000g/1000g_eur.snpidfix.{pgen,pvar.zst,psam}
#   - merge_1000g/common.snplist
#   - merge_1000g/imputed_all.concat.dose.vcf.gz
#
# Run from the project root:
#   bash scripts/merge_imputed_with_1000g.sh

cd "$(dirname "$0")/.."

module load cellgen/plink/2.00

# 1) Subset the EUR reference to the shared SNP set.
plink2 \
  --pfile merge_1000g/1000g_eur.snpidfix vzs \
  --extract merge_1000g/common.snplist \
  --make-pgen vzs \
  --out merge_1000g/1000g_eur.common

# 2) Import the concatenated imputed VCF with normalized CHR:POS:REF:ALT IDs,
#    then subset to the same shared SNP set.
plink2 \
  --vcf merge_1000g/imputed_all.concat.dose.vcf.gz dosage=DS \
  --double-id \
  --snps-only just-acgt \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 200 missing \
  --extract merge_1000g/common.snplist \
  --make-pgen vzs \
  --out merge_1000g/imputed_all.final

# 3) Export both sides to bgzipped VCF because this PLINK build does not support
#    sample-wise pmerge for these inputs.
plink2 \
  --pfile merge_1000g/1000g_eur.common vzs \
  --export vcf bgz \
  --out merge_1000g/1000g_eur.common

plink2 \
  --pfile merge_1000g/imputed_all.final vzs \
  --export vcf bgz \
  --out merge_1000g/imputed_all.final

module load cellgen/bcftools/1.21

# 4) Index the VCFs and merge samples across the shared variant set.
bcftools index -f merge_1000g/1000g_eur.common.vcf.gz
bcftools index -f merge_1000g/imputed_all.final.vcf.gz
bcftools merge \
  -m none \
  -Oz \
  --threads 8 \
  -o merge_1000g/eur_plus_imputed.vcf.gz \
  merge_1000g/1000g_eur.common.vcf.gz \
  merge_1000g/imputed_all.final.vcf.gz

tabix -f -p vcf merge_1000g/eur_plus_imputed.vcf.gz

# 5) Quick sanity check: should be 503 EUR + 11 imputed = 514 samples.
bcftools query -l merge_1000g/eur_plus_imputed.vcf.gz | wc -l
# another sanity check, number of variants should be ~30.16M
bcftools view -HG merge_1000g/eur_plus_imputed.vcf.gz | wc -l
