#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash scripts/run_prs_pipeline.sh <merged_genotypes_prefix> <pgs_id_file>
# Example:
#   bash scripts/run_prs_pipeline.sh data/merged_eur_11samples config/pgs_ids.tsv

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <merged_plink_prefix> <pgs_id_file>"
  exit 1
fi

MERGED_PREFIX="$1"
PGS_ID_FILE="$2"

mkdir -p scores/raw scores/plink2 results logs

# Step 1: fetch scores and convert to plink2 format
python3 scripts/download_pgs_weights.py \
  --ids "${PGS_ID_FILE}" \
  --outdir scores

# Step 2: run PLINK2 score for each PGS
while IFS=$'\t' read -r PGS_ID _; do
  [[ -z "${PGS_ID}" ]] && continue
  [[ "${PGS_ID}" =~ ^# ]] && continue

  SCORE_FILE="scores/plink2/${PGS_ID}.score.tsv"
  if [[ ! -s "${SCORE_FILE}" ]]; then
    echo "[WARN] Missing score file for ${PGS_ID}, skipping" | tee -a logs/prs.log
    continue
  fi

  OUT="results/${PGS_ID}"
  echo "[INFO] Scoring ${PGS_ID}" | tee -a logs/prs.log
  plink2 \
    --pfile "${MERGED_PREFIX}" \
    --score "${SCORE_FILE}" 1 2 3 header-read cols=+scoresums \
    --out "${OUT}" | tee -a logs/prs.log

done < "${PGS_ID_FILE}"

echo "[DONE] PRS outputs are in results/*.sscore"
