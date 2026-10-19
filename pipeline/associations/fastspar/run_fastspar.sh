#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 2 ]]; then
  echo "Usage: bash run_fastspar.sh COUNT_MATRIX OUTPUT_PREFIX" >&2
  exit 2
fi
PROFILE="$(realpath -- "$1")"
PREFIX="$2"
T="${MACKDB_THREADS:-1}"
FASTSPAR_BIN="${FASTSPAR_BIN:-fastspar}"
FASTSPAR_BOOTSTRAP_BIN="${FASTSPAR_BOOTSTRAP_BIN:-fastspar_bootstrap}"
FASTSPAR_PVALUES_BIN="${FASTSPAR_PVALUES_BIN:-fastspar_pvalues}"
PARALLEL_BIN="${PARALLEL_BIN:-parallel}"
for exe in "$FASTSPAR_BIN" "$FASTSPAR_BOOTSTRAP_BIN" "$FASTSPAR_PVALUES_BIN" "$PARALLEL_BIN"; do
  command -v "$exe" >/dev/null || { echo "Missing executable: $exe" >&2; exit 127; }
done
mkdir -p -- "$(dirname -- "$PREFIX")"
CORR="${PREFIX}_correlation.tsv"
COVAR="${PREFIX}_covariance.tsv"
BOOSTSTARP_COUNTS="${PREFIX}_bootstrap_counts"
BOOSTSTARP_CORR="${PREFIX}_bootstrap_correlation"
PVALUE="${PREFIX}_pvalues.tsv"
"$FASTSPAR_BIN" --iterations 50 --exclude_iterations 20 --threshold 0.2 --otu_table "$PROFILE" --correlation "$CORR" --covariance "$COVAR" --threads "$T"
mkdir -- "$BOOSTSTARP_COUNTS"
"$FASTSPAR_BOOTSTRAP_BIN" --otu_table "$PROFILE" --number 10000 --prefix "$BOOSTSTARP_COUNTS/fake_data" --threads "$T"
mkdir -- "$BOOSTSTARP_CORR"
"$PARALLEL_BIN" --jobs "$T" "$FASTSPAR_BIN" --otu_table {} --correlation "$BOOSTSTARP_CORR/cor_{/}" --covariance "$BOOSTSTARP_CORR/cov_{/}" -i 5 ::: "$BOOSTSTARP_COUNTS"/*
"$FASTSPAR_PVALUES_BIN" --otu_table "$PROFILE" --correlation "$CORR" --prefix "$BOOSTSTARP_CORR/cor_fake_data_" --permutations 10000 --outfile "$PVALUE" --threads "$T"
rm -rf -- "$BOOSTSTARP_COUNTS" "$BOOSTSTARP_CORR"
