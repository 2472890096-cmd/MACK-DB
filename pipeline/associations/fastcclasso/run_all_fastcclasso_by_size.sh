#!/usr/bin/env bash
set -Eeuo pipefail

BASE="${BASE:-/mnt/new_data/all_fastccl_species_count_GTDB226}"
ABUNDANCE_DIR="${ABUNDANCE_DIR:-/home/wangmengyang/project/GDTB226_species_abundance_result/abundance}"
SPECIES_DIR="${SPECIES_DIR:-/home/wangmengyang/project/GDTB226_species_abundance_result/species}"
RUNNER="${RUNNER:-/home/wangmengyang/project/fastcclasso_test/run_fastcclasso_parallel.R}"
FAST_CORE="${FAST_CORE:-/home/wangmengyang/project/fastcclasso_test/official_fastCCLasso/RCode/fastCCLasso.R}"

WORK_DIR="${WORK_DIR:-/home/wangmengyang/project/fastcclasso_test/all_project_lists}"
LOG_ROOT="${LOG_ROOT:-${BASE}/batch_logs}"

NBOOT="${NBOOT:-100}"
KMAX="${KMAX:-20}"
PVALUE="${PVALUE:-0.0001}"
COR_CUTOFF="${COR_CUTOFF:-0.1}"
WRITE_FULL_MATRICES="${WRITE_FULL_MATRICES:-T}"
WRITE_PAIR_CSV="${WRITE_PAIR_CSV:-T}"
SEED="${SEED:-20260804}"
DRY_RUN="${DRY_RUN:-F}"

WORKERS_LT500="${WORKERS_LT500:-24}"
WORKERS_500_999="${WORKERS_500_999:-24}"
WORKERS_1000_1499="${WORKERS_1000_1499:-12}"
WORKERS_1500_1999="${WORKERS_1500_1999:-6}"
WORKERS_GE2000="${WORKERS_GE2000:-2}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

log() {
  echo "[$(date '+%F %T')] $*"
}

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

for path in "$ABUNDANCE_DIR" "$SPECIES_DIR"; do
  [[ -d "$path" ]] || die "Directory does not exist: $path"
done
for path in "$RUNNER" "$FAST_CORE"; do
  [[ -r "$path" ]] || die "File does not exist or is not readable: $path"
done

mkdir -p "$BASE" "$WORK_DIR" "$LOG_ROOT"

exec 9>"$BASE/run_all_fastcclasso.lock"
if command -v flock >/dev/null 2>&1; then
  flock -n 9 || die "Another fastCCLasso batch appears to be running: $BASE"
fi
echo $$ > "$BASE/run_all_fastcclasso.pid"
trap 'rm -f "$BASE/run_all_fastcclasso.pid"' EXIT
mkdir -p \
  "$BASE/matrix/cor" \
  "$BASE/matrix/pv" \
  "$BASE/cor-pv_tsvformat" \
  "$BASE/cor-pv_csvformat" \
  "$BASE/proj_links" \
  "$BASE/project_logs" \
  "$BASE/project_status" \
  "$BASE/group_summaries" \
  "$BASE/fastCCLasso"

ALL_LIST="$WORK_DIR/all_abundance_projects.list"
COUNT_TABLE="$WORK_DIR/all_project_species_counts.tsv"
MISSING_SPECIES="$WORK_DIR/missing_species_files.tsv"
EMPTY_SPECIES="$WORK_DIR/empty_species_files.tsv"
TOO_FEW_TAXA="$WORK_DIR/projects_fewer_than_2_taxa.tsv"

LT500_LIST="$WORK_DIR/projects_lt500.list"
P500_999_LIST="$WORK_DIR/projects_500_999.list"
P1000_1499_LIST="$WORK_DIR/projects_1000_1499.list"
P1500_1999_LIST="$WORK_DIR/projects_1500_1999.list"
GE2000_LIST="$WORK_DIR/projects_ge2000.list"

find "$ABUNDANCE_DIR" -maxdepth 1 -type f -name '*_abundance.matrix' \
  -printf '%f\n' \
| sed 's/_abundance\.matrix$//' \
| LC_ALL=C sort -u > "$ALL_LIST"

: > "$COUNT_TABLE"
: > "$MISSING_SPECIES"
: > "$EMPTY_SPECIES"
: > "$TOO_FEW_TAXA"
: > "$LT500_LIST"
: > "$P500_999_LIST"
: > "$P1000_1499_LIST"
: > "$P1500_1999_LIST"
: > "$GE2000_LIST"

while IFS= read -r project; do
  [[ -n "$project" ]] || continue
  species_file="$SPECIES_DIR/${project}_species.matrix"

  if [[ ! -f "$species_file" ]]; then
    printf '%s\t%s\n' "$project" "$species_file" >> "$MISSING_SPECIES"
    continue
  fi

  taxa=$(awk 'NF {n++} END {print n+0}' "$species_file")
  printf '%s\t%s\n' "$project" "$taxa" >> "$COUNT_TABLE"

  if (( taxa == 0 )); then
    printf '%s\t%s\n' "$project" "$species_file" >> "$EMPTY_SPECIES"
    continue
  elif (( taxa < 2 )); then
    printf '%s\t%s\n' "$project" "$species_file" >> "$TOO_FEW_TAXA"
    continue
  elif (( taxa < 500 )); then
    echo "$project" >> "$LT500_LIST"
  elif (( taxa < 1000 )); then
    echo "$project" >> "$P500_999_LIST"
  elif (( taxa < 1500 )); then
    echo "$project" >> "$P1000_1499_LIST"
  elif (( taxa < 2000 )); then
    echo "$project" >> "$P1500_1999_LIST"
  else
    echo "$project" >> "$GE2000_LIST"
  fi
done < "$ALL_LIST"

# Run larger projects first within each size group to reduce the long tail.
sort_group_by_taxa_desc() {
  local list_file="$1"
  local tmp_file="${list_file}.tmp"
  awk 'NR==FNR {taxa[$1]=$2; next} NF {print taxa[$1] "\t" $1}' \
    "$COUNT_TABLE" "$list_file" \
  | sort -t $'\t' -k1,1nr -k2,2 \
  | cut -f2- > "$tmp_file"
  mv "$tmp_file" "$list_file"
}

for list_file in \
  "$LT500_LIST" \
  "$P500_999_LIST" \
  "$P1000_1499_LIST" \
  "$P1500_1999_LIST" \
  "$GE2000_LIST"
do
  sort_group_by_taxa_desc "$list_file"
done

count_lines() {
  awk 'NF {n++} END {print n+0}' "$1"
}

log "fastCCLasso batch configuration"
echo "  Output base       : $BASE"
echo "  Total abundance   : $(count_lines "$ALL_LIST")"
echo "  Missing species   : $(count_lines "$MISSING_SPECIES")"
echo "  Empty species     : $(count_lines "$EMPTY_SPECIES")"
echo "  Fewer than 2 taxa : $(count_lines "$TOO_FEW_TAXA")"
echo "  <500 taxa         : $(count_lines "$LT500_LIST") workers=$WORKERS_LT500"
echo "  500-999 taxa      : $(count_lines "$P500_999_LIST") workers=$WORKERS_500_999"
echo "  1000-1499 taxa    : $(count_lines "$P1000_1499_LIST") workers=$WORKERS_1000_1499"
echo "  1500-1999 taxa    : $(count_lines "$P1500_1999_LIST") workers=$WORKERS_1500_1999"
echo "  >=2000 taxa       : $(count_lines "$GE2000_LIST") workers=$WORKERS_GE2000"
echo "  n_boot            : $NBOOT"
echo "  k_max             : $KMAX"
echo "  p-value cutoff    : $PVALUE"
echo "  |cor| cutoff      : $COR_CUTOFF"
echo "  Full matrices     : $WRITE_FULL_MATRICES"
echo "  Pair CSV          : $WRITE_PAIR_CSV"
echo
df -h "$BASE" || true

if [[ "${DRY_RUN^^}" =~ ^(T|TRUE|1|Y|YES)$ ]]; then
  log "DRY_RUN enabled: lists/configuration generated; no project was executed."
  exit 0
fi

if pgrep -af 'run_remaining_cclasso_by_size[.]sh|run_cclasso_original_parallel[.]R' \
  > "$BASE/original_cclasso_processes_before_start.txt"; then
  cat "$BASE/original_cclasso_processes_before_start.txt" >&2
  die "Original CCLasso batch is still running. Stop it before starting fastCCLasso."
fi
rm -f "$BASE/original_cclasso_processes_before_start.txt"

echo -e "group\tprojects\tworkers\tstart_time\tend_time\texit_code" \
  > "$BASE/batch_group_status.tsv"

run_group() {
  local label="$1"
  local list_file="$2"
  local workers="$3"
  local n_projects
  local start_time
  local end_time
  local rc

  n_projects=$(count_lines "$list_file")
  if (( n_projects == 0 )); then
    log "Skip empty group: $label"
    return 0
  fi

  start_time=$(date '+%F %T')
  log "===== Running $label: projects=$n_projects workers=$workers ====="
  log "Main log: $LOG_ROOT/${label}.main.log"

  set +e
  Rscript "$RUNNER" \
    --file "$list_file" \
    --abundance_dir "$ABUNDANCE_DIR" \
    --species_dir "$SPECIES_DIR" \
    --output "$BASE" \
    --fast_core "$FAST_CORE" \
    --run_label "$label" \
    --nboot "$NBOOT" \
    --kmax "$KMAX" \
    --pvalue "$PVALUE" \
    --cor_cutoff "$COR_CUTOFF" \
    --workers "$workers" \
    --seed "$SEED" \
    --skip_existing T \
    --write_full_matrices "$WRITE_FULL_MATRICES" \
    --write_pair_csv "$WRITE_PAIR_CSV" \
    > "$LOG_ROOT/${label}.main.log" 2>&1
  rc=$?
  set -e

  end_time=$(date '+%F %T')
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$label" "$n_projects" "$workers" "$start_time" "$end_time" "$rc" \
    >> "$BASE/batch_group_status.tsv"

  if (( rc != 0 )); then
    log "[ERROR] group $label exited with code $rc; stopping batch"
    return "$rc"
  else
    log "[finished group] $label"
  fi
}

run_group "lt500_w${WORKERS_LT500}" "$LT500_LIST" "$WORKERS_LT500"
run_group "500_999_w${WORKERS_500_999}" "$P500_999_LIST" "$WORKERS_500_999"
run_group "1000_1499_w${WORKERS_1000_1499}" "$P1000_1499_LIST" "$WORKERS_1000_1499"
run_group "1500_1999_w${WORKERS_1500_1999}" "$P1500_1999_LIST" "$WORKERS_1500_1999"
run_group "ge2000_w${WORKERS_GE2000}" "$GE2000_LIST" "$WORKERS_GE2000"

summary_files=("$BASE"/group_summaries/project_timing_and_status_*.tsv)
if (( ${#summary_files[@]} > 0 )) && [[ -e "${summary_files[0]}" ]]; then
  awk 'FNR==1 && NR!=1 {next} {print}' "${summary_files[@]}" \
    > "$BASE/project_timing_and_status_all.tsv"
fi

find "$BASE/cor-pv_tsvformat" -maxdepth 1 -type f \
  -name '*_fastCCLasso_freq_sp_pairs.tsv' \
  -printf '%f\n' \
| sed 's/_fastCCLasso_freq_sp_pairs\.tsv$//' \
| LC_ALL=C sort -u > "$BASE/completed_projects.list"

comm -23 "$ALL_LIST" "$BASE/completed_projects.list" \
  > "$BASE/not_completed_projects.list"

log "All groups finished."
echo "Completed projects : $(count_lines "$BASE/completed_projects.list")"
echo "Not completed      : $(count_lines "$BASE/not_completed_projects.list")"
echo "Summary            : $BASE/project_timing_and_status_all.tsv"
echo "Group status       : $BASE/batch_group_status.tsv"
echo "Main output        : $BASE"
