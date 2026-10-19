#!/usr/bin/env bash
set -Eeuo pipefail

MODE="${MODE:-lt1500}"

BASE="${BASE:-/mnt/new_data/all_fastccl_genus_count_GTDB226}"

ABUNDANCE_DIR="${ABUNDANCE_DIR:-/home/wangmengyang/project/GDTB226_genus_abundance_result/abundance}"
GENUS_DIR="${GENUS_DIR:-/home/wangmengyang/project/GDTB226_genus_abundance_result/genus}"

RUNNER="${RUNNER:-/home/wangmengyang/project/fastcclasso_test/run_fastcclasso_parallel_genus.R}"
FAST_CORE="${FAST_CORE:-/home/wangmengyang/project/fastcclasso_test/official_fastCCLasso/RCode/fastCCLasso.R}"

WORK_DIR="${WORK_DIR:-/home/wangmengyang/project/fastcclasso_test/genus_project_lists}"
LOG_ROOT="${LOG_ROOT:-${BASE}/batch_logs}"

NBOOT="${NBOOT:-100}"
KMAX="${KMAX:-20}"
PVALUE="${PVALUE:-0.0001}"
COR_CUTOFF="${COR_CUTOFF:-0.1}"
WRITE_FULL_MATRICES="${WRITE_FULL_MATRICES:-T}"
WRITE_PAIR_CSV="${WRITE_PAIR_CSV:-T}"
SEED="${SEED:-20260806}"
WORKERS="${WORKERS:-10}"
DRY_RUN="${DRY_RUN:-F}"

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

[[ -d "$ABUNDANCE_DIR" ]] || die "Missing directory: $ABUNDANCE_DIR"
[[ -d "$GENUS_DIR" ]] || die "Missing directory: $GENUS_DIR"
[[ -r "$RUNNER" ]] || die "Missing runner: $RUNNER"
[[ -r "$FAST_CORE" ]] || die "Missing fastCCLasso core: $FAST_CORE"

mkdir -p \
    "$BASE" \
    "$WORK_DIR" \
    "$LOG_ROOT" \
    "$BASE/matrix/cor" \
    "$BASE/matrix/pv" \
    "$BASE/cor-pv_tsvformat" \
    "$BASE/cor-pv_csvformat" \
    "$BASE/proj_links" \
    "$BASE/project_logs" \
    "$BASE/project_status" \
    "$BASE/group_summaries" \
    "$BASE/fastCCLasso"

LOCK_FILE="$BASE/run_genus_fastcclasso.lock"

exec 9>"$LOCK_FILE"

if command -v flock >/dev/null 2>&1; then
    flock -n 9 || die "Another genus fastCCLasso batch is already running."
fi

PID_FILE="$BASE/run_genus_${MODE}.pid"
echo $$ > "$PID_FILE"
trap 'rm -f "$PID_FILE"' EXIT

ALL_LIST="$WORK_DIR/all_genus_abundance_projects.list"
COUNT_TABLE="$WORK_DIR/all_genus_project_counts.tsv"

MISSING_GENUS="$WORK_DIR/missing_genus_files.tsv"
EMPTY_GENUS="$WORK_DIR/empty_genus_files.tsv"
TOO_FEW_TAXA="$WORK_DIR/projects_fewer_than_2_genus.tsv"
TOO_FEW_SAMPLES="$WORK_DIR/projects_fewer_than_6_samples.tsv"

LT1500_LIST="$WORK_DIR/projects_genus_lt1500_samples_lt1500.list"
LT1500_SAMPLE_GE1500="$WORK_DIR/projects_genus_lt1500_samples_ge1500.list"
P1500_1999_LIST="$WORK_DIR/projects_genus_1500_1999.list"
GE2000_LIST="$WORK_DIR/projects_genus_ge2000.list"

find "$ABUNDANCE_DIR" \
    -maxdepth 1 -type f \
    -name '*_abundance.matrix' \
    -printf '%f\n' \
| sed 's/_abundance\.matrix$//' \
| LC_ALL=C sort -u > "$ALL_LIST"

: > "$COUNT_TABLE"
: > "$MISSING_GENUS"
: > "$EMPTY_GENUS"
: > "$TOO_FEW_TAXA"
: > "$TOO_FEW_SAMPLES"
: > "$LT1500_LIST"
: > "$LT1500_SAMPLE_GE1500"
: > "$P1500_1999_LIST"
: > "$GE2000_LIST"

while IFS= read -r project; do
    [[ -n "$project" ]] || continue

    abundance_file="$ABUNDANCE_DIR/${project}_abundance.matrix"
    genus_file="$GENUS_DIR/${project}_genus.matrix"

    if [[ ! -f "$genus_file" ]]; then
        printf '%s\t%s\n' \
            "$project" "$genus_file" >> "$MISSING_GENUS"
        continue
    fi

    taxa=$(
        awk 'NF {n++} END {print n+0}' \
        "$genus_file"
    )

    samples=$(
        awk -F'\t' '
            NR == 1 {
                print NF - 1
                exit
            }
        ' "$abundance_file"
    )

    printf '%s\t%s\t%s\n' \
        "$project" "$taxa" "$samples" >> "$COUNT_TABLE"

    if (( taxa == 0 )); then
        printf '%s\t%s\n' \
            "$project" "$genus_file" >> "$EMPTY_GENUS"

    elif (( taxa < 2 )); then
        printf '%s\t%s\t%s\n' \
            "$project" "$taxa" "$samples" >> "$TOO_FEW_TAXA"

    elif (( samples < 6 )); then
        printf '%s\t%s\t%s\n' \
            "$project" "$taxa" "$samples" >> "$TOO_FEW_SAMPLES"

    elif (( taxa < 1500 && samples < 1500 )); then
        echo "$project" >> "$LT1500_LIST"

    elif (( taxa < 1500 && samples >= 1500 )); then
        echo "$project" >> "$LT1500_SAMPLE_GE1500"

    elif (( taxa < 2000 )); then
        echo "$project" >> "$P1500_1999_LIST"

    else
        echo "$project" >> "$GE2000_LIST"
    fi

done < "$ALL_LIST"

sort_list() {
    local list_file="$1"
    local tmp_file="${list_file}.tmp"

    [[ -s "$list_file" ]] || return 0

    awk '
        NR == FNR {
            taxa[$1] = $2
            samples[$1] = $3
            next
        }

        NF {
            print taxa[$1] "\t" samples[$1] "\t" $1
        }
    ' "$COUNT_TABLE" "$list_file" \
    | sort -t $'\t' -k1,1nr -k2,2nr -k3,3 \
    | cut -f3- > "$tmp_file"

    mv "$tmp_file" "$list_file"
}

sort_list "$LT1500_LIST"
sort_list "$LT1500_SAMPLE_GE1500"
sort_list "$P1500_1999_LIST"
sort_list "$GE2000_LIST"

count_lines() {
    awk 'NF {n++} END {print n+0}' "$1"
}

log "Genus fastCCLasso configuration"
echo "  mode                         : $MODE"
echo "  output                       : $BASE"
echo "  total abundance projects     : $(count_lines "$ALL_LIST")"
echo "  missing genus files          : $(count_lines "$MISSING_GENUS")"
echo "  empty genus files            : $(count_lines "$EMPTY_GENUS")"
echo "  fewer than 2 genus           : $(count_lines "$TOO_FEW_TAXA")"
echo "  fewer than 6 samples         : $(count_lines "$TOO_FEW_SAMPLES")"
echo "  genus<1500, samples<1500     : $(count_lines "$LT1500_LIST")"
echo "  genus<1500, samples>=1500    : $(count_lines "$LT1500_SAMPLE_GE1500")"
echo "  genus 1500-1999              : $(count_lines "$P1500_1999_LIST")"
echo "  genus >=2000                 : $(count_lines "$GE2000_LIST")"
echo "  workers                      : $WORKERS"
echo "  n_boot                       : $NBOOT"
echo "  k_max                        : $KMAX"
echo "  p-value cutoff               : $PVALUE"
echo "  |cor| cutoff                 : $COR_CUTOFF"
echo

case "$MODE" in
    lt1500)
        LIST_FILE="$LT1500_LIST"
        LABEL="genus_lt1500_w${WORKERS}"
        ;;

    1500_1999)
        LIST_FILE="$P1500_1999_LIST"
        LABEL="genus_1500_1999_w${WORKERS}"
        ;;

    ge2000)
        LIST_FILE="$GE2000_LIST"
        LABEL="genus_ge2000_w${WORKERS}"
        ;;

    *)
        die "Unsupported MODE: $MODE"
        ;;
esac

PROJECT_COUNT=$(count_lines "$LIST_FILE")

log "Selected list: $LIST_FILE"
log "Selected projects: $PROJECT_COUNT"

if [[ "${DRY_RUN^^}" =~ ^(T|TRUE|1|Y|YES)$ ]]; then
    log "DRY_RUN enabled; no project was executed."
    exit 0
fi

if (( PROJECT_COUNT == 0 )); then
    log "No project in selected group."
    exit 0
fi

log "===== Running $LABEL: projects=$PROJECT_COUNT workers=$WORKERS ====="
log "Main log: $LOG_ROOT/${LABEL}.main.log"

Rscript "$RUNNER" \
    --file "$LIST_FILE" \
    --abundance_dir "$ABUNDANCE_DIR" \
    --species_dir "$GENUS_DIR" \
    --output "$BASE" \
    --fast_core "$FAST_CORE" \
    --run_label "$LABEL" \
    --nboot "$NBOOT" \
    --kmax "$KMAX" \
    --pvalue "$PVALUE" \
    --cor_cutoff "$COR_CUTOFF" \
    --workers "$WORKERS" \
    --seed "$SEED" \
    --skip_existing T \
    --write_full_matrices "$WRITE_FULL_MATRICES" \
    --write_pair_csv "$WRITE_PAIR_CSV" \
    > "$LOG_ROOT/${LABEL}.main.log" 2>&1

log "[finished group] $LABEL"

find "$BASE/cor-pv_tsvformat" \
    -maxdepth 1 -type f \
    -name '*_fastCCLasso_freq_genus_pairs.tsv' \
    -printf '%f\n' \
| sed 's/_fastCCLasso_freq_genus_pairs\.tsv$//' \
| LC_ALL=C sort -u \
> "$BASE/completed_genus_projects.list"

log "Completed genus projects: $(count_lines "$BASE/completed_genus_projects.list")"
log "Output: $BASE"
