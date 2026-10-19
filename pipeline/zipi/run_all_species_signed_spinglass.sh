#!/usr/bin/env bash

# ============================================================
# 批量运行全部 species-level project 的有符号网络分析
#
# 自动从 FASTSPAR_DIR 中识别 *_correlation.tsv 对应的 project。
# 对每个 project 检查：
#   1) correlation.tsv
#   2) pvalues.tsv
#   3) abundance.matrix
#   4) species.matrix
#
# 默认顺序运行，避免同时运行多个 signed spin-glass 占用过多内存。
# 已成功项目通过 done/*.done 自动跳过，可安全断点续跑。
# 若确实需要全部重跑：FORCE=1 ./run_all_species_signed_spinglass.sh
# ============================================================

set -u
set -o pipefail

# -------------------- 固定路径 --------------------
WORK_DIR="/home/wangmengyang/project/keygenus/no_isolate_node/two_pairs_test"
R_SCRIPT="${WORK_DIR}/FastSpar_signed_spinglass_network.R"

FASTSPAR_DIR="/home/wangmengyang/project/r226_fastspar_species_result"
ABUNDANCE_DIR="/home/wangmengyang/project/GDTB226_species_abundance_result/abundance"
SPECIES_DIR="/home/wangmengyang/project/GDTB226_species_abundance_result/species"

# 全部批量结果单独放在这个目录，避免与前面的少量测试结果混在一起。
OUT_ROOT="${WORK_DIR}/species_all_signed"

TOPOMAP_DIR="${OUT_ROOT}/topomap"
NETWORK_DIR="${OUT_ROOT}/network_html"
COABUND_DIR="${OUT_ROOT}/coabund_html"
GEXF_DIR="${OUT_ROOT}/gexf"
NODE_DIR="${OUT_ROOT}/node_tables"
LOG_DIR="${OUT_ROOT}/logs"
DONE_DIR="${OUT_ROOT}/done"
STATUS_DIR="${OUT_ROOT}/status"

# FORCE=0：已有 done 标记的项目自动跳过；FORCE=1：全部重跑。
FORCE="${FORCE:-0}"

# -------------------- 分析参数 --------------------
COR_THRESHOLD="0.3"
SPINS="25"
GAMMA="1"
GAMMA_MINUS="1"
SEED="12345"
PLOT_DEGREE_THRESHOLD="12"
PLOT_TOP_K="3"
POSITIVE_EDGE_COLOR="#D55E00"
NEGATIVE_EDGE_COLOR="#0072B2"

mkdir -p \
    "${TOPOMAP_DIR}" \
    "${NETWORK_DIR}" \
    "${COABUND_DIR}" \
    "${GEXF_DIR}" \
    "${NODE_DIR}" \
    "${LOG_DIR}" \
    "${DONE_DIR}" \
    "${STATUS_DIR}"

# -------------------- 前置检查 --------------------
if ! command -v Rscript >/dev/null 2>&1; then
    echo "错误：当前环境找不到 Rscript。请先执行 conda activate r_network。" >&2
    exit 1
fi

if [[ ! -s "${R_SCRIPT}" ]]; then
    echo "错误：找不到R脚本或脚本为空：${R_SCRIPT}" >&2
    exit 1
fi

for required_dir in "${FASTSPAR_DIR}" "${ABUNDANCE_DIR}" "${SPECIES_DIR}"; do
    if [[ ! -d "${required_dir}" ]]; then
        echo "错误：目录不存在：${required_dir}" >&2
        exit 1
    fi
done

# -------------------- 自动获得全部project --------------------
mapfile -t PROJECTS < <(
    find "${FASTSPAR_DIR}" -maxdepth 1 -type f -name '*_correlation.tsv' -printf '%f\n' \
        | sed 's/_correlation\.tsv$//' \
        | sort -V
)

TOTAL="${#PROJECTS[@]}"
if [[ "${TOTAL}" -eq 0 ]]; then
    echo "错误：${FASTSPAR_DIR} 中没有找到 *_correlation.tsv。" >&2
    exit 1
fi

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
STATUS_FILE="${STATUS_DIR}/status_${RUN_ID}.tsv"
LATEST_STATUS="${STATUS_DIR}/latest_status.tsv"

printf 'index\ttotal\tproject\tstatus\texit_code\tstart_time\tend_time\tduration_seconds\tmessage\n' \
    > "${STATUS_FILE}"

cleanup_project_outputs() {
    local project="$1"

    rm -f "${TOPOMAP_DIR}/${project}_"* 2>/dev/null || true

    # networkD3/chorddiag 可能生成与HTML同名的 *_files 依赖目录，必须一并清理。
    rm -rf "${NETWORK_DIR}/${project}_"* 2>/dev/null || true
    rm -rf "${COABUND_DIR}/${project}_"* 2>/dev/null || true

    rm -f "${GEXF_DIR}/${project}.gexf" \
          "${GEXF_DIR}/${project}.graphml" 2>/dev/null || true
    rm -f "${NODE_DIR}/${project}_"* 2>/dev/null || true
}

write_status() {
    local index="$1"
    local project="$2"
    local status="$3"
    local exit_code="$4"
    local start_time="$5"
    local end_time="$6"
    local duration="$7"
    local message="$8"

    # 避免message中的制表符或换行破坏TSV格式。
    message="${message//$'\t'/ }"
    message="${message//$'\n'/ }"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${index}" "${TOTAL}" "${project}" "${status}" "${exit_code}" \
        "${start_time}" "${end_time}" "${duration}" "${message}" \
        >> "${STATUS_FILE}"
}

SUCCESS=0
FAILED=0
MISSING=0
SKIPPED=0

BATCH_START_EPOCH="$(date +%s)"

echo "============================================================"
echo "species有符号网络批量分析开始"
echo "project总数：${TOTAL}"
echo "R脚本：${R_SCRIPT}"
echo "结果目录：${OUT_ROOT}"
echo "断点续跑：FORCE=${FORCE}（0=跳过已完成，1=全部重跑）"
echo "============================================================"

for i in "${!PROJECTS[@]}"; do
    index=$((i + 1))
    project="${PROJECTS[$i]}"

    COR_FILE="${FASTSPAR_DIR}/${project}_correlation.tsv"
    P_FILE="${FASTSPAR_DIR}/${project}_pvalues.tsv"
    ABUNDANCE_FILE="${ABUNDANCE_DIR}/${project}_abundance.matrix"
    SPECIES_FILE="${SPECIES_DIR}/${project}_species.matrix"
    LOG_FILE="${LOG_DIR}/${project}.log"
    DONE_FILE="${DONE_DIR}/${project}.done"

    start_time="$(date '+%F %T')"
    start_epoch="$(date +%s)"

    echo
    echo "============================================================"
    echo "[${index}/${TOTAL}] 开始处理：${project}"
    echo "============================================================"

    if [[ "${FORCE}" != "1" && -s "${DONE_FILE}" ]]; then
        end_time="$(date '+%F %T')"
        end_epoch="$(date +%s)"
        duration=$((end_epoch - start_epoch))
        write_status \
            "${index}" "${project}" "SKIPPED_DONE" "0" \
            "${start_time}" "${end_time}" "${duration}" \
            "已有成功完成标记，断点续跑时跳过"
        SKIPPED=$((SKIPPED + 1))
        echo "${project}：已完成，跳过。"
        continue
    fi

    missing_files=()
    for input_file in \
        "${COR_FILE}" \
        "${P_FILE}" \
        "${ABUNDANCE_FILE}" \
        "${SPECIES_FILE}"
    do
        if [[ ! -s "${input_file}" ]]; then
            missing_files+=("${input_file}")
        fi
    done

    if [[ "${#missing_files[@]}" -gt 0 ]]; then
        end_time="$(date '+%F %T')"
        end_epoch="$(date +%s)"
        duration=$((end_epoch - start_epoch))
        missing_message="缺少或为空：$(IFS=';'; echo "${missing_files[*]}")"
        printf '%s\n' "${missing_message}" > "${LOG_FILE}"
        write_status \
            "${index}" "${project}" "MISSING_INPUT" "NA" \
            "${start_time}" "${end_time}" "${duration}" \
            "${missing_message}"
        MISSING=$((MISSING + 1))
        echo "${project}：输入文件不完整，跳过。"
        continue
    fi

    # 对失败后重新运行的项目，先清除旧的残缺输出，避免新旧结果混杂。
    cleanup_project_outputs "${project}"
    rm -f "${DONE_FILE}"

    echo "相关性矩阵：${COR_FILE}"
    echo "P值矩阵：${P_FILE}"
    echo "丰度矩阵：${ABUNDANCE_FILE}"
    echo "species文件：${SPECIES_FILE}"
    echo "日志：${LOG_FILE}"

    Rscript "${R_SCRIPT}" \
        -p "${project}" \
        -l species \
        -f "${COR_THRESHOLD}" \
        -C "${COR_FILE}" \
        -P "${P_FILE}" \
        -a "${TOPOMAP_DIR}" \
        -b "${NETWORK_DIR}" \
        -c "${COABUND_DIR}" \
        -g "${GEXF_DIR}" \
        -d "${NODE_DIR}" \
        --spins "${SPINS}" \
        --gamma "${GAMMA}" \
        --gamma-minus "${GAMMA_MINUS}" \
        --seed "${SEED}" \
        --plot-degree-threshold "${PLOT_DEGREE_THRESHOLD}" \
        --plot-top-k "${PLOT_TOP_K}" \
        --positive-edge-color "${POSITIVE_EDGE_COLOR}" \
        --negative-edge-color "${NEGATIVE_EDGE_COLOR}" \
        > "${LOG_FILE}" 2>&1

    exit_code=$?
    end_time="$(date '+%F %T')"
    end_epoch="$(date +%s)"
    duration=$((end_epoch - start_epoch))

    if [[ "${exit_code}" -eq 0 ]]; then
        {
            echo -e "project\t${project}"
            echo -e "completed_at\t${end_time}"
            echo -e "script\t${R_SCRIPT}"
            echo -e "correlation_threshold\t${COR_THRESHOLD}"
            echo -e "spins\t${SPINS}"
            echo -e "gamma\t${GAMMA}"
            echo -e "gamma_minus\t${GAMMA_MINUS}"
            echo -e "seed\t${SEED}"
        } > "${DONE_FILE}"

        write_status \
            "${index}" "${project}" "SUCCESS" "${exit_code}" \
            "${start_time}" "${end_time}" "${duration}" \
            "运行成功"
        SUCCESS=$((SUCCESS + 1))
        echo "${project}：运行成功，用时 ${duration} 秒。"
    else
        last_error="$(tail -n 5 "${LOG_FILE}" 2>/dev/null | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
        write_status \
            "${index}" "${project}" "FAILED" "${exit_code}" \
            "${start_time}" "${end_time}" "${duration}" \
            "${last_error}"
        FAILED=$((FAILED + 1))
        echo "${project}：运行失败，退出码 ${exit_code}。"
        echo "查看日志：${LOG_FILE}"
    fi

done

BATCH_END_EPOCH="$(date +%s)"
BATCH_DURATION=$((BATCH_END_EPOCH - BATCH_START_EPOCH))

cp -f "${STATUS_FILE}" "${LATEST_STATUS}"

SUMMARY_FILE="${STATUS_DIR}/summary_${RUN_ID}.txt"
{
    echo "run_id=${RUN_ID}"
    echo "total=${TOTAL}"
    echo "success=${SUCCESS}"
    echo "failed=${FAILED}"
    echo "missing_input=${MISSING}"
    echo "skipped_done=${SKIPPED}"
    echo "duration_seconds=${BATCH_DURATION}"
    echo "status_file=${STATUS_FILE}"
} > "${SUMMARY_FILE}"

echo
echo "============================================================"
echo "全部species project处理结束"
echo "总数：${TOTAL}"
echo "成功：${SUCCESS}"
echo "失败：${FAILED}"
echo "缺少输入：${MISSING}"
echo "已完成跳过：${SKIPPED}"
echo "总用时：${BATCH_DURATION} 秒"
echo "状态表：${STATUS_FILE}"
echo "最新状态：${LATEST_STATUS}"
echo "汇总：${SUMMARY_FILE}"
echo "============================================================"

# 仅在存在失败或缺失输入时返回非零退出码，便于nohup后自动判断。
if [[ "${FAILED}" -gt 0 || "${MISSING}" -gt 0 ]]; then
    exit 2
fi

exit 0
