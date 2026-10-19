#!/usr/bin/env bash

# ============================================================
# 批量运行全部 genus-level project 的有符号网络分析
#
# 输入目录：
#   abundance : /home/wangmengyang/project/GDTB226_genus_abundance_result/abundance
#   genus     : /home/wangmengyang/project/GDTB226_genus_abundance_result/genus
#   cor       : /home/wangmengyang/project/GDTB226_genus_abundance_result/fastspar_result/cor
#   pv        : /home/wangmengyang/project/GDTB226_genus_abundance_result/fastspar_result/pv
#
# 特点：
#   1) 从 *_abundance.matrix 自动识别全部 project；
#   2) 自动兼容常见 correlation / pvalue 文件名；
#   3) 顺序运行，避免多个 signed spin-glass 同时占用大量内存；
#   4) done/*.done 支持断点续跑；
#   5) 单个 project 失败后继续处理后续 project；
#   6) FORCE=1 可强制全部重跑。
# ============================================================

set -u
set -o pipefail

# -------------------- 固定路径 --------------------
WORK_DIR="/home/wangmengyang/project/keygenus/no_isolate_node/two_pairs_test"
R_SCRIPT="${WORK_DIR}/FastSpar_signed_spinglass_network.R"

BASE_DIR="/home/wangmengyang/project/GDTB226_genus_abundance_result"
ABUNDANCE_DIR="${BASE_DIR}/abundance"
GENUS_DIR="${BASE_DIR}/genus"
COR_DIR="${BASE_DIR}/fastspar_result/cor"
PV_DIR="${BASE_DIR}/fastspar_result/pv"

# genus全部结果单独保存，不与species结果混合。
OUT_ROOT="${WORK_DIR}/genus_all_signed"

TOPOMAP_DIR="${OUT_ROOT}/topomap"
NETWORK_DIR="${OUT_ROOT}/network_html"
COABUND_DIR="${OUT_ROOT}/coabund_html"
GEXF_DIR="${OUT_ROOT}/gexf"
NODE_DIR="${OUT_ROOT}/node_tables"
LOG_DIR="${OUT_ROOT}/logs"
DONE_DIR="${OUT_ROOT}/done"
STATUS_DIR="${OUT_ROOT}/status"

# FORCE=0：跳过已有done标记；FORCE=1：全部重跑。
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

for required_dir in \
    "${ABUNDANCE_DIR}" \
    "${GENUS_DIR}" \
    "${COR_DIR}" \
    "${PV_DIR}"
do
    if [[ ! -d "${required_dir}" ]]; then
        echo "错误：目录不存在：${required_dir}" >&2
        exit 1
    fi
done

# -------------------- 文件名兼容函数 --------------------
# 优先使用标准FastSpar命名；若目录中某project只有一个候选文件，则自动采用。
resolve_matrix_file() {
    local directory="$1"
    local project="$2"
    local kind="$3"
    local candidate
    local -a candidates=()
    local -a fallback_matches=()

    if [[ "${kind}" == "cor" ]]; then
        candidates=(
            "${directory}/${project}_correlation.tsv"
            "${directory}/${project}_cor.tsv"
            "${directory}/${project}.correlation.tsv"
            "${directory}/${project}.cor.tsv"
            "${directory}/${project}_correlation.txt"
            "${directory}/${project}_cor.txt"
            "${directory}/${project}_correlation.matrix"
            "${directory}/${project}.tsv"
            "${directory}/${project}.txt"
            "${directory}/${project}"
        )
    else
        candidates=(
            "${directory}/${project}_pvalues.tsv"
            "${directory}/${project}_pvalue.tsv"
            "${directory}/${project}_p_values.tsv"
            "${directory}/${project}_p-value.tsv"
            "${directory}/${project}_pv.tsv"
            "${directory}/${project}_p.tsv"
            "${directory}/${project}.pvalues.tsv"
            "${directory}/${project}.pvalue.tsv"
            "${directory}/${project}_pvalues.txt"
            "${directory}/${project}_pvalue.txt"
            "${directory}/${project}_pvalues.matrix"
            "${directory}/${project}.tsv"
            "${directory}/${project}.txt"
            "${directory}/${project}"
        )
    fi

    for candidate in "${candidates[@]}"; do
        if [[ -s "${candidate}" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done

    # 兜底：若该目录中以project开头的非空文件只有一个，则使用它。
    mapfile -t fallback_matches < <(
        find "${directory}" -maxdepth 1 -type f -name "${project}*" -size +0c -print \
            | sort -V
    )

    if [[ "${#fallback_matches[@]}" -eq 1 ]]; then
        printf '%s\n' "${fallback_matches[0]}"
        return 0
    fi

    return 1
}

# -------------------- 自动获得全部project --------------------
mapfile -t PROJECTS < <(
    find "${ABUNDANCE_DIR}" -maxdepth 1 -type f -name '*_abundance.matrix' -printf '%f\n' \
        | sed 's/_abundance\.matrix$//' \
        | sort -V
)

TOTAL="${#PROJECTS[@]}"
if [[ "${TOTAL}" -eq 0 ]]; then
    echo "错误：${ABUNDANCE_DIR} 中没有找到 *_abundance.matrix。" >&2
    exit 1
fi

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
STATUS_FILE="${STATUS_DIR}/status_${RUN_ID}.tsv"
LATEST_STATUS="${STATUS_DIR}/latest_status.tsv"

printf 'index\ttotal\tproject\tstatus\texit_code\tstart_time\tend_time\tduration_seconds\tcor_file\tpvalue_file\tmessage\n' \
    > "${STATUS_FILE}"

cleanup_project_outputs() {
    local project="$1"

    rm -f "${TOPOMAP_DIR}/${project}_"* 2>/dev/null || true
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
    local cor_file="$8"
    local pvalue_file="$9"
    local message="${10}"

    message="${message//$'\t'/ }"
    message="${message//$'\n'/ }"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${index}" "${TOTAL}" "${project}" "${status}" "${exit_code}" \
        "${start_time}" "${end_time}" "${duration}" \
        "${cor_file}" "${pvalue_file}" "${message}" \
        >> "${STATUS_FILE}"
}

SUCCESS=0
FAILED=0
MISSING=0
SKIPPED=0

BATCH_START_EPOCH="$(date +%s)"

echo "============================================================"
echo "genus有符号网络批量分析开始"
echo "project总数：${TOTAL}"
echo "R脚本：${R_SCRIPT}"
echo "相关性目录：${COR_DIR}"
echo "P值目录：${PV_DIR}"
echo "结果目录：${OUT_ROOT}"
echo "断点续跑：FORCE=${FORCE}（0=跳过已完成，1=全部重跑）"
echo "============================================================"

for i in "${!PROJECTS[@]}"; do
    index=$((i + 1))
    project="${PROJECTS[$i]}"

    ABUNDANCE_FILE="${ABUNDANCE_DIR}/${project}_abundance.matrix"
    GENUS_FILE="${GENUS_DIR}/${project}_genus.matrix"
    COR_FILE="$(resolve_matrix_file "${COR_DIR}" "${project}" "cor" || true)"
    P_FILE="$(resolve_matrix_file "${PV_DIR}" "${project}" "pv" || true)"
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
            "${COR_FILE}" "${P_FILE}" \
            "已有成功完成标记，断点续跑时跳过"
        SKIPPED=$((SKIPPED + 1))
        echo "${project}：已完成，跳过。"
        continue
    fi

    missing_files=()
    [[ -s "${ABUNDANCE_FILE}" ]] || missing_files+=("${ABUNDANCE_FILE}")
    [[ -s "${GENUS_FILE}" ]] || missing_files+=("${GENUS_FILE}")
    [[ -n "${COR_FILE}" && -s "${COR_FILE}" ]] || missing_files+=("${COR_DIR}/${project} 对应相关性文件")
    [[ -n "${P_FILE}" && -s "${P_FILE}" ]] || missing_files+=("${PV_DIR}/${project} 对应P值文件")

    if [[ "${#missing_files[@]}" -gt 0 ]]; then
        end_time="$(date '+%F %T')"
        end_epoch="$(date +%s)"
        duration=$((end_epoch - start_epoch))
        missing_message="缺少或无法唯一识别：$(IFS=';'; echo "${missing_files[*]}")"
        printf '%s\n' "${missing_message}" > "${LOG_FILE}"
        write_status \
            "${index}" "${project}" "MISSING_INPUT" "NA" \
            "${start_time}" "${end_time}" "${duration}" \
            "${COR_FILE}" "${P_FILE}" "${missing_message}"
        MISSING=$((MISSING + 1))
        echo "${project}：输入文件不完整，跳过。"
        continue
    fi

    cleanup_project_outputs "${project}"
    rm -f "${DONE_FILE}"

    echo "相关性矩阵：${COR_FILE}"
    echo "P值矩阵：${P_FILE}"
    echo "丰度矩阵：${ABUNDANCE_FILE}"
    echo "genus文件：${GENUS_FILE}"
    echo "日志：${LOG_FILE}"

    Rscript "${R_SCRIPT}" \
        -p "${project}" \
        -l genus \
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
            echo -e "level\tgenus"
            echo -e "completed_at\t${end_time}"
            echo -e "script\t${R_SCRIPT}"
            echo -e "correlation_file\t${COR_FILE}"
            echo -e "pvalue_file\t${P_FILE}"
            echo -e "correlation_threshold\t${COR_THRESHOLD}"
            echo -e "spins\t${SPINS}"
            echo -e "gamma\t${GAMMA}"
            echo -e "gamma_minus\t${GAMMA_MINUS}"
            echo -e "seed\t${SEED}"
        } > "${DONE_FILE}"

        write_status \
            "${index}" "${project}" "SUCCESS" "${exit_code}" \
            "${start_time}" "${end_time}" "${duration}" \
            "${COR_FILE}" "${P_FILE}" "运行成功"
        SUCCESS=$((SUCCESS + 1))
        echo "${project}：运行成功，用时 ${duration} 秒。"
    else
        last_error="$(tail -n 8 "${LOG_FILE}" 2>/dev/null | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
        write_status \
            "${index}" "${project}" "FAILED" "${exit_code}" \
            "${start_time}" "${end_time}" "${duration}" \
            "${COR_FILE}" "${P_FILE}" "${last_error}"
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
    echo "level=genus"
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
echo "全部genus project处理结束"
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

if [[ "${FAILED}" -gt 0 || "${MISSING}" -gt 0 ]]; then
    exit 2
fi

exit 0
