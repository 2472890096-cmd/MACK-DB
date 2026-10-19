#!/bin/bash

# --- 1. 配置基础路径 ---
if [[ $# != 4 ]]; then
    echo "Usage: bash batch_matrix.sh PROJECT_LIST SAMPLE_DIR PROFILE_DIR OUTPUT_DIR" >&2
    exit 2
fi
PROJECT_LIST="$(realpath -- "$1")"
SAMPLE_DIR="$(realpath -- "$2")"
PROFILE_DIR="$(realpath -- "$3")"
mkdir -p -- "$4"
BASE_DIR="$(realpath -- "$4")"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/../../abundance/shape_sparcc_format.py"
TASK_TMP="$(mktemp -d "${TMPDIR:-/tmp}/mackdb-matrix.XXXXXX")"
trap 'rm -rf -- "$TASK_TMP"' EXIT

# --- 2. 环境检查 ---
if [ ! -f "$PROJECT_LIST" ]; then
    echo "错误: 找不到文件 $PROJECT_LIST，请先运行校验脚本生成它！"
    exit 1
fi

# --- 3. 循环批量运行 (临时缓存 + 动态分类) ---
echo "开始动态阈值批量生成 Project 丰度矩阵..."
count=0

while read -r project_id || [[ -n "$project_id" ]]; do
    project_id=$(echo "$project_id" | tr -d '[:space:]')
    if [ -z "$project_id" ]; then continue; fi

    list_file="${SAMPLE_DIR}/${project_id}_uniq.txt"
    if [ ! -f "$list_file" ]; then
        echo "警告: 找不到 ${project_id} 的样本列表 ($list_file)，已跳过。"
        continue
    fi

    ((count++))
    echo "--------------------------------------------------"
    echo "[${count}] 正在处理: ${project_id} ..."

    # 初始阈值
    f_val=0.1
    
    # 巧妙利用系统的 /tmp 目录存放测试文件，避免产生废弃的空文件夹
    tmp_profile="${TASK_TMP}/${project_id}_tmp_abund.matrix"
    tmp_species="${TASK_TMP}/${project_id}_tmp_spec.matrix"

    while true; do
        # 格式化保留一位小数 (比如 0.1, 0.2)
        f_val_str=$(printf "%.1f" $f_val)
        
        # 1. 在临时目录生成矩阵
        "${MACKDB_PYTHON:-python3}" "$PY_SCRIPT" -l "$list_file" -i "$PROFILE_DIR" -f "$f_val_str" -p "$tmp_profile" -s "$tmp_species" || exit "$?"
        
        # 2. 统计生成的物种数
        species_count=$(wc -l < "$tmp_species")

        # 3. 判断是否达标 (或是否已经达到最高阈值 1.0)
        if [ "$species_count" -le 2000 ] || (( $(awk "BEGIN {print ($f_val >= 1.0)}") )); then
            
            if [ "$species_count" -le 2000 ]; then
                echo "  -> 成功达标！(采用阈值: $f_val_str, 最终物种数: $species_count)"
            else
                echo "  -> [警告] 阈值已达 1.0，物种数($species_count)仍大于 2000，强制保留结果。"
            fi

            # --- 核心修改：只有成功时，才创建对应的阈值文件夹 ---
            FINAL_DIR="${BASE_DIR}/filtered_${f_val_str}"
            mkdir -p "${FINAL_DIR}/abundance" "${FINAL_DIR}/species"

            # 把临时文件移动到最终归宿
            mv "$tmp_profile" "${FINAL_DIR}/abundance/${project_id}_abundance.matrix"
            mv "$tmp_species" "${FINAL_DIR}/species/${project_id}_species.matrix"
            
            # 打断死循环，进入下一个 project
            break
        else
            echo "  -> 超标 (阈值: $f_val_str, 物种数: $species_count > 2000)。正在提升阈值重试..."
            # 删除超标的临时文件，清理现场
            rm -f "$tmp_profile" "$tmp_species"
            
            # 阈值加 0.1，继续循环
            f_val=$(awk "BEGIN {print $f_val + 0.1}")
        fi
    done

done < "$PROJECT_LIST"

echo "=================================================="
echo "恭喜！动态降维与自动分类处理完成，共处理了 $count 个 Projects。"
