#!/usr/bin/env bash
set -uo pipefail

BASE=/mnt/new_data/3_18_DKI_genus
NEWBASE=$BASE/new_DKI_rk4_fast

SCRIPT=$NEWBASE/script/DKI_formula_fixed.py
LIST=$NEWBASE/all_project.txt

LOGDIR=$NEWBASE/logs
REPORTDIR=$NEWBASE/reports
REPORT=$REPORTDIR/genus_dki_formula_fixed_status.tsv

mkdir -p "$LOGDIR" "$REPORTDIR"

if [ ! -s "$SCRIPT" ]; then
    echo "ERROR: 找不到 DKI 脚本：$SCRIPT"
    exit 1
fi

if [ ! -s "$LIST" ]; then
    echo "ERROR: 找不到 project 列表：$LIST"
    exit 1
fi

TOTAL=$(grep -cve '^[[:space:]]*$' "$LIST")

echo "============================================================"
echo "Genus DKI formula-fixed batch"
echo "Original data base: $BASE"
echo "New output base:    $NEWBASE"
echo "Project list:       $LIST"
echo "Projects:           $TOTAL"
echo "Script:             $SCRIPT"
echo "============================================================"

echo -e \
"project\tstatus\tstart_time\tend_time\truntime_seconds\tN_species\tM_train\tM_test\tmessage" \
> "$REPORT"

CURRENT=0

while IFS= read -r PROJECT_RAW; do
    PROJECT=$(printf '%s' "$PROJECT_RAW" | xargs)
    [ -z "$PROJECT" ] && continue

    CURRENT=$((CURRENT + 1))

    if [[ ! "$PROJECT" =~ ^PRJ[A-Z]{2}[0-9]+$ ]]; then
        echo "[$CURRENT/$TOTAL] [SKIP_INVALID_PROJECT] $PROJECT"
        continue
    fi

    DATA_DIR=$BASE/$PROJECT/data
    TRAIN=$DATA_DIR/Ptrain.csv
    TEST=$DATA_DIR/Ptest.csv
    SPECIES_ID=$DATA_DIR/species_id.csv
    SAMPLE_ID=$DATA_DIR/sample_id.csv

    OUT=$NEWBASE/$PROJECT/results
    LOG=$LOGDIR/${PROJECT}.dki_formula_fixed.log

    START_TIME=$(date '+%F %T')
    START_SEC=$(date +%s)

    echo
    echo "============================================================"
    echo "[$CURRENT/$TOTAL] Project: $PROJECT"
    echo "Start: $START_TIME"
    echo "Train: $TRAIN"
    echo "Test:  $TEST"
    echo "Out:   $OUT"
    echo "============================================================"

    # --------------------------------------------------------
    # 检查原始 DKI 输入
    # --------------------------------------------------------
    MISSING=""

    [ -s "$TRAIN" ]      || MISSING="${MISSING}Ptrain.csv,"
    [ -s "$TEST" ]       || MISSING="${MISSING}Ptest.csv,"
    [ -s "$SPECIES_ID" ] || MISSING="${MISSING}species_id.csv,"
    [ -s "$SAMPLE_ID" ]  || MISSING="${MISSING}sample_id.csv,"

    if [ -n "$MISSING" ]; then
        END_TIME=$(date '+%F %T')
        END_SEC=$(date +%s)
        RUNTIME=$((END_SEC - START_SEC))

        echo "[MISSING_INPUT] $PROJECT：${MISSING%,}"

        echo -e \
"${PROJECT}\tMISSING_INPUT\t${START_TIME}\t${END_TIME}\t${RUNTIME}\tNA\tNA\tNA\t${MISSING%,}" \
        >> "$REPORT"

        continue
    fi

    # --------------------------------------------------------
    # 如果新版完整结果已经存在，则跳过
    # --------------------------------------------------------
    if [ -s "$OUT/qtst.csv" ] \
       && [ -s "$OUT/qtrn.csv" ] \
       && [ -s "$OUT/model_metrics.json" ]; then

        END_TIME=$(date '+%F %T')
        END_SEC=$(date +%s)
        RUNTIME=$((END_SEC - START_SEC))

        echo "[SKIP] 新版结果已经存在：$PROJECT"

        echo -e \
"${PROJECT}\tSKIP\t${START_TIME}\t${END_TIME}\t${RUNTIME}\tNA\tNA\tNA\texisting complete outputs" \
        >> "$REPORT"

        continue
    fi

    # 如果只有部分输出，删除不完整结果后重跑
    if [ -d "$OUT" ]; then
        echo "[INFO] 删除不完整的新版结果目录：$OUT"
        rm -rf "$OUT"
    fi

    mkdir -p "$OUT"

    # 记录输入维度，便于排查超大项目
    N_SPECIES=$(wc -l < "$TRAIN")
    M_TRAIN=$(awk -F',' 'NR==1{print NF; exit}' "$TRAIN")
    M_TEST=$(awk -F',' 'NR==1{print NF; exit}' "$TEST")

    echo "Input dimensions:"
    echo "  N_genus=$N_SPECIES"
    echo "  M_train_samples=$M_TRAIN"
    echo "  M_test_perturbations=$M_TEST"

    # --------------------------------------------------------
    # 运行修正公式后的 DKI
    # --------------------------------------------------------
    python3 "$SCRIPT" \
      --train-file "$TRAIN" \
      --test-file "$TEST" \
      --output-dir "$OUT" \
      --project "$PROJECT" \
      --max-epochs 1000 \
      --early-stop-patience 80 \
      --lr 0.005 \
      --batch-size 32 \
      --activation relu \
      --bottleneck auto \
      --normalize-z true \
      --time-mode endpoint \
      --ode-method rk4 \
      --ode-step-size 0.05 \
      --device cpu \
      --seed 20260629 \
      > "$LOG" 2>&1

    EXIT_CODE=$?

    END_TIME=$(date '+%F %T')
    END_SEC=$(date +%s)
    RUNTIME=$((END_SEC - START_SEC))

    if [ "$EXIT_CODE" -ne 0 ]; then
        echo "[FAILED_DKI] $PROJECT，exit code=$EXIT_CODE"
        echo "详细日志：$LOG"

        echo -e \
"${PROJECT}\tFAILED_DKI\t${START_TIME}\t${END_TIME}\t${RUNTIME}\t${N_SPECIES}\t${M_TRAIN}\t${M_TEST}\texit_code_${EXIT_CODE}; see ${LOG}" \
        >> "$REPORT"

        continue
    fi

    # --------------------------------------------------------
    # 检查必要输出
    # --------------------------------------------------------
    if [ -s "$OUT/qtst.csv" ] \
       && [ -s "$OUT/qtrn.csv" ] \
       && [ -s "$OUT/model_metrics.json" ]; then

        echo "[DONE] $PROJECT，运行时间：$RUNTIME 秒"

        echo -e \
"${PROJECT}\tDONE\t${START_TIME}\t${END_TIME}\t${RUNTIME}\t${N_SPECIES}\t${M_TRAIN}\t${M_TEST}\tok" \
        >> "$REPORT"
    else
        echo "[FAILED_OUTPUT_CHECK] $PROJECT"
        echo "输出目录：$OUT"

        echo -e \
"${PROJECT}\tFAILED_OUTPUT_CHECK\t${START_TIME}\t${END_TIME}\t${RUNTIME}\t${N_SPECIES}\t${M_TRAIN}\t${M_TEST}\toutputs incomplete" \
        >> "$REPORT"
    fi

done < "$LIST"

echo
echo "============================================================"
echo "Genus DKI batch finished"
echo "Status report: $REPORT"
echo "============================================================"

awk -F'\t' '
NR > 1 {
    count[$2]++
}
END {
    for (status in count) {
        print status, count[status]
    }
}' "$REPORT"
