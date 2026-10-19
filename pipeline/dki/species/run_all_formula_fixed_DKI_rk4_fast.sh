#!/usr/bin/env bash
set -u

BASE=/mnt/new_data/3_18_DKI_species
NEWBASE=${BASE}/new_DKI_rk4_fast
SCRIPT=${NEWBASE}/scripts/DKI_formula_fixed.py
LIST=${NEWBASE}/all_project.txt
LOGDIR=${NEWBASE}/logs
REPORT=${NEWBASE}/reports/fixed_dki_rk4_fast_all_status.tsv

mkdir -p "$LOGDIR" "${NEWBASE}/reports"

echo -e "project\tstatus\tstart_time\tend_time\truntime_seconds\ttrain_file\ttest_file\toutput_dir\tlog_file\tmessage" > "$REPORT"

while IFS= read -r PROJECT_RAW; do
  PROJECT=$(echo "$PROJECT_RAW" | xargs)
  [ -z "$PROJECT" ] && continue

  TRAIN=${BASE}/${PROJECT}/data/Ptrain.csv
  TEST=${BASE}/${PROJECT}/data/Ptest.csv
  OUT=${NEWBASE}/${PROJECT}/results
  LOG=${LOGDIR}/${PROJECT}.formula_fixed_rk4_fast.log

  START_TIME=$(date "+%F %T")
  START_SEC=$(date +%s)

  echo "===== Running $PROJECT ====="
  echo "TRAIN=$TRAIN"
  echo "TEST=$TEST"
  echo "OUT=$OUT"

  if [ ! -s "$TRAIN" ] || [ ! -s "$TEST" ]; then
    END_TIME=$(date "+%F %T")
    END_SEC=$(date +%s)
    RUNTIME=$((END_SEC - START_SEC))
    echo "[MISSING_INPUT] $PROJECT"
    echo -e "${PROJECT}\tMISSING_INPUT\t${START_TIME}\t${END_TIME}\t${RUNTIME}\t${TRAIN}\t${TEST}\t${OUT}\t${LOG}\tmissing Ptrain or Ptest" >> "$REPORT"
    continue
  fi

  if [ -s "${OUT}/qtst.csv" ] && [ -s "${OUT}/qtrn.csv" ] && [ -s "${OUT}/model_metrics.json" ]; then
    END_TIME=$(date "+%F %T")
    END_SEC=$(date +%s)
    RUNTIME=$((END_SEC - START_SEC))
    echo "[SKIP] Existing outputs found: $PROJECT"
    echo -e "${PROJECT}\tSKIP\t${START_TIME}\t${END_TIME}\t${RUNTIME}\t${TRAIN}\t${TEST}\t${OUT}\t${LOG}\texisting outputs" >> "$REPORT"
    continue
  fi

  mkdir -p "$OUT"

  python "$SCRIPT" \
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
  END_TIME=$(date "+%F %T")
  END_SEC=$(date +%s)
  RUNTIME=$((END_SEC - START_SEC))

  if [ $EXIT_CODE -ne 0 ]; then
    echo "[FAILED] $PROJECT, see $LOG"
    echo -e "${PROJECT}\tFAILED\t${START_TIME}\t${END_TIME}\t${RUNTIME}\t${TRAIN}\t${TEST}\t${OUT}\t${LOG}\texit_code_${EXIT_CODE}" >> "$REPORT"
  else
    echo "[DONE] $PROJECT"
    echo -e "${PROJECT}\tDONE\t${START_TIME}\t${END_TIME}\t${RUNTIME}\t${TRAIN}\t${TEST}\t${OUT}\t${LOG}\tok" >> "$REPORT"
  fi

done < "$LIST"

echo "All projects finished or skipped."
echo "Status report: $REPORT"
