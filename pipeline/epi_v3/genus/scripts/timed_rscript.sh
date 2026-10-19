#!/usr/bin/env bash
set -euo pipefail
RESOURCE_LOG="${MACKDB_EPI_RESOURCE_LOG:-${TMPDIR:?}/resource_usage.txt}"
if [[ "$TMPDIR" =~ [[:space:]] ]]; then
  R_TMP_WORK="$(mktemp -d /tmp/mackdb_genus_r.XXXXXX)"
  trap 'rm -rf -- "$R_TMP_WORK"' EXIT
  TMPDIR="$R_TMP_WORK" TMP="$R_TMP_WORK" TEMP="$R_TMP_WORK" \
    /usr/bin/time -v -o "$RESOURCE_LOG" "${MACKDB_RSCRIPT:-Rscript}" "$@"
else
  exec /usr/bin/time -v -o "$RESOURCE_LOG" "${MACKDB_RSCRIPT:-Rscript}" "$@"
fi
