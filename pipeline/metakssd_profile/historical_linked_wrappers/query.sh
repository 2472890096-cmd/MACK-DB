#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 3 ]]; then
  echo "Usage: bash query.sh REFERENCE_MARKER_DB SAMPLE_SKETCH OUTPUT_TXT" >&2
  exit 2
fi
META_KSSD="${METAKSSD_BIN:-metakssd}"
command -v "$META_KSSD" >/dev/null || { echo "Missing MetaKSSD executable: $META_KSSD" >&2; exit 127; }
mkdir -p -- "$(dirname -- "$3")"
"$META_KSSD" composite -r "$1" -q "$2" > "$3"
