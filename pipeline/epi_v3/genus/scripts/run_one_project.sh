#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${MACKDB_EPI_RUN_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd)}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export HOME="$ROOT/tmp/runtime_home" XDG_CACHE_HOME="$ROOT/tmp/cache" TMPDIR="$ROOT/tmp"
export TMP="$TMPDIR" TEMP="$TMPDIR" LC_ALL=C.UTF-8
cd -- "$ROOT"
exec "${MACKDB_PYTHON:-python3}" -B "$SCRIPT_DIR/species_production_workflow.py" worker "$@"
