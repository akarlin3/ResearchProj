#!/usr/bin/env bash
# One-command re-validation for Echo. Mirrors Minos/future/reproduce.sh.
#
#   bash reproduce.sh          # default
#   FULL=1 bash reproduce.sh   # full-N harness self-test
#
# CP1 method self-test always runs (SOLID, synthetic). CP2 data check + CP3 validation run
# only if real test-retest data is present (download-on-demand; not committed). CP4
# consistency runs only if the manuscript exists. Stops at the first hard failure.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python}"
N_SELFTEST=$([ "${FULL:-0}" = "1" ] && echo 50000 || echo 20000)

run_stage() {
  local name="$1"; shift
  echo ">>> ${name}: running"
  if "$@"; then echo ">>> ${name}: PASS"; else
    local rc=$?
    if [ "$rc" = "3" ]; then echo ">>> ${name}: SKIP (data gate not satisfied)"; return 0; fi
    echo ">>> ${name}: FAIL (rc=$rc)"; exit "$rc"
  fi
}

echo "=== Echo re-validation (CP1 -> CP4) ==="
run_stage "CP1 method self-test (SOLID)" "$PY" "$HERE/scripts/run_harness.py" --n "$N_SELFTEST"

echo ">>> CP2 data check"
"$PY" "$HERE/scripts/fetch_invivo.py" --check || true

run_stage "CP3 real-data validation (PROVISIONAL)" "$PY" "$HERE/scripts/run_validation.py"

if [ -f "$HERE/paper/consistency.py" ] && [ -f "$HERE/paper/echo.tex" ]; then
  run_stage "CP4 manuscript consistency" "$PY" "$HERE/paper/consistency.py"
else
  echo ">>> CP4 manuscript consistency: SKIP (manuscript not built — PASS-only)"
fi

echo "=== done ==="
