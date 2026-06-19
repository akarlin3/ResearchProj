#!/usr/bin/env bash
# One-command re-validation of the Minos future/ build.
#
# Runs, in order: CP1 (theory reproduces via read-only import) -> CP2 (applied gap) ->
# CP3 (applied monitor) -> CP4 (manuscript consistency). Each stage is a gate; the script
# stops at the first failure. Stages not yet built are reported as PENDING (not failures).
#
# This is the re-validation contract from ASSUMPTIONS.md: after Fashion/Gauge publish, bump
# the pinned versions and run this once to learn what still holds.
#
# Usage:  bash Minos/future/reproduce.sh            # FAST (smoke; default)
#         FULL=1 bash Minos/future/reproduce.sh     # full-N gates
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
FULLFLAG=""; [ "${FULL:-0}" = "1" ] && FULLFLAG="--full"

if [ ! -x "$PROT" ]; then
  echo "FATAL: python interpreter not found at PROT=$PROT (set PROT=... to override)"; exit 2
fi

rc=0
run_stage() {  # name  script  [args...]
  local name="$1"; local script="$2"; shift 2
  echo; echo "######## $name ########"
  if [ ! -f "$script" ]; then
    echo ">>> PENDING: $script not built yet"; return 0
  fi
  if "$PROT" "$script" "$@"; then
    echo ">>> $name: PASS"
  else
    echo ">>> $name: FAIL"; rc=1
  fi
}

run_stage "CP1 theory reproduces (read-only import)" "$HERE/verify_cp1.py" $FULLFLAG
run_stage "CP2 applied decision-calibration gap"     "$HERE/applied/gap_applied.py" $FULLFLAG
run_stage "CP3 applied validity monitor"             "$HERE/applied/monitor_applied.py" $FULLFLAG
run_stage "CP4 manuscript consistency"               "$HERE/paper/consistency.py"

echo
if [ "$rc" -eq 0 ]; then
  echo "================ reproduce.sh: all built stages GREEN ================"
else
  echo "================ reproduce.sh: a stage FAILED (rc=$rc) ================"
fi
exit "$rc"
