#!/usr/bin/env bash
# Augur one-command re-validation.
#
# Augur runs no data and computes no new numbers — it is a synthesis. "Re-validation" therefore
# means: (1) confirm the submission block reflects the current anchor states, and (2) run the test
# suite (block engaged + citations/anchors well-formed). When an anchor publishes, update
# check_anchors.py + ASSUMPTIONS.md and re-run this; the block lifts only per ASSUMPTIONS.md §3.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python3}"

echo "### Augur re-validation ($(basename "$HERE"))"
echo

echo "## 1. Submission-block gate (check_anchors.py)"
"$PY" "$HERE/check_anchors.py"
BLOCK_RC=$?
echo

echo "## 2. Tests (block engaged + citations/anchors well-formed)"
if "$PY" -m pytest "$HERE/tests" -q; then
  echo "tests: PASS"
else
  echo "tests: FAIL"
  exit 1
fi
echo

if [ "$BLOCK_RC" -ne 0 ]; then
  echo ">>> Augur status: PROVISIONAL, SUBMISSION BLOCKED (expected until Fashion + Minos publish)."
else
  echo ">>> Augur status: anchors published — re-verify CITATIONS.md Tier B, then submit."
fi
# The script itself succeeds (the block being engaged is the correct state, not a failure).
exit 0
