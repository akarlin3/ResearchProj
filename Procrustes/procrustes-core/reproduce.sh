#!/usr/bin/env bash
# Procrustes -- one-touch reproduction. Regenerates EVERY gate number from seeds and
# rebuilds the manuscript, then reports the (author-armed) release posture.
#
# Substrate: Procrustes is fully synthetic -- ground truth is the Lattice digital
# reference object (seed-generated, no data files; no clinical substrate by design).
# "Both substrates of the experiments" are therefore the two evidence layers:
#   (1) the encoded refute-gate tests (pytest: R1/R2/R3 + boundaries + diagnostic),
#   (2) the seeded gate drivers (phase 1/2/3 -> results/*.json -> manuscript numbers).
# Both are regenerated and checked here.
#
# Uses the proteus env (system python3 lacks numpy). Lattice is resolved via the
# initialised submodule, or set LATTICE_PATH.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
if [ ! -x "$PROT" ]; then PROT="$(command -v python3)"; fi
# resolve Lattice (submodule sibling of the project, or env override)
if [ -z "${LATTICE_PATH:-}" ]; then
  for c in "$HERE/../Lattice" "$HERE/../../Lattice" "$HOME/researchProj/Lattice"; do
    [ -f "$c/lattice/__init__.py" ] && export LATTICE_PATH="$(cd "$c" && pwd)" && break
  done
fi
echo "Procrustes reproduce :: python=$PROT  LATTICE_PATH=${LATTICE_PATH:-<import-resolved>}"
cd "$HERE"

echo "== [0/5] install (editable) =="
"$PROT" -m pip install -e . -q || { echo "pip install FAILED"; exit 1; }

echo "== [1/5] refute-gate tests (R1/R2/R3 + boundaries + diagnostic) =="
"$PROT" -m pytest -q || { echo "pytest FAILED"; exit 1; }

echo "== [2/5] GATE B -- apples-to-apples Gauge separation (8 seeds) =="
"$PROT" experiments/run_phase1.py 8 || { echo "phase1 FAILED (GATE B refute fired?)"; exit 1; }

echo "== [3/5] GATE C -- diagnostic reach / honest scope (8 seeds) =="
"$PROT" experiments/run_phase2.py 8 || { echo "phase2 FAILED"; exit 1; }

echo "== [4/5] GATE D -- robustness envelope (8 seeds + 16-seed headline) =="
"$PROT" experiments/run_phase3.py 8 || { echo "phase3 FAILED"; exit 1; }

echo "== [5/5] manuscript: consistency gate -> numbers.tex -> figures -> compile =="
bash paper/build.sh || { echo "manuscript build FAILED"; exit 1; }
echo "== citation verification =="
"$PROT" paper/verify_citations.py || { echo "citation gate FAILED"; exit 1; }

echo "== release posture (author-armed; reporting only) =="
"$PROT" release_gate.py --dry-run || true   # HELD is the expected default; non-zero is informational

echo ""
echo "REPRODUCE: GREEN -- every gate number regenerated; manuscript compiled; citations verified."
