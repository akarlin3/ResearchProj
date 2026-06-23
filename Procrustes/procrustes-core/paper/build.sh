#!/usr/bin/env bash
# Build the Procrustes manuscript. Runs the consistency gate (which regenerates numbers.tex from
# the seeded gate results and checks macro coverage + locked-scope invariants), then compiles
# with tectonic (fallback pdflatex x2 + bibtex). Mirrors Gauge/Augur/Matrix paper/build.sh.
#
# Assumes the gate results exist (results/phase{1,2,3}_*.json). To regenerate them too, run the
# repo-root reproduce.sh instead. Uses the proteus env (system python3 lacks numpy).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORE="$(cd "$HERE/.." && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
if [ ! -x "$PROT" ]; then PROT="$(command -v python3)"; fi

echo "== consistency gate (regenerates numbers.tex; checks macros + locked-scope invariants) =="
"$PROT" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

echo "== figures (regenerate from seeded results; non-fatal -- committed PDFs are the fallback) =="
"$PROT" "$HERE/figures/make_figures.py" 2>/dev/null \
  || echo "  [note] figure regeneration skipped/failed; using committed figures/*.pdf"

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling procrustes.tex with tectonic =="
  tectonic procrustes.tex && { echo "built: $HERE/procrustes.pdf"; exit 0; }
  echo "tectonic failed"; exit 1
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling with pdflatex (x2 + bibtex) =="
  pdflatex -interaction=nonstopmode -halt-on-error procrustes.tex >/dev/null 2>&1
  bibtex procrustes >/dev/null 2>&1 || true
  pdflatex -interaction=nonstopmode -halt-on-error procrustes.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error procrustes.tex && { echo "built procrustes.pdf"; exit 0; }
  echo "pdflatex failed (see procrustes.log)"; exit 1
else
  echo "== no LaTeX engine found; consistency gate passed, skipping PDF compile =="
  exit 0
fi
