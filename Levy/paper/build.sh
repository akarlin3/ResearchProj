#!/usr/bin/env bash
# Build the Levy manuscript: run the traceability gate, then compile.
#
#   1. python consistency.py   regenerate numbers.tex from seeded results + verify
#   2. tectonic levy.tex       compile to levy.pdf (fallback: pdflatex x2)
#
# The consistency gate MUST pass before compiling; numbers.tex is auto-generated and never
# hand-edited.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
PY="$PROT"; [ -x "$PY" ] || PY="${PYTHON:-python3}"

echo "== consistency gate (numbers trace to seeded results) =="
"$PY" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling levy.tex with tectonic =="
  tectonic levy.tex && { echo "built levy.pdf (tectonic)"; exit 0; }
  echo "tectonic failed"; exit 1
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling levy.tex with pdflatex (x2) =="
  pdflatex -interaction=nonstopmode -halt-on-error levy.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error levy.tex && { echo "built levy.pdf (pdflatex)"; exit 0; }
  echo "pdflatex failed (see levy.log)"; exit 1
else
  echo "== no LaTeX engine found; consistency gate passed, skipping PDF compile =="
  exit 0
fi
