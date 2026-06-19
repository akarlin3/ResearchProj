#!/usr/bin/env bash
# Build the Minos future/ manuscript. Runs the CP4 consistency gate first (regenerates
# numbers.tex from the seeded results and checks minos.tex traceability), then compiles.
#
# Prefers `tectonic` (self-contained, Overleaf-like). Falls back to pdflatex x2. If neither
# is installed, the consistency gate still runs (the CP4 numbers gate does not need LaTeX).
#
# Usage: bash build.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"

echo "== CP4 consistency gate =="
"$PROT" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling with tectonic =="
  tectonic minos.tex && { echo "built minos.pdf (tectonic)"; exit 0; }
  echo "tectonic failed"; exit 1
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling with pdflatex (x2) =="
  pdflatex -interaction=nonstopmode -halt-on-error minos.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error minos.tex && { echo "built minos.pdf (pdflatex)"; exit 0; }
  echo "pdflatex failed (see minos.log)"; exit 1
else
  echo "== no LaTeX engine found; consistency gate passed, skipping PDF compile =="
  echo "   (on Overleaf, the ebgaramond+microtype preamble builds as-is)"
  exit 0
fi
