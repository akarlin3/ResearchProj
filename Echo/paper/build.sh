#!/usr/bin/env bash
# Build the Echo manuscript with the CP4 consistency gate first. Mirrors
# Minos/future/paper/build.sh. No-op-with-message until echo.tex exists (PASS-only, CP4).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python}"

if [ ! -f "$HERE/echo.tex" ]; then
  echo "echo.tex not present — manuscript is built only on a CP3 PASS (see README.md)."
  exit 0
fi

echo ">>> CP4 consistency gate"
"$PY" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  tectonic echo.tex
else
  pdflatex -interaction=nonstopmode echo.tex && pdflatex -interaction=nonstopmode echo.tex
fi
echo ">>> built echo.pdf"
