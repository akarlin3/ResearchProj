#!/usr/bin/env bash
# Build the Gauge manuscript and run the GATE 3 consistency check.
# Requires `tectonic` (self-contained LaTeX; fetches packages on first run).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

echo "[1/2] GATE 3 consistency check (numbers trace to gated CP printouts)"
python consistency.py

echo "[2/2] compiling gauge_v3_revised.tex with tectonic"
tectonic gauge_v3_revised.tex
# Canonical committed artifact keeps its established name.
cp -f gauge_v3_revised.pdf gauge_v3_revised_R2.pdf
echo "built: $HERE/gauge_v3_revised_R2.pdf"
