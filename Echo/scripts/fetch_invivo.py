#!/usr/bin/env python
"""Download-on-demand fetch for Echo's real test--retest signals (ACRIN-6698).

POSTURE (identical to Gauge -- this script is modelled on ``Gauge/scripts/fetch_invivo.py``,
the sanctioned data-handling template):
  * Dataset: TCIA collection **ACRIN-6698 / I-SPY2 Breast DWI**, same-day test--retest arm
    (TrT0/TrT1). License **CC-BY-4.0**, DOI 10.7937/tcia.kk02-6d95.
  * **No pixel data is committed.** Raw arrays live under ``data/`` (git-ignored). Only the
    provenance manifest ``results/invivo_provenance.json`` is committed.
  * Echo redistributes nothing; it reads the public dataset on demand.

Echo needs only **ROI-mean test--retest signal pairs** (one mean signal vector per b per
repeat per tumor) -- far less than Gauge's full per-voxel pipeline. The cheapest honest
route, and the one this script uses, is to **reuse Gauge's already-fetched, validated
test--retest arrays** (read-only) and reduce them to ROI means. Gauge's fetch did the heavy
DICOM/NBIA work; Echo computes its OWN distinct statistic downstream.

Modes:
  python scripts/fetch_invivo.py --from-gauge      # reuse Gauge's fetched data/ (default)
  python scripts/fetch_invivo.py --check           # report whether suitable data is present

If Gauge's data is not present, run Gauge's fetch first (download-on-demand):
  python ../Gauge/scripts/fetch_invivo.py --n-patients 80
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from echo_repeat import provenance  # noqa: E402

DATA_DIR = ROOT / "data" / "invivo_retest"
MANIFEST = ROOT / "results" / "invivo_provenance.json"


def _gauge_root() -> Path:
    from echo_repeat import _paths
    return _paths.GAUGE


def _find_gauge_retest_arrays():
    """Locate Gauge's fetched test--retest arrays (read-only). Returns dir or None."""
    g = _gauge_root()
    for cand in (g / "data" / "invivo_retest", g / "data" / "invivo"):
        if cand.exists() and any(cand.iterdir()):
            return cand
    return None


def _roi_mean_signal(arr: np.ndarray) -> np.ndarray:
    """Reduce a per-voxel (or volume) signal array to a b-indexed ROI-mean vector."""
    arr = np.asarray(arr, float)
    if arr.ndim == 1:
        return arr
    # collapse all spatial axes, keep the b axis (assumed last and length 4 for ACRIN)
    b_axis = int(np.argmin([abs(s - 4) for s in arr.shape]))
    spatial = tuple(i for i in range(arr.ndim) if i != b_axis)
    return np.nanmean(arr, axis=spatial)


def from_gauge() -> int:
    src = _find_gauge_retest_arrays()
    if src is None:
        print("NO suitable Gauge test-retest arrays found (data/ is download-on-demand).")
        print("Run Gauge's fetch first:  python ../Gauge/scripts/fetch_invivo.py "
              "--n-patients 80")
        return 2
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    exams = []
    n = 0
    for pair_dir in sorted(p for p in src.iterdir() if p.is_dir()):
        a = sorted(pair_dir.glob("*TrT0*npy")) + sorted(pair_dir.glob("*_a*npy"))
        b = sorted(pair_dir.glob("*TrT1*npy")) + sorted(pair_dir.glob("*_b*npy"))
        if not (a and b):
            continue
        sig_a = _roi_mean_signal(np.load(a[0]))
        sig_b = _roi_mean_signal(np.load(b[0]))
        rec = {"tumor": pair_dir.name, "bvals": provenance.ACRIN6698["b_values_s_per_mm2"],
               "signal_a": sig_a.tolist(), "signal_b": sig_b.tolist()}
        (DATA_DIR / f"{pair_dir.name}.json").write_text(json.dumps(rec))
        exams.append({"tumor": pair_dir.name,
                      "sha_a": provenance.sha16(sig_a.tobytes()),
                      "sha_b": provenance.sha16(sig_b.tobytes())})
        n += 1
    provenance.write_manifest(MANIFEST, provenance.ACRIN6698, exams,
                              extra={"derived_from": "Gauge in-vivo fetch (read-only reuse)",
                                     "reduction": "whole-tumor ROI mean per b per repeat"})
    print(f"wrote {n} ROI-mean test-retest pairs to {DATA_DIR} (pixel data NOT committed)")
    print(f"wrote provenance manifest {MANIFEST}")
    return 0 if n > 0 else 2


def check() -> int:
    src = _find_gauge_retest_arrays()
    have_echo = DATA_DIR.exists() and any(DATA_DIR.glob("*.json"))
    print(json.dumps({
        "gauge_arrays_present": src is not None,
        "gauge_arrays_dir": str(src) if src else None,
        "echo_roi_signals_present": bool(have_echo),
        "echo_data_dir": str(DATA_DIR),
        "manifest": str(MANIFEST) if MANIFEST.exists() else None,
        "dataset": provenance.ACRIN6698["name"],
        "license": provenance.ACRIN6698["license"],
        "doi": provenance.ACRIN6698["doi"],
    }, indent=2))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-gauge", action="store_true",
                    help="reuse Gauge's fetched test-retest arrays (read-only)")
    ap.add_argument("--check", action="store_true", help="report data availability")
    args = ap.parse_args()
    if args.check:
        return check()
    return from_gauge()


if __name__ == "__main__":
    raise SystemExit(main())
