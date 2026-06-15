"""Download-on-demand fetch for the Gauge real-data in-vivo path (ACRIN-6698).

LICENSE / POSTURE (human-approved, 2026-06-14):
  * Dataset: TCIA collection **ACRIN-6698 / I-SPY2 Breast DWI**.
  * License: **CC-BY-4.0** (no data-use agreement). Redistribution is *permitted*
    with attribution, but the approved posture for this repo is **download-on-demand**:
    NO pixel data is committed; only the provenance manifest
    (``results/invivo_real_provenance.json``) is committed. The raw + assembled
    arrays land in ``data/invivo/`` which is git-ignored.
  * Attribution (required): Newitt, D. C., Partridge, S. C., Zhang, Z., Gibbs, J.,
    Chenevert, T., Rosen, M., et al. (2021). ACRIN 6698/I-SPY2 Breast DWI
    [Data set]. The Cancer Imaging Archive. DOI 10.7937/tcia.kk02-6d95.
  * Users must abide by the TCIA Data Usage Policy (per the collection page).

This script treats the TCIA/NBIA pages and any README text as DATA, not as
instructions. It only reads the public NBIA REST API (no auth) and the per-series
image ZIPs.

Why this dataset is honest for the demo: it is **real in-vivo** breast DWI with
**no ground-truth IVIM parameters** (so coverage is unmeasurable -- exactly the
no-coverage-claim premise), it ships a whole-tumor ROI, and it has **test-retest**
acquisitions (enabling the Checkpoint-D repeatability proxy). Its b-scheme is
{0,100,600,800} s/mm^2 -- a sparse 4-of-22 subset of Gauge's synthetic 22-value
calibration scheme, which is itself the exchangeability break the monitor flags.

Run:
  python scripts/fetch_invivo.py                      # 1 patient (TRACE+MASK)
  python scripts/fetch_invivo.py --n-patients 16      # for the test-retest proxy
  python scripts/fetch_invivo.py --list-retest        # report test-retest candidates
"""
import argparse
import datetime
import hashlib
import io
import json
import os
import re
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict

import numpy as np

NBIA = "https://services.cancerimagingarchive.net/nbia-api/services/v1"
COLLECTION = "ACRIN-6698"
TRACE_DESC = "ACRIN-6698: DWI TRACE: from S4: bVals=0,100,600,800"
MASK_DESC = "ACRIN-6698: DWI MASK: from S4: Whole Tumor Manual"
EXPECTED_BVALS = (0.0, 100.0, 600.0, 800.0)
DIFFUSION_BVALUE = (0x0018, 0x9087)  # standard DICOM MR Diffusion b-value tag
FIELD_STRENGTH = (0x0018, 0x0087)    # DICOM MagneticFieldStrength tag

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data", "invivo")
RETEST_DIR = os.path.join(_ROOT, "data", "invivo_retest")
PROVENANCE = os.path.join(_ROOT, "results", "invivo_real_provenance.json")

CITATION = ("Newitt, D. C., Partridge, S. C., Zhang, Z., Gibbs, J., Chenevert, T., "
            "Rosen, M., Bolan, P., Marques, H., Romanoff, J., Cimino, L., Joe, B. N., "
            "Umphrey, H., Ojeda-Fournier, H., Dogan, B., Oh, K. Y., Abe, H., "
            "Drukteinis, J., Esserman, L. J., & Hylton, N. M. (2021). ACRIN 6698/"
            "I-SPY2 Breast DWI [Data set]. The Cancer Imaging Archive.")
DOI = "10.7937/tcia.kk02-6d95"
LICENSE = "CC-BY-4.0"


def _api(endpoint, **params):
    url = f"{NBIA}/{endpoint}?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "gauge-invivo/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def _get_series_zip(series_uid):
    url = f"{NBIA}/getImage?" + urllib.parse.urlencode(
        {"SeriesInstanceUID": series_uid})
    req = urllib.request.Request(url, headers={"User-Agent": "gauge-invivo/0.1"})
    with urllib.request.urlopen(req, timeout=600) as r:
        return r.read()


def _read_dicoms(zip_bytes):
    """Yield pydicom datasets that carry pixel data from a series ZIP."""
    import pydicom
    z = zipfile.ZipFile(io.BytesIO(zip_bytes))
    for name in z.namelist():
        ds = pydicom.dcmread(io.BytesIO(z.read(name)), force=True)
        if hasattr(ds, "Rows") and hasattr(ds, "PixelData"):
            yield ds


def _bvalue(ds):
    if DIFFUSION_BVALUE in ds:
        try:
            return float(ds[DIFFUSION_BVALUE].value)
        except (TypeError, ValueError):
            return None
    return None


def _rescaled(ds):
    arr = ds.pixel_array.astype(float)
    slope = float(getattr(ds, "RescaleSlope", 1.0) or 1.0)
    inter = float(getattr(ds, "RescaleIntercept", 0.0) or 0.0)
    return arr * slope + inter


def assemble_trace(zip_bytes):
    """Assemble a 4-b TRACE series into (signals4d (X,Y,Z,4), z_locs, bvals)."""
    by_b_z = {}
    rows = cols = None
    for ds in _read_dicoms(zip_bytes):
        bv = _bvalue(ds)
        if bv is None:
            continue
        z = round(float(getattr(ds, "SliceLocation",
                               getattr(ds, "InstanceNumber", 0))), 3)
        rows, cols = int(ds.Rows), int(ds.Columns)
        by_b_z[(bv, z)] = _rescaled(ds)
    bvals = sorted({b for (b, _) in by_b_z})
    zlocs = sorted({z for (_, z) in by_b_z})
    if tuple(bvals) != EXPECTED_BVALS:
        raise ValueError(f"b-values {bvals} != expected {EXPECTED_BVALS}")
    vol = np.full((rows, cols, len(zlocs), len(bvals)), np.nan, float)
    for (b, z), img in by_b_z.items():
        vol[:, :, zlocs.index(z), bvals.index(b)] = img
    if np.isnan(vol).any():
        raise ValueError("incomplete (b, slice) grid -- missing images")
    return vol, np.asarray(zlocs, float), np.asarray(bvals, float)


def assemble_mask(zip_bytes, trace_zlocs):
    """Map a tumor-MASK series onto the TRACE slice grid -> (X,Y,Z) bool."""
    rows = cols = None
    by_z = {}
    for ds in _read_dicoms(zip_bytes):
        z = round(float(getattr(ds, "SliceLocation",
                               getattr(ds, "InstanceNumber", 0))), 3)
        rows, cols = int(ds.Rows), int(ds.Columns)
        by_z[z] = ds.pixel_array.astype(float) > 0
    mask = np.zeros((rows, cols, len(trace_zlocs)), bool)
    for z, m in by_z.items():
        k = int(np.argmin(np.abs(trace_zlocs - z)))
        mask[:, :, k] |= m
    return mask


def field_strength_from_zip(zip_bytes):
    """Read DICOM MagneticFieldStrength (0018,0087) from a series ZIP, or None.

    The standardized TRACE series carry the field-strength tag per pixel-bearing
    instance; we take the first non-empty value (consistent within a series).
    """
    for ds in _read_dicoms(zip_bytes):
        if FIELD_STRENGTH in ds:
            try:
                return float(ds[FIELD_STRENGTH].value)
            except (TypeError, ValueError):
                continue
    return None


def _sha256(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()[:16]


def list_trace_exams():
    """Return {patient: [exam dicts]} for every TRACE series, with its MASK."""
    series = _api("getSeries", Collection=COLLECTION, Modality="MR")
    traces = [s for s in series if s.get("SeriesDescription") == TRACE_DESC]
    masks = [s for s in series if s.get("SeriesDescription") == MASK_DESC]
    mask_by_study = {s["StudyInstanceUID"]: s for s in masks}
    by_patient = defaultdict(list)
    for t in traces:
        by_patient[t["PatientID"]].append({
            "patient": t["PatientID"], "study": t["StudyInstanceUID"],
            "date": t.get("StudyDate", ""), "trace_uid": t["SeriesInstanceUID"],
            "mask_uid": (mask_by_study.get(t["StudyInstanceUID"]) or {}).get(
                "SeriesInstanceUID"),
        })
    return by_patient


def fetch_exam(exam, out_root=DATA_DIR):
    """Download + assemble one exam into ``out_root/<patient>/<study8>/``."""
    sub = os.path.join(out_root, exam["patient"], exam["study"][-8:])
    os.makedirs(sub, exist_ok=True)
    vol, zlocs, bvals = assemble_trace(_get_series_zip(exam["trace_uid"]))
    np.save(os.path.join(sub, "signals_4d.npy"), vol.astype(np.float32))
    np.savetxt(os.path.join(sub, "bvals.txt"), bvals, fmt="%g")
    np.save(os.path.join(sub, "slice_locs.npy"), zlocs)
    rec = {**exam, "shape": list(vol.shape), "bvals": bvals.tolist(),
           "sha_signals": _sha256(vol.astype(np.float32))}
    if exam.get("mask_uid"):
        mask = assemble_mask(_get_series_zip(exam["mask_uid"]), zlocs)
        np.save(os.path.join(sub, "tumor_mask.npy"), mask)
        rec["n_tumor_voxels"] = int(mask.sum())
    rec["path"] = os.path.relpath(sub, _ROOT)
    json.dump(rec, open(os.path.join(sub, "meta.json"), "w"), indent=2)
    return rec


# =========================================================================== #
# GENUINE SAME-DAY TEST-RETEST ARM (ACRIN-6698 "coffee-break" TrT0/TrT1).
#
# ACRIN-6698 includes a same-visit scan-RESCAN arm: the patient is scanned, gets
# off the table, is repositioned, and re-scanned in the SAME session. Both
# repeats share one StudyDate / StudyInstanceUID and are tagged "TrT0:" / "TrT1:"
# in the standardized DWI SeriesDescription, each with its OWN whole-tumor mask.
# This is the genuine repeatability arm -- it supersedes the earlier longitudinal
# pairing (timepoints months apart, treatment-confounded), which is NOT same-day.
#
# Pattern (per repeat, same study):
#   TRACE: "ACRIN-6698: TrT{0,1}: DWI TRACE: from S...: bVals=0,100,600,800"
#   MASK:  "ACRIN-6698: TrT{0,1}: DWI MASK: from S...: Whole Tumor Manual"
# StudyDesc encodes the treatment timepoint: ..._T0 (pretreatment) / ..._T1.
# =========================================================================== #
_TRT_RE = re.compile(r"TrT([01])")


def _trt_label(desc):
    m = _TRT_RE.search(desc or "")
    return f"TrT{m.group(1)}" if m else None


def _is_trt_trace(desc):
    return "DWI TRACE" in desc and "bVals=0,100,600,800" in desc


def _is_trt_mask(desc):
    return "DWI MASK" in desc and "Whole Tumor Manual" in desc


def list_test_retest_cells():
    """Return same-day TrT0/TrT1 scan-rescan cells with both TRACE+MASK present.

    A "cell" is one (patient, study) with BOTH TrT0 and TrT1 standardized DWI
    TRACE (b=0,100,600,800) AND a TrT0 and TrT1 whole-tumor MASK. Same StudyDate
    by construction (one StudyInstanceUID). Returns a list of dicts, one per cell.
    """
    series = _api("getSeries", Collection=COLLECTION, Modality="MR")
    cells = defaultdict(lambda: {"trace": {}, "mask": {}, "dates": set(),
                                 "studydesc": set(), "manufacturer": set()})
    for s in series:
        desc = s.get("SeriesDescription", "")
        lab = _trt_label(desc)
        if lab is None:
            continue
        key = (s["PatientID"], s["StudyInstanceUID"])
        c = cells[key]
        c["dates"].add(s.get("StudyDate", ""))
        c["studydesc"].add(s.get("StudyDesc", ""))
        c["manufacturer"].add(s.get("Manufacturer", ""))
        if _is_trt_trace(desc):
            c["trace"][lab] = s["SeriesInstanceUID"]
        elif _is_trt_mask(desc):
            c["mask"][lab] = s["SeriesInstanceUID"]
    out = []
    for (pid, study), c in cells.items():
        if not ({"TrT0", "TrT1"} <= set(c["trace"])
                and {"TrT0", "TrT1"} <= set(c["mask"])):
            continue
        studydesc = sorted(c["studydesc"])[0] if c["studydesc"] else ""
        timepoint = "T0" if studydesc.endswith("_T0") else (
            "T1" if studydesc.endswith("_T1") else "unknown")
        out.append({
            "patient": pid, "study": study,
            "date": sorted(c["dates"])[0] if c["dates"] else "",
            "study_desc": studydesc,
            "timepoint": timepoint,
            "pretreatment": timepoint == "T0",
            "manufacturer": sorted(c["manufacturer"])[0] if c["manufacturer"] else "",
            "trace": dict(c["trace"]), "mask": dict(c["mask"]),
        })
    return sorted(out, key=lambda d: (d["patient"], d["study"]))


def fetch_retest_repeat(cell, label, out_root=RETEST_DIR):
    """Download + assemble ONE repeat (TrT0 or TrT1) of a same-day cell.

    Writes ``out_root/<patient>/<TrT0|TrT1>/`` with signals_4d.npy, tumor_mask.npy,
    bvals.txt, slice_locs.npy and a per-repeat meta.json carrying full provenance
    (incl. field strength + treatment timepoint). Returns the provenance record.
    Raises on a failed bvals assertion / empty mask so the caller can skip+log.
    """
    sub = os.path.join(out_root, cell["patient"], label)
    os.makedirs(sub, exist_ok=True)
    trace_zip = _get_series_zip(cell["trace"][label])
    vol, zlocs, bvals = assemble_trace(trace_zip)            # asserts bvals
    field_t = field_strength_from_zip(trace_zip)
    mask = assemble_mask(_get_series_zip(cell["mask"][label]), zlocs)
    if int(mask.sum()) == 0:
        raise ValueError(f"empty tumor mask for {cell['patient']} {label}")
    np.save(os.path.join(sub, "signals_4d.npy"), vol.astype(np.float32))
    np.save(os.path.join(sub, "tumor_mask.npy"), mask)
    np.save(os.path.join(sub, "slice_locs.npy"), zlocs)
    np.savetxt(os.path.join(sub, "bvals.txt"), bvals, fmt="%g")
    rec = {
        "patient": cell["patient"], "study": cell["study"],
        "date": cell["date"], "study_desc": cell["study_desc"],
        "timepoint": cell["timepoint"], "pretreatment": cell["pretreatment"],
        "repeat": label,
        "trace_uid": cell["trace"][label], "mask_uid": cell["mask"][label],
        "manufacturer": cell["manufacturer"],
        "field_strength_T": field_t,
        "shape": list(vol.shape), "bvals": bvals.tolist(),
        "sha_signals": _sha256(vol.astype(np.float32)),
        "n_tumor_voxels": int(mask.sum()),
        "path": os.path.relpath(sub, _ROOT),
    }
    json.dump(rec, open(os.path.join(sub, "meta.json"), "w"), indent=2)
    return rec


def fetch_test_retest_arm(out_root=RETEST_DIR, max_patients=None):
    """Fetch the genuine same-day TrT0/TrT1 scan-rescan arm of ACRIN-6698.

    For every complete same-day cell, downloads BOTH repeats (TRACE+MASK) into
    separate exam dirs. Skips a patient if EITHER repeat fails assembly (bvals
    assertion / empty mask) so each kept patient has a clean pair. Returns
    ``(records, skips)`` where records is a flat list of per-repeat provenance.
    """
    cells = list_test_retest_cells()
    # one cell per patient (one patient has 2 complete cells; take the first)
    by_patient = {}
    for c in cells:
        by_patient.setdefault(c["patient"], c)
    patients = sorted(by_patient)
    if max_patients:
        patients = patients[:max_patients]
    print(f"[retest] same-day TrT0/TrT1 cells: {len(cells)} "
          f"across {len(by_patient)} patients; fetching {len(patients)}")
    records, skips = [], []
    for i, pid in enumerate(patients, 1):
        cell = by_patient[pid]
        print(f"[retest {i}/{len(patients)}] {pid} study ...{cell['study'][-8:]} "
              f"({cell['timepoint']}, {cell['manufacturer']})")
        pair = []
        ok = True
        for label in ("TrT0", "TrT1"):
            try:
                pair.append(fetch_retest_repeat(cell, label, out_root))
            except Exception as e:                          # noqa: BLE001
                skips.append({"patient": pid, "repeat": label, "reason": str(e)})
                print(f"    SKIP {pid} {label}: {e}")
                ok = False
                break
        if ok:
            records.extend(pair)
    return records, skips


def write_retest_manifest(records, skips, out_root=RETEST_DIR):
    """Write a standalone manifest for the fetched same-day arm (data/-side).

    The committed provenance ``retest`` block is written by gauge.invivo; this
    JSON is a fetch-side ledger that also gets merged into the dataset block.
    """
    manifest = os.path.join(out_root, "manifest.json")
    os.makedirs(out_root, exist_ok=True)
    n_pairs = len(records) // 2
    payload = {
        "arm": "ACRIN-6698 same-day test-retest (TrT0/TrT1)",
        "description": ("same-visit scan-rescan ('coffee-break'): both repeats "
                        "share one StudyDate/StudyInstanceUID, each with its own "
                        "whole-tumor mask. Supersedes the earlier longitudinal "
                        "n=11 pairing (timepoints months apart)."),
        "n_pairs": n_pairs, "n_repeats": len(records),
        "b_values_s_per_mm2": list(EXPECTED_BVALS),
        "license": LICENSE, "doi": DOI, "citation": CITATION,
        "nbia_api": NBIA, "collection": COLLECTION,
        "access_date_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "repeats": records, "skips": skips,
    }
    json.dump(payload, open(manifest, "w"), indent=2)
    return manifest


def write_provenance(records):
    os.makedirs(os.path.dirname(PROVENANCE), exist_ok=True)
    prov = {
        "dataset": {
            "name": "ACRIN-6698 / I-SPY2 Breast DWI", "source": "TCIA",
            "collection_page": "https://www.cancerimagingarchive.net/collection/acrin-6698/",
            "doi": DOI, "license": LICENSE, "citation": CITATION,
            "modality": "MR DWI (real in-vivo, breast)",
            "ground_truth_ivim": False,
            "b_values_s_per_mm2": list(EXPECTED_BVALS),
            "posture": "download-on-demand; pixel data NOT committed (data/ is git-ignored)",
            "fetched_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "nbia_api": NBIA,
            "exams": records,
        }
    }
    # preserve any existing "run" block written by gauge.invivo
    if os.path.exists(PROVENANCE):
        try:
            old = json.load(open(PROVENANCE))
            if "run" in old:
                prov["run"] = old["run"]
        except (json.JSONDecodeError, OSError):
            pass
    json.dump(prov, open(PROVENANCE, "w"), indent=2)
    return PROVENANCE


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-patients", type=int, default=1)
    ap.add_argument("--out", default=DATA_DIR)
    ap.add_argument("--list-retest", action="store_true",
                    help="report test-retest candidates (>=2 TRACE exams) and exit")
    ap.add_argument("--retest-arm", action="store_true",
                    help="fetch the genuine same-day TrT0/TrT1 scan-rescan arm "
                         "into data/invivo_retest/ (the repeatability arm)")
    ap.add_argument("--list-retest-arm", action="store_true",
                    help="report same-day TrT0/TrT1 cells (count + split) and exit")
    ap.add_argument("--max-patients", type=int, default=None,
                    help="cap the number of retest patients (debug)")
    args = ap.parse_args()

    if args.list_retest_arm:
        cells = list_test_retest_cells()
        pats = sorted({c["patient"] for c in cells})
        t0 = sum(1 for c in cells if c["timepoint"] == "T0")
        t1 = sum(1 for c in cells if c["timepoint"] == "T1")
        mfg = defaultdict(int)
        for c in cells:
            mfg[c["manufacturer"]] += 1
        print(f"same-day TrT0/TrT1 complete cells: {len(cells)} "
              f"across {len(pats)} patients")
        print(f"  timepoint split: T0 (pretreatment)={t0}  T1={t1}")
        print(f"  manufacturer split: {dict(mfg)}")
        return 0

    if args.retest_arm:
        out_root = args.out if args.out != DATA_DIR else RETEST_DIR
        records, skips = fetch_test_retest_arm(out_root, max_patients=args.max_patients)
        manifest = write_retest_manifest(records, skips, out_root)
        n_pairs = len(records) // 2
        print("\n" + "=" * 76)
        print(f"RETEST FETCH SUMMARY -- {COLLECTION} same-day TrT0/TrT1 "
              f"(license {LICENSE}, DOI {DOI})")
        print("=" * 76)
        print(f"  usable same-day pairs: {n_pairs}   skipped: {len(skips)}")
        for s in skips:
            print(f"    SKIP {s['patient']} {s['repeat']}: {s['reason']}")
        print(f"  pixel data -> {os.path.relpath(out_root, _ROOT)}/ (git-ignored)")
        print(f"  manifest   -> {os.path.relpath(manifest, _ROOT)} (git-ignored)")
        print("=" * 76)
        return 0

    by_patient = list_trace_exams()
    if args.list_retest:
        rt = {p: ex for p, ex in by_patient.items() if len(ex) >= 2}
        print(f"patients with >=2 TRACE exams (test-retest candidates): {len(rt)}")
        for p, ex in list(rt.items())[:10]:
            print(f"  {p}: {len(ex)} exams  dates={[e['date'] for e in ex]}")
        return 0

    patients = sorted(by_patient)[:args.n_patients]
    records = []
    for p in patients:
        for exam in by_patient[p][:2]:  # up to test+retest per patient
            print(f"[fetch] {p} study ...{exam['study'][-8:]} "
                  f"(mask={'yes' if exam.get('mask_uid') else 'no'})")
            records.append(fetch_exam(exam, args.out))
    prov = write_provenance(records)

    print("\n" + "=" * 76)
    print(f"FETCH SUMMARY -- {COLLECTION} (license {LICENSE}, DOI {DOI})")
    print("=" * 76)
    for r in records:
        print(f"  {r['patient']} | shape {tuple(r['shape'])} | "
              f"b={r['bvals']} | tumor voxels={r.get('n_tumor_voxels','-')}")
    n_rt = sum(1 for ex in by_patient.values() if len(ex) >= 2)
    print(f"  test-retest available in collection: YES ({n_rt} patients with >=2 exams)")
    print(f"  pixel data -> {os.path.relpath(args.out, _ROOT)}/ (git-ignored)")
    print(f"  provenance (committed) -> {os.path.relpath(prov, _ROOT)}")
    print("=" * 76)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
