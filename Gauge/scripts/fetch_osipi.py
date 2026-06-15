#!/usr/bin/env python
"""Download-on-demand fetcher for the OSIPI TF2.4 IVIM reference phantom (DRO).

Task 1 (external-phantom replication) validates Gauge's conformal/CQR intervals
against an EXTERNAL, community-standard **synthetic** reference object whose
ground-truth (D, D*, f) we did NOT generate -- the OSIPI ISMRM Open Science
Initiative for Perfusion Imaging IVIM code collection. This neutralises the
"validated only on the authors' own forward model" critique: the parameter maps
(the digital reference object, DRO) come from OSIPI, not from us.

Data posture (matches the ACRIN-6698 pattern in scripts/fetch_invivo.py):
  * download-on-demand; raw phantom data lives under data/ (git-ignored).
  * commit ONLY this provenance manifest (URLs, DOI, versions, checksums, access
    date, license) -- never the raw arrays.
  * licenses: OSIPI code repo is Apache-2.0; the Zenodo DATA record is CC-BY-4.0.
    Both are recorded below; attribution preserved.

The DRO itself: ``Utilities/DRO.npy`` inside the OSIPI Zenodo data record is an
object array of 5000 voxels, each a dict ``{D, f, Dp, S0, bvals, signals}`` with
D, Dp(=D*) in mm^2/s (absolute D*, exact-match to our convention) and a fixed
sparse 7-b acquisition, pre-noised at ~SNR 80. This is the external ground truth.

Honesty: the OSIPI phantom is an *independent synthetic* reference (a DRO), NOT
in vivo. No in-vivo coverage claim is made anywhere in Task 1.

Run:  python scripts/fetch_osipi.py
"""
import argparse
import datetime
import hashlib
import json
import os
import urllib.request
import zipfile

import numpy as np

# --------------------------------------------------------------------------- #
# Pinned sources (immutable): OSIPI code repo tag v0.1.0 + Zenodo data record.
# --------------------------------------------------------------------------- #
REPO = "https://github.com/OSIPI/TF2.4_IVIM-MRI_CodeCollection"
REPO_TAG = "v0.1.0"
REPO_COMMIT = "23577722f9849e5229acaad8b61abc89de6f9542"
ZENODO_RECORD = "14605039"
ZENODO_DOI = "10.5281/zenodo.14605039"
ZIP_URL = (f"https://zenodo.org/records/{ZENODO_RECORD}/files/"
           "OSIPI_TF24_data_phantoms.zip?download=1")
ZIP_NAME = "OSIPI_TF24_data_phantoms.zip"
ZIP_MD5 = "e7b3fe1d811a7a45c5aaf6c604c82793"          # from Zenodo file API
DRO_MEMBER = "Utilities/DRO.npy"                        # the DRO inside the zip

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DATA_DIR = os.path.join(_ROOT, "data", "osipi")
_EXTRACT_DIR = os.path.join(_DATA_DIR, "extracted")
_RESULTS_DIR = os.path.join(_ROOT, "results")
_PROVENANCE = os.path.join(_RESULTS_DIR, "osipi_provenance.json")
DRO_PATH = os.path.join(_EXTRACT_DIR, "DRO.npy")


def _md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def _sha256_arr(arr):
    """Truncated SHA-256 of an array's raw bytes (provenance, matches fetch_invivo)."""
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()[:16]


def download_zip(force=False):
    os.makedirs(_DATA_DIR, exist_ok=True)
    zip_path = os.path.join(_DATA_DIR, ZIP_NAME)
    if force or not os.path.exists(zip_path):
        print(f"[osipi] downloading {ZIP_URL} -> {zip_path}")
        urllib.request.urlretrieve(ZIP_URL, zip_path)
    got = _md5(zip_path)
    if got != ZIP_MD5:
        raise RuntimeError(
            f"OSIPI zip md5 mismatch: got {got}, expected {ZIP_MD5}. "
            "Delete data/osipi and re-run.")
    print(f"[osipi] zip md5 OK ({got})")
    return zip_path


def extract_dro(zip_path):
    os.makedirs(_EXTRACT_DIR, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(DRO_MEMBER) as src, open(DRO_PATH, "wb") as dst:
            dst.write(src.read())
    print(f"[osipi] extracted {DRO_MEMBER} -> {DRO_PATH}")
    return DRO_PATH


def verify_and_summarize(dro_path):
    """Load the DRO, sanity-check shapes/units/ranges, return a provenance summary."""
    dro = np.load(dro_path, allow_pickle=True)
    assert dro.ndim == 1 and dro.size > 0, "DRO.npy is not a 1-D object array"
    keys = set(dro[0].keys())
    assert {"D", "f", "Dp", "S0", "bvals", "signals"} <= keys, \
        f"unexpected DRO entry keys: {keys}"
    D = np.array([float(e["D"]) for e in dro])
    Dstar = np.array([float(e["Dp"]) for e in dro])
    f = np.array([float(e["f"]) for e in dro])
    bvals = np.asarray(dro[0]["bvals"], float)
    # all entries share the fixed b-scheme
    same_b = all(np.array_equal(np.asarray(e["bvals"], float), bvals) for e in dro)
    # unit sanity: D, D* in mm^2/s; D* physically >> D and within fit bounds.
    assert 1e-4 < D.min() and D.max() < 5e-3, f"D out of range: {D.min()},{D.max()}"
    assert 1e-2 < Dstar.min() and Dstar.max() < 3e-1, \
        f"D* out of range: {Dstar.min()},{Dstar.max()}"
    assert 0.0 <= f.min() and f.max() <= 1.0, f"f out of range: {f.min()},{f.max()}"
    summary = {
        "n_voxels": int(dro.size),
        "entry_keys": sorted(keys),
        "sha256_16_D": _sha256_arr(D),
        "sha256_16_Dstar": _sha256_arr(Dstar),
        "sha256_16_f": _sha256_arr(f),
        "native_b_values_s_per_mm2": bvals.tolist(),
        "native_b_scheme_constant_across_voxels": bool(same_b),
        "units": "D, D* in mm^2/s (absolute D*); f dimensionless -- matches Gauge",
        "ground_truth_ranges": {
            "D_mm2_s": [float(D.min()), float(D.max())],
            "Dstar_mm2_s": [float(Dstar.min()), float(Dstar.max())],
            "f": [float(f.min()), float(f.max())],
        },
        "note_distribution_shift_vs_gauge": (
            "OSIPI true D* spans ~0.05-0.20 mm^2/s, shifted HIGHER than Gauge's "
            "calibration prior (0.010-0.100); a genuine out-of-distribution test."),
    }
    print(f"[osipi] DRO verified: {summary['n_voxels']} voxels, "
          f"D*∈[{Dstar.min():.3f},{Dstar.max():.3f}] mm^2/s")
    return summary


def write_provenance(summary):
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    prov = {
        "dataset": {
            "name": "OSIPI TF2.4 IVIM-MRI Code Collection -- reference phantom (DRO)",
            "kind": "EXTERNAL community-standard SYNTHETIC digital reference object "
                    "(DRO) -- NOT in vivo; ground-truth (D, D*, f) NOT generated by us",
            "code_repo": REPO,
            "code_repo_tag": REPO_TAG,
            "code_repo_commit": REPO_COMMIT,
            "code_license": "Apache-2.0",
            "data_record_doi": ZENODO_DOI,
            "data_record_url": f"https://zenodo.org/records/{ZENODO_RECORD}",
            "data_license": "CC-BY-4.0",
            "zip_name": ZIP_NAME,
            "zip_md5": ZIP_MD5,
            "dro_member": DRO_MEMBER,
            "forward_model": "bi-exponential S/S0 = (1-f)exp(-bD) + f exp(-bD*) "
                             "(OSIPI consensus model; identical to Gauge's)",
            "citation_software": (
                "Jalnefjord O, Rashid IA, Kuppens D, van der Thiel MM, van Houdt PJ, "
                "Voorter PHM, Peterson ET, Gurney-Champion OJ. OSIPI TF2.4 IVIM-MRI "
                "Code Collection (v0.1.0); 2026."),
            "citation_data": (
                "Gurney-Champion O, Rashid I, van der Thiel M, Kuppens D, Voorter P, "
                "van Houdt P, Peterson E, Jalnefjord O. Data to "
                "github.com/OSIPI/TF2.4_IVIM-MRI_CodeCollection. Zenodo; "
                f"doi:{ZENODO_DOI}."),
            "citation_paper_pending": (
                "Jalnefjord et al., 'An open-source code repository for intravoxel "
                "incoherent motion analysis: OSIPI', Magn Reson Med (2026, in press; "
                "DOI pending)."),
            "posture": "download-on-demand; raw phantom NOT committed (data/ git-ignored)",
            "no_in_vivo_coverage_claim": True,
            "fetched_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            **summary,
        }
    }
    # read-modify-write: preserve any analysis-written blocks (e.g. "run").
    if os.path.exists(_PROVENANCE):
        try:
            old = json.load(open(_PROVENANCE))
            old.update(prov)
            prov = old
        except Exception:
            pass
    with open(_PROVENANCE, "w") as fh:
        json.dump(prov, fh, indent=1)
    print(f"[osipi] wrote provenance -> {_PROVENANCE}")
    return _PROVENANCE


def fetch(force=False):
    zip_path = download_zip(force=force)
    dro_path = extract_dro(zip_path)
    summary = verify_and_summarize(dro_path)
    write_provenance(summary)
    return dro_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    args = ap.parse_args()
    fetch(force=args.force)
