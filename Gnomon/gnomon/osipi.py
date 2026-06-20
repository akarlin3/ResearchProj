"""OSIPI open abdomen scan: download-on-demand + provenance.  [CP2 — spec below]

Target T1 (the 54.7% NLLS D* railing rate) is a REAL-data number. The data is the
OSIPI TF2.4 abdomen acquisition in the open Zenodo archive **14605039**
(CC-licensed). Mirroring Lattice's ``osipi.py`` / Lethe's ``fetch_invivo.py``:

* fetched **on demand** (``zenodo_get`` / HTTP) into a **gitignored** ``download/``
  dir -- never redistributed in-tree;
* a provenance manifest (record id, file hashes, retrieval date) is written;
* importing this module touches no network.

CP2 deliverables: ``fetch_abdomen()`` (download + verify), ``load_high_snr_roi()``
(b-values, signals, the ROI/SNR selection that defines the 1618-voxel set -- a
selection Gnomon documents explicitly, since Fashion's prose under-specifies it).
Needs the optional ``[data]`` extra (requests / zenodo_get).
"""
from __future__ import annotations

ZENODO_RECORD = "14605039"
LICENSE = "CC (OSIPI TF2.4 open data)"


def fetch_abdomen(*args, **kwargs):  # [CP2]
    raise NotImplementedError(
        "Gnomon CP2: download-on-demand OSIPI abdomen scan (gitignored + provenance)."
    )
