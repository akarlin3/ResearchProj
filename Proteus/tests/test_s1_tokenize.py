"""S1 tokenize-to-3Di — positive-output test (Foldseek-native ProstT5 backend).

Asserts a POSITIVE artifact: every S0 representative gets a NON-EMPTY 3Di string
whose length equals its amino-acid length, and the output is a valid Foldseek
query DB (the ``_ss`` 3Di records + index) that S2 can consume directly.

Skips cleanly if Foldseek is absent, if this Foldseek build lacks
``--prostt5-model``, or if ProstT5 weights are not locally available (the test
does NOT download the ~2.4 GB weights — point it at them via paths.prostt5_weights
or PROTEUS_PROSTT5_MODEL).
"""
from __future__ import annotations

import os
import shutil

import pytest

from proteus.s1_tokenize import (
    foldseek_supports_prostt5,
    parse_fasta,
    resolve_prostt5_weights,
    tokenize,
)
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MINI = os.path.join(REPO, "tests", "data", "mini_corpus.fasta")


def _representatives(tmp_path):
    """Build S0 representatives to feed S1; skip if mmseqs is unavailable."""
    if shutil.which("mmseqs") is None:
        pytest.skip("mmseqs2 not installed — S1 needs S0 representatives")
    from proteus.s0_dereplicate import dereplicate  # noqa: PLC0415
    cfg = load_config()
    rep = tmp_path / "s0_representatives.fasta"
    tsv = tmp_path / "s0_clusters.tsv"
    dereplicate(MINI, str(rep), str(tsv), cfg, tmp_dir=str(tmp_path / "s0work"))
    return rep


def _weights(cfg):
    w = resolve_prostt5_weights(cfg, allow_download=False)
    if w is None:
        pytest.skip("ProstT5 weights not local — set paths.prostt5_weights or "
                    "PROTEUS_PROSTT5_MODEL (test does not download 2.4 GB)")
    return w


def _guard():
    if shutil.which("foldseek") is None:
        pytest.skip("foldseek not installed")
    if not foldseek_supports_prostt5():
        pytest.skip("this Foldseek build lacks --prostt5-model (no native ProstT5)")
    if not os.path.exists(MINI):
        pytest.skip("mini_corpus.fasta missing — run tests/data/make_mini_corpus.py")


def test_s1_every_representative_gets_length_matched_3di(tmp_path):
    _guard()
    cfg = load_config()
    _weights(cfg)  # skip early if weights absent
    rep = _representatives(tmp_path)
    out = tmp_path / "s1_3di"

    summary = tokenize(str(rep), str(out), cfg)

    assert summary["backend"] == "foldseek-native"
    # every representative produced a non-empty 3Di string
    assert summary["n_records"] > 0
    assert summary["n_nonempty"] == summary["n_records"], "some 3Di strings are empty"
    # 3Di length == amino-acid length for EVERY record (one token per residue)
    assert summary["all_lengths_match"] is True
    for r in summary["records"]:
        assert r["len_3di"] == r["length"] > 0, f"{r['id']}: 3Di/seq length mismatch"

    # the inspectable 3Di FASTA exists and has one record per representative,
    # over the 3Di alphabet (uppercase letters)
    threedi = dict(parse_fasta(summary["threedi_fasta"]))
    rep_ids = [rid for rid, _ in parse_fasta(str(rep))]
    assert set(threedi) == set(rep_ids)
    for rid, td in threedi.items():
        assert td and td.upper() == td and td.isalpha(), f"{rid}: not a 3Di string"


def test_s1_emits_valid_foldseek_query_db_for_s2(tmp_path):
    """S2 consumes a Foldseek query DB directly; assert createdb produced one
    (the 3Di ``_ss`` sub-DB plus its index)."""
    _guard()
    cfg = load_config()
    _weights(cfg)
    rep = _representatives(tmp_path)
    out = tmp_path / "s1_3di"

    summary = tokenize(str(rep), str(out), cfg)
    qdb = summary["querydb"]
    assert qdb is not None, "Foldseek-native backend must emit a query DB"
    # The amino-acid DB and the 3Di structural DB + indices S2 will search against.
    for suffix in ("", ".index", ".dbtype", "_ss", "_ss.index"):
        path = qdb + suffix
        assert os.path.exists(path), f"missing Foldseek query DB file: {path}"
