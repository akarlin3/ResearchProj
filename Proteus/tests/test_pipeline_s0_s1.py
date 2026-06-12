"""CP4 — wire S0 -> S1 on the mini corpus and confirm the artifacts S2 needs.

This is the end-to-end front-of-pipeline check: dereplicate the mini corpus (S0),
then tokenize the survivors to 3Di (S1), and assert the hand-off contract for S2
holds — a representative FASTA and a Foldseek query DB (with 3Di records) whose
record count matches the S0 representative count. S2 itself is NOT implemented here.

Skips cleanly if MMseqs2/Foldseek/ProstT5 weights are unavailable.
"""
from __future__ import annotations

import os
import shutil

import pytest

from proteus.s0_dereplicate import dereplicate, parse_fasta
from proteus.s1_tokenize import (
    foldseek_supports_prostt5,
    resolve_prostt5_weights,
    tokenize,
)
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MINI = os.path.join(REPO, "tests", "data", "mini_corpus.fasta")

EXPECTED_REPS = 8  # documented known answer for the mini corpus


def _guard(cfg):
    if shutil.which("mmseqs") is None:
        pytest.skip("mmseqs2 not installed")
    if shutil.which("foldseek") is None:
        pytest.skip("foldseek not installed")
    if not foldseek_supports_prostt5():
        pytest.skip("Foldseek build lacks --prostt5-model")
    if not os.path.exists(MINI):
        pytest.skip("mini_corpus.fasta missing — run tests/data/make_mini_corpus.py")
    if resolve_prostt5_weights(cfg, allow_download=False) is None:
        pytest.skip("ProstT5 weights not local (set paths.prostt5_weights / "
                    "PROTEUS_PROSTT5_MODEL)")


def test_s0_to_s1_chain_produces_s2_inputs(tmp_path):
    cfg = load_config()
    _guard(cfg)

    # --- S0 -------------------------------------------------------------- #
    rep = tmp_path / "s0_representatives.fasta"
    tsv = tmp_path / "s0_clusters.tsv"
    s0 = dereplicate(MINI, str(rep), str(tsv), cfg, tmp_dir=str(tmp_path / "s0work"))
    assert s0["n_representatives"] == EXPECTED_REPS
    assert s0["homology_gate"] is False

    # --- S1 -------------------------------------------------------------- #
    out = tmp_path / "s1_3di"
    s1 = tokenize(str(rep), str(out), cfg)

    # Hand-off contract for S2: one 3Di record per S0 representative.
    rep_ids = [rid for rid, _ in parse_fasta(str(rep))]
    assert len(rep_ids) == s0["n_representatives"]
    assert s1["n_records"] == s0["n_representatives"], (
        "S1 must tokenize exactly the S0 representatives")
    assert s1["n_nonempty"] == s1["n_records"]
    assert s1["all_lengths_match"] is True

    # The two concrete artifacts S2 will read: the representative FASTA and the
    # Foldseek query DB (with its 3Di _ss sub-DB).
    assert rep.exists() and rep.stat().st_size > 0
    assert s1["querydb"] is not None
    assert os.path.exists(s1["querydb"] + "_ss"), "Foldseek query DB has no 3Di records"
    assert os.path.exists(s1["threedi_fasta"]), "inspectable 3Di artifact missing"
