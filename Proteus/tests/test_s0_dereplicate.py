"""S0 dereplication — positive-output test on the known-answer mini corpus.

The fixture (tests/data/mini_corpus.fasta, built by make_mini_corpus.py) has a
DOCUMENTED expected outcome: 10 input sequences, of which IsPETase + its two
near-duplicate variants collapse into one cluster, leaving 8 representatives.
This test asserts that known answer — and that NO homology gate was applied
(every distinct sequence, negatives and decoys included, survives).

Skips cleanly if MMseqs2 is not installed (so the suite still runs on a thin host).
"""
from __future__ import annotations

import os
import shutil

import pytest

from proteus.s0_dereplicate import dereplicate, parse_clusters_tsv, parse_fasta
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MINI = os.path.join(REPO, "tests", "data", "mini_corpus.fasta")

# Documented known answer (see make_mini_corpus.py).
EXPECTED_INPUTS = 10
EXPECTED_REPS = 8
ISPETASE_CLUSTER = {"IsPETase", "IsPETase_var1", "IsPETase_var2"}
SURVIVING_SINGLETONS = {"LCC_WT", "CalB", "AChE", "CRL", "Est2",
                        "decoy_allalpha", "decoy_random"}


def _guard():
    if shutil.which("mmseqs") is None:
        pytest.skip("mmseqs2 not installed")
    if not os.path.exists(MINI):
        pytest.skip("mini_corpus.fasta missing — run tests/data/make_mini_corpus.py")


def test_mini_corpus_has_expected_inputs():
    if not os.path.exists(MINI):
        pytest.skip("mini_corpus.fasta missing — run tests/data/make_mini_corpus.py")
    ids = [rid for rid, _ in parse_fasta(MINI)]
    assert len(ids) == EXPECTED_INPUTS, f"fixture should have {EXPECTED_INPUTS} records"
    assert ISPETASE_CLUSTER | SURVIVING_SINGLETONS == set(ids)


def test_s0_collapses_near_duplicates_to_documented_rep_count(tmp_path):
    _guard()
    cfg = load_config()
    rep = tmp_path / "s0_representatives.fasta"
    tsv = tmp_path / "s0_clusters.tsv"
    summary = dereplicate(MINI, str(rep), str(tsv), cfg, tmp_dir=str(tmp_path / "work"))

    # The documented known answer: 10 inputs -> 8 representatives.
    assert summary["n_input"] == EXPECTED_INPUTS
    assert summary["n_representatives"] == EXPECTED_REPS, (
        f"expected {EXPECTED_REPS} reps, got {summary['n_representatives']} "
        f"({summary['representatives']})")
    assert rep.exists() and tsv.exists()

    # The two IsPETase variants must collapse INTO the IsPETase cluster (one rep,
    # three members), and every distinct sequence must survive as a rep.
    clusters = parse_clusters_tsv(str(tsv))
    member_to_rep = {m: r for r, members in clusters.items() for m in members}
    ip_rep = member_to_rep["IsPETase"]
    assert member_to_rep["IsPETase_var1"] == ip_rep
    assert member_to_rep["IsPETase_var2"] == ip_rep
    assert set(clusters[ip_rep]) == ISPETASE_CLUSTER

    reps = {rid for rid, _ in parse_fasta(str(rep))}
    # the 7 distinct non-IsPETase sequences each survive as their own rep
    assert SURVIVING_SINGLETONS <= reps
    # the variants are NOT representatives (they were collapsed)
    assert "IsPETase_var1" not in reps and "IsPETase_var2" not in reps


def test_s0_applies_no_homology_gate(tmp_path):
    """The negatives (CalB/AChE/CRL/Est2) and decoys are NOT PETases, yet they
    must all survive S0: dereplication never filters by similarity to a reference."""
    _guard()
    cfg = load_config()
    rep = tmp_path / "rep.fasta"
    tsv = tmp_path / "clu.tsv"
    summary = dereplicate(MINI, str(rep), str(tsv), cfg, tmp_dir=str(tmp_path / "w"))
    assert summary["homology_gate"] is False
    reps = {rid for rid, _ in parse_fasta(str(rep))}
    for non_petase in ("CalB", "AChE", "CRL", "Est2", "decoy_allalpha", "decoy_random"):
        assert non_petase in reps, f"{non_petase} wrongly dropped — looks like a homology gate"
