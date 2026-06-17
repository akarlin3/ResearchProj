"""S4 control-recovery tests (Checkpoint 2).

These are the *supervised correctness* checks for an otherwise unsupervised method:
the blind triad detector must (a) recover the documented catalytic serine in each
positive control, and (b) correctly find NO catalytic triad at the mutated site of
the S165A inactivated control (6THS) — the deliberate trap.

Structures must be fetched first (controls/fetch_controls.py); tests skip if absent.
"""
from __future__ import annotations

import os

import pytest

from proteus.s4_geometry import analyze_model
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(REPO, "structures")


def _pdb(pdb_id: str) -> str:
    return os.path.join(STRUCT, f"{pdb_id}.pdb")


def _require(pdb_id: str) -> str:
    p = _pdb(pdb_id)
    if not os.path.exists(p):
        pytest.skip(f"{pdb_id}.pdb not fetched — run controls/fetch_controls.py")
    return p


def _triad_with_ser(result: dict, res_id: int):
    return [t for t in result["triads"] if t["ser"]["res_id"] == res_id]


def test_s4_recovers_ispetase_ser160():
    """6EQE (IsPETase): blind detection finds the Ser160/His237/Asp206 triad."""
    cfg = load_config()
    res = analyze_model(_require("6EQE"), cfg)
    hits = _triad_with_ser(res, 160)
    assert hits, "S4 failed to recover IsPETase catalytic Ser160"
    t = hits[0]
    assert t["his"]["res_id"] == 237, f"expected His237 partner, got {t['his']}"
    assert t["acid"]["res_id"] == 206, f"expected Asp206 partner, got {t['acid']}"
    assert t["passes"], "Ser160 triad present but failed full gate (oxyanion hole)"
    assert t["ser_og_his_ne2"] <= cfg["s4_geometry"]["ser_og_his_ne2_max"]
    assert t["his_nd1_acid"] <= cfg["s4_geometry"]["his_nd1_acid_max"]


def test_s4_recovers_lcc_wt_ser165():
    """4EB0 (LCC wild-type): blind detection finds the intact Ser165 triad."""
    cfg = load_config()
    res = analyze_model(_require("4EB0"), cfg)
    hits = _triad_with_ser(res, 165)
    assert hits, "S4 failed to recover LCC-WT catalytic Ser165"
    assert hits[0]["passes"], "Ser165 triad present but failed full gate"


def test_s4_6ths_s165a_trap_reports_no_catalytic_triad():
    """6THS (LCC-ICCG S165A): the catalytic serine is mutated to Ala, so there must
    be NO triad anchored at residue 165. This is EXPECTED, not a detector failure."""
    cfg = load_config()
    res = analyze_model(_require("6THS"), cfg)
    assert not _triad_with_ser(res, 165), (
        "6THS residue 165 is ALA (S165A) — no triad should be anchored there")
