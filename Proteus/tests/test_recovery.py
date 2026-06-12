"""Divergent-positive recovery test (held out of the calibration anchor).

GuaPA and MG8 (the recovery PETases in references.csv) are unresolvable here — no
reachable sequence/structure — so the archaeal PET46 structure (8B4U) stands in as
a real sequence-divergent positive. It is scored against the FINISHED IsPETase/LCC
anchor (never used to build it), so this measures GENERALIZATION, not fit.

The known answer (per the "held out + widened line" decision):
  * PET46 recovers a catalytic triad + pocket (S4/S5 fire on an archaeal PETase).
  * It scores ABOVE every negative (the fold+cleft signal generalizes) ...
  * ... but BELOW the IsPETase/LCC production operating point (that N=2 line is too
    strict for divergent PETases).
  * A widened operating point that also keeps PET46 holds precision 1.0 (no negative
    sneaks in) — the recommended next-gen threshold.
  * Adding PET46 must NOT change the production anchor/operating point (it is held out).

Requires fpocket + fetched controls incl. 8B4U; skips otherwise.
"""
from __future__ import annotations

import copy
import os
import shutil

import pytest

from proteus.calibrate import analyze_controls, recovery_screen, score_analysis
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(REPO, "structures")


def _guard():
    if shutil.which("fpocket") is None:
        pytest.skip("fpocket not installed")
    needed = ["6EQE", "4EB0", "6THS", "1TCA", "1EA5", "1CRL", "1EVQ", "8B4U"]
    missing = [p for p in needed if not os.path.exists(os.path.join(STRUCT, f"{p}.pdb"))]
    if missing:
        pytest.skip(f"controls not fetched: {missing} — run controls/fetch_controls.py")


def _calibrate(cfg):
    return score_analysis(analyze_controls(cfg, STRUCT), cfg)


def test_recovery_excluded_from_anchor_keeps_production_calibration():
    """The recovery structure must NOT be in the scored controls / anchor."""
    _guard()
    cfg = load_config()
    cal = _calibrate(cfg)
    # PET46 is never one of the scored controls (it's held out of the anchor path)
    assert "PET46" not in cal["per_control"], "recovery control leaked into the anchor"
    assert cal["verdict"]["separated"] is True  # production still separates


def test_divergent_positive_recovers_above_negatives_below_line():
    _guard()
    cfg = load_config()
    cal = _calibrate(cfg)
    rec = recovery_screen(cfg, STRUCT, cal)

    pet46 = next(r for r in rec["recovery"] if r["id"] == "PET46")
    assert pet46["present"] and pet46["triad_found"] and pet46["pocket_ok"], (
        "S4/S5 should fire on the archaeal PETase 8B4U")
    assert pet46["composite"] is not None
    # generalization: above every negative, below the mesophilic production line
    assert pet46["above_all_negatives"] is True
    assert pet46["above_production_line"] is False
    assert pet46["status"] == "above_negatives_below_line"
    assert pet46["composite"] > rec["max_negative"]
    assert pet46["composite"] < rec["production_threshold"]


def test_widened_operating_point_recovers_divergent_positive_at_precision_1():
    _guard()
    cfg = load_config()
    cal = _calibrate(cfg)
    rec = recovery_screen(cfg, STRUCT, cal)

    w = rec["widened_operating_point"]
    assert w is not None, "a widened line should be proposed once PET46 clears the negatives"
    assert "PET46" in w["includes_recovery"]
    # lowering the line to keep PET46 must not let any negative in
    assert w["false_positives"] == 0
    assert w["precision"] == 1.0
    assert w["threshold"] < rec["production_threshold"]  # it is a WIDER (lower) line


def test_recovery_screen_does_not_mutate_calibration():
    """recovery_screen reads the anchor; it must not change the calibration result."""
    _guard()
    cfg = load_config()
    cal = _calibrate(cfg)
    before = copy.deepcopy(cal["operating_point"])
    _ = recovery_screen(cfg, STRUCT, cal)
    assert cal["operating_point"] == before
