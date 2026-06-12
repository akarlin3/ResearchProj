"""Live AutoDock Vina docking smoke — the REAL docking path end-to-end.

Unlike test_docking.py (which drives the orchestration with a fake scorer), this
runs actual Vina: it preps a control receptor (Open Babel), docks the committed
PET-mimic ligand (controls/ligands/bhet.pdbqt) into the catalytic-Ser-OG box found
by S4, and asserts a finite binding affinity. This exercises proteus.docking's real
`vina_scorer` (receptor auto-prep + Vina maps + dock).

Skips cleanly unless the LOCAL docking toolchain is present: the `vina` Python
bindings, Open Babel on PATH, the control structure, and the ligand fixture. On a
host without them (this is a LOCAL/M4 step) the test is a clean skip, not a failure.
"""
from __future__ import annotations

import math
import os
import shutil

import pytest

from proteus import docking
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(REPO, "structures")
BHET = os.path.join(REPO, "controls", "ligands", "bhet.pdbqt")


def _guard():
    pytest.importorskip("vina", reason="AutoDock Vina python bindings not installed")
    if shutil.which("obabel") is None:
        pytest.skip("Open Babel (obabel) not on PATH — needed for receptor prep")
    if not os.path.exists(os.path.join(STRUCT, "6EQE.pdb")):
        pytest.skip("6EQE.pdb not fetched — run controls/fetch_controls.py")
    if not os.path.exists(BHET):
        pytest.skip("BHET ligand fixture missing")


def test_prepare_receptor_pdbqt_from_pdb(tmp_path):
    _guard()
    out = docking.prepare_receptor_pdbqt(
        os.path.join(STRUCT, "6EQE.pdb"), str(tmp_path / "rec.pdbqt"))
    assert os.path.exists(out) and os.path.getsize(out) > 0
    text = open(out).read()
    assert "ATOM" in text or "HETATM" in text, "receptor PDBQT has no atoms"


def test_live_dock_bhet_into_ispetase_active_site():
    """Dock BHET into IsPETase's Ser160 box; expect a finite (favourable) affinity."""
    _guard()
    cfg = load_config()
    cfg["docking"]["ligand_pdbqt"] = BHET
    cfg["docking"]["exhaustiveness"] = 4  # keep the smoke fast

    rec = docking.dock_model(os.path.join(STRUCT, "6EQE.pdb"), cfg,
                             docking.vina_scorer(), cand_id="6EQE")
    assert rec["docked"] is True, f"docking failed: {rec.get('stage_failed')}"
    assert rec["catalytic_ser"] == 160          # box centred on the S4 catalytic Ser
    assert rec["n_poses"] >= 1
    assert math.isfinite(rec["affinity"]), "Vina affinity is not finite"
    # a real ligand in a real active site scores favourably (negative kcal/mol)
    assert rec["affinity"] < 0.0
