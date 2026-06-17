"""S5 cleft-filter tests (Checkpoint 4).

Positive-output: assert a catalytic pocket is selected near the Ser OG and that the
metrics (including the primary exposure discriminator) come back finite. Requires
fpocket on PATH and fetched controls; skips otherwise.

Size-invariance: scale a control's coordinates x2 and assert that the size-invariant
peripherality forms (`percentile`, `rg_norm`) are unchanged while `absolute` doubles.
This proves the invariance directly rather than inferring it from the calibration.
"""
from __future__ import annotations

import math
import os
import shutil

import pytest

from proteus.s5_cleft_filter import _catalytic_og, _peripherality_modes, analyze_cleft
from proteus.utils import load_config, load_structure, protein_atoms

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (pdb_id, catalytic Ser res id) — IsPETase, a clearly surface-peripheral positive.
_IS_PETASE = ("6EQE", 160)


def _require(pdb_id: str) -> str:
    if shutil.which("fpocket") is None:
        pytest.skip("fpocket not installed")
    p = os.path.join(REPO, "structures", f"{pdb_id}.pdb")
    if not os.path.exists(p):
        pytest.skip(f"{pdb_id}.pdb not fetched — run controls/fetch_controls.py")
    return p


def _structure_only(pdb_id: str):
    """Load a control's protein atoms without needing fpocket (geometry-only tests)."""
    p = os.path.join(REPO, "structures", f"{pdb_id}.pdb")
    if not os.path.exists(p):
        pytest.skip(f"{pdb_id}.pdb not fetched — run controls/fetch_controls.py")
    return protein_atoms(load_structure(p))


def test_s5_selects_catalytic_pocket_and_returns_finite_metrics():
    cfg = load_config()
    r = analyze_cleft(_require("6EQE"), 160, cfg)
    assert r["pocket_id"] is not None, "no catalytic pocket selected for IsPETase"
    assert r["dist_og_pocket"] <= cfg["s5_cleft_filter"]["catalytic_pocket_max_dist"]
    for k in ("exposure", "volume", "druggability", "depth", "aromatics", "polarity",
              "hydrophobicity"):
        assert k in r["metrics"], f"missing metric {k}"
        assert math.isfinite(r["metrics"][k]), f"metric {k} is not finite"
    # IsPETase catalytic Ser is surface-peripheral -> positive exposure value in every
    # mode (percentile/rg_norm/absolute are all > 0 for a surface Ser).
    assert r["metrics"]["exposure"] > 0
    assert r["peripherality_mode"] == cfg["s5_cleft_filter"]["peripherality_mode"]
    assert set(r["peripherality"]) == {"absolute", "rg_norm", "percentile"}


def test_peripherality_size_invariance_under_coordinate_scaling():
    """Scaling ALL coordinates x2 must leave the size-invariant forms unchanged and
    double the absolute distance. This is the direct proof of invariance (Checkpoint 4).
    Geometry-only: no fpocket required."""
    pdb_id, ser = _IS_PETASE
    protein = _structure_only(pdb_id)
    og = _catalytic_og(protein, ser)
    assert og is not None, f"no Ser OG at residue {ser} in {pdb_id}"

    base = _peripherality_modes(protein, og, ser)

    # Scale every atomic coordinate by 2 (uniform isotropic scaling about the origin).
    scaled = protein.copy()
    scaled.coord = scaled.coord * 2.0
    og2 = _catalytic_og(scaled, ser)
    after = _peripherality_modes(scaled, og2, ser)

    # Size-invariant forms: unchanged within float tolerance.
    assert after["percentile"] == pytest.approx(base["percentile"], abs=1e-9)
    assert after["rg_norm"] == pytest.approx(base["rg_norm"], rel=1e-9)
    # Size-DEPENDENT form: doubles.
    assert after["absolute"] == pytest.approx(2.0 * base["absolute"], rel=1e-9)
    # And the absolute distance actually changed (guards against a no-op structure).
    assert after["absolute"] > base["absolute"]


def test_percentile_is_a_fraction():
    """percentile is, by construction, a fraction of CA atoms -> always in [0, 1].

    Note: this is the *CA*-based residue rank (per spec), which for an alpha/beta-
    hydrolase nucleophile-elbow Ser is fairly central even though the OG side-chain tip
    is solvent-exposed. Whether percentile *separates* PETases from decoys is a relative
    question answered by the calibration, not by an absolute threshold here."""
    pdb_id, ser = _IS_PETASE
    protein = _structure_only(pdb_id)
    og = _catalytic_og(protein, ser)
    p = _peripherality_modes(protein, og, ser)
    assert 0.0 <= p["percentile"] <= 1.0
