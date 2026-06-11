"""S5 cleft-filter tests (Checkpoint 4).

Positive-output: assert a catalytic pocket is selected near the Ser OG and that the
metrics (including the primary exposure discriminator) come back finite. Requires
fpocket on PATH and fetched controls; skips otherwise.
"""
from __future__ import annotations

import math
import os
import shutil

import pytest

from proteus.s5_cleft_filter import analyze_cleft
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _require(pdb_id: str) -> str:
    if shutil.which("fpocket") is None:
        pytest.skip("fpocket not installed")
    p = os.path.join(REPO, "structures", f"{pdb_id}.pdb")
    if not os.path.exists(p):
        pytest.skip(f"{pdb_id}.pdb not fetched — run controls/fetch_controls.py")
    return p


def test_s5_selects_catalytic_pocket_and_returns_finite_metrics():
    cfg = load_config()
    r = analyze_cleft(_require("6EQE"), 160, cfg)
    assert r["pocket_id"] is not None, "no catalytic pocket selected for IsPETase"
    assert r["dist_og_pocket"] <= cfg["s5_cleft_filter"]["catalytic_pocket_max_dist"]
    for k in ("exposure", "volume", "druggability", "depth", "aromatics", "polarity",
              "hydrophobicity"):
        assert k in r["metrics"], f"missing metric {k}"
        assert math.isfinite(r["metrics"][k]), f"metric {k} is not finite"
    # IsPETase catalytic Ser is surface-peripheral -> positive exposure value
    assert r["metrics"]["exposure"] > 0
