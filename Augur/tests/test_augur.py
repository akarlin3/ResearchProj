"""Augur gate tests: the submission block is engaged and the manifests are well-formed.

These are the CP2 gate assertions (argument coherent / citations verified / submission block in
place), in machine-checkable form. They intentionally assert the block is ENGAGED — that is the
correct state until Fashion + Minos publish.
"""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import check_anchors as ca  # noqa: E402


def test_submission_block_is_engaged():
    """Load-bearing anchors are unpublished -> block engaged -> gate exits non-zero."""
    assert ca.block_engaged(), "submission block must be engaged while anchors are unpublished"
    assert ca.main() == 1, "check_anchors.main() must signal BLOCKED (exit 1) today"


def test_load_bearing_anchors_present_and_unpublished():
    by_name = {a.name: a for a in ca.ANCHORS}
    for name in ("Fashion", "Minos"):
        assert by_name[name].load_bearing, f"{name} must be load-bearing"
        assert not by_name[name].published, f"{name} must be pinned unpublished (no DOI yet)"
        assert by_name[name].doi is None


def test_block_would_lift_only_with_both_load_bearing_published():
    """Flipping a single load-bearing anchor is not enough; both must publish."""
    from dataclasses import replace

    one = tuple(replace(a, published=True, doi="10.x/x") if a.name == "Fashion" else a
                for a in ca.ANCHORS)
    assert ca.block_engaged(one), "block must stay engaged until BOTH Fashion and Minos publish"

    both = tuple(replace(a, published=True, doi="10.x/x")
                 if a.name in ("Fashion", "Minos") else a for a in ca.ANCHORS)
    assert not ca.block_engaged(both), "block must lift once both load-bearing anchors publish"


def test_citations_have_verified_tier_a_with_dois():
    text = (ROOT / "CITATIONS.md").read_text()
    # Tier A primary sources verified this build, with their DOIs.
    assert "10.1016/j.acra.2018.08.012" in text, "Sun 2019 DOI must be present"
    assert "10.1177/0284185118791201" in text, "Yang 2019 DOI must be present"
    assert "r = 0.389" in text, "the verified D*-Ktrans value must be quoted"
    assert "Tier A" in text and "Tier B" in text, "citations must separate verified vs inherited"


def test_synthesis_and_assumptions_present():
    for fname in ("synthesis.md", "ASSUMPTIONS.md", "SUBMISSION_BLOCK.md", "README.md"):
        p = ROOT / fname
        assert p.exists() and p.stat().st_size > 0, f"{fname} must exist and be non-empty"
    # PROVISIONAL discipline is declared where it matters.
    assert "PROVISIONAL" in (ROOT / "ASSUMPTIONS.md").read_text()
    assert "SUBMISSION-BLOCKED" in (ROOT / "synthesis.md").read_text() \
        or "SUBMISSION BLOCK" in (ROOT / "synthesis.md").read_text().upper()
