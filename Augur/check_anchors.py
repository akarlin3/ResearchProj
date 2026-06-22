#!/usr/bin/env python3
"""Augur submission-block gate.

Augur is the end-stage *synthesis* of the IVIM-UQ program. It is NOT submittable until its
load-bearing anchors publish. This script is the machine-checkable form of that rule: it pins each
anchor's published state and exits non-zero (block ENGAGED) while the load-bearing anchors are
unpublished.

Re-validation contract (ASSUMPTIONS.md §3): when an anchor publishes, set published=True with its
DOI here, then run `bash reproduce.sh`. The block lifts only when Fashion AND Minos are both
published (Lethe strongly recommended).

Exit code: 0 = clear to submit; 1 = submission BLOCKED (the expected state today).
"""
from __future__ import annotations

import sys
from dataclasses import dataclass


@dataclass(frozen=True)
class Anchor:
    name: str
    role: str
    load_bearing: bool   # must publish before submission
    published: bool
    doi: str | None


# Pinned 2026-06-22 (ASSUMPTIONS.md). Flip published/doi here when an anchor lands.
ANCHORS = (
    Anchor("Fashion", "trust (ruler)",                 load_bearing=True,  published=False, doi=None),
    Anchor("Minos",   "value-of-information (decision)", load_bearing=True,  published=False, doi=None),
    Anchor("Lethe",   "action (wrong-size limit)",      load_bearing=False, published=False, doi=None),
    Anchor("Gauge",   "action (identifiability wall)",  load_bearing=False, published=False, doi=None),
)


def block_engaged(anchors=ANCHORS) -> bool:
    """The block is engaged while ANY load-bearing anchor is unpublished."""
    return any(a.load_bearing and not a.published for a in anchors)


def main() -> int:
    print("Augur submission-block gate")
    print("=" * 60)
    for a in ANCHORS:
        flag = "load-bearing" if a.load_bearing else "recommended"
        state = f"PUBLISHED ({a.doi})" if a.published else "UNPUBLISHED (no DOI)"
        print(f"  {a.name:8s} [{flag:12s}] {a.role:32s} -> {state}")
    print("=" * 60)
    if block_engaged():
        pending = [a.name for a in ANCHORS if a.load_bearing and not a.published]
        print(f"SUBMISSION BLOCKED: load-bearing anchors unpublished: {', '.join(pending)}.")
        print("Augur is PROVISIONAL. Do NOT submit. See SUBMISSION_BLOCK.md / ASSUMPTIONS.md.")
        return 1
    print("All load-bearing anchors published. Block CLEAR — re-run the citation re-verification")
    print("(CITATIONS.md Tier B) before submitting.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
