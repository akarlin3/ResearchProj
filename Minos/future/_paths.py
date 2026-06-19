"""Read-only path wiring for ``Minos/future/``.

``future/`` depends on the validated theory half (``Minos/theory``, ``Minos/minos-core``)
and on the two upstream papers (``Fashion``, ``Gauge``) **by import only**. This module
inserts their locations on ``sys.path``; it never copies, edits, or shadows them. Every
path is resolved relative to *this file*, so imports work regardless of the current working
directory.

This is the single chokepoint for the dependency graph. When ``future/`` is promoted into
the real Minos paper (see ``PROMOTION.md``) these paths move in exactly one place.

Dependency provenance and pinned versions live in ``ASSUMPTIONS.md``. The Fashion/Gauge
wiring is **PROVISIONAL**: both papers are in review, and anything imported through
``add_fashion()`` / ``add_gauge()`` feeds results that must be re-validated if those papers
change in revision.
"""
from __future__ import annotations

import sys
from pathlib import Path

_FUTURE = Path(__file__).resolve().parent      # Minos/future
MINOS = _FUTURE.parent                          # Minos
REPO = MINOS.parent                             # researchProj (repo root)

THEORY = MINOS / "theory"
MINOS_CORE = MINOS / "minos-core"
FASHION = REPO / "Fashion"                      # PROVISIONAL dependency (in review)
GAUGE = REPO / "Gauge"                          # PROVISIONAL dependency (in review)


def _prepend(p: Path) -> str:
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)
    return s


def add_minos_core() -> str:
    """Make ``import minos`` resolve to the validated v1-v3 package (read-only)."""
    if not (MINOS_CORE / "minos" / "__init__.py").exists():
        raise FileNotFoundError(f"minos-core package not found at {MINOS_CORE}")
    return _prepend(MINOS_CORE)


def add_theory() -> str:
    """Expose the ``theory/`` script directory (read-only)."""
    if not (THEORY / "confirm.py").exists():
        raise FileNotFoundError(f"theory/ not found at {THEORY}")
    return _prepend(THEORY)


def add_fashion() -> str:
    """Make ``import uq`` resolve to Fashion's calibration package (read-only, PROVISIONAL)."""
    if not (FASHION / "uq" / "__init__.py").exists():
        raise FileNotFoundError(f"Fashion uq package not found at {FASHION}")
    return _prepend(FASHION)


def add_gauge() -> str:
    """Make ``import gauge`` resolve to Gauge's coverage/monitor package (read-only, PROVISIONAL)."""
    if not (GAUGE / "gauge" / "__init__.py").exists():
        raise FileNotFoundError(f"Gauge package not found at {GAUGE}")
    return _prepend(GAUGE)


def add_all() -> dict[str, str]:
    """Wire everything; return the resolved paths (useful for provenance printing)."""
    return {
        "minos_core": add_minos_core(),
        "theory": add_theory(),
        "fashion": add_fashion(),
        "gauge": add_gauge(),
    }


if __name__ == "__main__":  # pragma: no cover - manual provenance dump
    import json

    print(json.dumps(add_all(), indent=2))
