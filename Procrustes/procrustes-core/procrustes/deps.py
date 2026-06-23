"""Locate the Lattice digital reference object (DRO).

Procrustes draws its NON-bi-exponential ground truth from Lattice's clean-room,
seed-generated generators -- it never ships data files of its own.  In the
monorepo Lattice is a (private) submodule; in a standalone checkout it may be
``pip install``-ed or pointed to via ``LATTICE_PATH``.  This module resolves it
once, at import time, so the rest of the package can simply ``import lattice``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

__all__ = ["ensure_lattice"]


def _has_lattice(path: str | os.PathLike) -> bool:
    return (Path(path) / "lattice" / "__init__.py").is_file()


def ensure_lattice() -> str:
    """Make ``lattice`` importable; return the directory that provides it.

    Resolution order (first hit wins):
      1. already importable (installed / on PYTHONPATH),
      2. ``$LATTICE_PATH``,
      3. a sibling ``Lattice/`` dir walking up from this file (monorepo layout),
      4. the monorepo dev fallback ``~/researchProj/Lattice``.
    Raises ``ImportError`` with a clear remedy if none are found.
    """
    try:
        import lattice  # noqa: F401
        return os.path.dirname(os.path.dirname(lattice.__file__))
    except ImportError:
        pass

    candidates: list[str] = []
    env = os.environ.get("LATTICE_PATH")
    if env:
        candidates.append(env)
    for parent in Path(__file__).resolve().parents:
        candidates.append(str(parent / "Lattice"))
    candidates.append(os.path.expanduser("~/researchProj/Lattice"))

    for cand in candidates:
        if _has_lattice(cand):
            sys.path.insert(0, cand)
            return cand

    raise ImportError(
        "Lattice DRO not found. Install the `lattice` package, set LATTICE_PATH "
        "to the dir containing `lattice/`, or initialise the Lattice submodule "
        "(`git submodule update --init Lattice`)."
    )
