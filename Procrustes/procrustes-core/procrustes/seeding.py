"""Deterministic seeding -- every random draw in Procrustes is reproducible.

Mirrors the Minos house convention: a single GLOBAL_SEED and an explicit
``np.random.Generator`` factory, never bare ``np.random``.  The non-bi-exp
ground truth itself is seeded inside Lattice; these seeds drive only the
calibration/test splits and any Procrustes-side resampling.
"""
from __future__ import annotations

import numpy as np

GLOBAL_SEED = 20240517


def make_rng(seed: int | None = None) -> np.random.Generator:
    return np.random.default_rng(GLOBAL_SEED if seed is None else seed)
