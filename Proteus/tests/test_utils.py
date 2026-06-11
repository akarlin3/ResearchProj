"""Unit tests for the shared utils layer (Checkpoint 1).

Positive-output: assert real returned values, not merely "no exception".
"""
from __future__ import annotations

import math

from proteus.utils import euclidean, get_seed, load_config


def test_load_config_returns_seed_and_stage_blocks():
    cfg = load_config()
    assert isinstance(cfg, dict)
    # the single global seed every stochastic stage reads
    assert cfg["random_seed"] == 1729
    # S4/S5 blocks must be present — calibration cannot run without them
    assert "s4_geometry" in cfg, "missing s4_geometry block in config"
    assert "s5_cleft_filter" in cfg, "missing s5_cleft_filter block in config"


def test_get_seed_reads_global_seed():
    cfg = load_config()
    assert get_seed(cfg) == 1729


def test_euclidean_distance_3_4_5():
    # classic 3-4-5 right triangle -> distance 5.0 (a positive, known artifact)
    d = euclidean([0.0, 0.0, 0.0], [3.0, 4.0, 0.0])
    assert math.isclose(d, 5.0, rel_tol=1e-9)
    # symmetry + non-negativity
    assert math.isclose(euclidean([3.0, 4.0, 0.0], [0.0, 0.0, 0.0]), 5.0, rel_tol=1e-9)
    assert euclidean([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]) == 0.0
