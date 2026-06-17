"""Tests for the split-conformal / CQR wrapper (numpy-only)."""
from statistics import NormalDist

import numpy as np
import pytest

from caliper import conformal as C
from caliper import metrics as M

LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
_Z = np.array([NormalDist().inv_cdf(p) for p in LEVELS])


def _overconfident(n, rng, spread=0.5):
    mu = rng.normal(size=(n, 1))
    y = mu + rng.normal(size=(n, 1))               # true sigma = 1
    q = mu[:, :, None] + spread * _Z[None, None, :]  # under-dispersed
    return y, q


def test_conformity_scores_definition():
    s = C.conformity_scores(np.array([0.0]), np.array([2.0]), np.array([3.0]))
    assert s[0] == pytest.approx(1.0)  # y above hi by 1
    s = C.conformity_scores(np.array([0.0]), np.array([2.0]), np.array([-1.0]))
    assert s[0] == pytest.approx(1.0)  # y below lo by 1
    s = C.conformity_scores(np.array([0.0]), np.array([2.0]), np.array([1.0]))
    assert s[0] == pytest.approx(-1.0)  # inside -> negative


def test_offset_increases_with_lower_alpha():
    rng = np.random.default_rng(0)
    scores = rng.normal(size=1000)
    q90 = C.conformal_offset(scores, alpha=0.1)
    q50 = C.conformal_offset(scores, alpha=0.5)
    assert q90 > q50


def test_conformal_restores_marginal_coverage():
    rng = np.random.default_rng(0)
    y_cal, q_cal = _overconfident(3000, rng)
    y_te, q_te = _overconfident(5000, rng)
    raw = M.score_quantiles(y_te, q_te, LEVELS, alpha=0.1)[0]
    cq = C.SplitConformalQuantile(LEVELS).calibrate(q_cal, y_cal)
    q_corr = cq.apply(q_te)
    cor = M.score_quantiles(y_te, q_corr, LEVELS, alpha=0.1)[0]
    assert raw.coverage < 0.75              # raw is over-confident
    assert abs(cor.coverage_gap) < 0.03     # conformal restores nominal


def test_apply_preserves_monotonic_quantiles():
    rng = np.random.default_rng(1)
    y_cal, q_cal = _overconfident(2000, rng)
    _, q_te = _overconfident(1000, rng)
    cq = C.SplitConformalQuantile(LEVELS).calibrate(q_cal, y_cal)
    q_corr = cq.apply(q_te)
    assert np.all(np.diff(q_corr, axis=2) >= -1e-9)


def test_apply_before_calibrate_raises():
    cq = C.SplitConformalQuantile(LEVELS)
    with pytest.raises(RuntimeError):
        cq.apply(np.zeros((4, 1, 5)))


def test_calibrate_shape_validation():
    cq = C.SplitConformalQuantile(LEVELS)
    with pytest.raises(ValueError):
        cq.calibrate(np.zeros((4, 2, 5)), np.zeros((4, 3)))
