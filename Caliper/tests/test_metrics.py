"""Pytest port of the metrics ruler invariants (numpy-only)."""
from statistics import NormalDist

import numpy as np
import pytest

from caliper import metrics as M


def test_pinball_known_values():
    assert M.pinball_loss(1.0, 0.0, 0.5) == pytest.approx(0.5)
    assert M.pinball_loss(0.0, 1.0, 0.5) == pytest.approx(0.5)
    assert M.pinball_loss(2.0, 0.0, 0.9) == pytest.approx(1.8)
    assert M.pinball_loss(0.0, 2.0, 0.9) == pytest.approx(0.2)


def test_interval_score_inside_is_width():
    assert M.interval_score(0.0, 2.0, 1.0, alpha=0.1) == pytest.approx(2.0)


def test_interval_score_penalizes_misses_symmetrically():
    assert M.interval_score(0.0, 2.0, 3.0, alpha=0.1) == pytest.approx(22.0)
    assert M.interval_score(0.0, 2.0, -1.0, alpha=0.1) == pytest.approx(22.0)


def test_interval_score_vectorized():
    lo = np.zeros(3)
    hi = np.full(3, 2.0)
    y = np.array([1.0, 3.0, -1.0])
    np.testing.assert_allclose(M.interval_score(lo, hi, y, alpha=0.1),
                               [2.0, 22.0, 22.0])


def test_empirical_coverage_bounds():
    y = np.array([0.0, 1.0, 2.0])
    assert M.empirical_coverage(y, np.full(3, -1.0), np.full(3, 3.0)) == 1.0
    assert M.empirical_coverage(y, np.full(3, 5.0), np.full(3, 6.0)) == 0.0


def test_central_interval_selects_levels():
    qlev = np.array([0.05, 0.5, 0.95])
    qp = np.array([[1.0, 2.0, 3.0]])
    lo, hi = M.central_interval(qp, qlev, alpha=0.1)
    assert lo[0] == pytest.approx(1.0)
    assert hi[0] == pytest.approx(3.0)


def test_tercile_groups_balanced():
    rng = np.random.default_rng(0)
    g = M.tercile_groups(rng.normal(size=9000))
    counts = np.bincount(g)
    assert counts.shape[0] == 3
    assert np.all(np.abs(counts - 3000) < 250)


def _well_specified_quantiles(n, levels, rng):
    mu = rng.normal(size=(n, 1))
    y_true = mu + rng.normal(size=(n, 1))
    z = np.array([NormalDist().inv_cdf(p) for p in levels])
    q_pred = mu[:, :, None] + z[None, None, :]
    return mu, y_true, q_pred


def test_well_specified_coverage_and_ece():
    rng = np.random.default_rng(0)
    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    mu, y_true, q_pred = _well_specified_quantiles(20000, levels, rng)
    scores = M.score_quantiles(y_true, q_pred, levels, alpha=0.1,
                               param_names=["g"], conditioning=mu)
    s = scores[0]
    assert abs(s.coverage - 0.9) < 0.02
    assert s.ece < 0.02
    assert s.coverage_gap == pytest.approx(s.coverage - 0.9)


def test_well_specified_conditional_on_covariate_is_nominal():
    rng = np.random.default_rng(1)
    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    mu, y_true, q_pred = _well_specified_quantiles(20000, levels, rng)
    scores = M.score_quantiles(y_true, q_pred, levels, alpha=0.1,
                               conditioning=mu)
    for g, v in scores[0].conditional.items():
        assert abs(v - 0.9) < 0.05


def test_overconfident_undercovers():
    rng = np.random.default_rng(2)
    levels = np.array([0.05, 0.5, 0.95])
    n = 5000
    mu = rng.normal(size=(n, 1))
    y_true = mu + rng.normal(size=(n, 1))
    z = np.array([NormalDist().inv_cdf(p) for p in levels])
    q_pred = mu[:, :, None] + 0.5 * z[None, None, :]  # half the true spread
    scores = M.score_quantiles(y_true, q_pred, levels, alpha=0.1)
    assert scores[0].coverage < 0.8  # clearly undercovers
    assert scores[0].coverage_gap < 0


def test_score_quantiles_validates_shapes():
    levels = np.array([0.1, 0.9])
    with pytest.raises(ValueError):
        M.score_quantiles(np.zeros((4,)), np.zeros((4, 1, 2)), levels)
    with pytest.raises(ValueError):
        M.score_quantiles(np.zeros((4, 1)), np.zeros((4, 1, 3)), levels)
    with pytest.raises(ValueError):  # non-ascending levels
        M.score_quantiles(np.zeros((4, 1)), np.zeros((4, 1, 2)),
                          np.array([0.9, 0.1]))


def test_format_scorecard_runs():
    rng = np.random.default_rng(3)
    levels = np.array([0.05, 0.5, 0.95])
    mu, y_true, q_pred = _well_specified_quantiles(500, levels, rng)
    scores = M.score_quantiles(y_true, q_pred, levels, param_names=["g"])
    txt = M.format_scorecard(scores, title="my-card")
    assert "my-card" in txt and "param" in txt and "conditional" in txt
