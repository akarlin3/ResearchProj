"""Estimator tests -- skipped automatically when torch is unavailable.

Kept small/fast (tiny flow, few epochs) so it can run in CI when the
``[estimator]`` extra is present. The numpy-only core is tested elsewhere.
"""
import numpy as np
import pytest

torch = pytest.importorskip("torch")  # noqa: F841  skip whole module w/o torch

from caliper import metrics as M  # noqa: E402
from caliper.conformal import SplitConformalQuantile  # noqa: E402
from caliper.estimator_maf import MAFPosterior  # noqa: E402
from caliper.forward import synthetic_cohort  # noqa: E402

LEVELS = np.array([0.05, 0.5, 0.95])


def _small_est():
    return MAFPosterior(n_bvalues=synthetic_cohort(n=1).bvalues.size,
                        hidden=32, n_layers=3, epochs=15, n_posterior=200,
                        seed=0)


def test_predict_quantiles_shape_and_monotone():
    train = synthetic_cohort(n=1500, snr=50.0, seed=0)
    est = _small_est().fit(train.signals, train.params)
    test = synthetic_cohort(n=200, snr=50.0, seed=2)
    q = est.predict_quantiles(test.signals, LEVELS)
    assert q.shape == (200, 3, 3)
    assert np.all(np.diff(q, axis=2) >= -1e-6)  # monotone quantiles
    assert np.isfinite(q).all()


def test_training_nll_decreases():
    train = synthetic_cohort(n=1500, snr=50.0, seed=0)
    est = _small_est().fit(train.signals, train.params)
    assert est.history_[-1] < est.history_[0]  # learned something
    assert np.isfinite(est.history_[-1])


def test_predict_before_fit_raises():
    est = _small_est()
    with pytest.raises(RuntimeError):
        est.predict_quantiles(np.zeros((2, est.n_bvalues)), LEVELS)


def test_conformal_improves_marginal_coverage_under_shift():
    # train high SNR, deploy low SNR -> raw over-confident -> conformal fixes it
    train = synthetic_cohort(n=2000, snr=60.0, seed=0)
    cal = synthetic_cohort(n=800, snr=25.0, seed=1)
    test = synthetic_cohort(n=800, snr=25.0, seed=2)
    est = _small_est().fit(train.signals, train.params)
    q_cal = est.predict_quantiles(cal.signals, LEVELS)
    q_test = est.predict_quantiles(test.signals, LEVELS)
    raw = M.score_quantiles(test.params, q_test, LEVELS, alpha=0.1)
    cq = SplitConformalQuantile(LEVELS).calibrate(q_cal, cal.params)
    cor = M.score_quantiles(test.params, cq.apply(q_test), LEVELS, alpha=0.1)
    # raw under-covers; conformal moves every parameter closer to nominal
    for r, c in zip(raw, cor):
        assert abs(c.coverage_gap) <= abs(r.coverage_gap) + 0.02
    # at least one parameter was clearly under-covered raw and is fixed.
    # (This is a small/fast flow trained for only 15 epochs, so the
    #  over-confidence is milder than the full demo's ~0.53 coverage.)
    assert min(r.coverage for r in raw) < 0.82
    assert min(abs(c.coverage_gap) for c in cor) < 0.06
