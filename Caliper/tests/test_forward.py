"""Tests for the IVIM forward model and synthetic cohort (numpy-only)."""
import numpy as np

from caliper import forward as F


def test_ivim_signal_b0_is_s0():
    s = F.ivim_signal(np.array([0.0]), D=1.0, f=0.2, Dstar=20.0, s0=1.0)
    np.testing.assert_allclose(s.ravel(), [1.0])


def test_ivim_signal_monotone_decay():
    b = F.DEFAULT_BVALUES
    s = F.ivim_signal(b, D=1.5, f=0.2, Dstar=30.0)
    assert np.all(np.diff(s.ravel()) <= 1e-12)


def test_ivim_signal_vectorized_shape():
    b = F.DEFAULT_BVALUES
    D = np.array([1.0, 2.0])
    f = np.array([0.1, 0.3])
    Ds = np.array([20.0, 50.0])
    s = F.ivim_signal(b, D, f, Ds)
    assert s.shape == (2, b.size)


def test_perfusion_fraction_limits():
    # f=0 -> pure tissue exp(-bD); f=1 -> pure perfusion exp(-bD*)
    b = np.array([0.0, 200.0])
    s_tissue = F.ivim_signal(b, D=1.0, f=0.0, Dstar=50.0)
    np.testing.assert_allclose(s_tissue.ravel(),
                               np.exp(-b * 1.0 * 1e-3))
    s_perf = F.ivim_signal(b, D=1.0, f=1.0, Dstar=50.0)
    np.testing.assert_allclose(s_perf.ravel(),
                               np.exp(-b * 50.0 * 1e-3))


def test_synthetic_cohort_reproducible():
    a = F.synthetic_cohort(n=500, seed=7)
    b = F.synthetic_cohort(n=500, seed=7)
    np.testing.assert_array_equal(a.params, b.params)
    np.testing.assert_array_equal(a.signals, b.signals)


def test_synthetic_cohort_seed_changes_data():
    a = F.synthetic_cohort(n=500, seed=1)
    b = F.synthetic_cohort(n=500, seed=2)
    assert not np.allclose(a.params, b.params)


def test_param_priors_within_bounds():
    c = F.synthetic_cohort(n=3000, seed=0)
    D, f, Ds = c.params[:, 0], c.params[:, 1], c.params[:, 2]
    assert D.min() >= 0.5 and D.max() <= 2.5
    assert f.min() >= 0.05 and f.max() <= 0.40
    assert Ds.min() >= 10.0 - 1e-6 and Ds.max() <= 100.0 + 1e-6


def test_snr_controls_residual_scale():
    lo = F.synthetic_cohort(n=2000, snr=20.0, noise="gaussian", seed=0)
    hi = F.synthetic_cohort(n=2000, snr=100.0, noise="gaussian", seed=0)
    res_lo = (lo.signals - lo.signals_clean).std()
    res_hi = (hi.signals - hi.signals_clean).std()
    assert res_lo > res_hi  # lower SNR -> noisier
