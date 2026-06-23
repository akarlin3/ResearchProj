"""Mittag-Leffler E_alpha(-t) numerics -- the kernel of the joint CTRW forward model.

Validated against the two closed forms that exist:
  * E_1(-t)     = exp(-t)
  * E_{1/2}(-t) = exp(t^2) * erfc(t)     (t >= 0)
plus the structural properties (E_alpha(0)=1, complete monotonicity, the small-t series).
"""
import numpy as np
from scipy import special

from levy import mittag_leffler as ml


def test_alpha_one_is_exponential():
    t = np.linspace(0.0, 30.0, 60)
    assert np.allclose(ml.mlf_neg(1.0, t), np.exp(-t), rtol=1e-9, atol=1e-12)


def test_alpha_half_closed_form():
    # E_{1/2}(-t) = exp(t^2) erfc(t); erfcx(t) = exp(t^2) erfc(t) is the stable form.
    t = np.linspace(0.0, 8.0, 40)
    expected = special.erfcx(t)
    assert np.allclose(ml.mlf_neg(0.5, t), expected, rtol=1e-6, atol=1e-9)


def test_at_zero_is_one():
    for a in (0.3, 0.5, 0.7, 0.85, 0.95, 1.0):
        assert np.isclose(ml.mlf_neg(a, 0.0), 1.0, rtol=1e-10)


def test_completely_monotone_decreasing():
    # E_alpha(-t) is positive, in (0,1], and strictly decreasing in t for t>0.
    t = np.linspace(0.0, 50.0, 200)
    for a in (0.4, 0.6, 0.8, 0.95):
        y = ml.mlf_neg(a, t)
        assert np.all(y > 0.0)
        assert np.all(y <= 1.0 + 1e-12)
        assert np.all(np.diff(y) < 1e-9)  # non-increasing (numerically)


def test_small_t_series_leading_order():
    # E_alpha(-t) = 1 - t/Gamma(1+alpha) + O(t^2).
    t = 1e-3
    for a in (0.5, 0.7, 0.9):
        approx = 1.0 - t / special.gamma(1.0 + a)
        assert np.isclose(ml.mlf_neg(a, t), approx, rtol=1e-4, atol=1e-7)


def test_vectorized_matches_scalar():
    a = 0.75
    t = np.array([0.0, 0.1, 1.0, 5.0, 20.0])
    vec = ml.mlf_neg(a, t)
    sca = np.array([ml.mlf_neg(a, float(x)) for x in t])
    assert np.allclose(vec, sca, rtol=1e-10, atol=1e-12)


def test_alpha_at_or_above_threshold_is_exact_exp():
    # alpha >= 0.99 is returned as exp(-t) exactly (documented behaviour: the near-delta
    # kernel is unresolvable there and alpha=1 IS the CP0 exponential lead lane).
    t = np.array([0.0, 0.5, 2.0, 10.0])
    for a in (0.99, 0.995, 1.0):
        assert np.allclose(ml.mlf_neg(a, t), np.exp(-t), rtol=1e-12, atol=1e-15)


def test_approaches_exp_at_small_t_in_spectral_region():
    # In the spectral region (alpha < 0.99), as alpha->1 the function -> exp(-t) where the
    # algebraic tail t^{-1}/Gamma(1-alpha) is negligible, i.e. at small t.
    t = np.array([0.05, 0.1, 0.3])
    near = ml.mlf_neg(0.97, t)
    assert np.allclose(near, np.exp(-t), rtol=1e-2, atol=2e-3)


def _quad_reference(alpha, t):
    """Independent adaptive-quadrature reference for E_alpha(-t) (the Mainardi integral)."""
    from scipy import integrate
    s = t ** (1.0 / alpha)
    def K(r):
        ra = r ** alpha
        return (np.sin(alpha * np.pi) / np.pi) * r ** (alpha - 1.0) / (ra * ra + 2 * ra * np.cos(alpha * np.pi) + 1.0)
    val, _ = integrate.quad(lambda r: np.exp(-r * s) * K(r), 0.0, np.inf, limit=400)
    return val


def test_matches_adaptive_quadrature():
    # The fast fixed-grid integral must match an independent ADAPTIVE quadrature across the
    # physiological CTRW regime (alpha in [0.4,0.95]) and a wide t range, including the tail.
    for a in (0.4, 0.6, 0.7, 0.85, 0.95):
        for t in (1.5, 2.0, 3.0, 5.0, 10.0, 20.0):
            got = float(ml.mlf_neg(a, np.array([t]))[0])
            ref = _quad_reference(a, t)
            # 5e-6 relative is ~3 orders below diffusion-MRI measurement noise (SNR 20-60);
            # the small-alpha corner has a ~2e-7 grid-truncation floor on the r^{alpha-1} tail.
            assert abs(got - ref) <= 5e-6 * abs(ref) + 1e-9, f"alpha={a} t={t}: {got} vs {ref}"
