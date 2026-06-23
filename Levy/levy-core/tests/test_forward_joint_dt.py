"""Two-diffusion-time CTRW forward model. At fixed Delta the time/space orders are confounded
(degenerate); exposing a SECOND diffusion time separates them, because the Delta-dependence
carries the exponent (alpha - beta/2):

    S(b, Delta) = S0 * E_alpha( -(b D)^{beta/2} * (Delta/Delta0)^{alpha - beta/2} ).

At Delta = Delta0 (ratio 1) this is exactly the single-Delta joint model.
"""
import numpy as np

from levy import forward


def test_reduces_to_single_dt_at_ratio_one():
    b = np.array([0.0, 300.0, 800.0, 1600.0, 2500.0])
    theta = np.array([1.0, 1.5e-3, 0.8, 1.7])
    s_dt = forward.forward_joint_dt(b, 1.0, theta)
    s_single = forward.forward_joint(b, theta)
    assert np.allclose(s_dt, s_single, rtol=1e-12, atol=1e-14)


def test_b0_returns_S0_any_dt():
    for ratio in (0.5, 1.0, 2.0, 3.0):
        s = forward.forward_joint_dt([0.0], ratio, np.array([1.3, 1.5e-3, 0.7, 1.8]))
        assert np.isclose(s[0], 1.3)


def test_dt_ratio_changes_signal_when_alpha_ne_beta_half():
    # Unless alpha == beta/2 the Delta-exponent (alpha - beta/2) is nonzero, so a different
    # diffusion time genuinely changes the attenuation (this is what breaks the degeneracy).
    b = np.array([1000.0, 2000.0])
    theta = np.array([1.0, 1.5e-3, 0.8, 1.7])   # alpha - beta/2 = 0.8 - 0.85 = -0.05 != 0
    s1 = forward.forward_joint_dt(b, 1.0, theta)
    s2 = forward.forward_joint_dt(b, 2.0, theta)
    assert not np.allclose(s1, s2, rtol=1e-3)


def test_dt_ratio_irrelevant_when_alpha_eq_beta_half():
    # When alpha == beta/2 the Delta-exponent vanishes -> diffusion time has no effect (the
    # confound is exact). This is the ridge direction of the degeneracy.
    b = np.array([1000.0, 2000.0])
    theta = np.array([1.0, 1.5e-3, 0.85, 1.7])  # alpha = beta/2 = 0.85
    s1 = forward.forward_joint_dt(b, 1.0, theta)
    s2 = forward.forward_joint_dt(b, 2.5, theta)
    assert np.allclose(s1, s2, rtol=1e-9)
