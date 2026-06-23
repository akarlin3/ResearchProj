"""Joint CTRW / fractional Bloch-Torrey forward model S(b; S0, D, alpha, beta) and its
finite-difference Jacobian. The model is

    S(b) = S0 * E_alpha( -(b D)^{beta/2} )                     (single diffusion time)

with E_alpha the one-parameter Mittag-Leffler function (mittag_leffler.mlf_neg). It reduces
EXACTLY to the CP0 stretched-exponential lead lane at the time-order alpha = 1:

    E_1(-x) = e^{-x}  =>  S(b) = S0 * exp(-(b D)^{beta/2}),

i.e. the CP0 heterogeneity exponent equals beta/2 with the time-order pinned to 1.
"""
import numpy as np

from levy import forward


def test_param_order_convention():
    assert forward.PARAM_NAMES_JOINT == ("S0", "D", "alpha", "beta")
    assert forward.IDX_JOINT == {"S0": 0, "D": 1, "alpha": 2, "beta": 3}


def test_b0_returns_S0():
    s = forward.forward_joint([0.0], np.array([2.3, 1.5e-3, 0.7, 1.8]))
    assert np.isclose(s[0], 2.3)


def test_reduces_to_cp0_stretched_exp_at_alpha_one():
    # At alpha=1, S(b) = S0 exp(-(bD)^{beta/2}) == CP0 stretched-exp with exponent beta/2.
    b = np.array([0.0, 200.0, 500.0, 1000.0, 2000.0, 3000.0])
    S0, D, beta = 1.0, 1.5e-3, 1.7
    joint = forward.forward_joint(b, np.array([S0, D, 1.0, beta]))
    cp0 = forward.signal(b, np.array([S0, D, beta / 2.0]))  # CP0 alpha == beta/2
    assert np.allclose(joint, cp0, rtol=1e-9, atol=1e-12)


def test_reduces_to_monoexponential_at_alpha1_beta2():
    b = np.array([0.0, 500.0, 1000.0, 2000.0])
    s = forward.forward_joint(b, np.array([1.0, 1.5e-3, 1.0, 2.0]))
    assert np.allclose(s, np.exp(-b * 1.5e-3), rtol=1e-9)


def test_signal_monotone_decreasing_and_bounded():
    b = np.linspace(0, 3000, 80)
    for alpha, beta in ((0.6, 1.4), (0.8, 1.8), (0.95, 1.95)):
        s = forward.forward_joint(b, np.array([1.0, 1.5e-3, alpha, beta]))
        assert np.all(s > 0.0)
        assert np.all(s <= 1.0 + 1e-12)
        assert np.all(np.diff(s) <= 1e-12)


def test_joint_jacobian_matches_finite_difference():
    b = np.array([0.0, 200.0, 600.0, 1200.0, 2000.0, 3000.0])
    theta = np.array([1.0, 1.5e-3, 0.75, 1.7])
    J = forward.jacobian_joint(b, theta)
    assert J.shape == (len(b), 4)
    steps = np.array([1e-6, 1e-9, 1e-5, 1e-5])
    for k in range(4):
        tp = theta.copy(); tp[k] += steps[k]
        tm = theta.copy(); tm[k] -= steps[k]
        fd = (forward.forward_joint(b, tp) - forward.forward_joint(b, tm)) / (2 * steps[k])
        assert np.allclose(J[:, k], fd, rtol=1e-3, atol=1e-6), f"param {k}"


def test_joint_jacobian_reduces_to_cp0_at_alpha1_beta2():
    # At alpha=1, beta=2: dS/dS0, dS/dD must match the analytic CP0 mono-exponential Jacobian
    # (alpha_cp0 = 1). The alpha/beta columns are the joint-model extras.
    b = np.array([0.0, 500.0, 1000.0, 2000.0, 3000.0])
    theta_j = np.array([1.0, 1.5e-3, 1.0, 2.0])
    Jj = forward.jacobian_joint(b, theta_j)
    Jc = forward.jacobian(b, np.array([1.0, 1.5e-3, 1.0]))  # CP0 (S0, D, alpha=1)
    assert np.allclose(Jj[:, forward.IDX_JOINT["S0"]], Jc[:, forward.IDX["S0"]], rtol=1e-4, atol=1e-7)
    assert np.allclose(Jj[:, forward.IDX_JOINT["D"]], Jc[:, forward.IDX["D"]], rtol=1e-3, atol=1e-6)
