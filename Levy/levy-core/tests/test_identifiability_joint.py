"""Joint Rician MLE of (S0, D, alpha, beta) and its parametric bootstrap.

Two regimes, two behaviours:
  * a TWO-diffusion-time design separates the time/space orders -> the MLE recovers (alpha, beta);
  * a SINGLE-diffusion-time design is degenerate -> the bootstrap (alpha_hat, beta_hat) collapse
    onto a ridge (wide marginals, strong anticorrelation). This IS the Phase-3 finding.
"""
import numpy as np

from levy import identifiability_joint as ij, seeding


TRUTH = np.array([1.0, 1.5e-3, 0.80, 1.7])


def test_two_dt_design_recovers_truth_at_high_snr():
    # Two diffusion times (ratios 1.0 and 2.5), 10 b-values each, low noise -> identifiable.
    b, dt = ij.two_dt_design(b_max=2500.0, n_b=10, ratios=(1.0, 2.5))
    rng = seeding.make_rng(7)
    from levy import forward, noise
    sigma = noise.sigma_from_snr(1.0, 300.0)
    M = noise.rician_sample(forward.signal_multidt(b, dt, TRUTH), sigma, rng)
    theta0 = np.array([1.0, 1.2e-3, 0.7, 1.6])
    theta_hat, _, ok = ij.mle_joint(b, dt, M, sigma, theta0)
    assert ok
    assert abs(theta_hat[2] - 0.80) < 0.05, f"alpha_hat={theta_hat[2]}"
    assert abs(theta_hat[3] - 1.70) < 0.10, f"beta_hat={theta_hat[3]}"


def test_single_dt_bootstrap_shows_ridge():
    # One diffusion time: the (alpha_hat, beta_hat) cloud is a degenerate ridge.
    b, dt = ij.two_dt_design(b_max=2500.0, n_b=8, ratios=(1.0,))  # single dt
    res = ij.parametric_bootstrap_joint(TRUTH, b, dt, snr=40.0, n_boot=40, rng=seeding.make_rng(3))
    assert res.alpha_hats.shape == (40,)
    assert res.beta_hats.shape == (40,)
    # strong anticorrelation between alpha_hat and beta_hat (the ridge)
    assert res.corr_alpha_beta < -0.5, f"corr={res.corr_alpha_beta}"
    # marginal CIs are wide (degeneracy): alpha CI width is a large fraction of alpha
    assert res.alpha_ci[1] - res.alpha_ci[0] > 0.15


def test_bootstrap_ci_brackets_and_orders():
    b, dt = ij.two_dt_design(b_max=2500.0, n_b=8, ratios=(1.0,))
    res = ij.parametric_bootstrap_joint(TRUTH, b, dt, snr=40.0, n_boot=30, rng=seeding.make_rng(5))
    assert res.alpha_ci[0] <= res.alpha_ci[1]
    assert res.beta_ci[0] <= res.beta_ci[1]
    assert -1.0 - 1e-9 <= res.corr_alpha_beta <= 1.0 + 1e-9
