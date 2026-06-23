"""Joint Rician MLE of (S0, D, alpha, beta) and its parametric bootstrap -- the empirical,
finite-sample companion to the analytic joint CRLB (fisher_joint.py).

Why both, again: the joint FIM gives the *information-geometry* degeneracy (rho_alpha_beta, the
condition number) at the truth; the parametric bootstrap shows what actually happens to the
joint MLE under noise. At a single diffusion time the (alpha_hat, beta_hat) bootstrap cloud
collapses onto a RIDGE -- strong anticorrelation, wide marginals -- which is the honest
finite-sample face of the structural degeneracy. A two-diffusion-time design breaks it.

Everything is seeded (explicit Generators) so the CIs reproduce.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize

from . import forward, noise

# Optimisation works in a transformed space keeping S0>0, D>0, alpha in (0, ALPHA_MAX],
# beta in (BETA_LO, BETA_HI):  p = (log S0, log D, logit(alpha/ALPHA_MAX), logit((beta-1)/(2-1))).
ALPHA_MAX = forward.ALPHA_JOINT_MAX   # 0.98 -- stay in the Mittag-Leffler-resolved region
BETA_LO, BETA_HI = 1.0, 2.0


def two_dt_design(b_max: float = 2500.0, n_b: int = 8, ratios=(1.0,)):
    """Build a (b, dt) acquisition: the same 0..b_max ramp repeated at each diffusion-time ratio.

    ratios=(1.0,) is a single-diffusion-time clinical design; (1.0, 2.5) adds a second Delta.
    Returns (b_arr, dt_arr) of equal length n_b*len(ratios).
    """
    from .wall import default_b_design
    bs, dts = [], []
    for r in ratios:
        bb = default_b_design(b_max=b_max, n_b=n_b)
        bs.append(bb)
        dts.append(np.full(bb.shape, float(r)))
    return np.concatenate(bs), np.concatenate(dts)


def _pack(theta):
    S0, D, alpha, beta = (float(theta[0]), float(theta[1]), float(theta[2]), float(theta[3]))
    a = np.clip(alpha / ALPHA_MAX, 1e-6, 1 - 1e-6)
    bb = np.clip((beta - BETA_LO) / (BETA_HI - BETA_LO), 1e-6, 1 - 1e-6)
    return np.array([np.log(S0), np.log(D), np.log(a / (1 - a)), np.log(bb / (1 - bb))])


def _unpack(p):
    S0 = np.exp(np.clip(p[0], -6.9, 6.9))
    D = np.exp(np.clip(p[1], -13.8, -2.3))
    a = 1.0 / (1.0 + np.exp(-np.clip(p[2], -30.0, 30.0)))
    bb = 1.0 / (1.0 + np.exp(-np.clip(p[3], -30.0, 30.0)))
    return np.array([S0, D, a * ALPHA_MAX, BETA_LO + bb * (BETA_HI - BETA_LO)])


def nll(theta, b, dt, M, sigma):
    """Negative Rician log-likelihood of magnitude data ``M`` for design (b, dt)."""
    nu = forward.signal_multidt(b, dt, theta)
    return -float(np.sum(noise.rician_logpdf(M, nu, sigma)))


def mle_joint(b, dt, M, sigma, theta0):
    """Joint Rician MLE of (S0, D, alpha, beta). Returns (theta_hat, nll_min, success)."""
    def obj(p):
        return nll(_unpack(p), b, dt, M, sigma)

    res = optimize.minimize(obj, _pack(theta0), method="Nelder-Mead",
                            options=dict(maxiter=4000, xatol=1e-7, fatol=1e-9))
    return _unpack(res.x), float(res.fun), bool(res.success)


@dataclass(frozen=True)
class BootstrapJointResult:
    truth: np.ndarray
    alpha_hats: np.ndarray
    beta_hats: np.ndarray
    alpha_ci: tuple          # (lo, hi) percentile CI on alpha_hat
    beta_ci: tuple
    corr_alpha_beta: float   # empirical Pearson correlation of (alpha_hat, beta_hat) -> the ridge
    alpha_bias: float
    beta_bias: float

    @property
    def alpha_rel_width(self) -> float:
        return float((self.alpha_ci[1] - self.alpha_ci[0]) / self.truth[2])

    @property
    def beta_rel_width(self) -> float:
        return float((self.beta_ci[1] - self.beta_ci[0]) / self.truth[3])


def parametric_bootstrap_joint(truth, b, dt, snr, n_boot, rng, level=0.95):
    """Refit (alpha, beta) on ``n_boot`` Rician replicates at ``truth`` for design (b, dt).

    Returns the (alpha_hat, beta_hat) clouds, marginal percentile CIs, and their empirical
    correlation -- the finite-sample signature of the (alpha, beta) degeneracy.
    """
    truth = np.asarray(truth, dtype=float)
    b = np.asarray(b, dtype=float)
    dt = np.asarray(dt, dtype=float)
    sigma = noise.sigma_from_snr(float(truth[0]), snr)
    nu = forward.signal_multidt(b, dt, truth)

    alpha_hats = np.empty(n_boot)
    beta_hats = np.empty(n_boot)
    for k in range(n_boot):
        M = noise.rician_sample(nu, sigma, rng)
        theta_hat, _, _ = mle_joint(b, dt, M, sigma, truth)
        alpha_hats[k] = theta_hat[forward.IDX_JOINT["alpha"]]
        beta_hats[k] = theta_hat[forward.IDX_JOINT["beta"]]

    lo_q, hi_q = (1 - level) / 2, 1 - (1 - level) / 2
    a_ci = tuple(np.quantile(alpha_hats, [lo_q, hi_q]))
    b_ci = tuple(np.quantile(beta_hats, [lo_q, hi_q]))
    if np.std(alpha_hats) > 0 and np.std(beta_hats) > 0:
        corr = float(np.corrcoef(alpha_hats, beta_hats)[0, 1])
    else:
        corr = float("nan")
    return BootstrapJointResult(
        truth=truth, alpha_hats=alpha_hats, beta_hats=beta_hats,
        alpha_ci=(float(a_ci[0]), float(a_ci[1])), beta_ci=(float(b_ci[0]), float(b_ci[1])),
        corr_alpha_beta=corr,
        alpha_bias=float(np.mean(alpha_hats) - truth[forward.IDX_JOINT["alpha"]]),
        beta_bias=float(np.mean(beta_hats) - truth[forward.IDX_JOINT["beta"]]),
    )
