"""Joint 4-parameter Fisher information / Cramer-Rao bounds for the CTRW forward model
theta = (S0, D, alpha, beta) -- the Phase-3 structural-identifiability layer.

Same statistical experiment as the CP0 lead lane (the parameters enter the likelihood only
through nu_i = S(b_i; theta) under Rician magnitude noise), now with BOTH the time-fractional
order alpha and the space-fractional order beta free:

    FIM_R(theta) = sum_i f(nu_i/sigma)/sigma^2 * g_i g_i^T,   g_i = dS(b_i)/dtheta  (4-vector)

with g_i from the finite-difference Jacobian ``forward.jacobian_joint`` and f the Rician info
factor (``noise.rician_info_factor``). The deliverable is the DEGENERACY structure of (alpha,
beta): the correlation rho_alpha_beta implied by FIM^{-1} (-> +-1 == the time and space orders
trade off and cannot be separated at a single diffusion time), the alpha-D / beta-D trade-offs,
and the FIM condition number. CRLB = identifiability, scoped to the regime; never impossibility.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import forward, noise


def fisher_matrix_joint(b, theta_joint, snr, model: str = "rician", dt=None):
    """4x4 Fisher information matrix for theta = (S0, D, alpha, beta).

    model="gaussian": (1/sigma^2) J^T J   (high-SNR reference)
    model="rician":   J^T diag(f(nu_i/sigma)/sigma^2) J   (honest finite-SNR)

    ``dt`` is the per-measurement diffusion-time ratio array (None -> single diffusion time,
    i.e. ones). A multi-diffusion-time design breaks the (alpha, beta) degeneracy.
    """
    b = np.asarray(b, dtype=float)
    theta = np.asarray(theta_joint, dtype=float)
    S0 = float(theta[0])
    sigma = noise.sigma_from_snr(S0, snr)
    dt = np.ones_like(b) if dt is None else np.asarray(dt, dtype=float)

    J = forward.jacobian_multidt(b, dt, theta)      # (m, 4)
    nu = forward.signal_multidt(b, dt, theta)       # (m,)

    if model == "gaussian":
        weights = np.full(nu.shape, 1.0 / (sigma * sigma))
    elif model == "rician":
        a = nu / sigma
        weights = noise.rician_info_factor(a) / (sigma * sigma)
    else:
        raise ValueError(f"unknown model {model!r} (use 'gaussian' or 'rician')")

    return J.T @ (weights[:, None] * J)


@dataclass(frozen=True)
class CRLBJointResult:
    """Joint-CRLB diagnostics at one (truth, b-design, SNR) cell."""

    theta: np.ndarray          # (S0, D, alpha, beta) ground truth
    snr: float
    sigma: float
    fim: np.ndarray            # 4x4 Fisher information matrix
    cov: np.ndarray            # FIM^{-1} = CRLB covariance
    model: str

    def _i(self, name: str) -> int:
        return forward.IDX_JOINT[name]

    @property
    def crlb(self) -> np.ndarray:
        return np.diag(self.cov)

    @property
    def crlb_alpha(self) -> float:
        return float(self.cov[self._i("alpha"), self._i("alpha")])

    @property
    def crlb_beta(self) -> float:
        return float(self.cov[self._i("beta"), self._i("beta")])

    @property
    def crlb_D(self) -> float:
        return float(self.cov[self._i("D"), self._i("D")])

    @property
    def se_alpha(self) -> float:
        return float(np.sqrt(max(self.crlb_alpha, 0.0)))

    @property
    def se_beta(self) -> float:
        return float(np.sqrt(max(self.crlb_beta, 0.0)))

    @property
    def cv_alpha(self) -> float:
        """Relative CRLB on alpha: sqrt(CRLB_alpha)/alpha."""
        return self.se_alpha / float(self.theta[self._i("alpha")])

    @property
    def cv_beta(self) -> float:
        """Relative CRLB on beta: sqrt(CRLB_beta)/beta."""
        return self.se_beta / float(self.theta[self._i("beta")])

    def _rho(self, a: str, c: str) -> float:
        ia, ic = self._i(a), self._i(c)
        denom = np.sqrt(self.cov[ia, ia] * self.cov[ic, ic])
        if denom <= 0:
            return float("nan")
        return float(self.cov[ia, ic] / denom)

    @property
    def rho_alpha_beta(self) -> float:
        """(alpha, beta) correlation implied by FIM^{-1}; +-1 == time/space orders degenerate."""
        return self._rho("alpha", "beta")

    @property
    def rho_alpha_D(self) -> float:
        return self._rho("alpha", "D")

    @property
    def rho_beta_D(self) -> float:
        return self._rho("beta", "D")

    @property
    def cond(self) -> float:
        """FIM condition number (2-norm). Large => the design cannot separate the parameters."""
        return float(np.linalg.cond(self.fim))


def crlb_joint(b, theta_joint, snr, model: str = "rician", dt=None) -> CRLBJointResult:
    """Compute the joint-CRLB diagnostics at one cell (pinv fallback if FIM singular).

    ``dt`` is the per-measurement diffusion-time ratio array (None -> single diffusion time).
    """
    b = np.asarray(b, dtype=float)
    theta = np.asarray(theta_joint, dtype=float)
    sigma = noise.sigma_from_snr(float(theta[0]), snr)
    fim = fisher_matrix_joint(b, theta, snr, model=model, dt=dt)
    try:
        cov = np.linalg.inv(fim)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(fim)
    return CRLBJointResult(theta=theta, snr=float(snr), sigma=float(sigma),
                           fim=fim, cov=cov, model=model)
