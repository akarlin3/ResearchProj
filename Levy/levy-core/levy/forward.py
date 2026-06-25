"""The diffusion-MRI signal-decay forward model S(b; theta) and its Jacobian.

LEAD LANE -- stretched-exponential (single fractional order alpha):

    S(b; S0, D, alpha) = S0 * exp( -(b * D)^alpha )                          (Bennett 2003)

with alpha in (0, 1] the stretching / intra-voxel-heterogeneity exponent, D the
diffusion coefficient (units mm^2/s), b the diffusion weighting (units s/mm^2), S0
the non-diffusion-weighted amplitude. alpha = 1 recovers the mono-exponential
S0 * exp(-b D).

This is the object whose identifiability Levy bounds. Crucially, the parameter enters
the likelihood ONLY through the b-indexed *signal attenuation* S(b; theta) -- NOT through
a trajectory / mean-squared-displacement / increment process. That is the forward-model
difference from the fBm/Hurst CRB of Coeurjolly-Istas (2001): same "fractional exponent"
word, different statistical experiment.

The joint CTRW / fractional Bloch-Torrey two-exponent model (Phase 3) is built in
``forward_joint``. At a single diffusion time the CTRW signal (Magin, Ingo, Colon-Perez,
Triplett & Mareci, *Micropor. Mesopor. Mater.* 178:39-43, 2013, Eq. 6) reduces -- via the
narrow-pulse relation b = q^2 * (Delta - delta/3) at fixed Delta, so |q|^beta proportional
to b^{beta/2} -- to

    S(b; S0, D, alpha, beta) = S0 * E_alpha( -(b * D)^{beta/2} )

with E_alpha the one-parameter Mittag-Leffler function (``mittag_leffler.mlf_neg``), the
time-fractional (Caputo) order alpha in (0, 1] and the space-fractional (Riesz) order
beta in (1, 2]. It reduces EXACTLY to the stretched-exponential lead lane at alpha = 1
(E_1(-x) = e^{-x}): the CP0 heterogeneity exponent equals beta/2 with the time-order
pinned to 1. Phase 3 prices whether alpha (time) and beta (space) are separately
recoverable at a single clinical diffusion time, or degenerate.

PARAMETER ORDER conventions (used by the Jacobians, Fisher matrices, and CRLB indexing):
    lead lane:  theta = (S0, D, alpha)          indices 0, 1, 2
    joint:      theta = (S0, D, alpha, beta)     indices 0, 1, 2, 3
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import mittag_leffler

PARAM_NAMES = ("S0", "D", "alpha")
IDX = {name: i for i, name in enumerate(PARAM_NAMES)}

# Joint CTRW / fractional Bloch-Torrey parameter order (Phase 3): theta = (S0, D, alpha, beta)
PARAM_NAMES_JOINT = ("S0", "D", "alpha", "beta")
IDX_JOINT = {name: i for i, name in enumerate(PARAM_NAMES_JOINT)}

# The joint model's time-order alpha is evaluated through the Mittag-Leffler kernel, which is
# resolved only up to ~0.99 (see mittag_leffler). Cap the WORKING alpha range safely below
# that; physiological CTRW alpha <= 0.95 (Magin/Ingo 2013), so 0.98 is a generous ceiling.
ALPHA_JOINT_MAX = 0.98


@dataclass(frozen=True)
class StretchedExp:
    """Stretched-exponential ground-truth parameters and a fixed b-value design.

    Parameters
    ----------
    S0, D, alpha:
        Ground-truth parameters. D in mm^2/s (e.g. 1.5e-3 for tissue), alpha in (0, 1].
    """

    S0: float = 1.0
    D: float = 1.5e-3
    alpha: float = 0.75

    @property
    def theta(self) -> np.ndarray:
        return np.array([self.S0, self.D, self.alpha], dtype=float)


def signal(b, theta):
    """S(b; theta) = S0 * exp(-(b D)^alpha).  ``b`` array-like, ``theta`` = (S0, D, alpha)."""
    b = np.asarray(b, dtype=float)
    S0, D, alpha = float(theta[0]), float(theta[1]), float(theta[2])
    u = np.power(b * D, alpha)  # (b D)^alpha; (b=0) -> 0
    return S0 * np.exp(-u)


def jacobian(b, theta):
    """Closed-form dS/dtheta at each b. Returns array of shape (len(b), 3).

    With u = (b D)^alpha and S = S0 exp(-u):
        dS/dS0    =  S / S0
        dS/dD     = -S * alpha * u / D
        dS/dalpha = -S * u * ln(b D)
    At b = 0: u = 0 so dS/dD = 0 and dS/dalpha = 0 (only S0 is informed).
    """
    b = np.asarray(b, dtype=float)
    S0, D, alpha = float(theta[0]), float(theta[1]), float(theta[2])
    bD = b * D
    u = np.power(bD, alpha)
    S = S0 * np.exp(-u)

    dS0 = S / S0
    dD = -S * alpha * u / D
    # ln(bD) is -inf at b=0, but u=0 there so u*ln(bD) -> 0. Evaluate log only where bD>0
    # to avoid 0*(-inf)=nan; the derivative is exactly 0 at b=0 (only S0 is informed).
    safe_bD = np.where(bD > 0.0, bD, 1.0)
    ln_bD = np.log(safe_bD)
    dalpha = np.where(bD > 0.0, -S * u * ln_bD, 0.0)

    return np.stack([dS0, dD, dalpha], axis=-1)


def forward_joint(b, theta_joint):
    """Joint CTRW / fractional Bloch-Torrey forward model (Phase 3).

        S(b) = S0 * E_alpha( -(b D)^{beta/2} )

    theta_joint = (S0, D, alpha, beta), with the time-fractional order alpha in (0, 1] and
    the space-fractional order beta in (1, 2]. E_alpha is the Mittag-Leffler kernel
    (``mittag_leffler.mlf_neg``). At alpha = 1 this is exactly the CP0 stretched-exponential
    with heterogeneity exponent beta/2 (E_1(-x) = e^{-x}).
    """
    b = np.asarray(b, dtype=float)
    S0, D, alpha, beta = (float(theta_joint[0]), float(theta_joint[1]),
                          float(theta_joint[2]), float(theta_joint[3]))
    arg = np.power(b * D, beta / 2.0)            # (b D)^{beta/2}; (b=0) -> 0 -> E_alpha(0)=1
    return S0 * mittag_leffler.mlf_neg(alpha, arg)


def forward_joint_dt(b, dt_ratio, theta_joint):
    """Two-diffusion-time CTRW signal at diffusion time Delta = dt_ratio * Delta0:

        S(b, Delta) = S0 * E_alpha( -(b D)^{beta/2} * (Delta/Delta0)^{alpha - beta/2} ).

    ``dt_ratio`` = Delta/Delta0 (1.0 -> the single-Delta ``forward_joint``); a scalar, or an
    array broadcasting against ``b`` for a multi-diffusion-time design. The Delta-exponent
    (alpha - beta/2) is exactly what separates the time order alpha from the space order beta:
    at a single Delta they are confounded, and only when alpha != beta/2 does a different Delta
    change the attenuation. (Magin/Ingo separated alpha and beta precisely by acquiring at more
    than one experiment family.)
    """
    b = np.asarray(b, dtype=float)
    dt_ratio = np.asarray(dt_ratio, dtype=float)
    S0, D, alpha, beta = (float(theta_joint[0]), float(theta_joint[1]),
                          float(theta_joint[2]), float(theta_joint[3]))
    arg = np.power(b * D, beta / 2.0) * (dt_ratio ** (alpha - beta / 2.0))
    return S0 * mittag_leffler.mlf_neg(alpha, arg)


def signal_multidt(b, dt, theta_joint):
    """Multi-measurement CTRW signal: per-measurement b-values ``b`` and diffusion-time ratios
    ``dt`` (same shape). Each S_i = S0 E_alpha(-(b_i D)^{beta/2} dt_i^{alpha-beta/2}).

    A single-diffusion-time design is ``dt = ones_like(b)`` (== ``forward_joint``).
    """
    return forward_joint_dt(b, dt, theta_joint)


def jacobian_multidt(b, dt, theta_joint, rel_step: float = 1e-6):
    """Central finite-difference dS/dtheta_joint for a multi-diffusion-time design. (len(b), 4)."""
    b = np.asarray(b, dtype=float)
    dt = np.asarray(dt, dtype=float)
    theta = np.asarray(theta_joint, dtype=float)
    J = np.empty((b.shape[0], 4))
    floors = np.array([1e-6, 1e-9, 1e-5, 1e-5])
    for k in range(4):
        h = max(rel_step * abs(theta[k]), floors[k])
        tp = theta.copy(); tp[k] += h
        tm = theta.copy(); tm[k] -= h
        J[:, k] = (signal_multidt(b, dt, tp) - signal_multidt(b, dt, tm)) / (2.0 * h)
    return J


def jacobian_joint(b, theta_joint, rel_step: float = 1e-6):
    """dS/dtheta_joint at each b via central finite differences. Shape (len(b), 4).

    The Mittag-Leffler order-derivative (dE_alpha/dalpha) has no elementary closed form, so
    the whole 4-parameter Jacobian is taken by central differences with a per-parameter step
    scaled to each parameter's magnitude. The signal is a single smooth Mittag-Leffler call
    (no series/integral crossover), which keeps these differences clean. Validated against the
    analytic CP0 Jacobian at the alpha=1, beta=2 mono-exponential reduction.
    """
    b = np.asarray(b, dtype=float)
    theta = np.asarray(theta_joint, dtype=float)
    J = np.empty((b.shape[0], 4))
    # absolute steps: floor each so tiny parameters (D ~ 1e-3) still get a sensible step.
    floors = np.array([1e-6, 1e-9, 1e-5, 1e-5])
    for k in range(4):
        h = max(rel_step * abs(theta[k]), floors[k])
        tp = theta.copy(); tp[k] += h
        tm = theta.copy(); tm[k] -= h
        J[:, k] = (forward_joint(b, tp) - forward_joint(b, tm)) / (2.0 * h)
    return J
