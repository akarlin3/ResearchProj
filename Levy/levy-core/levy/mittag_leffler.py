r"""One-parameter Mittag-Leffler function E_alpha(-t) for real t >= 0, 0 < alpha <= 1.

This is the kernel of the joint CTRW / fractional Bloch-Torrey forward model
(``forward.forward_joint``): the diffusion signal at a single diffusion time is
``S(b)/S0 = E_alpha(-(b D)^{beta/2})`` (Magin, Ingo, Colon-Perez, Triplett & Mareci,
*Micropor. Mesopor. Mater.* 178:39-43, 2013 -- the CTRW signal interpolates a stretched
exponential and a power law). At alpha = 1 it reduces to ``exp(-t)`` (the CP0 lead lane).

We only ever need the NEGATIVE-REAL-ARGUMENT branch E_alpha(-t), t >= 0 -- the signal
decay. On that branch E_alpha(-t) is completely monotone (Pollard), which gives a stable
spectral (Mainardi) integral representation for 0 < alpha < 1:

    E_alpha(-t) = \int_0^\infty e^{-r t} K_alpha(r) dr,
    K_alpha(r)  = (1/pi) * r^{alpha-1} sin(alpha pi)
                  / (r^{2 alpha} + 2 r^alpha cos(alpha pi) + 1)   >= 0.

K_alpha integrates to 1 (so E_alpha(0)=1) and is non-negative for 0<alpha<1 (so the
function is positive and decreasing). For small t the alternating power series
``sum_k (-t)^k / Gamma(alpha k + 1)`` is exact and cheap; for larger t it loses precision
to cancellation, so we switch to the integral. alpha = 1 (and alpha within a hair of 1,
where the integral kernel degenerates to a delta at r=1) is special-cased to exp(-t).
"""
from __future__ import annotations

import numpy as np

# alpha within this of 1 -> exp(-t). Two reasons: (1) at alpha=1 the function IS the CP0
# exponential lead lane (exact reduction); (2) as alpha->1 the spectral kernel collapses to a
# near-delta at r=1 that the fixed grid cannot resolve (validated reliable only for
# alpha <= 0.99). Physiological CTRW alpha <= 0.95 and the joint model caps its working range
# at 0.98, so the spectral branch is always exercised in its validated-accurate region.
_ALPHA_ONE_TOL = 1e-2

# Fixed log-spaced quadrature grid in r = exp(u) for the Mainardi integral. The spectral
# integral is used for ALL t > 0 (no series/integral crossover -> the signal is a single
# smooth function of the parameters, which keeps the finite-difference Jacobian clean). Grid
# is dense/wide enough that the alpha=1/2 closed form (erfcx) and an independent adaptive
# quadrature both match to < 1e-6 across alpha in (0,1) and t in [0, 20].
# N=3201 over [-40,40] (du=0.025) matches the alpha=1/2 closed form and an independent
# adaptive quadrature to < 5e-6 across alpha in (0,0.99) and t in [0,20] -- i.e. ~3 orders
# below diffusion-MRI measurement noise -- while staying ~3x faster than a denser grid (which
# matters for the parametric-bootstrap MLE inner loop). The quadrature error is a SMOOTH
# function of the parameters, so it cancels in the finite-difference Jacobian.
_U = np.linspace(-40.0, 40.0, 3201)
_R = np.exp(_U)
_DU = _U[1] - _U[0]


def _spectral(alpha: float, t: np.ndarray) -> np.ndarray:
    r"""E_alpha(-t) via the Mainardi completely-monotone integral (0<alpha<1).

    The spectral (Pollard) representation is for the relaxation function
    ``E_alpha(-s^alpha) = \int_0^\infty e^{-r s} K_alpha(r) dr``. To get ``E_alpha(-t)``
    we substitute ``s = t^{1/alpha}`` so the exponential argument is ``r * t^{1/alpha}``.
    """
    r = _R
    # K_alpha(r) on the grid (independent of t).
    ra = r ** alpha
    K = (np.sin(alpha * np.pi) / np.pi) * (r ** (alpha - 1.0)) / (ra * ra + 2.0 * ra * np.cos(alpha * np.pi) + 1.0)
    s = np.asarray(t, dtype=float) ** (1.0 / alpha)          # E_alpha(-t) = E_alpha(-(s)^alpha), s=t^{1/alpha}
    # integrand in u: e^{-r s} K(r) * r  (because dr = r du); trapezoid over the fixed grid.
    # vectorize over t: shape (len(t), len(r))
    expo = np.exp(-np.outer(s, r))
    integ = expo * (K * r)[None, :]
    return np.trapezoid(integ, dx=_DU, axis=1)


def mlf_neg(alpha: float, t) -> np.ndarray | float:
    """E_alpha(-t) for t >= 0 and 0 < alpha <= 1. Scalar or array ``t``; returns same shape.

    Negative-argument branch only (the diffusion signal decay). At alpha within
    ``_ALPHA_ONE_TOL`` of 1 the exact exponential is returned.
    """
    scalar = np.isscalar(t)
    t = np.atleast_1d(np.asarray(t, dtype=float))
    if np.any(t < 0):
        raise ValueError("mlf_neg is the negative-argument branch: t must be >= 0")
    alpha = float(alpha)
    if not (0.0 < alpha <= 1.0 + _ALPHA_ONE_TOL):
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}")

    if alpha >= 1.0 - _ALPHA_ONE_TOL:
        out = np.exp(-t)
    else:
        out = np.empty_like(t)
        zero = t == 0.0
        out[zero] = 1.0          # E_alpha(0)=1 exactly (the integral's int K dr=1 up to ~1e-7)
        if np.any(~zero):
            out[~zero] = _spectral(alpha, t[~zero])
    return float(out[0]) if scalar else out
