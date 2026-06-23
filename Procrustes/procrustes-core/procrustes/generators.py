"""Ground-truth cohorts (from Lattice) + the bi-exponential estimator under test.

Procrustes never owns a generator: the non-bi-exponential truth is drawn from
Lattice's clean-room, seed-generated families.  The *estimator* is a bounded
bi-exponential NLLS fit -- deliberately MISSPECIFIED whenever the truth is not
bi-exponential.  We also extract, from the bi-exp fit residual, the observable
goodness-of-fit features the misspecification diagnostic is built on.

Key reuse property (relied on by the gate): for a fixed ``seed``, Lattice draws
the base ``(D, Dstar, f)`` identically across families and across ``extra`` knob
values -- so we observe the *same voxel population* under increasing departure,
with the bi-exp limit (e.g. beta=1, cv=0, g=0) as a built-in placebo.
"""
from __future__ import annotations

import numpy as np

from .deps import ensure_lattice

ensure_lattice()
from lattice import make_cohort                       # noqa: E402
from lattice.selfcheck import fit_biexp_nlls          # noqa: E402


def _biexp(b: np.ndarray, theta) -> np.ndarray:
    D, Dstar, f = theta
    return f * np.exp(-b * Dstar) + (1.0 - f) * np.exp(-b * D)


def cohort_at(family: str, knob: str, value: float, cfg) -> "object":
    """A Lattice cohort for ``family`` with departure ``knob=value``."""
    return make_cohort(family=family, n=cfg.n, snr=cfg.snr, seed=cfg.seed,
                       prior=cfg.prior, noise=cfg.noise, extra={knob: value})


def fit_cohort(cohort, snr: float) -> dict:
    """Bi-exp NLLS fit of every voxel; return point estimates + residual features.

    Residual features (observable, no ground truth):
      rss       -- residual sum of squares of the bi-exp fit,
      chi2_red  -- reduced chi-square at known sigma = S0/SNR,
      ac1       -- lag-1 autocorrelation of the b-ordered residual (structure),
      longrun   -- longest same-sign residual run (structure).
    Under correct specification residuals are ~white (ac1~0, short runs); under
    misspecification they carry coherent shape.
    """
    b = np.asarray(cohort.bvalues, float)
    sigma2 = (1.0 / snr) ** 2
    dof = max(len(b) - 3, 1)
    Dhat = np.empty(len(cohort))
    Dstarhat = np.empty(len(cohort))
    fhat = np.empty(len(cohort))
    rss = np.empty(len(cohort))
    chi2 = np.empty(len(cohort))
    ac1 = np.empty(len(cohort))
    longrun = np.empty(len(cohort))

    for i, s in enumerate(cohort.signals):
        theta = fit_biexp_nlls(b, s)
        r = _biexp(b, theta) - np.asarray(s, float)
        Dhat[i], Dstarhat[i], fhat[i] = theta
        rss[i] = float(np.sum(r * r))
        chi2[i] = rss[i] / (dof * sigma2)
        ac1[i] = float(np.corrcoef(r[:-1], r[1:])[0, 1]) if np.std(r) > 1e-12 else 0.0
        longrun[i] = _longest_sign_run(r)

    return dict(Dhat=Dhat, Dstarhat=Dstarhat, fhat=fhat,
                rss=rss, chi2_red=chi2, ac1=ac1, longrun=longrun)


def _longest_sign_run(r: np.ndarray) -> float:
    sign = np.sign(r)
    sign[sign == 0] = 1
    best = run = 1
    for k in range(1, len(sign)):
        run = run + 1 if sign[k] == sign[k - 1] else 1
        best = max(best, run)
    return float(best)
