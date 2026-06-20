"""Box-constrained NLLS fit + boundary-railing diagnostic.  [CP2 — spec below]

Reimplements the classical baseline whose failure mode Fashion's 54.7% number
reports: a four-parameter (S0, D, f, Dstar) box-constrained non-linear least-squares
fit (``scipy.optimize.least_squares``, trust-region reflective), with an asymptotic
Gaussian (Laplace/CRLB) covariance from the Jacobian -- over-confident by
construction, which is the point.

Boundary railing (pinned in ``manifest.RAILING``): a parameter is *railed* when
``|x_hat - bound| / (upper - lower) < rail_tol`` (rail_tol = 1e-3 primary, 1e-2
sensitivity). The railing **rate** for D* is the load-bearing diagnostic --
reproduced on synthetic cohorts (characterization) and on the real OSIPI abdomen
ROI (target T1: 54.7% +/- 5pp).

CP2 deliverables: ``NLLSEstimator`` (fit, sigma, railed mask, predict_quantiles),
``railing_rate(...)``. Box bounds, init, solver settings documented inline + in
METHODS.md (a Fashion completeness gap).
"""
from __future__ import annotations


class NLLSEstimator:  # [CP2]
    def __init__(self, *args, **kwargs):
        raise NotImplementedError("Gnomon CP2: clean-room box-constrained NLLS.")


def railing_rate(*args, **kwargs):  # [CP2]
    raise NotImplementedError("Gnomon CP2: D* boundary-railing rate.")
