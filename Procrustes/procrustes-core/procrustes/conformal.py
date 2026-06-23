"""Split-conformal coverage on the target parameter.

Distribution-free split conformal gives finite-sample *marginal* coverage for
any base predictor under exchangeability -- including a misspecified bi-exp fit
(Lei et al. 2018).  Procrustes is about what that guarantee does NOT buy:
*conditional* coverage along the (latent, unobservable) misspecification axis.
"""
from __future__ import annotations

import numpy as np


def conformal_radius(abs_residuals: np.ndarray, alpha: float) -> float:
    """Split-conformal absolute-residual quantile (finite-sample inflated)."""
    r = np.sort(np.asarray(abs_residuals, float))
    n = len(r)
    k = min(int(np.ceil((n + 1) * (1.0 - alpha))), n)
    return float(r[k - 1])


def coverage(abs_residuals: np.ndarray, radius: float) -> float:
    return float((np.asarray(abs_residuals, float) <= radius).mean())


def split_indices(n: int, rng) -> tuple[np.ndarray, np.ndarray]:
    idx = np.arange(n)
    rng.shuffle(idx)
    return idx[: n // 2], idx[n // 2:]
