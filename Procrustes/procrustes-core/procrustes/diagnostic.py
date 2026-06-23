"""The observable misspecification diagnostic (the risk-bearing half).

Can a deployer DETECT that the bi-exp fit is misspecified -- and thus that the
trusted D map's conditional coverage is at risk -- from observables alone, with
no ground truth?  Gauge's naive-transfer monitor was a *hidden channel* for pure
dispersion (AUC ~ 0.50).  Procrustes instead reads the bi-exp fit's residual
STRUCTURE (reduced-chi2 / lag-1 autocorrelation / longest same-sign run).

Two CP-native measures, distinct from model-order goodness-of-fit (AIC/BIC) and
from the misspecified-CRB variance test (Wang-Tamir-Bush 2026, ASL):
  - detection AUC : worst-departure vs bi-exp-limit voxels,
  - rho(obs, |D-error|) : does the observable RANK voxels by their D-error, i.e.
                          predict *where* conditional coverage fails?
"""
from __future__ import annotations

import numpy as np

STATS = ("chi2_red", "ac1", "longrun")


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUC = P(statistic larger on `pos` than on `neg`)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    allv = np.concatenate([pos, neg])
    ranks = allv.argsort().argsort() + 1.0
    rp = ranks[: len(pos)].sum()
    return float((rp - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg)))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = np.asarray(x).argsort().argsort().astype(float)
    ry = np.asarray(y).argsort().argsort().astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def diagnose(fits: dict, truth: np.ndarray, limit_value, worst_value,
             test_idx: np.ndarray) -> dict:
    """AUC + |D-error| rank-correlation for each residual-structure statistic.

    ``fits`` maps departure value -> the dict returned by generators.fit_cohort.
    """
    out = {}
    for stat in STATS:
        out[f"auc_{stat}"] = auc(fits[worst_value][stat][test_idx],
                                 fits[limit_value][stat][test_idx])
        derr = np.abs(fits[worst_value]["Dhat"][test_idx] - truth[test_idx])
        out[f"rho_{stat}"] = spearman(fits[worst_value][stat][test_idx], derr)
    out["auc_best"] = max(out[f"auc_{s}"] for s in STATS)
    return out
