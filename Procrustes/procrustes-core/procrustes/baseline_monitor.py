"""The naive-transfer drift-monitor baseline (the hidden-channel baseline-to-beat).

Procrustes' diagnostic claim is "the residual-STRUCTURE diagnostic detects the
heavy-tail misspecification channel that breaks D's conditional coverage, beating
the naive deployment monitor". To make that claim REPRODUCIBLE rather than an
asserted '~0.501', this module *computes* the baseline on the same cohorts.

The baseline is the standard label-free deployment-drift monitor (the recipe in
Gauge's ``gauge/monitor.py`` and projMinos: a Mahalanobis OOD score on observable
SIGNAL-SHAPE summaries, optionally augmented with the fit-residual NORM). It is
re-implemented here clean-room on Procrustes' own observables -- it imports
nothing from Gauge -- because the shape/Mahalanobis recipe is generic, not
Gauge-specific IP.

Detection task (identical to the Procrustes diagnostic): separate worst-departure
voxels (positives) from bi-exp-limit / placebo voxels (negatives) using observ-
ables only. AUC ~ 0.5 means the channel is HIDDEN to that monitor.

  * ``signal_shape_features``  -- low-b slope, high-b slope, curvature, early drop
                                  (no model fit needed; observable in vivo).
  * ``mahalanobis_monitor_auc``-- Gauge/Minos-style Mahalanobis drift AUC on the
                                  shape features (+ optional residual norm).
"""
from __future__ import annotations

import numpy as np

from .diagnostic import auc


def signal_shape_features(signals: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Label-free signal-shape features (cf. Gauge ``_signal_shape_features``).

    Columns (slope_low, slope_high, curvature, early_drop):
      slope_low/high  -- log-signal decay rate over low-b (<=60) / high-b (>=200),
      curvature       -- slope_low - slope_high (the bi-exp perfusion signature),
      early_drop      -- 1 - S(b~30)/S0.
    Computed from the magnitude signal alone (no ground truth, no fit).
    """
    signals = np.atleast_2d(np.asarray(signals, float))
    b = np.asarray(b, float)
    s0 = signals[:, int(np.argmin(b))]
    s0 = np.where(s0 > 0, s0, signals.max(1))
    eps = 1e-6
    logs = np.log(np.clip(signals, eps, None))

    def _slope(mask):
        bb = b[mask]
        A = np.column_stack([bb, np.ones_like(bb)])
        coef, *_ = np.linalg.lstsq(A, logs[:, mask].T, rcond=None)
        return -coef[0]

    low, high = b <= 60.0, b >= 200.0
    slope_low = _slope(low) if low.sum() >= 2 else np.zeros(signals.shape[0])
    slope_high = _slope(high) if high.sum() >= 2 else np.zeros(signals.shape[0])
    j30 = int(np.argmin(np.abs(b - 30.0)))
    early_drop = 1.0 - signals[:, j30] / np.clip(s0, eps, None)
    return np.column_stack([slope_low, slope_high,
                            slope_low - slope_high, early_drop])


def _mahalanobis(feat, mu, inv_cov, sd):
    Z = (np.atleast_2d(np.asarray(feat, float)) - mu) / sd
    d = np.einsum("ni,ij,nj->n", Z, inv_cov, Z)
    return np.sqrt(np.maximum(d, 0.0))


def mahalanobis_monitor_auc(feat_neg_cal, feat_neg_test, feat_pos_test,
                            resid_norm=None) -> float:
    """Gauge/Minos-style drift-monitor detection AUC (worst vs placebo).

    Fit the placebo (negative) calibration distribution of observable features,
    score placebo-test (neg) and worst-test (pos) voxels by Mahalanobis distance,
    and return the separation AUC. ~0.5 => the channel is hidden to the monitor.
    If ``resid_norm`` is given as a tuple (neg_cal, neg_test, pos_test) of per-
    voxel residual norms, it is appended as an extra observable column (the full
    monitor; otherwise the shape-only monitor).
    """
    fnc = np.atleast_2d(np.asarray(feat_neg_cal, float))
    fnt = np.atleast_2d(np.asarray(feat_neg_test, float))
    fpt = np.atleast_2d(np.asarray(feat_pos_test, float))
    if resid_norm is not None:
        rnc, rnt, rpt = (np.asarray(x, float).reshape(-1, 1) for x in resid_norm)
        fnc = np.hstack([fnc, rnc]); fnt = np.hstack([fnt, rnt]); fpt = np.hstack([fpt, rpt])
    mu = fnc.mean(0)
    sd = fnc.std(0) + 1e-12
    Zc = (fnc - mu) / sd
    cov = np.atleast_2d(np.cov(Zc, rowvar=False)) + np.eye(Zc.shape[1]) * 1e-6
    inv_cov = np.linalg.inv(cov)
    score_neg = _mahalanobis(fnt, mu, inv_cov, sd)
    score_pos = _mahalanobis(fpt, mu, inv_cov, sd)
    return auc(score_pos, score_neg)
