"""caliper.conformal -- split-conformal / CQR coverage correction.

Numpy-only and estimator-agnostic. Works on the quantile arrays produced by any
estimator exposing ``predict_quantiles(signals, q_levels) -> (n, n_params,
n_levels)``.

The method is conformalized quantile regression (CQR; Romano, Patterson &
Candes 2019) applied independently to each symmetric pair of quantile levels.
For a central pair (q_lo, q_hi) with miss-rate ``a = 2 * q_lo``, the calibration
conformity score is

    E_i = max(q_lo(x_i) - y_i,  y_i - q_hi(x_i)),

and the correction is the finite-sample-adjusted (1 - a) empirical quantile of
the {E_i}. The corrected interval is [q_lo - Q, q_hi + Q]. This restores
*marginal* coverage to nominal under exchangeability; it does not by itself fix
*conditional* coverage.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "conformity_scores",
    "conformal_offset",
    "SplitConformalQuantile",
]


def conformity_scores(q_lo, q_hi, y) -> np.ndarray:
    """CQR conformity scores E_i = max(q_lo - y, y - q_hi)."""
    q_lo = np.asarray(q_lo, dtype=float)
    q_hi = np.asarray(q_hi, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.maximum(q_lo - y, y - q_hi)


def conformal_offset(scores, alpha: float) -> float:
    """Finite-sample (1 - alpha) quantile of the conformity scores.

    Uses the standard CQR level ceil((n+1)(1-alpha))/n, clipped to [0, 1].
    """
    scores = np.asarray(scores, dtype=float)
    n = scores.shape[0]
    if n == 0:
        return 0.0
    level = np.ceil((n + 1) * (1.0 - alpha)) / n
    level = min(max(level, 0.0), 1.0)
    return float(np.quantile(scores, level, method="higher"))


def _symmetric_pairs(n_levels: int):
    """Yield (j_lo, j_hi) index pairs symmetric about the centre."""
    for j in range(n_levels // 2):
        yield j, n_levels - 1 - j


@dataclass
class SplitConformalQuantile:
    """Split-conformal / CQR wrapper over an estimator's quantile output.

    Usage
    -----
    >>> cq = SplitConformalQuantile(q_levels)
    >>> cq.calibrate(q_pred_cal, y_cal)      # q_pred_cal: (n_cal, P, L)
    >>> q_corr = cq.apply(q_pred_test)       # (n_test, P, L)

    The correction is per-parameter and per symmetric level-pair. The median
    level (if present) is left unchanged. Output quantiles are re-sorted along
    the level axis to guarantee monotonicity.
    """

    q_levels: np.ndarray

    def __post_init__(self) -> None:
        self.q_levels = np.asarray(self.q_levels, dtype=float)
        if not np.all(np.diff(self.q_levels) > 0):
            raise ValueError("q_levels must be strictly ascending")
        # offsets_[p][(j_lo, j_hi)] = Q ; filled by calibrate()
        self.offsets_: list[dict[tuple[int, int], float]] = []
        self._n_params: int | None = None

    def calibrate(self, q_pred_cal, y_cal) -> "SplitConformalQuantile":
        q_pred_cal = np.asarray(q_pred_cal, dtype=float)
        y_cal = np.asarray(y_cal, dtype=float)
        n, P, L = q_pred_cal.shape
        if y_cal.shape != (n, P):
            raise ValueError("y_cal must be (n_cal, n_params)")
        if L != self.q_levels.shape[0]:
            raise ValueError("q_pred_cal level axis != len(q_levels)")
        self._n_params = P
        self.offsets_ = []
        for p in range(P):
            off: dict[tuple[int, int], float] = {}
            for j_lo, j_hi in _symmetric_pairs(L):
                a = 2.0 * self.q_levels[j_lo]  # nominal miss-rate of this pair
                scores = conformity_scores(
                    q_pred_cal[:, p, j_lo], q_pred_cal[:, p, j_hi], y_cal[:, p]
                )
                off[(j_lo, j_hi)] = conformal_offset(scores, a)
            self.offsets_.append(off)
        return self

    def apply(self, q_pred) -> np.ndarray:
        if not self.offsets_:
            raise RuntimeError("call calibrate() before apply()")
        q_pred = np.asarray(q_pred, dtype=float)
        n, P, L = q_pred.shape
        if P != self._n_params:
            raise ValueError("n_params mismatch with calibration")
        out = q_pred.copy()
        for p in range(P):
            for (j_lo, j_hi), Q in self.offsets_[p].items():
                out[:, p, j_lo] = q_pred[:, p, j_lo] - Q
                out[:, p, j_hi] = q_pred[:, p, j_hi] + Q
        # enforce monotonic, non-crossing quantiles along the level axis
        out = np.sort(out, axis=2)
        return out

    def calibrate_apply(self, q_pred_cal, y_cal, q_pred_test) -> np.ndarray:
        """Convenience: calibrate on one split, correct another."""
        return self.calibrate(q_pred_cal, y_cal).apply(q_pred_test)


if __name__ == "__main__":
    # Self-contained sanity demo: an over-confident Gaussian estimator that
    # conformal restores to nominal marginal coverage.
    from statistics import NormalDist

    from caliper import metrics as M

    rng = np.random.default_rng(0)
    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    z = np.array([NormalDist().inv_cdf(p) for p in levels])

    def make(n):
        mu = rng.normal(size=(n, 1))
        y = mu + rng.normal(size=(n, 1))           # true sigma = 1
        q = mu[:, :, None] + 0.5 * z[None, None, :]  # believes sigma = 0.5
        return y, q

    y_cal, q_cal = make(2000)
    y_te, q_te = make(4000)
    raw = M.score_quantiles(y_te, q_te, levels, alpha=0.1)[0]
    cq = SplitConformalQuantile(levels).calibrate(q_cal, y_cal)
    q_corr = cq.apply(q_te)
    cor = M.score_quantiles(y_te, q_corr, levels, alpha=0.1)[0]
    print(f"raw  coverage={raw.coverage:.3f} gap={raw.coverage_gap:+.3f}")
    print(f"conf coverage={cor.coverage:.3f} gap={cor.coverage_gap:+.3f}")
