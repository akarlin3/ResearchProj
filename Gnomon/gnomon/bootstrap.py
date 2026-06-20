"""Bootstrap confidence intervals for load-bearing numbers.  [CP2 — spec below]

Guardrail 6 requires bootstrap CIs on every load-bearing number (railing rate,
coverages, ECE/sharpness gaps). Clean-room, numpy-only resampling (percentile + BCa
where the statistic is biased), seeded from ``manifest.BOOTSTRAP``. Used by CP3 to
decide PASS/FAIL: a target passes if the claimed value lands inside the rebuild's CI
(or the point estimate is within the frozen tolerance), and directional gaps pass
only if their CI excludes 0.
"""
from __future__ import annotations


def bootstrap_ci(*args, **kwargs):  # [CP2]
    raise NotImplementedError("Gnomon CP2: seeded bootstrap CI (percentile/BCa).")
