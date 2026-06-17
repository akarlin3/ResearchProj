"""Caliper -- an IVIM uncertainty-quantification calibration toolkit.

The public surface is intentionally small:

* ``caliper.metrics``  -- the model-agnostic calibration ruler (numpy only).
* ``caliper.forward``  -- bi-exponential IVIM forward model + synthetic cohorts.
* ``caliper.conformal``-- split-conformal / CQR coverage correction.
* ``caliper.estimator_maf`` -- a masked-autoregressive-flow posterior (needs
  the optional ``[estimator]`` extra: torch).

Only ``metrics`` and ``forward`` and ``conformal`` are imported eagerly; the
torch-dependent estimator is imported lazily so the core stays numpy-only.
"""
from __future__ import annotations

from . import conformal, forward, metrics

__all__ = ["metrics", "forward", "conformal", "__version__"]
__version__ = "0.1.0"
