"""Bayesian IVIM posteriors: Laplace + per-voxel MCMC.  [CP2 — spec below]

Reproduces Fashion's headline coverage table (targets T3a/b/c) by constructing two
posteriors over (D, Dstar, f) and reading TWO interval kinds from each:

* **Laplace** -- Gaussian at the MAP with curvature covariance. Its symmetric
  SD interval under-covers skewed, bound-pinned D* (target T3a: cov ~ 0.30).
* **MCMC** -- per-voxel sampler (clean-room random-walk / adaptive Metropolis over
  a Rician/Gaussian likelihood with the manifest's priors). From the SAME chain:
    - a Gaussian SD interval (still symmetric -> overconfident; target T3b ~ 0.67);
    - a 2.5/97.5 **quantile** interval (shape-correct -> near-nominal; target T3c
      ~ 0.94 for D*, and ~ 0.94 for D, f).

The headline mechanism: the right *shape* (quantiles), not a bigger SD, fixes D*.
CP2 documents the full sampler spec (proposal, chain length, burn-in, thinning,
seed, convergence check) -- training/fitting completeness Fashion was flagged for.
"""
from __future__ import annotations


class LaplacePosterior:  # [CP2]
    def __init__(self, *args, **kwargs):
        raise NotImplementedError("Gnomon CP2: Laplace posterior + SD interval.")


class MCMCPosterior:  # [CP2]
    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "Gnomon CP2: per-voxel MCMC; SD vs quantile intervals."
        )
