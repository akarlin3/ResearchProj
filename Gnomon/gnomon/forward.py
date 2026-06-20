"""Clean-room IVIM bi-exponential forward model.  [CP2 — spec below]

Implements the standard intravoxel-incoherent-motion signal equation (Le Bihan
et al., Radiology 1988), reimplemented from the physics, not copied from Fashion or
Caliper:

    S(b)/S0 = (1 - f) * exp(-b * D)  +  f * exp(-b * (D + Dstar))

with D the tissue diffusivity, Dstar the pseudo-diffusion coefficient, f the
perfusion fraction, b the diffusion weighting (s/mm^2). Dstar >> D, which is exactly
why the fast compartment is weakly identified at clinical-sparse sampling -- the
mechanism behind the boundary-railing and the under-coverage Gnomon must reproduce.

CP2 deliverables (each documented inline when implemented):
* ``ivim(b, D, Dstar, f, s0=1.0)`` -- vectorized forward signal.
* ``jacobian(...)`` -- analytic Jacobian wrt (S0, D, f, Dstar) for NLLS + the
  Gaussian/Laplace CRLB covariance (the approximation whose D*-skew weakness CP4
  documents).
* self-consistency: clean-signal round-trip recovers truth (CP2 gate).
"""
from __future__ import annotations


def ivim(*args, **kwargs):  # noqa: D401  [CP2]
    raise NotImplementedError("Gnomon CP2: clean-room IVIM forward model.")
