"""Synthetic cohort substrate via read-only Lattice.  [CP2 — spec below]

Gnomon's synthetic data is the **Lattice** DRO, imported read-only through
:mod:`gnomon._paths` (``lattice`` only; Caliper forbidden). This module is a thin
adapter that:

* calls ``lattice.make_cohort(family, n, snr, seed, prior, noise, bvalues)`` to draw
  reproducible (D, Dstar, f) ground truth + multi-b signals;
* exposes the headline 9-cell design (3 truths x SNR{10,20,40} x 200 noise) and the
  clinical-sparse / dense b-schemes pinned in :mod:`gnomon.manifest`;
* never reimplements the generators -- the forward model in :mod:`gnomon.forward` is
  Gnomon's *own* (for the estimators/CRLB); the *data* comes from Lattice so results
  are comparable across the IVIM program.

CP2 gate: a clean (noise-free) Lattice cohort round-trips through Gnomon's NLLS to
recover truth within tolerance (continuity / sanity check).
"""
from __future__ import annotations

from . import _paths


def load_lattice():  # [CP2]
    """Resolve and return the read-only ``lattice`` module (Caliper stays out)."""
    _paths.ensure_lattice(strict=True)
    import lattice  # noqa: WPS433  (resolved read-only via _paths)
    return lattice


def headline_cohorts(*args, **kwargs):  # [CP2]
    raise NotImplementedError("Gnomon CP2: Lattice-backed 9-cell headline design.")
