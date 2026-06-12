"""Synthetic latent field + measurement model, with common-random-number (CRN)
realisation so every ``(tau, delta)`` reuses one set of base variates.

The latent severity ``theta`` and the measurement map live behind a small seam so a
real IVIM parameter map + Fashion posterior can later replace the synthetic source
without touching the decision / VoI / gate core. See the marked region below.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import MinosConfig


@dataclass(frozen=True)
class BaseDraws:
    """CRN base variates, sampled once and reused across the whole sweep.

    ``theta`` is the true severity scored by the utility (independent of ``tau``/``delta``);
    ``z_eta`` and ``z_w`` are standard-normal base draws transformed per shift in
    :func:`realise`. ``z_skew`` is the second standard normal that skews the error when
    ``cfg.kappa>0`` (DESIGN_B Section 2); at ``kappa=0`` it is unused and the error is
    exactly v1's. ``report_center`` is the reported posterior mean ``mu`` in the
    *posterior-centric* model (``None`` in the v1 forward-measurement model).
    """

    theta: np.ndarray
    z_eta: np.ndarray
    z_w: np.ndarray
    z_skew: np.ndarray
    report_center: np.ndarray | None = None


# ----------------------------------------------------------------------------------
# IVIM seam — Fashion integration point (deferred).
#
# ``sample_latent`` is the ONLY place the synthetic latent source is defined. To wire
# in real data later, replace its body so it returns a per-voxel IVIM parameter map
# (e.g. pseudo-diffusion fraction f) as ``theta``, and have :func:`realise` consume a
# Fashion posterior ``(mu, sigma_rep)`` instead of the synthetic Gaussian measurement.
# Nothing in decision.py / voi.py / gate.py reads the latent source directly.
# ----------------------------------------------------------------------------------
def sample_latent(cfg: MinosConfig, rng: np.random.Generator, n: int) -> np.ndarray:
    """Sample ``n`` latent severities from the configured prior.

    ``latent_mode="mixture"`` (v1): the symmetric 3-component Gaussian mixture.
    ``latent_mode="gaussian"`` (v2 gap sweep): a single ``N(theta_mean, theta_std^2)``.
    """
    if cfg.latent_mode == "gaussian":
        return cfg.theta_mean + cfg.theta_std * rng.standard_normal(n)
    weights = np.asarray(cfg.mix_weights, dtype=float)
    means = np.asarray(cfg.mix_means, dtype=float)
    stds = np.asarray(cfg.mix_stds, dtype=float)
    comp = rng.choice(len(weights), size=n, p=weights)
    z = rng.standard_normal(n)
    return means[comp] + stds[comp] * z


def _unit_skew_error(z_eta: np.ndarray, z_skew: np.ndarray, cfg: MinosConfig) -> np.ndarray:
    """Zero-mean, unit-variance error with skew-normal shape ``cfg.kappa`` (DESIGN_B §2).

    With ``d(kappa)=kappa/sqrt(1+kappa^2)``, ``raw = d|z_skew| + sqrt(1-d^2) z_eta`` is a
    skew-normal(shape ``+kappa``) variate — a heavy *right* tail for ``kappa>0`` —
    standardised to mean 0, sd 1. In the posterior-centric model the truth is
    ``theta = mu + s*u``, so a right-tailed ``u`` is the heavy *under-treatment* tail
    (truth occasionally far more severe than the report); the symmetric report cannot
    represent it, which is what makes the decision hedge up (``tau*>1``). At ``kappa=0``
    (or ``family="gaussian"``) it returns ``z_eta`` unchanged, so the error reduces to
    v1's ``N(0, s^2)``.
    """
    if cfg.family == "gaussian" or cfg.kappa == 0.0:
        return z_eta
    d = cfg.kappa / np.sqrt(1.0 + cfg.kappa ** 2)
    raw = d * np.abs(z_skew) + np.sqrt(1.0 - d ** 2) * z_eta
    mean = d * np.sqrt(2.0 / np.pi)
    sd = np.sqrt(1.0 - 2.0 * d ** 2 / np.pi)
    return (raw - mean) / sd


def make_population(cfg: MinosConfig, rng: np.random.Generator) -> BaseDraws:
    """Draw the CRN base variates for ``cfg.n_voxels`` voxels.

    ``z_skew`` is drawn *after* ``(theta, z_eta, z_w)`` so the v1 variate stream is
    byte-for-byte unchanged: at ``kappa=0`` every v1 number reproduces exactly.

    Two generative directions (DESIGN_B §2):

    * v1 forward-measurement (``posterior_centric=False``): ``theta`` is the latent
      severity (prior draw); the report ``mu = theta + error`` is formed in
      :func:`realise`. The Bayes posterior ``theta|mu`` carries prior shrinkage.
    * posterior-centric (``posterior_centric=True``): the latent draw is the *report
      centre* ``mu``; the truth is ``theta = mu + s * u`` with ``u`` the unit skew error.
      Then ``E[theta|mu] = mu`` exactly, so at ``kappa=0`` the report ``N(mu, s^2)`` IS the
      true posterior and ``tau* = 1`` for any latent placement / utility asymmetry —
      isolating the gap as a pure *misspecification* (skew) effect.
    """
    n = cfg.n_voxels
    latent = sample_latent(cfg, rng, n)
    z_eta = rng.standard_normal(n)
    z_w = rng.standard_normal(n)
    z_skew = rng.standard_normal(n)
    if cfg.posterior_centric:
        report_center = latent
        theta = report_center + cfg.s * _unit_skew_error(z_eta, z_skew, cfg)
        return BaseDraws(theta=theta, z_eta=z_eta, z_w=z_w, z_skew=z_skew,
                         report_center=report_center)
    return BaseDraws(theta=latent, z_eta=z_eta, z_w=z_w, z_skew=z_skew)


def _as_shift_mask(shift, n: int) -> np.ndarray:
    if np.isscalar(shift) or (isinstance(shift, np.ndarray) and shift.ndim == 0):
        return np.full(n, bool(shift))
    shift = np.asarray(shift, dtype=bool)
    assert shift.shape == (n,), "shift mask must be scalar or length n_voxels"
    return shift


def standardized_error(base: BaseDraws, cfg: MinosConfig) -> np.ndarray:
    """Unit-variance, zero-mean base error with skew-normal shape ``cfg.kappa``.

    Wraps :func:`_unit_skew_error` on this population's base draws (DESIGN_B §2).
    """
    return _unit_skew_error(base.z_eta, base.z_skew, cfg)


def realise(base: BaseDraws, cfg: MinosConfig, *, delta: float, shift):
    """Map base variates to observed estimate ``mu`` and acquisition feature ``w``.

    ``shift`` is a bool scalar (applies to all voxels) or a per-voxel mask. Shifted
    voxels get inflated true noise ``s*(1+alpha*delta)``, a downward bias
    ``-beta*s*delta``, and feature mean ``delta``; unshifted voxels are nominal.

    In the v1 forward model ``mu = theta + bias + sigma_true * u`` (``u`` the standardised
    skew error). In the posterior-centric model the report centre is primary, so
    ``mu = report_center + bias`` (the truth already carries the skew, set at sampling).
    """
    n = base.theta.shape[0]
    mask = _as_shift_mask(shift, n)
    sigma_true = np.where(mask, cfg.s * (1.0 + cfg.alpha * delta), cfg.s)
    bias = np.where(mask, -cfg.beta * cfg.s * delta, 0.0)
    feat_mean = np.where(mask, delta, 0.0)
    if cfg.posterior_centric:
        mu = base.report_center + bias
    else:
        mu = base.theta + bias + sigma_true * standardized_error(base, cfg)
    w = cfg.w_train_mean + feat_mean + cfg.w_train_std * base.z_w
    return mu, w
