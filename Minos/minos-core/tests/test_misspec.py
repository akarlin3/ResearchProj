"""Misspecification generator: skew-normal error, posterior-centric truth, v1 limit.

The skew machinery (DESIGN_B Section 2) must (a) reduce to v1's Gaussian error exactly at
kappa=0, (b) be a standardised (mean 0, sd 1) right-skewed error for kappa>0, and (c) in the
posterior-centric model leave the report unbiased as a posterior mean.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import skew

from minos.config import MinosConfig, gaussian_latent_config
from minos.generative import _unit_skew_error, make_population, realise
from minos.seeding import make_rng

CFG = MinosConfig(n_voxels=400_000)


def _base(cfg=CFG):
    return make_population(cfg, make_rng(cfg.seed))


def test_unit_error_is_z_eta_at_kappa_zero():
    base = _base()
    u = _unit_skew_error(base.z_eta, base.z_skew, CFG)  # kappa=0 default
    assert np.array_equal(u, base.z_eta)


def test_unit_error_gaussian_family_ignores_kappa():
    base = _base()
    cfg = CFG.replace(family="gaussian", kappa=3.0)
    u = _unit_skew_error(base.z_eta, base.z_skew, cfg)
    assert np.array_equal(u, base.z_eta)


def test_unit_error_is_standardised_and_right_skewed():
    base = make_population(MinosConfig(n_voxels=2_000_000), make_rng(0))
    for k in (0.5, 1.0, 2.0, 4.0):
        cfg = CFG.replace(kappa=k)
        u = _unit_skew_error(base.z_eta, base.z_skew, cfg)
        assert abs(u.mean()) < 5e-3
        assert abs(u.std() - 1.0) < 5e-3
        assert skew(u) > 0.0  # heavy RIGHT (under-treatment) tail for kappa>0
    # skew grows with kappa
    sk = [skew(_unit_skew_error(base.z_eta, base.z_skew, CFG.replace(kappa=k)))
          for k in (0.5, 1.0, 2.0, 4.0)]
    assert np.all(np.diff(sk) > 0)


def test_forward_model_kappa0_residual_is_gaussian():
    # v1 forward model (posterior_centric=False): mu - theta ~ N(0, s^2) exactly at kappa=0.
    base = _base()
    mu, _ = realise(base, CFG, delta=0.0, shift=False)
    resid = mu - base.theta
    assert abs(resid.mean()) < 5e-3
    assert abs(resid.std() - CFG.s) < 5e-3
    assert abs(skew(resid)) < 1e-2


def test_posterior_centric_report_is_unbiased_posterior_mean():
    # theta = mu + s*u with E[u]=0  =>  E[theta - mu] = 0 (report centred on the truth's mean).
    cfg = gaussian_latent_config(rho=0.5, kappa=3.0, lam=3.0,
                                 base=MinosConfig(n_voxels=2_000_000))
    base = make_population(cfg, make_rng(cfg.seed))
    mu, _ = realise(base, cfg, delta=0.0, shift=False)
    resid = base.theta - mu  # truth minus report
    assert np.allclose(mu, base.report_center)          # in-distribution report == centre
    assert abs(resid.mean()) < 5e-3                       # unbiased
    assert abs(resid.std() - cfg.s) < 5e-3               # correct sd (moment-matched)
    assert skew(resid) > 0.3                              # right-skewed under-treatment tail


def test_posterior_centric_kappa0_is_well_specified():
    # kappa=0 -> theta - mu ~ N(0, s^2): the report IS the true posterior.
    cfg = gaussian_latent_config(rho=0.5, kappa=0.0, lam=3.0,
                                 base=MinosConfig(n_voxels=1_000_000))
    base = make_population(cfg, make_rng(cfg.seed))
    mu, _ = realise(base, cfg, delta=0.0, shift=False)
    resid = base.theta - mu
    assert abs(resid.std() - cfg.s) < 5e-3
    assert abs(skew(resid)) < 1e-2


def test_gaussian_latent_centres_at_target_distance_from_threshold():
    cfg = gaussian_latent_config(rho=1.0, kappa=2.0, lam=3.0, base=CFG)
    assert np.isclose(cfg.theta_mean, cfg.t2 - 1.0 * cfg.s)
    assert cfg.t1 < cfg.t2 - 1.0           # spare threshold pushed far below t2
    base = make_population(cfg, make_rng(cfg.seed))
    # report centres straddle the active threshold t2 (both sides populated).
    assert (base.report_center < cfg.t2).mean() > 0.1
    assert (base.report_center > cfg.t2).mean() > 0.1
