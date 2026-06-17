"""The deployment shift split: observable (``delta_obs``) vs hidden (``delta_hid``).

The honesty constraint (DESIGN_C Section 1) requires two orthogonal channels that perturb the
truth<->report relationship by the SAME amount, one moving the observable reported point ``mu`` and
one leaving it byte-identical. These tests pin that structure: ``delta_hid`` must not move ``mu`` at
all (so any label-free monitor is provably blind), and matched ``delta`` must induce the same
``theta - mu`` law (so the regret is the same either way).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import skew

from minos.config import MinosConfig, gaussian_latent_config
from minos.generative import make_population, realise_deploy
from minos.seeding import make_rng

N = 500_000


def _cell(n=N):
    cfg = gaussian_latent_config(rho=0.5, kappa=3.0, lam=3.0, base=MinosConfig(n_voxels=n))
    return cfg, make_population(cfg, make_rng(cfg.seed))


def test_requires_posterior_centric():
    # v1 forward model has no reported centre -> realise_deploy must refuse.
    cfg = MinosConfig(n_voxels=10_000)  # posterior_centric=False (v1 default)
    base = make_population(cfg, make_rng(cfg.seed))
    with pytest.raises(AssertionError):
        realise_deploy(base, cfg, delta_obs=0.1)


def test_delta_hid_leaves_mu_byte_identical():
    # The hidden channel must NOT move any observable summary: mu is unchanged for ANY delta_hid.
    cfg, base = _cell()
    mu0, _ = realise_deploy(base, cfg, delta_obs=0.0, delta_hid=0.0)
    for dh in (0.05, 0.2, 0.5, 2.0):
        mu_h, _ = realise_deploy(base, cfg, delta_obs=0.0, delta_hid=dh)
        assert np.array_equal(mu_h, mu0), f"delta_hid={dh} moved mu (monitor would not be blind)"


def test_delta_obs_translates_mu_down():
    # The observable channel biases mu down by exactly beta*s*delta_obs.
    cfg, base = _cell()
    mu0, _ = realise_deploy(base, cfg, delta_obs=0.0)
    d = 0.1
    mu_o, _ = realise_deploy(base, cfg, delta_obs=d)
    assert np.allclose(mu_o, mu0 - cfg.beta * cfg.s * d)
    assert mu_o.mean() < mu0.mean()


def test_delta_hid_biases_truth_up_only():
    # The hidden channel biases theta up by exactly beta*s*delta_hid; mu fixed.
    cfg, base = _cell()
    _, th0 = realise_deploy(base, cfg, delta_hid=0.0)
    d = 0.1
    mu_h, th_h = realise_deploy(base, cfg, delta_hid=d)
    assert np.allclose(th_h, th0 + cfg.beta * cfg.s * d)
    assert th_h.mean() > th0.mean()


def test_zero_shift_is_report_centre_and_skew_truth():
    # At zero shift: mu == report centre; theta = report_center + s*u with a right-skewed u.
    cfg, base = _cell()
    mu, theta = realise_deploy(base, cfg)
    assert np.array_equal(mu, base.report_center)
    resid = theta - mu  # = s*u
    assert abs(resid.mean()) < 5e-3          # unbiased (E[u]=0)
    assert abs(resid.std() - cfg.s) < 5e-3   # moment-matched sd
    assert skew(resid) > 0.3                 # heavy right (under-treatment) tail


def test_matched_delta_gives_same_theta_minus_mu_law():
    # The crux of the honesty construction: delta_obs(d) and delta_hid(d) induce the SAME theta-mu
    # distribution (hence the same deployment decision problem and the same regret), differing ONLY
    # in whether mu moved.
    cfg, base = _cell()
    d = 0.15
    mu_o, th_o = realise_deploy(base, cfg, delta_obs=d, delta_hid=0.0)
    mu_h, th_h = realise_deploy(base, cfg, delta_obs=0.0, delta_hid=d)
    gap_o = th_o - mu_o
    gap_h = th_h - mu_h
    assert np.allclose(gap_o, gap_h)                       # identical truth<->report gap
    assert not np.allclose(mu_o, mu_h)                     # but the observable channel differs
    assert np.allclose(mu_h, base.report_center)           # hidden leaves mu put
