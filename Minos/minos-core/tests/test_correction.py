"""The loss-calibration baseline and stale-correction regret (GATE 1).

Encodes: the baseline IS v2's decision-calibrated scale (reuse identity); consistency
``tau_hat_cal -> tau*_oracle`` as ``N_cal`` grows at zero shift; ``R(0,0) ~ 0``; and ``R`` increases
with each shift knob separately. The v2 gap still reproduces on this cell.
"""
from __future__ import annotations

import numpy as np

from minos.calibration import gap, tau_star
from minos.config import MinosConfig, gaussian_latent_config
from minos.correction import (
    fit_loss_calibration,
    oracle_deploy_scale,
    stale_regret,
)
from minos.generative import make_population
from minos.seeding import make_rng


def _cell(n, kappa=3.0, lam=3.0, rho=0.5):
    cfg = gaussian_latent_config(rho=rho, kappa=kappa, lam=lam, base=MinosConfig(n_voxels=n))
    return cfg, make_population(cfg, make_rng(cfg.seed))


def test_baseline_is_v2_decision_calibrated_scale():
    # The cited loss-calibration baseline == v2's tau_star on D_cal (pure reuse).
    cfg, base = _cell(200_000)
    assert fit_loss_calibration(base, cfg) == tau_star(base, cfg)[0]


def test_consistency_baseline_to_oracle_at_zero_shift():
    # Estimator consistency: as N_cal grows the fitted scale's variance shrinks (textbook
    # consistency, robust to the shallow optimum), and at large N it agrees with the
    # (zero-shift) deployment-optimal scale.
    seeds = (11, 12, 13, 14, 15)

    def fits(n):
        cfg = gaussian_latent_config(rho=0.5, kappa=3.0, lam=3.0, base=MinosConfig(n_voxels=n))
        return cfg, [fit_loss_calibration(make_population(cfg, make_rng(sd)), cfg) for sd in seeds]

    _, small = fits(80_000)
    cfg_l, large = fits(500_000)
    assert np.std(large) < np.std(small)          # variance shrinks with N -> consistent estimator
    # the fitted scale agrees with the deployment-optimal scale at large N (an independent draw)
    tau_or, _ = oracle_deploy_scale(make_population(cfg_l, make_rng(999)), cfg_l)
    assert abs(float(np.mean(large)) - tau_or) < 0.02


def test_regret_is_zero_at_zero_shift():
    cfg, base = _cell(250_000)
    tau_hat = fit_loss_calibration(base, cfg)
    base_dep = make_population(cfg, make_rng(cfg.seed + 777))
    R0 = stale_regret(base_dep, cfg, tau_hat, delta_obs=0.0, delta_hid=0.0)
    assert R0 >= 0.0             # regret is non-negative by construction
    assert R0 < 3e-3            # ~0 at the regret floor


def test_regret_increases_with_delta_obs():
    cfg, base = _cell(250_000)
    tau_hat = fit_loss_calibration(base, cfg)
    base_dep = make_population(cfg, make_rng(cfg.seed + 777))
    Rs = [stale_regret(base_dep, cfg, tau_hat, delta_obs=d, delta_hid=0.0)
          for d in (0.0, 0.05, 0.1, 0.2)]
    assert all(r >= 0.0 for r in Rs)
    assert np.all(np.diff(Rs) > 0)
    assert Rs[-1] - Rs[0] > 0.05


def test_regret_increases_with_delta_hid():
    cfg, base = _cell(250_000)
    tau_hat = fit_loss_calibration(base, cfg)
    base_dep = make_population(cfg, make_rng(cfg.seed + 777))
    Rs = [stale_regret(base_dep, cfg, tau_hat, delta_obs=0.0, delta_hid=d)
          for d in (0.0, 0.05, 0.1, 0.2)]
    assert all(r >= 0.0 for r in Rs)
    assert np.all(np.diff(Rs) > 0)
    assert Rs[-1] - Rs[0] > 0.05


def test_regret_non_negative_every_shift():
    # R >= 0 for every shift -> the stale correction never beats the deployment optimum.
    cfg, base = _cell(200_000)
    tau_hat = fit_loss_calibration(base, cfg)
    base_dep = make_population(cfg, make_rng(cfg.seed + 777))
    for do, dh in [(0.0, 0.0), (0.1, 0.0), (0.0, 0.1), (0.2, 0.0)]:
        assert stale_regret(base_dep, cfg, tau_hat, delta_obs=do, delta_hid=dh) >= 0.0


def test_v2_gap_still_reproduces_on_cell():
    # The motivating v2 result is intact on this misspecified cell: a positive gap, tau*>tau_stat.
    cfg, base = _cell(400_000)
    g = gap(base, cfg)
    assert g.gap > 0.03 and g.tau_star > g.tau_stat
