"""The decision-calibration gap: tau_stat, tau_star, G, and the break-even shift.

Encodes GATE 1 (kappa=0 -> G~0 and v1 reproduction; misspecified -> G>0 with tau*>tau_stat;
independence of the two code paths) and GATE 3 (break-even shift brackets a VoTG sign change).
"""
from __future__ import annotations

import numpy as np

from minos.config import BASELINE_V1, MinosConfig, gaussian_latent_config
from minos.calibration import (
    break_even_shift,
    gap,
    tau_star,
    tau_stat,
)
from minos.gate import votg
from minos.generative import make_population
from minos.seeding import make_rng
from minos.voi import posterior_eu_curve

N = 1_000_000


def _wellspec(kappa=0.0, lam=3.0, rho=0.5, n=N):
    cfg = gaussian_latent_config(rho=rho, kappa=kappa, lam=lam,
                                 base=MinosConfig(n_voxels=n))
    return cfg, make_population(cfg, make_rng(cfg.seed))


# ---- GATE 1: well-specified limit -----------------------------------------------
def test_tau_stat_is_one_at_kappa_zero_every_level():
    cfg, base = _wellspec(kappa=0.0)
    for L in (0.5, 0.8, 0.9, 0.95):
        assert abs(tau_stat(base, cfg, level=L) - 1.0) < 0.01


def test_tau_star_is_one_at_kappa_zero():
    cfg, base = _wellspec(kappa=0.0)
    t, _ = tau_star(base, cfg)
    assert abs(t - 1.0) < 0.03


def test_gap_vanishes_when_well_specified():
    cfg, base = _wellspec(kappa=0.0)
    g = gap(base, cfg)
    assert abs(g.gap) < 0.025


def test_v1_mixture_voc_argmin_at_calibration():
    # Reproduce v1: the symmetric mixture posterior EU peaks at tau=1 (VoC argmin = tau*).
    cfg = MinosConfig(n_voxels=N)  # BASELINE_V1 shape
    base = make_population(cfg, make_rng(cfg.seed))
    taus = np.round(np.arange(0.7, 1.31, 0.05), 3)
    eus = posterior_eu_curve(base, cfg, taus)
    assert abs(taus[int(np.argmax(eus))] - 1.0) < 1e-9


# ---- GATE 1: misspecified -> gap with the documented sign ------------------------
def test_gap_positive_and_underconfident_favored_under_misspecification():
    cfg, base = _wellspec(kappa=3.0, lam=3.0)
    g = gap(base, cfg)
    assert g.gap > 0.03
    assert g.tau_star > g.tau_stat        # decision wants WIDER than coverage
    assert g.tau_stat < 1.0 < g.tau_star + 1e-9  # calibration shrinks, decision holds/widens
    assert g.ratio > 1.0


def test_gap_increases_with_skew():
    gaps = []
    for k in (0.0, 1.0, 2.0, 3.0, 4.0):
        cfg, base = _wellspec(kappa=k, lam=3.0)
        gaps.append(gap(base, cfg).gap)
    # clearly resolved, strictly increasing over kappa>=1
    assert np.all(np.diff(gaps[1:]) > 0)
    assert gaps[-1] - gaps[0] > 0.08


def test_gap_increases_with_utility_asymmetry():
    gaps = []
    for lam in (1.0, 2.0, 4.0):
        cfg, base = _wellspec(kappa=3.0, lam=lam)
        gaps.append(gap(base, cfg).gap)
    assert gaps[0] > 0.0          # lambda=1 is the pure coverage offset (>0 under skew)
    assert np.all(np.diff(gaps) > 0)


# ---- GATE 0/1: the two paths are independent -------------------------------------
def test_tau_stat_ignores_the_utility():
    # tau_stat reads only coverage; doubling the cost asymmetry must not move it at all.
    cfg, base = _wellspec(kappa=3.0, lam=3.0)
    cfg2 = cfg.replace(k_under=cfg.k_under * 3.0)  # same base draws, different utility
    base2 = make_population(cfg2, make_rng(cfg2.seed))
    assert tau_stat(base, cfg) == tau_stat(base2, cfg2)


def test_tau_star_flat_eu_guard_at_symmetric_utility():
    # lambda=1: the escalate boundary is tau-independent -> EU flat -> tau* defaults to 1.
    cfg, base = _wellspec(kappa=3.0, lam=1.0)
    t, _ = tau_star(base, cfg)
    assert t == 1.0


# ---- GATE 3: break-even shift for the trust-gate ---------------------------------
def test_break_even_shift_brackets_a_votg_sign_change():
    cfg = BASELINE_V1.replace(n_voxels=600_000)
    base = make_population(cfg, make_rng(cfg.seed))
    d_be = break_even_shift(base, cfg)
    assert 0.0 < d_be < 1.5
    assert votg(base, cfg, delta=d_be - 0.1) < 0.0
    assert votg(base, cfg, delta=d_be + 0.2) > 0.0
