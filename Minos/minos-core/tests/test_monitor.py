"""The deployment-validity monitor M (GATE 2) and gated recovery (GATE 3).

Encodes the contribution's claims and — centrally — its honesty:
  * M is label-free: a pure function of the reported points; scrambling theta cannot move it.
  * Under delta_obs: M tracks regret (corr>0) and detects {R>tol} with AUC well above chance.
  * Under delta_hid: detection is AT CHANCE (AUC~0.5) — documented limitation, asserted as a band.
  * Zero-shift false alarms are controlled by m*.
  * Gated recovery beats the stale correction under delta_obs, does no harm at zero shift, and does
    NOT spuriously help under delta_hid (it cannot see it) — the leakage check.
The swappable interface is exercised for both monitor kinds.
"""
from __future__ import annotations

import numpy as np

from minos.config import MinosConfig, gaussian_latent_config
from minos.correction import (
    deploy_expected_utility,
    fit_loss_calibration,
    oracle_deploy_scale,
    stale_regret,
)
from minos.generative import make_population, realise_deploy
from minos.monitor import (
    build_reference,
    calibrate_threshold,
    detection_auc,
    gated_recovery_actions,
    monitor,
)
from minos.seeding import make_rng
from minos.voi import realised_utility

# Unit tests verify correctness, not statistical precision (the heavy, high-N runs live in
# experiments/run_c.py), so they use modest populations and batch counts to stay fast.
NCAL = 300_000
NDEP = 80_000
NULL_SEEDS = 30


def _setup(n_cal=NCAL):
    cfg = gaussian_latent_config(rho=0.5, kappa=3.0, lam=3.0, base=MinosConfig(n_voxels=n_cal))
    base_cal = make_population(cfg, make_rng(cfg.seed))
    tau_hat = fit_loss_calibration(base_cal, cfg)
    ref = build_reference(base_cal, cfg, tau_hat)
    return cfg, base_cal, tau_hat, ref


def _detection_batches(cfg, ref, tau_hat, *, channel, deltas, seeds, kind="utility_divergence"):
    """Per-batch (M_label_free, R_oracle) over seeds x deltas for one shift channel."""
    cfg_b = cfg.replace(n_voxels=NDEP)
    Ms, Rs = [], []
    for sd in seeds:
        base = make_population(cfg_b, make_rng(sd))
        for d in deltas:
            do, dh = (d, 0.0) if channel == "obs" else (0.0, d)
            mu, _ = realise_deploy(base, cfg_b, delta_obs=do, delta_hid=dh)
            Ms.append(monitor(mu, cfg_b, ref, kind=kind))
            Rs.append(stale_regret(base, cfg_b, tau_hat, delta_obs=do, delta_hid=dh))
    return np.array(Ms), np.array(Rs)


# ---- GATE 2: the no-label code path --------------------------------------------------
def test_monitor_is_label_free_under_theta_scramble():
    cfg, base_cal, tau_hat, ref = _setup(200_000)
    base_dep = make_population(cfg.replace(n_voxels=200_000), make_rng(cfg.seed + 777))
    mu, _ = realise_deploy(base_dep, cfg, delta_obs=0.1)
    m_true = monitor(mu, cfg, ref)
    # corrupt the labels arbitrarily; the monitor (a function of mu only) must be bit-identical.
    base_dep.theta[:] = make_rng(123).standard_normal(base_dep.theta.shape) * 50.0
    mu2, _ = realise_deploy(base_dep, cfg, delta_obs=0.1)
    assert monitor(mu2, cfg, ref) == m_true


# ---- GATE 2: detectable under obs, blind under hid -----------------------------------
def test_monitor_grows_under_obs_flat_under_hid():
    cfg, base_cal, tau_hat, ref = _setup()
    base_dep = make_population(cfg.replace(n_voxels=NDEP), make_rng(cfg.seed + 777))
    M_obs = [monitor(realise_deploy(base_dep, cfg.replace(n_voxels=NDEP), delta_obs=d)[0],
                     cfg.replace(n_voxels=NDEP), ref) for d in (0.0, 0.05, 0.1, 0.2)]
    M_hid = [monitor(realise_deploy(base_dep, cfg.replace(n_voxels=NDEP), delta_hid=d)[0],
                     cfg.replace(n_voxels=NDEP), ref) for d in (0.0, 0.05, 0.1, 0.2)]
    assert np.all(np.diff(M_obs) > 0)                       # rises with observable shift
    assert max(M_hid) - min(M_hid) < 1e-9                   # exactly flat under hidden shift
    assert M_hid[0] == M_hid[-1]


def test_corr_M_R_positive_under_obs():
    cfg, base_cal, tau_hat, ref = _setup()
    deltas = [0.0, 0.06, 0.12, 0.18, 0.24]
    M, R = _detection_batches(cfg, ref, tau_hat, channel="obs", deltas=deltas, seeds=[7001])
    # CRN sweep is deterministic and monotone -> strong positive Pearson correlation.
    assert np.corrcoef(M, R)[0, 1] > 0.85


def test_detection_auc_high_under_obs_chance_under_hid():
    cfg, base_cal, tau_hat, ref = _setup()
    seeds = [7001, 7002, 7003, 7004]
    deltas = [0.0, 0.06, 0.12, 0.2]
    M_o, R_o = _detection_batches(cfg, ref, tau_hat, channel="obs", deltas=deltas, seeds=seeds)
    M_h, R_h = _detection_batches(cfg, ref, tau_hat, channel="hid", deltas=deltas, seeds=seeds)
    # a tol that splits the regret range into detectable positives / negatives
    tol = 0.02
    auc_obs = detection_auc(M_o, R_o, tol)
    auc_hid = detection_auc(M_h, R_h, tol)
    assert auc_obs > 0.8                        # observable staleness is caught
    assert 0.35 <= auc_hid <= 0.65             # hidden staleness is at chance (documented limit)


# ---- GATE 3: false alarms, gated recovery, leakage check -----------------------------
def test_zero_shift_false_alarm_controlled():
    cfg, base_cal, tau_hat, ref = _setup()
    alpha = 0.05
    m_star = calibrate_threshold(cfg, ref, alpha=alpha, n_seeds=NULL_SEEDS, n_batch=NDEP)
    # fresh null seeds disjoint from the calibration of m*: empirical FPR ~ alpha, certainly small.
    cfg_b = cfg.replace(n_voxels=NDEP)
    fires = 0
    trials = 25
    for k in range(trials):
        base = make_population(cfg_b, make_rng(50_000 + k))
        mu, _ = realise_deploy(base, cfg_b)
        fires += monitor(mu, cfg_b, ref) > m_star
    assert fires / trials <= 0.2                # comfortably bounded near the 5% target


def _policy_regrets(cfg, ref, tau_hat, m_star, *, delta_obs, delta_hid, seed, kind="utility_divergence"):
    cfg_b = cfg.replace(n_voxels=NDEP)
    base = make_population(cfg_b, make_rng(seed))
    mu, theta = realise_deploy(base, cfg_b, delta_obs=delta_obs, delta_hid=delta_hid)
    tau_or, eu_or = oracle_deploy_scale(base, cfg_b, delta_obs=delta_obs, delta_hid=delta_hid)
    eu_stale = deploy_expected_utility(base, cfg_b, tau_hat, delta_obs=delta_obs, delta_hid=delta_hid)
    a_gate = gated_recovery_actions(mu, cfg_b, ref, m_star, kind=kind)
    eu_gated = float(np.mean(realised_utility(a_gate, theta, cfg_b)))
    return eu_or - eu_stale, eu_or - eu_gated     # (R_stale, R_gated)


def test_gated_recovery_beats_stale_under_obs():
    cfg, base_cal, tau_hat, ref = _setup()
    m_star = calibrate_threshold(cfg, ref, alpha=0.05, n_seeds=NULL_SEEDS, n_batch=NDEP)
    R_stale, R_gated = _policy_regrets(cfg, ref, tau_hat, m_star,
                                       delta_obs=0.2, delta_hid=0.0, seed=8001)
    assert R_gated < R_stale - 0.02              # clear recovery under observable shift


def test_gated_recovery_no_harm_zero_shift():
    cfg, base_cal, tau_hat, ref = _setup()
    m_star = calibrate_threshold(cfg, ref, alpha=0.05, n_seeds=NULL_SEEDS, n_batch=NDEP)
    R_stale, R_gated = _policy_regrets(cfg, ref, tau_hat, m_star,
                                       delta_obs=0.0, delta_hid=0.0, seed=8001)
    assert abs(R_gated - R_stale) < 1e-6         # gate doesn't fire -> identical to stale


def test_gated_recovery_no_spurious_help_under_hid():
    # The leakage check: the monitor cannot see the hidden shift, so the gate must NOT recover here.
    cfg, base_cal, tau_hat, ref = _setup()
    m_star = calibrate_threshold(cfg, ref, alpha=0.05, n_seeds=NULL_SEEDS, n_batch=NDEP)
    R_stale, R_gated = _policy_regrets(cfg, ref, tau_hat, m_star,
                                       delta_obs=0.0, delta_hid=0.2, seed=8001)
    assert abs(R_gated - R_stale) < 1e-6         # gate stays silent; gated == stale


# ---- the swappable interface is real -------------------------------------------------
def test_action_divergence_monitor_swaps_in():
    cfg, base_cal, tau_hat, ref = _setup()
    cfg_b = cfg.replace(n_voxels=NDEP)
    base = make_population(cfg_b, make_rng(cfg.seed + 777))
    m_obs = [monitor(realise_deploy(base, cfg_b, delta_obs=d)[0], cfg_b, ref, kind="action_divergence")
             for d in (0.0, 0.1, 0.2)]
    m_hid = [monitor(realise_deploy(base, cfg_b, delta_hid=d)[0], cfg_b, ref, kind="action_divergence")
             for d in (0.0, 0.1, 0.2)]
    assert np.all(np.diff(m_obs) > 0)            # detects observable shift
    assert max(m_hid) - min(m_hid) < 1e-9        # blind to hidden shift (same honesty property)
