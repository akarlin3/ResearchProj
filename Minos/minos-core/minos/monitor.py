"""The deployment-validity monitor ``M`` — the v3 contribution.

A **label-free** statistic of the unlabeled deployment reported points ``{mu_i}`` plus the **known**
utility, that predicts the stale-correction regret ``R`` (``correction.py``) *without* deployment
labels — so it flags when a learned loss-calibration correction has gone stale under shift.

Two hard boundaries (DESIGN_C §0, POSITIONING.md):

1. **Not a calibration method.** ``M`` does not propose a better scale; it detects staleness of the
   *cited* loss-calibration baseline (``correction.fit_loss_calibration``).
2. **Not a generic OOD detector.** ``M`` is **regret-targeted**: the deployment/calibration
   divergence is weighted by the *utility stakes near the decision threshold*, not by input density.
   A density-ratio OOD score would fire on far-from-threshold changes that cannot move the optimal
   scale; ``M`` is built to ignore exactly those.

THE NO-LABEL DISCIPLINE (asserted in tests). Every label-free function here — ``monitor``,
``calibrate_threshold``, ``gated_recovery_actions`` — takes the reported points ``mu_dep`` (and
``cfg``/``ref``), **never** ``theta``. ``theta`` enters this build only to *score* utility
(``correction.deploy_expected_utility``, the regret), never to compute ``M`` or to choose an action.
The honesty constraint follows structurally: ``M`` is a function of ``{mu}`` only, and the hidden
shift (``delta_hid``) does not move ``{mu}`` — so ``M`` is provably blind to it (detection at chance),
the documented limitation that motivates labeled repeatability spot-checks.

Math: DESIGN_C Sections 3.3-3.4.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import MinosConfig
from .decision import bayes_action
from .gate import _auc
from .generative import make_population, realise_deploy
from .seeding import make_rng
from .utility import Action

# Fixed reporting grid for the reported-coordinate histogram (in units of reported sd around t2).
Z_RANGE = (-5.0, 5.0)
Z_BINS = 50
# Default conservative-override half-width: escalate reported points within this many reported sd
# of the active threshold when the monitor fires (DESIGN_C §3.4).
Z_GUARD = 1.0

MONITOR_KINDS = ("utility_divergence", "action_divergence")


# --------------------------------------------------------------------------------------
# frozen calibration-side reference (built once, at calibration time)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class MonitorRef:
    """Everything the monitor compares deployment against, frozen at calibration time.

    ``tau_hat_cal`` is the fitted loss-calibration scale; ``z_edges``/``p_cal`` are the reference
    histogram of the decision-relevant reported coordinate ``z=(mu-t2)/s`` on ``D_cal``; ``a_cal`` is
    the reference action distribution the correction induces on ``D_cal`` (for the action-divergence
    monitor). No deployment information and no labels live here.
    """

    tau_hat_cal: float
    z_edges: np.ndarray
    p_cal: np.ndarray
    a_cal: np.ndarray


def _z_coord(mu: np.ndarray, cfg: MinosConfig) -> np.ndarray:
    """Signed reported distance to the active threshold ``t2``, in units of reported sd ``s``."""
    return (np.asarray(mu, dtype=float) - cfg.t2) / cfg.s


def _z_density(mu: np.ndarray, cfg: MinosConfig, z_edges: np.ndarray) -> np.ndarray:
    h, _ = np.histogram(_z_coord(mu, cfg), bins=z_edges, density=True)
    return h


def _stakes_weights(z_edges: np.ndarray, cfg: MinosConfig) -> np.ndarray:
    """Utility-stakes kernel at bin centres: ``omega(z) = k_under * phi(z)`` (DESIGN_C §3.3).

    ``phi`` (unit-normal pdf) localises to within ~1 reported sd of the threshold — the only region
    where the scale ``tau`` can flip the treat/escalate action — and ``k_under`` scales by the cost
    at stake. This weighting (utility stakes, not input density) is what makes ``M`` regret-targeted.
    """
    centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    return cfg.k_under * np.exp(-0.5 * centres ** 2) / np.sqrt(2.0 * np.pi)


def _action_distribution(mu: np.ndarray, cfg: MinosConfig, tau: float) -> np.ndarray:
    """Label-free action distribution the correction induces: fractions ``[spare, treat, escalate]``
    under ``a*(N(mu,(tau s)^2))``. Reads reported points only — no labels."""
    actions = np.asarray(bayes_action(np.asarray(mu, dtype=float), tau * cfg.s, cfg))
    return np.array([float(np.mean(actions == int(a))) for a in Action])


def build_reference(base_cal, cfg: MinosConfig, tau_hat_cal: float) -> MonitorRef:
    """Freeze the calibration-side reference (reported-coordinate histogram + induced action
    distribution). Label-free w.r.t. the *monitor*: it reads ``base_cal.report_center`` (the reported
    points), not ``theta``."""
    z_edges = np.linspace(Z_RANGE[0], Z_RANGE[1], Z_BINS + 1)
    mu_cal, _ = realise_deploy(base_cal, cfg, delta_obs=0.0, delta_hid=0.0)
    p_cal = _z_density(mu_cal, cfg, z_edges)
    a_cal = _action_distribution(mu_cal, cfg, tau_hat_cal)
    return MonitorRef(tau_hat_cal=float(tau_hat_cal), z_edges=z_edges,
                      p_cal=p_cal, a_cal=a_cal)


# --------------------------------------------------------------------------------------
# the monitor (one-function swappable interface; DEFAULT = utility_divergence)
# --------------------------------------------------------------------------------------
def _monitor_utility_divergence(mu_dep: np.ndarray, cfg: MinosConfig, ref: MonitorRef) -> float:
    """DEFAULT ``M``: utility-weighted L1 divergence of the reported-coordinate distribution.

    ``M = sum_b omega(z_b) * |p_dep(z_b) - p_cal(z_b)| * dz`` (DESIGN_C §3.3).
    """
    p_dep = _z_density(mu_dep, cfg, ref.z_edges)
    omega = _stakes_weights(ref.z_edges, cfg)
    dz = ref.z_edges[1] - ref.z_edges[0]
    return float(np.sum(omega * np.abs(p_dep - ref.p_cal)) * dz)


def _monitor_action_divergence(mu_dep: np.ndarray, cfg: MinosConfig, ref: MonitorRef) -> float:
    """Alternative ``M``: utility-weighted divergence of the *induced action* distribution.

    ``M = k_under * 0.5 * sum_a |p_dep(a) - p_cal(a)|`` — the (cost-scaled) total-variation between
    the action profile the stale correction induces on deployment and on calibration. A coarser,
    decision-grounded statistic than the DEFAULT (3 action bins vs the full reported-coordinate
    histogram), behind the same interface. Same honesty property: the actions are a function of
    ``mu`` only, so a hidden shift (``mu`` fixed) leaves the profile — and ``M`` — unmoved.

    Naturally regret-targeted: the action profile only moves when reported mass crosses a decision
    boundary, and the ``k_under`` scale ties the alarm to the under-treatment cost the scale hedges.
    """
    a_dep = _action_distribution(mu_dep, cfg, ref.tau_hat_cal)
    return float(cfg.k_under * 0.5 * np.sum(np.abs(a_dep - ref.a_cal)))


_MONITORS = {
    "utility_divergence": _monitor_utility_divergence,
    "action_divergence": _monitor_action_divergence,
}


def monitor(mu_dep: np.ndarray, cfg: MinosConfig, ref: MonitorRef,
            kind: str = "utility_divergence") -> float:
    """The deployment-validity monitor ``M(D_dep)`` — a pure function of the reported points.

    Takes ``mu_dep`` (reported points), ``cfg`` (known utility + reported ``s``), and the frozen
    calibration ``ref``. **Never reads ``theta``** — label-free by construction. ``kind`` selects the
    statistic behind this one interface (DEFAULT ``"utility_divergence"``; ``"selfconsistency"`` is a
    real alternative). Returns a non-negative staleness score.
    """
    if kind not in _MONITORS:
        raise ValueError(f"unknown monitor kind {kind!r}; choose from {MONITOR_KINDS}")
    return _MONITORS[kind](np.asarray(mu_dep, dtype=float), cfg, ref)


# --------------------------------------------------------------------------------------
# threshold m* (zero-shift false-alarm control) and the gated-recovery policy
# --------------------------------------------------------------------------------------
def calibrate_threshold(cfg: MinosConfig, ref: MonitorRef, *, alpha: float = 0.05,
                        n_seeds: int = 60, n_batch: int = 150_000,
                        kind: str = "utility_divergence", seed0: int = 9_000) -> float:
    """``m* = (1-alpha)`` quantile of ``M`` under the zero-shift null (DESIGN_C §3.4).

    Draws ``n_seeds`` fresh deployment-size populations from the calibration distribution (no shift),
    scores each with the monitor, and returns the upper-``alpha`` quantile — controlling the
    zero-shift false-alarm rate at ``alpha``. Label-free (only reported points enter).
    """
    cfg_b = cfg.replace(n_voxels=n_batch)
    scores = np.empty(n_seeds)
    for k in range(n_seeds):
        base = make_population(cfg_b, make_rng(seed0 + k))
        mu, _ = realise_deploy(base, cfg_b, delta_obs=0.0, delta_hid=0.0)
        scores[k] = monitor(mu, cfg_b, ref, kind=kind)
    return float(np.quantile(scores, 1.0 - alpha))


def gated_recovery_actions(mu_dep: np.ndarray, cfg: MinosConfig, ref: MonitorRef,
                           m_star: float, *, kind: str = "utility_divergence",
                           z_guard: float = Z_GUARD) -> np.ndarray:
    """Label-free gated-recovery policy (DESIGN_C §3.4).

    Apply the stale correction's actions ``a*(N(mu,(tau_hat_cal s)^2))`` everywhere; **if the monitor
    fires** (``M(mu_dep) > m_star``) override the decision-fragile reported points (``|z| < z_guard``)
    to the conservative ESCALATE arm. The override direction is fixed by the *known* cost asymmetry
    (``k_under > k_over``), not inferred from labels. Reads ``mu_dep`` and the scalar ``M`` only —
    never ``theta``.
    """
    mu_dep = np.asarray(mu_dep, dtype=float)
    actions = np.asarray(bayes_action(mu_dep, ref.tau_hat_cal * cfg.s, cfg)).copy()
    if monitor(mu_dep, cfg, ref, kind=kind) > m_star:
        fragile = np.abs(_z_coord(mu_dep, cfg)) < z_guard
        actions[fragile] = int(Action.ESCALATE)
    return actions


# --------------------------------------------------------------------------------------
# detection AUC helper (reuses gate._auc) — validation only (consumes the oracle regret)
# --------------------------------------------------------------------------------------
def detection_auc(monitor_scores: np.ndarray, regrets: np.ndarray, tol: float) -> float:
    """AUC of the (label-free) monitor scores against the oracle label ``{R > tol}``.

    Validation-only: the labels come from the simulation oracle ``R``; the *scores* are the
    label-free monitor. Reuses ``gate._auc`` (Mann-Whitney U).
    """
    labels = np.asarray(regrets, dtype=float) > tol
    return _auc(np.asarray(monitor_scores, dtype=float), labels)
