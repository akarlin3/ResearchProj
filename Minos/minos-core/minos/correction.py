"""The loss-calibration correction (a **cited baseline**) and its stale-correction regret.

v3 does **not** introduce a calibration method. Tuning the reported error bar so the *decision*
is optimal is an existing field — **loss-calibrated Bayesian inference** (Lacoste-Julien, Huszar
& Ghahramani, *Approximate inference for the loss-calibrated Bayesian*, AISTATS 2011; Cobb et al.
2018; Kusmierczyk et al. 2019; **post-hoc** form, Vadera et al., *Post-hoc loss-calibration*, UAI
2021) and **decision calibration** (Zhao, Kim, Sahoo, Ma & Ermon, NeurIPS 2021). Our scalar
``tau_hat_cal`` is the one-dimensional instance of that line (== v2's ``tau_star`` on a labeled
calibration set). See ``POSITIONING.md``.

What v3 adds (``monitor.py``) is the thing that line does not address: those methods assume the
calibration set represents deployment. This module supplies the **baseline** they prescribe and the
**stale-correction regret** ``R`` that arises when deployment shifts away from calibration — the
quantity the validity monitor is built to predict label-free. Math: DESIGN_C Sections 3.1-3.2.

Code paths: ``fit_loss_calibration`` reuses ``calibration.tau_star`` verbatim (the cited baseline);
the deploy evaluator reuses the ``utility``/``decision``/``voi`` primitives on the v3
``realise_deploy`` realisation. The oracle scale and ``R`` read deployment labels ``theta`` and so
are **simulation-only validation quantities** — never consumed by the monitor or the policy.
"""
from __future__ import annotations

import numpy as np

from .config import MinosConfig
from .generative import BaseDraws, realise_deploy
from .utility import Action, expected_utility_under_q, utility

# Chunk size for the vectorised tau-grid EU evaluation: caps the (3, n_tau, chunk) work array.
_EU_CHUNK = 300_000

# Oracle-scale search window. Wider at the top than v2's ``tau_star`` bounds (0.5, 2.5): a stale
# correction under a strong report->truth bias is optimally repaired by a much wider error bar, so
# the deployment optimum can sit well above 2.5; we bracket it generously and refine. The grid is
# coarse (step 0.1) because the parabola-fit vertex recovers sub-grid precision — the dense scan only
# has to localise the peak, and each grid point costs a scipy ``norm`` evaluation over all voxels.
TAU_DEPLOY_BOUNDS = (0.5, 4.5)
TAU_DEPLOY_STEP = 0.1
TAU_DEPLOY_FIT_HALFWIDTH = 0.4
FLAT_EU_TOL = 1e-4


# --------------------------------------------------------------------------------------
# the cited loss-calibration baseline
# --------------------------------------------------------------------------------------
def fit_loss_calibration(base_cal: BaseDraws, cfg: MinosConfig) -> float:
    """``tau_hat_cal = argmax_tau Ehat[U]`` on the labeled calibration set ``D_cal``.

    The scalar post-hoc loss-calibration instance (Lacoste-Julien 2011 / Vadera 2021), realised as
    v2's decision-calibrated scale on ``D_cal`` at zero shift. Reads the labels ``theta`` of
    ``D_cal`` — that is what makes ``D_cal`` a *labeled* calibration set. Reuses ``tau_star``.
    """
    from .calibration import tau_star

    tau, _ = tau_star(base_cal, cfg)
    return float(tau)


# --------------------------------------------------------------------------------------
# deployment expected utility, oracle scale, and stale-correction regret
# --------------------------------------------------------------------------------------
def _eu_curve(mu: np.ndarray, theta: np.ndarray, cfg: MinosConfig, taus) -> np.ndarray:
    """``Ehat[U(a*(N(mu,(tau s)^2)), theta)]`` for a whole grid of ``tau`` at once.

    The action rule and its realised utility are evaluated for every ``tau`` in one vectorised pass
    (the reported sd ``tau*s`` broadcasts against the voxels), collapsing what would be a Python loop
    of hundreds of small ``scipy.norm`` calls into a handful of large ones — the dominant cost. The
    chunking over voxels caps the ``(3, n_tau, chunk)`` work array. Returns ``EU`` per ``tau``.

    The decision rule (argmax of the reported-posterior expected utility) is identical to
    ``decision.bayes_action``; this is purely a batched evaluation of the same quantity.
    """
    taus = np.asarray(taus, dtype=float)
    sig = taus[:, None] * cfg.s                      # (T, 1)
    n = mu.shape[0]
    totals = np.zeros(taus.size)
    for s0 in range(0, n, _EU_CHUNK):
        m = mu[s0:s0 + _EU_CHUNK][None, :]           # (1, c)
        th = theta[s0:s0 + _EU_CHUNK]                # (c,)
        eus = np.stack([expected_utility_under_q(a, m, sig, cfg) for a in Action], axis=0)  # (3,T,c)
        acts = np.argmax(eus, axis=0)                # (T, c) chosen action per (tau, voxel)
        ucols = np.stack([utility(a, th, cfg) for a in Action], axis=0)  # (3, c) realised utilities
        realised = ucols[acts, np.arange(th.shape[0])[None, :]]          # (T, c)
        totals += realised.sum(axis=1)
    return totals / n


def _eu_at(mu: np.ndarray, theta: np.ndarray, cfg: MinosConfig, tau: float) -> float:
    """``Ehat[U(a*(N(mu,(tau s)^2)), theta)]`` at a single ``tau`` (one-row :func:`_eu_curve`)."""
    return float(_eu_curve(mu, theta, cfg, [tau])[0])


def deploy_expected_utility(base: BaseDraws, cfg: MinosConfig, tau: float, *,
                            delta_obs: float = 0.0, delta_hid: float = 0.0) -> float:
    """``Ehat[U(a*(N(mu,(tau s)^2)), theta)]`` under the deployment shift ``(delta_obs, delta_hid)``.

    Realises the deployment ``(mu, theta)`` and evaluates the reported-posterior decision rule's
    realised utility at ``tau`` (via :func:`_eu_curve`, the batched form of ``decision.bayes_action``).
    """
    mu, theta = realise_deploy(base, cfg, delta_obs=delta_obs, delta_hid=delta_hid)
    return _eu_at(mu, theta, cfg, tau)


def _oracle_from_realisation(mu: np.ndarray, theta: np.ndarray, cfg: MinosConfig, *,
                             bounds=TAU_DEPLOY_BOUNDS, step: float = TAU_DEPLOY_STEP,
                             halfwidth: float = TAU_DEPLOY_FIT_HALFWIDTH):
    """``argmax_tau Ehat[U(tau)]`` for a fixed realisation ``(mu, theta)`` — grid + parabola refine.

    Realises nothing itself (the caller supplies ``mu, theta`` once); the whole grid is scored in one
    vectorised :func:`_eu_curve` pass, then a local parabola fit refines the peak to sub-grid
    precision. Returns ``(tau, achievable_eu)`` where the EU is the actual value at the chosen ``tau``.
    """
    lo, hi = bounds
    grid = np.round(np.arange(lo, hi + step / 2, step), 4)
    eus = _eu_curve(mu, theta, cfg, grid)
    if float(np.ptp(eus)) < FLAT_EU_TOL:
        return 1.0, _eu_at(mu, theta, cfg, 1.0)
    i = int(np.argmax(eus))
    tau_best, eu_best = float(grid[i]), float(eus[i])
    win = np.abs(grid - grid[i]) <= halfwidth + 1e-9
    x, y = grid[win], eus[win]
    if x.size >= 3:
        a, b, c = np.polyfit(x, y, 2)
        if a < 0:  # concave -> interior maximum at the vertex
            vtx = -b / (2.0 * a)
            if x.min() <= vtx <= x.max():
                # Score the vertex by its ACTUAL EU (not the parabola's prediction) and keep it only
                # if it genuinely beats the grid best -> the returned EU is always achievable.
                eu_vtx = _eu_at(mu, theta, cfg, float(vtx))
                if eu_vtx > eu_best:
                    tau_best, eu_best = float(vtx), float(eu_vtx)
    return tau_best, eu_best


def oracle_deploy_scale(base: BaseDraws, cfg: MinosConfig, *, delta_obs: float = 0.0,
                        delta_hid: float = 0.0, bounds=TAU_DEPLOY_BOUNDS,
                        step: float = TAU_DEPLOY_STEP,
                        halfwidth: float = TAU_DEPLOY_FIT_HALFWIDTH):
    """``argmax_tau deploy_EU(tau; delta)`` — the deployment-optimal scale (needs deployment labels).

    Dense grid + local parabola-fit refine, mirroring ``calibration.tau_star`` but over the
    deployment objective. Returns ``(tau, eu)``. Simulation-only: it reads deployment ``theta`` and
    exists purely to validate the monitor (it is never consumed by ``monitor.py`` or the policy).
    """
    mu, theta = realise_deploy(base, cfg, delta_obs=delta_obs, delta_hid=delta_hid)
    return _oracle_from_realisation(mu, theta, cfg, bounds=bounds, step=step, halfwidth=halfwidth)


def stale_regret(base: BaseDraws, cfg: MinosConfig, tau_hat_cal: float, *,
                 delta_obs: float = 0.0, delta_hid: float = 0.0) -> float:
    """``R(delta) = max_tau deploy_EU(tau) - deploy_EU(tau_hat_cal) >= 0`` — utility the stale
    correction strands versus the deployment-optimal scale (DESIGN_C §3.2).

    The stale scale ``tau_hat_cal`` is itself a feasible deployment scale, so the deployment optimum
    is never worse: ``R >= 0`` by construction. The ``max(0, .)`` clamps the residual grid/parabola
    optimiser slack at zero shift (where the optimum *is* ``tau_hat_cal`` and any coarse search can
    only under-shoot it) — it cannot mask real regret, which is strictly positive under shift.

    Realises the deployment ``(mu, theta)`` **once** and reuses it for both the oracle sweep and the
    stale evaluation.
    """
    mu, theta = realise_deploy(base, cfg, delta_obs=delta_obs, delta_hid=delta_hid)
    _, eu_or = _oracle_from_realisation(mu, theta, cfg)
    eu_stale = _eu_at(mu, theta, cfg, tau_hat_cal)
    return float(max(0.0, eu_or - eu_stale))
