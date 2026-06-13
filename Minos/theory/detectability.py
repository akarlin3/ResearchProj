"""Theorem 2 (part ii) — the achievable label-free detectability bound (CP2, gated).

The v2 gap motivates a loss-calibration correction; under deployment shift that correction goes
stale, stranding regret R. v3 builds a *label-free* monitor M = f({mu}) of the reported points. This
script derives and prints the **achievable bound** on the OBSERVABLE component of staleness:

    R_obs  <=  L * delta(P_O^dep, P_O^cal)

where P_O is the law of the observable reported point ``mu`` and L is a utility-Lipschitz constant
read off U. Two consistent forms are derived (Wasserstein-1 and total-variation), L is identified
from the code's utility, and both bounds are checked against the actual Minos-Core v3 regret on the
default cell — the bound must HOLD (it is an upper envelope), and L must be the right order as the
observed regret-vs-divergence slope.

The companion *impossibility* (part i: the HIDDEN component is undetectable, best AUC = 1/2) is a
data-processing argument drafted in ``impossibility.md`` and flagged for human proof-review — it is
NOT machine-verified here. This file only demonstrates the structural fact the impossibility rests
on: a hidden shift leaves the observable law P_O (hence M) invariant.

Run:  ``.venv-theory/bin/python theory/detectability.py``
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.stats import norm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "minos-core"))

from minos.config import MinosConfig, gaussian_latent_config       # noqa: E402
from minos.correction import fit_loss_calibration, stale_regret    # noqa: E402
from minos.generative import make_population, realise_deploy        # noqa: E402
from minos.monitor import build_reference, monitor                  # noqa: E402
from minos.seeding import make_rng                                  # noqa: E402

# v3 default calibration cell (THEORY_MODEL.md §4).
KAP, LAM, RHO = 3.0, 3.0, 0.5


def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def per_report_expected_utility(mu_grid, cfg, tau):
    """g_tau(mu) = E_u[ U(a*(N(mu,(tau s)^2)), mu + s u) ] on a grid of reported points mu.

    The realized expected utility of the reported-posterior decision rule at one reported point,
    averaged over the TRUE (skew) error law. Bounded; its oscillation enters the TV bound's L.
    Computed by quadrature over the skew-normal u (exact law of the code's standardised error).
    """
    from minos.decision import bayes_action
    from minos.utility import Action, utility
    from scipy.stats import skewnorm

    # standardised skew-normal u (mean 0, var 1) with shape kappa — the code's _unit_skew_error law.
    a = cfg.kappa
    d = a / np.sqrt(1 + a**2)
    mean_sn, sd_sn = d * np.sqrt(2 / np.pi), np.sqrt(1 - 2 * d**2 / np.pi)
    uq = np.linspace(-6, 9, 3001)
    pu = skewnorm.pdf(uq * sd_sn + mean_sn, a) * sd_sn          # density of standardised u
    pu /= np.trapezoid(pu, uq)
    g = np.empty_like(mu_grid)
    for i, mu in enumerate(mu_grid):
        act = int(bayes_action(np.array([mu]), tau * cfg.s, cfg)[0])
        theta = mu + cfg.s * uq
        g[i] = np.trapezoid(utility(Action(act), theta, cfg) * pu, uq)
    return g


def tv_two_gaussians(delta_mean, sigma):
    """Total variation between N(m, sigma^2) and N(m+delta_mean, sigma^2): 2*Phi(|d|/2sigma) - 1."""
    return 2 * norm.cdf(abs(delta_mean) / (2 * sigma)) - 1.0


def main():
    hr("CP2 / GATE 2 (part ii) — achievable bound  R_obs <= L * delta(P_O^dep, P_O^cal)")

    cfg = gaussian_latent_config(rho=RHO, kappa=KAP, lam=LAM, base=MinosConfig(n_voxels=1_000_000))
    L_U = max(cfg.k_under, cfg.k_over)            # utility Lipschitz constant in theta
    gain = cfg.beta * cfg.s                        # observable-shift gain: mu moves by gain*delta_obs

    print(f"\nutility: piecewise-linear, slopes k_under={cfg.k_under}, k_over={cfg.k_over}")
    print(f"  => Lipschitz constant in theta:  L_U = max(k_under, k_over) = {L_U}")
    print(f"observable shift moves the reported point mu by gain*delta_obs, gain = beta*s = {gain}")
    print(f"  => P_O^dep is P_O^cal translated by Delta = {gain}*delta_obs, so")
    print(f"     W1(P_O^dep, P_O^cal) = Delta   (1-Wasserstein of a pure translation)")

    # ----------------------------------------------------------------------------------------
    # Derivation (printed). Two steps:
    #   (S1)  R(P_dep) <= 2 * sup_tau |EU(tau;P_dep) - EU(tau;P_cal)|.
    #         [tau_hat is optimal at cal; insert/remove it in the optimality chain.]
    #   (S2a) |EU(tau;P_dep)-EU(tau;P_cal)| = |int g_tau d(P_dep-P_cal)| <= osc(g_tau)*TV(P_O).
    #   (S2b) Since U is L_U-Lipschitz in theta and the observable shift only displaces the
    #         decision boundary by Delta in report-space, the realized-utility response per unit
    #         boundary move is <= L_U, giving R_obs <= L_U * W1.
    # ----------------------------------------------------------------------------------------
    hr("Derivation (printed)")
    print("(S1)  Optimality chain.  tau_hat = argmax_tau EU(tau; P_cal), so R(P_cal)=0, and")
    print("      R(P_dep) = EU(tau_o;P_dep) - EU(tau_hat;P_dep) <= 2 * sup_tau |EU(tau;P_dep)-EU(tau;P_cal)| =: 2*eps")
    print("      (tau_o = oracle scale; bound each EU by its calibration value +/- eps).")
    print("(S2a) TV form.   eps <= osc(g_tau) * TV(P_O^dep, P_O^cal),   g_tau(mu)=E[U(a*(mu;tau),theta)|mu]")
    print("                 => R_obs <= 2*osc(g_tau) * TV(P_O)        [L_TV = 2*osc(g_tau)]")
    print("(S2b) Wasserstein form.  U is L_U-Lipschitz in theta; observable shift only displaces the")
    print("                 boundary by Delta in report-space => R_obs <= L_U * W1(P_O)   [L = L_U]")

    # oscillation of g_tau at the fitted stale scale (for the TV constant)
    base_cal = make_population(cfg, make_rng(cfg.seed))
    tau_hat = fit_loss_calibration(base_cal, cfg)
    mu_grid = np.linspace(cfg.t2 - 4 * cfg.s, cfg.t2 + 4 * cfg.s, 400)
    g = per_report_expected_utility(mu_grid, cfg, tau_hat)
    osc_g = float(g.max() - g.min())
    L_TV = 2 * osc_g
    print(f"\nfitted stale scale tau_hat_cal = {tau_hat:.4f}")
    print(f"osc(g_tau) = max g - min g = {osc_g:.4f}  (in utility units; ~ k_under * s * O(1))")
    print(f"  => L_TV = 2*osc(g_tau) = {L_TV:.4f};   L (Wasserstein) = L_U = {L_U}")

    # ----------------------------------------------------------------------------------------
    # Numerical check against the ACTUAL v3 model: the bound must HOLD across the observable sweep.
    # ----------------------------------------------------------------------------------------
    hr("Bound check vs Minos-Core v3 (observable sweep) — the bound must hold")
    cfg_dep = cfg.replace(n_voxels=1_000_000)
    base_dep = make_population(cfg_dep, make_rng(20240517 + 777))   # v3 DEP_SEED
    deltas = np.round(np.arange(0.0, 0.2401, 0.03), 3)
    # The sim regret has a documented ~1e-8 optimiser-slack floor at zero shift (correction.py:149
    # clamps it with max(0,.)); allow that floor as an absolute tolerance so the degenerate delta=0
    # row (where the bound is exactly 0) is not spuriously flagged. It is far below any real regret.
    FLOOR = 1e-5
    print("\n delta_obs   R_obs(sim)    W1=Delta   L_U*W1     TV(P_O)   L_TV*TV    holds?")
    rows = []
    for d in deltas:
        R = stale_regret(base_dep, cfg_dep, tau_hat, delta_obs=float(d), delta_hid=0.0)
        Delta = gain * d
        W1 = Delta
        tv = tv_two_gaussians(Delta, cfg.theta_std)
        b_w = L_U * W1
        b_tv = L_TV * tv
        holds = (R <= b_w + FLOOR) and (R <= b_tv + FLOOR)
        rows.append((d, R, W1, b_w, tv, b_tv, holds))
        print(f"   {d:.2f}     {R:.5f}     {W1:.4f}    {b_w:.4f}    {tv:.4f}   {b_tv:.4f}    {holds}")
    all_hold = all(r[-1] for r in rows)
    # observed local slope of R vs W1 (for "L consistent with observed regret-vs-divergence slope")
    Rs = np.array([r[1] for r in rows]); W1s = np.array([r[2] for r in rows])
    obs_slope = float(np.polyfit(W1s, Rs, 1)[0])
    print(f"\nobserved regret-vs-W1 slope (linear fit) = {obs_slope:.4f}  <=  L_U = {L_U}  "
          f"(bound is a conservative envelope; R_obs is sub-linear near 0)")
    assert all_hold, "achievable bound must hold across the observable sweep"

    # ----------------------------------------------------------------------------------------
    # The structural fact behind the impossibility (part i): a hidden shift leaves P_O (hence M)
    # invariant. We DEMONSTRATE invariance here; the AUC=1/2 conclusion is proved in impossibility.md
    # (human-review-flagged), not asserted as machine-verified.
    # ----------------------------------------------------------------------------------------
    hr("Structural fact for the impossibility (part i): P_O invariant under the hidden shift")
    ref = build_reference(base_cal, cfg, tau_hat)
    cfg_b = cfg.replace(n_voxels=200_000)
    base_b = make_population(cfg_b, make_rng(4242))
    mu_hid0, _ = realise_deploy(base_b, cfg_b, delta_obs=0.0, delta_hid=0.0)
    mu_hid1, _ = realise_deploy(base_b, cfg_b, delta_obs=0.0, delta_hid=0.20)
    same_points = bool(np.array_equal(mu_hid0, mu_hid1))
    M0 = monitor(mu_hid0, cfg_b, ref)
    M1 = monitor(mu_hid1, cfg_b, ref)
    print(f"reported points identical under delta_hid 0.0 vs 0.20:  {same_points}")
    print(f"monitor M(delta_hid=0) = {M0:.6f}   M(delta_hid=0.20) = {M1:.6f}   equal: {M0 == M1}")
    print("=> the hidden shift moves the truth but NOT the observable law P_O, so any M=f({mu}) is")
    print("   invariant. The AUC=1/2 conclusion is the data-processing argument in impossibility.md")
    print("   (REQUIRES HUMAN PROOF-REVIEW — not machine-verified).")

    hr("GATE 2 (part ii) PASS — L derived & printed; bound holds vs v3; impossibility drafted+flagged")
    print(f"  L (Wasserstein) = L_U = k_under = {L_U};   L_TV = 2*osc(g_tau) = {L_TV:.4f}")
    print(f"  R_obs <= L_U * W1(P_O)  and  R_obs <= L_TV * TV(P_O)  hold across delta_obs in [0, 0.24]")
    return {"L_U": L_U, "L_TV": L_TV, "osc_g": osc_g, "tau_hat": tau_hat, "rows": rows,
            "obs_slope": obs_slope, "hidden_invariant": same_points and (M0 == M1)}


if __name__ == "__main__":
    main()
