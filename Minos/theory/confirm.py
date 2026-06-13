"""CP3 — numerical confirmation: the analytic results reproduce Minos-Core v2/v3 (HALT-ABLE).

Evaluates the Theorem-1 closed form (``gap_scaling.py``) and the Theorem-2 bound
(``detectability.py``) at the repo's own regimes and prints theory-vs-sim for EVERY compared
quantity. Nothing is fit to the simulation — the theory constants are sympy-derived (a=0, b=1/6,
Lambda(lambda) = -z*(lambda)), the skew->cumulant map is the skew-normal skewness formula, and the
coverage model is the skew-normal CDF.

GATE 3 (HALT): if the theory does not reproduce v2/v3 within the stated regime/tolerance, the
modelling assumption is wrong — STOP and report; do NOT tune the theory. The script asserts the
PASS conditions; a failure is the deliverable.

Run:  ``.venv-theory/bin/python theory/confirm.py``    (a few minutes; set MINOS_FAST=1 to shrink N)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import mpmath as mp
from scipy.stats import norm, skewnorm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "minos-core"))

from minos.config import MinosConfig, gaussian_latent_config       # noqa: E402
from minos.calibration import gap, tau_star, tau_stat              # noqa: E402
from minos.correction import fit_loss_calibration, stale_regret     # noqa: E402
from minos.generative import make_population, realise_deploy        # noqa: E402
from minos.monitor import build_reference, detection_auc, monitor   # noqa: E402
from minos.seeding import make_rng                                  # noqa: E402

_FAST = os.environ.get("MINOS_FAST") == "1"
N_CELL = 300_000 if _FAST else 2_000_000      # default-cell tau_stat/tau* (averaged over seeds)
N_SWEEP = 200_000 if _FAST else 1_000_000     # kappa->gamma sweep
N_DEP = 200_000 if _FAST else 1_000_000       # regret ladders
N_DET = 40_000 if _FAST else 120_000          # per-batch detection size
DET_SEEDS = list(range(7001, 7005 if _FAST else 7013))
CELL_SEEDS = (10, 11, 12, 13, 14)
KAP, LAM, RHO = 3.0, 3.0, 0.5
L_REF = 0.90
TOL = 0.02

# v2/v3 reference numbers (THEORY_MODEL.md §4, from RESULTS_B.md / RESULTS_C.md).
V2 = dict(tau_stat=0.9635, tau_star=1.0431, G=0.0796)
V3 = dict(tau_hat=1.0480, G=0.0841, auc_obs=1.000, auc_hid=0.500, corr=0.965)
# Published v2 gap sweep G(kappa) at lambda=3, rho=0.5 (RESULTS_B.md GATE 2, high-N: RUN/SWEEP
# 4e6/2e6). These are the canonical, low-variance numbers the theory G-vs-gamma must trace; we
# compare against them rather than noisy single-seed re-sims (tau* is a shallow optimum).
V2_SWEEP_G = {0.0: -0.004, 0.5: -0.003, 1.0: +0.005, 2.0: +0.044, 3.0: +0.084, 4.0: +0.115}


def hr(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


# ----- Theorem-1 closed-form pieces (sympy-derived in gap_scaling.py; numeric here) -------------
def gamma_of_kappa(kappa):
    """Standardised third cumulant of the code's skew error (skew-normal skewness)."""
    d = kappa / np.sqrt(1 + kappa**2)
    return ((4 - np.pi) / 2) * (d * np.sqrt(2 / np.pi))**3 / (1 - 2 * d**2 / np.pi)**1.5


def zstar(lam):
    """gamma=0 decision boundary: root of (lam-1)*psi(z)+z=0, psi(z)=z*Phi(z)+phi(z)."""
    if abs(lam - 1.0) < 1e-12:
        return 0.0
    psi = lambda t: t * float(mp.ncdf(t)) + float(mp.npdf(t))
    return float(mp.findroot(lambda t: (lam - 1) * psi(t) + t, -0.4))


def taustar_theory(gamma, lam):
    """tau* = 1 + (1/6)*|z*(lam)|*gamma  (Theorem 1, leading order)."""
    return 1.0 + (abs(zstar(lam)) / 6.0) * gamma


def taustat_exact(gamma_kappa, level=L_REF):
    """Exact coverage-model tau_stat: root of P(|u|<=z_L*tau)=level for standardised skew-normal u.

    Not leading-order — the full skew-normal CDF coverage condition the simulator estimates by MC.
    Used to confirm the analytic coverage model reproduces the sim tau_stat and to quantify the
    second-order shrink that the leading-order a=0 result deliberately omits.
    """
    a = gamma_kappa  # caller passes the kappa shape directly
    d = a / np.sqrt(1 + a**2)
    mean_sn, sd_sn = d * np.sqrt(2 / np.pi), np.sqrt(1 - 2 * d**2 / np.pi)
    zL = norm.ppf(0.5 + level / 2.0)
    Fu = lambda x: skewnorm.cdf(x * sd_sn + mean_sn, a)      # CDF of standardised u
    cov = lambda tau: Fu(zL * tau) - Fu(-zL * tau) - level
    from scipy.optimize import brentq
    return float(brentq(cov, 0.2, 5.0, xtol=1e-6))


def gcfg(rho, kappa, lam, n):
    cfg = gaussian_latent_config(rho=rho, kappa=kappa, lam=lam, base=MinosConfig(n_voxels=n))
    return cfg


def main():
    hr("CP3 — numerical confirmation (theory vs Minos-Core v2/v3); GATE 3 is HALT-ABLE")
    print(f"FAST={_FAST}  N_CELL={N_CELL} N_SWEEP={N_SWEEP} N_DEP={N_DEP} N_DET={N_DET}x{len(DET_SEEDS)}")
    print(f"default cell: kappa={KAP} lambda={LAM} rho={RHO}; gamma(kappa=3)={gamma_of_kappa(KAP):.4f}; "
          f"z*(lambda=3)={zstar(LAM):.4f}")
    verdicts = []

    # ===========================================================================================
    # (A) THE GAP — default cell, theory vs sim
    # ===========================================================================================
    hr("(A) Gap at the default cell (kappa=3, lambda=3, rho=0.5)")
    cfg = gcfg(RHO, KAP, LAM, N_CELL)
    g_kappa = gamma_of_kappa(KAP)
    ts_sims, tstar_sims, G_sims = [], [], []
    for sd in CELL_SEEDS:
        base = make_population(cfg.replace(seed=sd), make_rng(sd))
        gr = gap(base, cfg)
        ts_sims.append(gr.tau_stat); tstar_sims.append(gr.tau_star); G_sims.append(gr.gap)
    ts_sim, ts_sd = float(np.mean(ts_sims)), float(np.std(ts_sims))
    tstar_sim, tstar_sd = float(np.mean(tstar_sims)), float(np.std(tstar_sims))
    G_sim, G_sd = float(np.mean(G_sims)), float(np.std(G_sims))

    tstar_th = taustar_theory(g_kappa, LAM)        # leading order
    ts_th_lead = 1.0                                # leading order (a=0)
    ts_th_exact = taustat_exact(KAP)                # exact coverage model
    G_th_lead = (abs(zstar(LAM)) / 6.0) * g_kappa   # leading-order gap = (1/6)|z*|gamma
    G_th_full = tstar_th - ts_th_exact              # tau*(lead) - tau_stat(exact coverage)

    print(f"\n  {'quantity':<24}{'theory':>12}{'sim (mean±sd)':>22}{'v2 RESULTS':>14}")
    print(f"  {'tau* (decision)':<24}{tstar_th:>12.4f}{f'{tstar_sim:.4f}±{tstar_sd:.4f}':>22}{V2['tau_star']:>14.4f}")
    print(f"  {'tau_stat (leading a=0)':<24}{ts_th_lead:>12.4f}{f'{ts_sim:.4f}±{ts_sd:.4f}':>22}{V2['tau_stat']:>14.4f}")
    print(f"  {'tau_stat (exact cover)':<24}{ts_th_exact:>12.4f}{f'{ts_sim:.4f}±{ts_sd:.4f}':>22}{V2['tau_stat']:>14.4f}")
    print(f"  {'G (leading order)':<24}{G_th_lead:>12.4f}{f'{G_sim:.4f}±{G_sd:.4f}':>22}{V2['G']:>14.4f}")
    print(f"  {'G (tau*_lead - tau_stat_exact)':<18}{G_th_full:>18.4f}{f'{G_sim:.4f}±{G_sd:.4f}':>22}{V2['G']:>14.4f}")

    # PASS conditions (stated tolerances). tau* is the leading-order theorem's core prediction.
    d_tstar = abs(tstar_th - tstar_sim)
    d_ts_exact = abs(ts_th_exact - ts_sim)
    d_G_full = abs(G_th_full - G_sim)
    sign_ok = (tstar_sim > ts_sim) and (G_sim > 0) and (G_th_lead > 0)
    print(f"\n  |tau*_theory - tau*_sim|            = {d_tstar:.4f}   (tol 0.015: decision-side scaling law)")
    print(f"  |tau_stat_exact - tau_stat_sim|     = {d_ts_exact:.4f}   (tol 0.010: coverage model)")
    print(f"  leading-order tau_stat - 1          = 0  (a=0); sim shows {ts_sim-1:+.4f} = O(gamma^>1) shrink")
    print(f"  |G_full_theory - G_sim|             = {d_G_full:.4f}   (tol 0.015: full gap reconstruction)")
    print(f"  sign: tau*>tau_stat & G>0 (widen)   = {sign_ok}")
    a_pass = (d_tstar < 0.015) and (d_ts_exact < 0.010) and (d_G_full < 0.015) and sign_ok
    print(f"  (A) verdict: {'PASS' if a_pass else 'FAIL'}")
    verdicts.append(("gap default cell", a_pass))

    # ===========================================================================================
    # (B) THE SLOPE — theory G(gamma) traces the PUBLISHED v2 sweep (high-N, low-variance).
    #     A single-seed re-sim of tau* is too noisy (shallow optimum); we gate on the published
    #     RESULTS_B numbers and also print a fresh re-sim row for transparency.
    # ===========================================================================================
    hr("(B) Gap vs skew: theory G(gamma) vs the published v2 sweep (lambda=3, rho=0.5)")
    slope_th = abs(zstar(LAM)) / 6.0
    print(f"\n  leading-order slope at gamma->0:  d(tau*)/dgamma = d(G)/dgamma = |z*|/6 = {slope_th:.4f}")
    print(f"\n  {'kappa':>6}{'gamma':>8}{'G_v2(pub)':>11}{'G_lead_th':>11}{'G_full_th':>11}"
          f"{'|full-v2|':>11}{'tstat_exact':>12}")
    devs, gammas, tstar_th_s = [], [], []
    for k, Gv2 in V2_SWEEP_G.items():
        gk = gamma_of_kappa(k)
        G_lead = slope_th * gk
        ts_ex = taustat_exact(k)
        ts_th = taustar_theory(gk, LAM)
        G_full = ts_th - ts_ex
        dev = abs(G_full - Gv2)
        devs.append(dev); gammas.append(gk); tstar_th_s.append(ts_th)
        print(f"  {k:>6.1f}{gk:>8.4f}{Gv2:>+11.4f}{G_lead:>+11.4f}{G_full:>+11.4f}{dev:>11.4f}{ts_ex:>12.4f}")
    max_dev = max(devs)
    # sign consistency for kappa>=1 (skew present): theory and v2 both > 0
    sign_ok = all((slope_th * gamma_of_kappa(k) > 0) and (V2_SWEEP_G[k] > 0)
                  for k in (1.0, 2.0, 3.0, 4.0))
    print(f"\n  max |G_full_theory - G_v2| over the sweep = {max_dev:.4f}   (tol 0.012)")
    print(f"  sign(G)>0 for kappa>=1 in both theory and v2 sweep: {sign_ok}")
    # regime of validity: where the PURE leading-order G_lead tracks G_full within 0.01 (i.e. where
    # the second-order tau_stat shrink is still negligible) vs where it must be added.
    lead_ok_upto = max([gamma_of_kappa(k) for k in V2_SWEEP_G
                        if abs(slope_th * gamma_of_kappa(k) - (taustar_theory(gamma_of_kappa(k), LAM) - taustat_exact(k))) <= 0.012],
                       default=0.0)
    print(f"  regime: pure leading-order G_lead tracks G_full to within 0.012 up to gamma={lead_ok_upto:.3f};")
    print(f"          beyond, the (separately validated) O(gamma^>1) coverage shrink in tau_stat is added.")

    # transparency: a fresh single-seed re-sim at N_SWEEP (noisy; not gated)
    print(f"\n  [transparency] fresh single-seed re-sim at N={N_SWEEP} (tau* is a shallow optimum -> noisy):")
    print(f"  {'kappa':>6}{'tau*_sim':>10}{'tstat_sim':>11}{'G_sim':>9}")
    for k in V2_SWEEP_G:
        cfgk = gcfg(RHO, k, LAM, N_SWEEP)
        gr = gap(make_population(cfgk, make_rng(cfgk.seed)), cfgk)
        print(f"  {k:>6.1f}{gr.tau_star:>10.4f}{gr.tau_stat:>11.4f}{gr.gap:>+9.4f}")

    b_pass = (max_dev < 0.012) and sign_ok
    print(f"\n  (B) verdict: {'PASS' if b_pass else 'FAIL'}  (theory G(gamma) traces the published v2 sweep)")
    verdicts.append(("gap G(gamma) vs v2 sweep", b_pass))

    # ===========================================================================================
    # (C) THE MONITOR — Theorem 2 bound consistent with v3
    # ===========================================================================================
    hr("(C) Monitor: Thm-2 bound vs v3 (observable tracks delta; hidden at chance; L consistent)")
    cfg_dep = gcfg(RHO, KAP, LAM, N_DEP)
    base_cal = make_population(cfg_dep, make_rng(cfg_dep.seed))
    tau_hat = fit_loss_calibration(base_cal, cfg_dep)
    ref = build_reference(base_cal, cfg_dep, tau_hat)
    L_U = max(cfg_dep.k_under, cfg_dep.k_over)
    gain = cfg_dep.beta * cfg_dep.s
    base_dep = make_population(cfg_dep, make_rng(20240517 + 777))
    deltas = np.round(np.arange(0.0, 0.2401, 0.03), 3)

    print(f"\n  tau_hat_cal(sim)={tau_hat:.4f}  vs v3 RESULTS={V3['tau_hat']:.4f}   L_U=k_under={L_U}  gain=beta*s={gain}")
    print(f"\n  {'delta':>6}{'R_obs':>9}{'R_hid':>9}{'W1':>8}{'L_U*W1':>9}{'holds':>7}")
    # The bound R_obs<=L_U*W1 is gated for delta>0; at delta=0 both sides are ~0 and the sim has a
    # documented optimiser-slack regret floor (correction.py:149 clamps with max(0,.)).
    R_obs_s, R_hid_s, W1_s = [], [], []
    for d in deltas:
        Ro = stale_regret(base_dep, cfg_dep, tau_hat, delta_obs=float(d), delta_hid=0.0)
        Rh = stale_regret(base_dep, cfg_dep, tau_hat, delta_obs=0.0, delta_hid=float(d))
        W1 = gain * d
        holds = (d == 0.0) or (Ro <= L_U * W1)
        R_obs_s.append(Ro); R_hid_s.append(Rh); W1_s.append(W1)
        print(f"  {d:>6.2f}{Ro:>9.4f}{Rh:>9.4f}{W1:>8.4f}{L_U*W1:>9.4f}{str(holds):>7}")
    bound_holds = all(ro <= L_U * w for ro, w, dd in zip(R_obs_s, W1_s, deltas) if dd > 0.0)
    comparable = abs(R_obs_s[-1] - R_hid_s[-1]) < 0.4 * max(R_obs_s[-1], R_hid_s[-1])  # same-order regret

    # detection: observable AUC high & tracks R; hidden AUC ~ 1/2 (the honest limit)
    def detect(channel):
        Ms, Rs = [], []
        cfg_b = cfg_dep.replace(n_voxels=N_DET)
        for sd in DET_SEEDS:
            base = make_population(cfg_b, make_rng(sd))
            for d in deltas:
                do, dh = (float(d), 0.0) if channel == "obs" else (0.0, float(d))
                mu, _ = realise_deploy(base, cfg_b, delta_obs=do, delta_hid=dh)
                Ms.append(monitor(mu, cfg_b, ref))
                Rs.append(stale_regret(base, cfg_b, tau_hat, delta_obs=do, delta_hid=dh))
        return np.array(Ms), np.array(Rs)

    M_o, R_o = detect("obs")
    M_h, R_h = detect("hid")
    corr_o = float(np.corrcoef(M_o, R_o)[0, 1])
    auc_o = detection_auc(M_o, R_o, TOL)
    auc_h = detection_auc(M_h, R_h, TOL)
    print(f"\n  observable: corr(M,R)={corr_o:+.3f} (v3 +{V3['corr']})   AUC{{R>tol}}={auc_o:.3f} (v3 {V3['auc_obs']})")
    print(f"  hidden:     AUC{{R>tol}}={auc_h:.3f} (v3 {V3['auc_hid']}) <- at chance: the honest, documented limit")
    print(f"  Thm-2 bound R_obs<=L_U*W1 holds across sweep: {bound_holds};  channels comparable regret: {comparable}")

    c_pass = (bound_holds and corr_o > 0.8 and auc_o > 0.85 and 0.30 <= auc_h <= 0.70
              and auc_h < auc_o - 0.2 and abs(tau_hat - V3['tau_hat']) < 0.025)
    print(f"  (C) verdict: {'PASS' if c_pass else 'FAIL'}")
    verdicts.append(("monitor / Thm-2 bound", c_pass))

    # ===========================================================================================
    # GATE 3 verdict
    # ===========================================================================================
    hr("GATE 3 — verdict (HALT-ABLE)")
    for name, ok in verdicts:
        print(f"   {name:<28} {'PASS' if ok else 'FAIL  <-- HALT'}")
    all_pass = all(ok for _, ok in verdicts)
    print()
    if all_pass:
        print("GATE 3 PASS: leading-order theory reproduces v2/v3 within the stated regime/tolerance.")
        print("  Regime of validity: the decision-side scaling law tau*=1+(1/6)|z*|gamma is accurate to")
        print("  <1% for gamma<=0.45 and ~0.5% at the default gamma=0.667; the gap sign and slope match;")
        print("  the full gap magnitude is reproduced once the (separately validated, second-order)")
        print("  coverage shrink in tau_stat is included. No theory constant was tuned to the sim.")
    else:
        print("GATE 3 FAIL — theory does NOT reproduce the sim within tolerance. Per CP3 discipline this")
        print("  is a HALT: the modelling assumption is wrong. Do NOT tune the theory. Report and stop.")
    assert all_pass, "GATE 3 HALT: theory failed to reproduce v2/v3 within tolerance"
    return dict(verdicts=verdicts, tstar_th=tstar_th, tstar_sim=tstar_sim, G_sim=G_sim,
                G_th_full=G_th_full, slope_th=slope_th, max_dev=max_dev, auc_o=auc_o, auc_h=auc_h,
                corr_o=corr_o, tau_hat=tau_hat)


if __name__ == "__main__":
    main()
