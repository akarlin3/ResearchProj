"""CP4 — publication figures for Plumbline (deterministic, vector PDF; no re-fabrication).

Two figures, each regenerating from seed:

  Fig 1  fig1_scaling_law.pdf   — the gap scaling law: G vs gamma, with
                                    * the leading-order line  G_lead = (1/6)|z*(lambda)|*gamma  (Thm 1),
                                    * the full-theory curve   G_full = tau*_lead - tau_stat_exact,
                                    * the published v2 empirical points (RESULTS_B.md), and
                                    * the CP1 error-barred gap (multi-seed mean +/- 95% CI, gap_ci.py),
                                  with the gamma<=0.14 leading-order regime shaded and the operating
                                  point gamma=0.667 marked.

  Fig 2  fig2_monitor_bound.pdf — the label-free detectability picture: R_obs vs W1 with the
                                  achievable envelope L_U*W1 (Thm 2(ii)); the observable and hidden
                                  channels' regret (recomputed live via the minos API, deterministic),
                                  annotating observable AUC 1.00 vs hidden AUC 0.50 (Thm 2(i), v3).

Analytic curves are computed live (sympy/scipy). Empirical points are GATED CONSTANTS reused from
this session's printouts (v2: RESULTS_B.md; CP1: gap_ci.py) — never re-fit. The monitor sweep is
recomputed deterministically at the v3 cell/seed (no stored numbers). Output is vector PDF.

Run:  ``.venv-theory/bin/python theory/figures/make_figures.py``
"""
from __future__ import annotations

import os
import sys

import numpy as np
import mpmath as mp
from scipy.stats import norm, skewnorm
from scipy.optimize import brentq

import matplotlib
matplotlib.use("pdf")                       # deterministic vector backend, no display
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "minos-core"))

from minos.config import MinosConfig, gaussian_latent_config       # noqa: E402
from minos.correction import fit_loss_calibration, stale_regret    # noqa: E402
from minos.generative import make_population                       # noqa: E402
from minos.seeding import make_rng                                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
KAP, LAM, RHO = 3.0, 3.0, 0.5
L_REF = 0.90

# ---- gated empirical constants (reused from this session's printouts; NOT re-fit here) ----------
# v2 published gap sweep G(kappa) at lambda=3, rho=0.5 (RESULTS_B.md GATE 2, high-N single seed).
PUB_SWEEP_G = {0.0: -0.004, 0.5: -0.003, 1.0: +0.005, 2.0: +0.044, 3.0: +0.084, 4.0: +0.115}
# CP1 error-barred gap (gap_ci.py this session: B=32 @ N=5e5/kappa; mean + 95% bootstrap CI).
CP1_SWEEP_G = {  # kappa: (mean, ci_lo, ci_hi)
    0.0: (0.0046, 0.0008, 0.0083), 0.5: (0.0064, 0.0028, 0.0100),
    1.0: (0.0148, 0.0111, 0.0184), 2.0: (0.0500, 0.0465, 0.0537),
    3.0: (0.0892, 0.0854, 0.0930), 4.0: (0.1196, 0.1158, 0.1235),
}
# CP1 default-cell headline (gap_ci.py: B=64 @ N=1e6) — for the annotation.
CP1_DEFAULT_G = (0.0876, 0.0855, 0.0897)
REGIME_GAMMA = 0.14          # pure leading-order tracks G_full within 0.012 up to ~here (confirm.py)


# ---- analytic pieces (sympy result of gap_scaling.py, numerically evaluated) --------------------
def zstar(lam):
    if abs(lam - 1.0) < 1e-12:
        return 0.0
    psi = lambda t: t * float(mp.ncdf(t)) + float(mp.npdf(t))
    return float(mp.findroot(lambda t: (lam - 1) * psi(t) + t, -0.4))


def gamma_of_kappa(kappa):
    d = kappa / np.sqrt(1 + kappa ** 2)
    return ((4 - np.pi) / 2) * (d * np.sqrt(2 / np.pi)) ** 3 / (1 - 2 * d ** 2 / np.pi) ** 1.5


def taustat_exact(kappa, level=L_REF):
    """Exact coverage-model tau_stat (the 2nd-order coverage shrink the a=0 leading term omits)."""
    a = kappa
    d = a / np.sqrt(1 + a ** 2)
    mean_sn, sd_sn = d * np.sqrt(2 / np.pi), np.sqrt(1 - 2 * d ** 2 / np.pi)
    zL = norm.ppf(0.5 + level / 2.0)
    Fu = lambda x: skewnorm.cdf(x * sd_sn + mean_sn, a)
    return float(brentq(lambda tau: Fu(zL * tau) - Fu(-zL * tau) - level, 0.2, 5.0, xtol=1e-8))


# =================================================================================================
# Fig 1 — the scaling law
# =================================================================================================
def fig1_scaling_law():
    slope = abs(zstar(LAM)) / 6.0
    g_op = gamma_of_kappa(KAP)

    # analytic curves over a fine kappa grid -> gamma
    kg = np.linspace(0.0, 4.3, 300)
    gg = np.array([gamma_of_kappa(k) for k in kg])
    G_lead = slope * gg                                          # leading-order line (in gamma)
    G_full = np.array([(1 + slope * gamma_of_kappa(k)) - taustat_exact(k) for k in kg])

    fig, ax = plt.subplots(figsize=(6.2, 4.4))

    # leading-order regime shading (gamma <= REGIME_GAMMA, where pure leading-order tracks full)
    ax.axvspan(0.0, REGIME_GAMMA, color="#dfeaf7", alpha=0.7, lw=0,
               label=fr"leading-order regime ($\gamma\lesssim{REGIME_GAMMA:g}$)")

    ax.plot(gg, G_full, "-", color="#1f4e79", lw=2.2, zorder=3,
            label=r"full theory  $\tau^*_{\rm lead}-\tau_{\rm stat}^{\rm exact}$")
    ax.plot(gg, G_lead, "--", color="#c0504d", lw=2.0, zorder=3,
            label=r"leading order  $\frac{1}{6}|z^*(\lambda)|\,\gamma$")

    # v2 published empirical points
    gk = [gamma_of_kappa(k) for k in PUB_SWEEP_G]
    ax.plot(gk, list(PUB_SWEEP_G.values()), "s", ms=6, color="#7f7f7f",
            mec="k", mew=0.5, zorder=4, label="v2 empirical (RESULTS_B)")

    # CP1 error-barred gap (mean +/- 95% CI)
    gk_cp1 = np.array([gamma_of_kappa(k) for k in CP1_SWEEP_G])
    means = np.array([v[0] for v in CP1_SWEEP_G.values()])
    lo = np.array([v[1] for v in CP1_SWEEP_G.values()])
    hi = np.array([v[2] for v in CP1_SWEEP_G.values()])
    ax.errorbar(gk_cp1, means, yerr=[means - lo, hi - means], fmt="o", ms=5, color="#2e7d32",
                ecolor="#2e7d32", elinewidth=1.4, capsize=3, zorder=5,
                label="CP1 error-barred (95% CI)")

    # operating point gamma=0.667
    ax.axvline(g_op, color="#555555", ls=":", lw=1.2, zorder=2)
    ax.annotate(fr"operating point $\gamma={g_op:.2f}$", xy=(g_op, 0.0), xytext=(g_op - 0.02, -0.018),
                ha="right", va="top", fontsize=8.5, color="#333333", rotation=90)

    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"posterior skewness  $\gamma$", fontsize=11)
    ax.set_ylabel(r"decision–calibration gap  $G=\tau^*-\tau_{\rm stat}$", fontsize=11)
    ax.set_title("Theorem 1 — the gap is a scaling law in posterior skewness", fontsize=11.5)
    ax.set_xlim(-0.01, 0.82)
    ax.set_ylim(-0.025, 0.135)
    ax.legend(fontsize=8.3, loc="upper left", framealpha=0.95)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    out = os.path.join(HERE, "fig1_scaling_law.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out}")
    print(f"    slope=|z*({LAM:g})|/6={slope:.4f}; operating gamma={g_op:.4f}; "
          f"G_lead(op)={slope*g_op:.4f}; G_full(op)={(1+slope*g_op)-taustat_exact(KAP):.4f}")


# =================================================================================================
# Fig 2 — the monitor bound (observable vs hidden), recomputed live & deterministically
# =================================================================================================
def fig2_monitor_bound():
    cfg = gaussian_latent_config(rho=RHO, kappa=KAP, lam=LAM, base=MinosConfig(n_voxels=1_000_000))
    L_U = max(cfg.k_under, cfg.k_over)
    gain = cfg.beta * cfg.s
    base_cal = make_population(cfg, make_rng(cfg.seed))
    tau_hat = fit_loss_calibration(base_cal, cfg)
    base_dep = make_population(cfg, make_rng(20240517 + 777))     # v3 DEP_SEED (matches confirm.py)
    deltas = np.round(np.arange(0.0, 0.2401, 0.03), 3)
    W1 = gain * deltas
    R_obs = np.array([stale_regret(base_dep, cfg, tau_hat, delta_obs=float(d), delta_hid=0.0) for d in deltas])
    R_hid = np.array([stale_regret(base_dep, cfg, tau_hat, delta_obs=0.0, delta_hid=float(d)) for d in deltas])
    env = L_U * W1

    obs_slope = float(np.polyfit(W1, R_obs, 1)[0])     # over all deltas, matching detectability.py
    fig, ax = plt.subplots(figsize=(6.4, 4.6))

    # --- main panel: the achievable envelope contains the realised regret ---
    ax.plot(W1, env, "-", color="#1f4e79", lw=2.0, zorder=2,
            label=fr"achievable bound  $R_{{\rm obs}}\leq L_U\,W_1$  ($L_U={L_U:g}$)")
    ax.fill_between(W1, env, R_obs, color="#dfeaf7", alpha=0.7, zorder=1)
    ax.plot(W1, R_obs, "o-", color="#c0504d", lw=1.8, ms=5, zorder=3,
            label=r"observable channel $R_{\rm obs}$  (detectable, AUC$=1.00$)")
    ax.plot(W1, R_hid, "^--", color="#2e7d32", lw=1.8, ms=5, zorder=3,
            label=r"hidden channel $R_{\rm hid}$  (undetectable, AUC$=0.50$)")
    ax.annotate(fr"conservative envelope" "\n" fr"(slope $L_U={L_U:g}$; observed $R$-slope ${obs_slope:.2f}$)",
                xy=(0.17, L_U * 0.17), xytext=(0.035, 0.46),
                fontsize=8.2, color="#1f4e79",
                arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0))

    # --- inset: zoom on the regret scale so the detectability contrast is legible ---
    axin = ax.inset_axes([0.50, 0.30, 0.46, 0.42])
    axin.plot(W1, R_obs, "o-", color="#c0504d", lw=1.6, ms=4)
    axin.plot(W1, R_hid, "^--", color="#2e7d32", lw=1.6, ms=4)
    axin.set_title("regret detail (same scale, two channels)", fontsize=7.8)
    axin.annotate("observable\nAUC = 1.00", xy=(W1[-1], R_obs[-1]),
                  xytext=(W1[-1] * 0.44, R_obs[-1] - 0.085), fontsize=7.6, color="#c0504d",
                  arrowprops=dict(arrowstyle="->", color="#c0504d", lw=0.9))
    axin.annotate("hidden: same-order\nregret, AUC = 0.50", xy=(W1[-1], R_hid[-1]),
                  xytext=(W1[1] * 0.2, R_hid[-1] + 0.02), fontsize=7.6, color="#2e7d32",
                  arrowprops=dict(arrowstyle="->", color="#2e7d32", lw=0.9))
    axin.set_xlim(0.0, W1[-1] * 1.02)
    axin.set_ylim(0.0, max(R_obs[-1], R_hid[-1]) * 1.30)
    axin.grid(True, alpha=0.25)
    axin.tick_params(labelsize=7)

    ax.set_xlabel(r"observable divergence  $W_1(P_O^{\rm dep},P_O^{\rm cal})=\beta s\,\delta$", fontsize=11)
    ax.set_ylabel(r"stale-correction regret  $R$", fontsize=11)
    ax.set_title("Theorem 2 — observable regret is bounded-detectable;\nhidden regret is invisible", fontsize=11.5)
    ax.set_xlim(0.0, W1[-1] * 1.02)
    ax.set_ylim(0.0, env[-1] * 1.08)
    ax.legend(fontsize=8.3, loc="upper left", framealpha=0.95)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    out = os.path.join(HERE, "fig2_monitor_bound.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out}")
    print(f"    L_U={L_U:g}; gain=beta*s={gain:g}; tau_hat_cal={tau_hat:.4f}; "
          f"R_obs(max)={R_obs[-1]:.4f} <= envelope {env[-1]:.4f}; R_hid(max)={R_hid[-1]:.4f}")
    # the bound must hold for delta>0 (degenerate delta=0 row is ~0 with optimiser-slack floor)
    assert all(ro <= e + 1e-5 for ro, e, d in zip(R_obs, env, deltas) if d > 0), \
        "Fig 2: observable regret must stay under the L_U*W1 envelope"


def main():
    print("CP4 — generating publication figures (deterministic, vector PDF):")
    fig1_scaling_law()
    fig2_monitor_bound()
    print("done.")


if __name__ == "__main__":
    main()
