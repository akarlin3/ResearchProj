"""Minos-Core v2 experiment driver — the decision-calibration gap (Headline B).

Seeded, config-driven. Reproduces the four checkpoint gates from a clean seed, prints every
number RESULTS_B.md cites, and writes the four publication figures as vector PDFs:

  (a) VoC(tau) well-specified vs misspecified, marking tau_stat and tau* — the gap made visible
  (b) headline: gap G vs kappa at several lambda
  (c) G surface over (kappa, rho)
  (d) VoTG(delta) with the break-even shift delta_be marked

Run from the project root:  ``python experiments/run_b.py``
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from minos.calibration import (  # noqa: E402
    break_even_shift,
    gap,
    posterior_eu,
    tau_stat,
    tau_star,
)
from minos.config import BASELINE_V1, MinosConfig, gaussian_latent_config  # noqa: E402
from minos.gate import votg  # noqa: E402
from minos.generative import make_population  # noqa: E402
from minos.seeding import make_rng  # noqa: E402
from minos.voi import voc  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(os.path.dirname(HERE), "figures")

# Estimation sizes (CRN; one base draw per config). Headline numbers use RUN_N; the sweep
# grids (many configs) use SWEEP_N to stay reproducible at reasonable cost.
RUN_N = 4_000_000
SWEEP_N = 2_000_000

# Sweep axes (DESIGN_B Section 4).
KAPPAS = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0]
LAMBDAS = [1.0, 2.0, 3.0, 4.0]
# rho is kept in the resolved regime [0, 2]: far from the threshold the decision is
# tau-insensitive and tau* is unidentified (DESIGN_B Section 5).
RHOS = [0.0, 0.5, 1.0, 1.5, 2.0]
KAPPAS_C = [0.0, 1.0, 2.0, 3.0, 4.0]   # kappa axis for the (kappa, rho) surface
RHO0 = 0.5                              # reference proximity for the kappa/lambda sweeps
LAM0 = 3.0                              # reference asymmetry for kappa/rho sweeps
KAP0 = 3.0                              # reference skew for lambda/rho slices
DELTAS = np.round(np.arange(0.0, 1.6001, 0.1), 3)
L_PROFILE = [0.5, 0.68, 0.8, 0.9, 0.95, 0.99]

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def hr(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def gcfg(rho, kappa, lam, n):
    cfg = gaussian_latent_config(rho=rho, kappa=kappa, lam=lam, base=MinosConfig(n_voxels=n))
    return cfg, make_population(cfg, make_rng(cfg.seed))


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    out = {}

    hr("CONFIG")
    ref = gaussian_latent_config(rho=RHO0, kappa=KAP0, lam=LAM0, base=MinosConfig())
    print(f"seed={ref.seed}  RUN_N={RUN_N}  SWEEP_N={SWEEP_N}")
    print(f"thresholds t1={ref.t1:.1f} t2={ref.t2:.1f} (single active threshold t2); s={ref.s}")
    print(f"gap latent: gaussian report centres, posterior_centric={ref.posterior_centric}, "
          f"theta_std={ref.theta_std}")
    print(f"axes: kappa={KAPPAS}  lambda={LAMBDAS}  rho={RHOS}")
    print(f"reference cell for headline: kappa={KAP0} lambda={LAM0} rho={RHO0}  "
          f"(tau_stat level L_ref=0.90)")

    # ---- GATE 1: gap estimator --------------------------------------------------
    hr("GATE 1 — gap estimator (definitions + well-specified limit + misspecified sign)")
    cfg0, base0 = gcfg(RHO0, 0.0, LAM0, RUN_N)
    g0 = gap(base0, cfg0)
    print(f"well-specified kappa=0 : tau_stat={g0.tau_stat:.4f}  tau*={g0.tau_star:.4f}  "
          f"G={g0.gap:+.4f}  ratio={g0.ratio:.4f}")
    assert abs(g0.gap) < 0.02, "well-specified gap should vanish"
    assert abs(g0.tau_stat - 1.0) < 0.01 and abs(g0.tau_star - 1.0) < 0.03

    # v1 mixture reproduction: VoC argmin (= tau*) at tau=1, and the v1 VoC numbers reproduce.
    cfgm = MinosConfig(n_voxels=RUN_N)
    basem = make_population(cfgm, make_rng(cfgm.seed))
    tstar_v1, _ = tau_star(basem, cfgm)
    # Exact reproduction at the v1 sample size (n=1e6, seed 20240517) -> v1 RESULTS.md numbers.
    cfgv1 = MinosConfig(n_voxels=1_000_000)
    basev1 = make_population(cfgv1, make_rng(cfgv1.seed))
    voc_05, voc_20 = voc(basev1, cfgv1, 0.5), voc(basev1, cfgv1, 2.0)
    print(f"v1 mixture reproduction: VoC argmin tau*={tstar_v1:.3f}; at n=1e6 "
          f"VoC(0.5)={voc_05:+.6f} VoC(2.0)={voc_20:+.6f}  (v1 RESULTS.md: +0.001169 / +0.004361)")
    assert abs(tstar_v1 - 1.0) < 0.05, "v1 symmetric baseline must peak at tau=1"
    assert abs(voc_05 - 0.001169) < 1e-5 and abs(voc_20 - 0.004361) < 1e-5, "v1 VoC must reproduce"

    cfgd, based = gcfg(RHO0, KAP0, LAM0, RUN_N)
    gd = gap(based, cfgd)
    print(f"default misspecified   : tau_stat={gd.tau_stat:.4f}  tau*={gd.tau_star:.4f}  "
          f"G={gd.gap:+.4f}  ratio={gd.ratio:.4f}")
    assert gd.gap > 0.03 and gd.tau_star > gd.tau_stat, "misspecified gap should be positive"
    assert gd.tau_stat < 1.0 < gd.tau_star + 1e-9, "calibration shrinks, decision holds/widens"

    # Independence of the two code paths.
    cfgd2 = cfgd.replace(k_under=cfgd.k_under * 3.0)
    based2 = make_population(cfgd2, make_rng(cfgd2.seed))
    ts_same = tau_stat(based, cfgd) == tau_stat(based2, cfgd2)
    print(f"independence: tau_stat unchanged when utility x3 -> {ts_same}; "
          f"tau* takes no coverage level -> True")
    assert ts_same, "tau_stat must be invariant to the utility (disjoint path)"

    # tau_stat level profile (robustness of the gap to the coverage criterion).
    prof = {L: tau_stat(based, cfgd, level=L) for L in L_PROFILE}
    print("tau_stat(L) profile @ default misspec:",
          {L: round(v, 3) for L, v in prof.items()}, f"  tau*={gd.tau_star:.3f}")
    print("GATE 1 PASS: G~0 well-specified & v1 reproduces; G>0 (tau*>tau_stat) misspecified; "
          "paths independent")
    out.update(g0=g0, gd=gd, prof=prof, tstar_v1=tstar_v1, voc_05=voc_05, voc_20=voc_20)

    # ---- GATE 2: gap map --------------------------------------------------------
    hr("GATE 2 — gap map: kappa x lambda and kappa x rho")
    grid_kl = np.zeros((len(KAPPAS), len(LAMBDAS)))
    for i, k in enumerate(KAPPAS):
        for j, lam in enumerate(LAMBDAS):
            cfg, base = gcfg(RHO0, k, lam, SWEEP_N)
            grid_kl[i, j] = gap(base, cfg).gap
    print("G(kappa, lambda) at rho=0.5  [rows kappa, cols lambda]:")
    print("       " + "  ".join(f"l={l:<4}" for l in LAMBDAS))
    for i, k in enumerate(KAPPAS):
        print(f"k={k:<4} " + "  ".join(f"{grid_kl[i, j]:+.3f}" for j in range(len(LAMBDAS))))

    grid_kr = np.zeros((len(KAPPAS_C), len(RHOS)))
    for i, k in enumerate(KAPPAS_C):
        for j, r in enumerate(RHOS):
            cfg, base = gcfg(r, k, LAM0, SWEEP_N)
            grid_kr[i, j] = gap(base, cfg).gap
    print("\nG(kappa, rho) at lambda=3  [rows kappa, cols rho]:")
    print("       " + "  ".join(f"r={r:<4}" for r in RHOS))
    for i, k in enumerate(KAPPAS_C):
        print(f"k={k:<4} " + "  ".join(f"{grid_kr[i, j]:+.3f}" for j in range(len(RHOS))))

    # Sanity checks.
    j3 = LAMBDAS.index(LAM0)
    kap_col = grid_kl[:, j3]                       # G vs kappa at lambda=3
    i3 = KAPPAS.index(KAP0)
    lam_row = grid_kl[i3, :]                       # G vs lambda at kappa=3
    ik3 = KAPPAS_C.index(KAP0)
    rho_row = grid_kr[ik3, :]                      # G vs rho at kappa=3
    rho_slope = float(np.polyfit(RHOS, rho_row, 1)[0])
    print(f"\nkappa=0 slice (should be ~0): lambda-row max |G| = {np.abs(grid_kl[0]).max():.3f}, "
          f"rho-row max |G| = {np.abs(grid_kr[0]).max():.3f}")
    print(f"monotone increasing in kappa (k>=1) at lambda=3: {[round(float(v),3) for v in kap_col[2:]]}")
    print(f"monotone increasing in lambda at kappa=3       : {[round(float(v),3) for v in lam_row]}")
    print(f"decreasing in rho at kappa=3 (largest at rho=0): {[round(float(v),3) for v in rho_row]}  "
          f"slope={rho_slope:+.4f}")
    assert np.abs(grid_kl[0]).max() < 0.025 and np.abs(grid_kr[0]).max() < 0.025
    assert np.all(np.diff(kap_col[2:]) > 0), "G not increasing in kappa"
    assert np.all(np.diff(lam_row) > 0), "G not increasing in lambda"
    # rho: largest at the threshold (rho=0), a clear downward trend (robust to MC wiggles).
    assert int(np.argmax(rho_row)) == 0, "G not largest at rho=0"
    assert rho_slope < 0 and (rho_row[0] - rho_row[-1]) > 0.01, "G not decreasing in rho"
    print("corner table G: "
          f"(k0,l1)={grid_kl[0,0]:+.3f} (k0,l4)={grid_kl[0,-1]:+.3f} "
          f"(k4,l1)={grid_kl[-1,0]:+.3f} (k4,l4)={grid_kl[-1,-1]:+.3f} | "
          f"(k0,r0)={grid_kr[0,0]:+.3f} (k4,r0)={grid_kr[-1,0]:+.3f} "
          f"(k4,r4)={grid_kr[-1,-1]:+.3f}")
    print("GATE 2 PASS: kappa=0 flat-zero; G up in kappa, up in lambda, down in rho")
    out.update(grid_kl=grid_kl, grid_kr=grid_kr)

    # ---- GATE 3: break-even shift -----------------------------------------------
    hr("GATE 3 — break-even shift for the trust-gate")
    cfgb = BASELINE_V1.replace(n_voxels=RUN_N)
    baseb = make_population(cfgb, make_rng(cfgb.seed))
    votg_curve = np.array([votg(baseb, cfgb, delta=d) for d in DELTAS])
    d_be = break_even_shift(baseb, cfgb)
    v_lo = votg(baseb, cfgb, delta=d_be - 0.1)
    v_hi = votg(baseb, cfgb, delta=d_be + 0.2)
    print(f"VoTG(delta=0)        = {votg_curve[0]:+.6f}")
    print(f"break-even delta_be  = {d_be:.4f}")
    print(f"VoTG(delta_be-0.1)   = {v_lo:+.6f}   (below break-even, should be < 0)")
    print(f"VoTG(delta_be+0.2)   = {v_hi:+.6f}   (above break-even, should be > 0)")
    assert v_lo < 0.0 < v_hi and 0.0 < d_be < 1.5
    print("GATE 3 PASS: VoTG<0 below delta_be, VoTG>0 above, delta_be finite")
    out.update(d_be=d_be, votg_curve=votg_curve)

    # ---- Figures ----------------------------------------------------------------
    hr("FIGURES")
    _fig_a(gd, g0)
    _fig_b(grid_kl)
    _fig_c(grid_kr)
    _fig_d(votg_curve, d_be)
    for name in ("fig_gap_a_voc", "fig_gap_b_G_vs_kappa", "fig_gap_c_surface",
                 "fig_gap_d_votg_breakeven"):
        print("wrote", os.path.join("figures", name + ".pdf"))

    hr("ALL GATES PASS")
    return out


def _voc_curve(rho, kappa, lam):
    cfg, base = gcfg(rho, kappa, lam, RUN_N)
    taus = np.round(np.arange(0.6, 1.8001, 0.025), 3)
    eu1 = posterior_eu(base, cfg, 1.0)
    voc = np.array([eu1 - posterior_eu(base, cfg, t) for t in taus])
    return taus, voc


def _fig_a(gd, g0):
    taus, voc_mis = _voc_curve(RHO0, KAP0, LAM0)
    _, voc_well = _voc_curve(RHO0, 0.0, LAM0)
    ymin = min(voc_mis.min(), voc_well.min())
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axhline(0.0, color="#404040", lw=0.8)
    ax.plot(taus, voc_well, lw=2, color="#9aa7b2", label="well-specified (kappa=0)")
    ax.plot(taus, voc_mis, lw=2, color="#1f4e79", label=f"misspecified (kappa={KAP0:.0f})")
    # the decision optimum: VoC dips BELOW zero -> the underconfident bar beats the calibrated one.
    ax.plot([gd.tau_star], [voc_mis.min()], "o", color="#2e8b57", ms=6, zorder=6)
    ax.axvline(gd.tau_stat, ls="--", color="#cc5500",
               label=f"tau_stat={gd.tau_stat:.3f} (coverage: shrink)")
    ax.axvline(gd.tau_star, ls=":", color="#2e8b57",
               label=f"tau*={gd.tau_star:.3f} (decision: widen)")
    yarrow = ymin - 0.0006
    ax.annotate("", xy=(gd.tau_star, yarrow), xytext=(gd.tau_stat, yarrow),
                arrowprops=dict(arrowstyle="<->", color="#000000", lw=1.3))
    ax.text(0.5 * (gd.tau_stat + gd.tau_star), yarrow - 0.0007,
            f"gap G={gd.gap:+.3f}", ha="center", fontsize=9)
    ax.set_ylim(ymin - 0.0018, 0.0075)
    ax.set_xlabel("reported-error-bar scale  tau  (1 = moment-matched)")
    ax.set_ylabel("value of calibration  VoC(tau)")
    ax.set_title("(a) Statistical vs decision calibration: the gap made visible")
    ax.legend(frameon=False, fontsize=8.5, loc="upper center", ncol=1)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_gap_a_voc.pdf"))
    plt.close(fig)


def _fig_b(grid_kl):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    colors = ["#9aa7b2", "#5b8cb0", "#1f4e79", "#0b2f52"]
    for j, lam in enumerate(LAMBDAS):
        ax.plot(KAPPAS, grid_kl[:, j], lw=2, marker="o", ms=3.5,
                color=colors[j % len(colors)], label=f"lambda={lam:.0f}")
    ax.axhline(0.0, color="#404040", lw=0.8)
    ax.set_xlabel("posterior skew  kappa")
    ax.set_ylabel("decision-calibration gap  G = tau* - tau_stat")
    ax.set_title("(b) The gap grows with skew and with cost asymmetry")
    ax.legend(frameon=False, fontsize=9, title="under:over cost")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_gap_b_G_vs_kappa.pdf"))
    plt.close(fig)


def _fig_c(grid_kr):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    K, R = np.meshgrid(KAPPAS_C, RHOS, indexing="ij")
    vmax = np.abs(grid_kr).max()
    cf = ax.contourf(K, R, grid_kr, levels=14, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    cs = ax.contour(K, R, grid_kr, levels=[0.0], colors="k", linewidths=1.0)
    ax.clabel(cs, fmt="G=0", fontsize=8)
    fig.colorbar(cf, ax=ax, label="gap  G")
    ax.set_xlabel("posterior skew  kappa")
    ax.set_ylabel("threshold proximity  rho  (distance to t2 / s)")
    ax.set_title("(c) Gap surface: largest near the threshold, strong skew")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_gap_c_surface.pdf"))
    plt.close(fig)


def _fig_d(votg_curve, d_be):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.axhline(0.0, color="#404040", lw=0.8)
    ax.plot(DELTAS, votg_curve, lw=2, marker="o", ms=3, color="#1f4e79", label="VoTG(delta)")
    ax.axvline(d_be, ls="--", color="#cc5500", label=f"break-even delta_be={d_be:.3f}")
    ax.fill_between(DELTAS, votg_curve, 0, where=votg_curve < 0, color="#cc5500", alpha=0.12)
    ax.set_xlabel("distribution shift  delta")
    ax.set_ylabel("value of the trust-gate  VoTG(delta)")
    ax.set_title("(d) Trust-gate break-even shift")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_gap_d_votg_breakeven.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    summary = main()
