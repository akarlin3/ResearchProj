"""Minos-Core v3 experiment driver — a deployment-validity monitor for a stale loss-calibration
correction under shift (the v3 contribution).

Seeded, config-driven. Reproduces the four checkpoint gates from a clean seed, prints every number
RESULTS_C.md cites, and writes the four publication figures as vector PDFs:

  (a) motivation: the v2 decision-calibration gap G vs skew kappa (why a loss-calibration correction)
  (b) stale-correction regret R vs shift, both knobs (observable and hidden) overlaid
  (c) the honesty figure: monitor M vs regret R scatter + detection ROC, observable vs hidden
  (d) regret ladder {stale, gated, oracle} across the observable sweep, operating point marked

Run from the project root:  ``python experiments/run_c.py``
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from minos.calibration import gap  # noqa: E402
from minos.config import MinosConfig, gaussian_latent_config  # noqa: E402
from minos.correction import (  # noqa: E402
    deploy_expected_utility,
    fit_loss_calibration,
    oracle_deploy_scale,
    stale_regret,
)
from minos.generative import make_population, realise_deploy  # noqa: E402
from minos.monitor import (  # noqa: E402
    build_reference,
    calibrate_threshold,
    detection_auc,
    gated_recovery_actions,
    monitor,
)
from minos.seeding import make_rng  # noqa: E402
from minos.voi import realised_utility  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(os.path.dirname(HERE), "figures")

# Estimation sizes (CRN; one base draw per config). The MINOS_FAST env knob shrinks every size for a
# quick end-to-end gate check; the committed RESULTS_C numbers come from the default (full) sizes.
_FAST = os.environ.get("MINOS_FAST") == "1"
CAL_N = 200_000 if _FAST else 2_000_000      # labeled calibration set D_cal (headline baseline fit)
DEP_N = 200_000 if _FAST else 1_000_000      # deployment headline regret curves
DETECT_N = 60_000 if _FAST else 150_000      # per-batch deployment size for the detection ROC
DETECT_SEEDS = list(range(7001, 7007 if _FAST else 7017))   # independent deployment seeds for the ROC
NULL_SEEDS = 20 if _FAST else 80         # zero-shift null batches for the m* threshold
ALPHA = 0.05             # target zero-shift false-alarm rate
TOL = 0.02               # regret tolerance for the detection label {R > tol} (documented choice)
DEP_SEED = 20240517 + 777

# Calibration cell = the v2 default misspecification (skewed, asymmetric, near the threshold).
KAP_CAL, LAM_CAL, RHO_CAL = 3.0, 3.0, 0.5
# Shift sweeps (matched grids so observable and hidden induce the same regret ladder).
DELTAS = np.round(np.arange(0.0, 0.2401, 0.03), 3)
# Motivation (a): the v2 gap over skew, at several cost asymmetries.
KAPPAS = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0]
LAMBDAS = [1.0, 2.0, 3.0, 4.0]
GAP_N = 200_000 if _FAST else 1_000_000

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


def cal_config(n=CAL_N):
    return gaussian_latent_config(rho=RHO_CAL, kappa=KAP_CAL, lam=LAM_CAL,
                                  base=MinosConfig(n_voxels=n))


def _roc_curve(scores, labels):
    """ROC (fpr, tpr) by sweeping the score threshold; labels is a bool array."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, bool)
    order = np.argsort(-scores, kind="mergesort")
    s = labels[order]
    tp = np.cumsum(s)
    fp = np.cumsum(~s)
    npos, nneg = max(int(labels.sum()), 1), max(int((~labels).sum()), 1)
    tpr = np.concatenate([[0.0], tp / npos])
    fpr = np.concatenate([[0.0], fp / nneg])
    return fpr, tpr


def detection_set(cfg, ref, tau_hat, channel):
    """(M, R) over DETECT_SEEDS x DELTAS for one shift channel — M label-free, R the sim oracle."""
    cfg_b = cfg.replace(n_voxels=DETECT_N)
    Ms, Rs = [], []
    for sd in DETECT_SEEDS:
        base = make_population(cfg_b, make_rng(sd))
        for d in DELTAS:
            do, dh = (d, 0.0) if channel == "obs" else (0.0, d)
            mu, _ = realise_deploy(base, cfg_b, delta_obs=do, delta_hid=dh)
            Ms.append(monitor(mu, cfg_b, ref))
            Rs.append(stale_regret(base, cfg_b, tau_hat, delta_obs=do, delta_hid=dh))
    return np.array(Ms), np.array(Rs)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    out = {}

    hr("CONFIG")
    cfg = cal_config()
    print(f"seed={cfg.seed}  CAL_N={CAL_N}  DEP_N={DEP_N}  DETECT_N={DETECT_N} "
          f"x {len(DETECT_SEEDS)} seeds")
    print(f"calibration cell: kappa={KAP_CAL} lambda={LAM_CAL} rho={RHO_CAL} "
          f"(posterior-centric, s={cfg.s}); single active threshold t2={cfg.t2}")
    print(f"shift split: delta_obs biases reported mu DOWN (observable); "
          f"delta_hid biases truth theta UP (hidden); gain beta*s={cfg.beta*cfg.s}")
    print(f"deltas={list(DELTAS)}  alpha={ALPHA}  regret tol={TOL}")

    # ---- GATE 1: loss-calibration baseline + stale-regret -----------------------
    hr("GATE 1 — loss-calibration baseline (cited) + stale-correction regret")
    base_cal = make_population(cfg, make_rng(cfg.seed))
    tau_hat = fit_loss_calibration(base_cal, cfg)
    print(f"loss-calibration baseline tau_hat_cal (= v2 decision-calibrated scale on D_cal) "
          f"= {tau_hat:.4f}")

    # Consistency: the zero-shift EU optimum is shallow, so a single-seed point estimate jitters at
    # the grid scale. The robust statement is the textbook one — the estimator's across-seed VARIANCE
    # shrinks with N — together with agreement to the deployment oracle at large N.
    print("\nconsistency: tau_hat_cal across seeds as N_cal grows (zero shift):")
    print("   N_cal     mean(tau_hat)   std(tau_hat)")
    cseeds = (11, 12, 13, 14, 15)
    cons_ns = (50_000, 150_000, 500_000) if _FAST else (100_000, 500_000, 2_000_000)
    cons = []
    for n in cons_ns:
        cfg_n = cal_config(n)
        ths = np.array([fit_loss_calibration(make_population(cfg_n, make_rng(sd)), cfg_n)
                        for sd in cseeds])
        cons.append((n, float(ths.mean()), float(ths.std())))
        print(f"  {n:>9}   {ths.mean():.4f}          {ths.std():.4f}")
    cfg_big = cal_config(cons_ns[-1])
    tau_or_big, _ = oracle_deploy_scale(make_population(cfg_big, make_rng(999)), cfg_big)
    print(f"tau*_oracle@dep (N=2e6, independent draw) = {tau_or_big:.4f}  "
          f"vs mean(tau_hat)@2e6 = {cons[-1][1]:.4f}")
    assert cons[-1][2] < cons[0][2], "estimator variance must shrink with N (consistency)"
    assert abs(cons[-1][1] - tau_or_big) < 0.06, "baseline must agree with the oracle at large N"

    # v2 reproduction: the motivating gap is intact on this cell.
    g_cell = gap(base_cal, cfg)
    print(f"\nv2 gap on the calibration cell: tau_stat={g_cell.tau_stat:.4f} "
          f"tau*={g_cell.tau_star:.4f}  G={g_cell.gap:+.4f} (reproduces the v2 headline sign)")
    assert g_cell.gap > 0.03 and g_cell.tau_star > g_cell.tau_stat

    # regret sweeps on the headline deployment sample
    cfg_dep = cfg.replace(n_voxels=DEP_N)
    base_dep = make_population(cfg_dep, make_rng(DEP_SEED))
    R_obs = np.array([stale_regret(base_dep, cfg_dep, tau_hat, delta_obs=d, delta_hid=0.0)
                      for d in DELTAS])
    R_hid = np.array([stale_regret(base_dep, cfg_dep, tau_hat, delta_obs=0.0, delta_hid=d)
                      for d in DELTAS])
    print(f"\nR(0,0) = {R_obs[0]:.6f}  (~0 at the regret floor)")
    print("R vs delta_obs:", "  ".join(f"{d:.2f}:{r:.4f}" for d, r in zip(DELTAS, R_obs)))
    print("R vs delta_hid:", "  ".join(f"{d:.2f}:{r:.4f}" for d, r in zip(DELTAS, R_hid)))
    assert R_obs[0] < 3e-3 and R_hid[0] < 3e-3, "regret must vanish at zero shift"
    # R rises with each knob: a clear net increase, monotone up to MC slack near the flat floor.
    assert np.all(np.diff(R_obs) > -1e-3) and (R_obs[-1] - R_obs[0]) > 0.05
    assert np.all(np.diff(R_hid) > -1e-3) and (R_hid[-1] - R_hid[0]) > 0.05
    print("GATE 1 PASS: baseline consistent (variance shrinks, agrees with oracle); R(0)~0; "
          "R grows with each shift knob")
    out.update(tau_hat=tau_hat, cons=cons, g_cell=g_cell, R_obs=R_obs, R_hid=R_hid)

    # ---- GATE 2: validity monitor M ---------------------------------------------
    hr("GATE 2 — validity monitor M (label-free): corr, detection AUC, honest hidden limit")
    ref = build_reference(base_cal, cfg, tau_hat)

    # no-label code path: scrambling the labels cannot move M.
    base_chk = make_population(cfg.replace(n_voxels=DETECT_N), make_rng(DEP_SEED))
    mu_chk, _ = realise_deploy(base_chk, cfg.replace(n_voxels=DETECT_N), delta_obs=0.1)
    m_before = monitor(mu_chk, cfg.replace(n_voxels=DETECT_N), ref)
    base_chk.theta[:] = make_rng(1).standard_normal(base_chk.theta.shape) * 99.0
    mu_chk2, _ = realise_deploy(base_chk, cfg.replace(n_voxels=DETECT_N), delta_obs=0.1)
    m_after = monitor(mu_chk2, cfg.replace(n_voxels=DETECT_N), ref)
    print(f"no-label path: M(scrambled theta)==M(real theta) -> {m_before == m_after} "
          f"(M={m_before:.6f})")
    assert m_before == m_after, "monitor must not depend on labels"

    M_obs, Rd_obs = detection_set(cfg, ref, tau_hat, "obs")
    M_hid, Rd_hid = detection_set(cfg, ref, tau_hat, "hid")
    corr_obs = float(np.corrcoef(M_obs, Rd_obs)[0, 1])
    auc_obs = detection_auc(M_obs, Rd_obs, TOL)
    auc_hid = detection_auc(M_hid, Rd_hid, TOL)
    pos_obs = int((Rd_obs > TOL).sum())
    pos_hid = int((Rd_hid > TOL).sum())
    print(f"under delta_obs: corr(M,R) = {corr_obs:+.3f}   detection AUC{{R>tol}} = {auc_obs:.3f}"
          f"   (positives {pos_obs}/{len(Rd_obs)})")
    print(f"under delta_hid: detection AUC{{R>tol}} = {auc_hid:.3f}   "
          f"(positives {pos_hid}/{len(Rd_hid)}) <- AT CHANCE: the honest, documented limit")

    m_star = calibrate_threshold(cfg, ref, alpha=ALPHA, n_seeds=NULL_SEEDS, n_batch=DETECT_N)
    cfg_b = cfg.replace(n_voxels=DETECT_N)
    fa = 0
    for k in range(NULL_SEEDS):
        base = make_population(cfg_b, make_rng(40_000 + k))   # disjoint from m*'s null seeds
        mu, _ = realise_deploy(base, cfg_b)
        fa += monitor(mu, cfg_b, ref) > m_star
    fpr0 = fa / NULL_SEEDS
    print(f"m* (alpha={ALPHA}) = {m_star:.5f}   zero-shift false-alarm rate = {fpr0:.3f}")
    assert corr_obs > 0.8 and auc_obs > 0.85
    # The hidden-shift AUC ~ 0.5 is a DOCUMENTED RESULT, not a red gate (DESIGN_C §5). We only sanity
    # bound it well away from a competent detector — it must NOT look like the observable channel.
    assert 0.30 <= auc_hid <= 0.70 and auc_hid < auc_obs - 0.2, "hidden staleness must be ~at chance"
    assert fpr0 <= 0.15
    print("GATE 2 PASS: M tracks R and detects under delta_obs; AT CHANCE under delta_hid "
          "(documented); zero-shift false alarms controlled")
    out.update(m_star=m_star, M_obs=M_obs, Rd_obs=Rd_obs, M_hid=M_hid, Rd_hid=Rd_hid,
               corr_obs=corr_obs, auc_obs=auc_obs, auc_hid=auc_hid, fpr0=fpr0)

    # ---- GATE 3: gated recovery + regret ladder ---------------------------------
    hr("GATE 3 — gated recovery + monitor calibration (regret ladder)")

    def ladder(delta_obs, delta_hid):
        mu, theta = realise_deploy(base_dep, cfg_dep, delta_obs=delta_obs, delta_hid=delta_hid)
        tau_or, eu_or = oracle_deploy_scale(base_dep, cfg_dep,
                                            delta_obs=delta_obs, delta_hid=delta_hid)
        eu_stale = deploy_expected_utility(base_dep, cfg_dep, tau_hat,
                                           delta_obs=delta_obs, delta_hid=delta_hid)
        a_gate = gated_recovery_actions(mu, cfg_dep, ref, m_star)
        eu_gate = float(np.mean(realised_utility(a_gate, theta, cfg_dep)))
        fires = monitor(mu, cfg_dep, ref) > m_star
        return eu_or, eu_stale, eu_gate, bool(fires)

    print("delta_obs:  EU_oracle   EU_gated   EU_stale  | R_stale  R_gated  fires")
    ladder_obs = []
    for d in DELTAS:
        eo, es, eg, fr = ladder(d, 0.0)
        ladder_obs.append((d, eo, es, eg, fr))
        print(f"   {d:.2f}:   {eo:+.5f}  {eg:+.5f}  {es:+.5f} |  {eo-es:+.5f}  {eo-eg:+.5f}  {fr}")
    # operating point: smallest observable shift at which the monitor fires
    fired = [d for d, *_rest, fr in ladder_obs if fr]
    op_delta = min(fired) if fired else None
    print(f"monitor operating point (first delta_obs that fires): {op_delta}")

    # zero-shift: no false-alarm harm; under delta_hid: no spurious help (leakage check)
    eo0, es0, eg0, fr0 = ladder(0.0, 0.0)
    print(f"\nzero shift:   R_stale={eo0-es0:+.5f}  R_gated={eo0-eg0:+.5f}  fires={fr0} "
          f"(gate silent -> no harm)")
    print("under delta_hid (gate is blind -> must NOT spuriously help):")
    ladder_hid = []
    for d in DELTAS:
        eo, es, eg, fr = ladder(0.0, d)
        ladder_hid.append((d, eo, es, eg, fr))
        print(f"   {d:.2f}:   R_stale={eo-es:+.5f}  R_gated={eo-eg:+.5f}  fires={fr}")

    # checks: recovery under obs at the strongest shift; silence elsewhere
    eo_m, es_m, eg_m, _ = ladder(DELTAS[-1], 0.0)
    assert (eo_m - eg_m) < (eo_m - es_m) - 0.02, "gated must beat stale under observable shift"
    assert abs((eo0 - eg0) - (eo0 - es0)) < 1e-6, "no false-alarm harm at zero shift"
    for d, eo, es, eg, fr in ladder_hid:
        assert abs((eo - eg) - (eo - es)) < 1e-6, "gate must not spuriously help under delta_hid"
    print("GATE 3 PASS: gated < stale under delta_obs; gated == stale at zero shift and under "
          "delta_hid (no spurious help)")
    out.update(ladder_obs=ladder_obs, ladder_hid=ladder_hid, op_delta=op_delta)

    # ---- FIGURES ----------------------------------------------------------------
    hr("FIGURES")
    _fig_a()
    _fig_b(R_obs, R_hid)
    _fig_c(M_obs, Rd_obs, M_hid, Rd_hid, m_star, auc_obs, auc_hid)
    _fig_d(ladder_obs, op_delta)
    for name in ("fig_c_a_motivation_gap", "fig_c_b_regret_vs_shift",
                 "fig_c_c_monitor_honesty", "fig_c_d_regret_ladder"):
        print("wrote", os.path.join("figures", name + ".pdf"))

    hr("ALL GATES PASS")
    return out


# --------------------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------------------
def _fig_a():
    """(a) motivation: the v2 decision-calibration gap G vs skew kappa at several lambda."""
    grid = np.zeros((len(KAPPAS), len(LAMBDAS)))
    for i, k in enumerate(KAPPAS):
        for j, lam in enumerate(LAMBDAS):
            cfg = gaussian_latent_config(rho=RHO_CAL, kappa=k, lam=lam,
                                         base=MinosConfig(n_voxels=GAP_N))
            grid[i, j] = gap(make_population(cfg, make_rng(cfg.seed)), cfg).gap
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    colors = ["#9aa7b2", "#5b8cb0", "#1f4e79", "#0b2f52"]
    for j, lam in enumerate(LAMBDAS):
        ax.plot(KAPPAS, grid[:, j], lw=2, marker="o", ms=3.5, color=colors[j],
                label=f"lambda={lam:.0f}")
    ax.axhline(0.0, color="#404040", lw=0.8)
    ax.set_xlabel("posterior skew  kappa")
    ax.set_ylabel("decision-calibration gap  G = tau* - tau_stat")
    ax.set_title("(a) Motivation: the v2 gap a loss-calibration correction repairs")
    ax.legend(frameon=False, fontsize=9, title="under:over cost")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_c_a_motivation_gap.pdf"))
    plt.close(fig)


def _fig_b(R_obs, R_hid):
    """(b) stale-correction regret R vs shift, both knobs overlaid."""
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(DELTAS, R_obs, lw=2, marker="o", ms=4, color="#1f4e79",
            label="delta_obs (observable)")
    ax.plot(DELTAS, R_hid, lw=2, marker="s", ms=4, color="#cc5500",
            label="delta_hid (hidden)")
    ax.set_xlabel("deployment shift  delta")
    ax.set_ylabel("stale-correction regret  R")
    ax.set_title("(b) The loss-calibration correction goes stale under shift")
    ax.legend(frameon=False, fontsize=9.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_c_b_regret_vs_shift.pdf"))
    plt.close(fig)


def _fig_c(M_obs, R_obs, M_hid, R_hid, m_star, auc_obs, auc_hid):
    """(c) the honesty figure: M vs R scatter + detection ROC, observable vs hidden."""
    fig, (axs, axr) = plt.subplots(1, 2, figsize=(10.4, 4.3))
    axs.scatter(R_obs, M_obs, s=14, color="#1f4e79", alpha=0.7, label="delta_obs (observable)")
    axs.scatter(R_hid, M_hid, s=14, color="#cc5500", alpha=0.7, label="delta_hid (hidden)")
    axs.axhline(m_star, ls="--", color="#404040", lw=1.0, label=f"m* (alpha quantile)")
    axs.set_xlabel("stale-correction regret  R  (oracle, sim only)")
    axs.set_ylabel("validity monitor  M  (label-free)")
    axs.set_title("(c) Monitor vs regret: observable rises, hidden is flat")
    axs.legend(frameon=False, fontsize=8.5, loc="upper left")

    fpr_o, tpr_o = _roc_curve(M_obs, R_obs > TOL)
    fpr_h, tpr_h = _roc_curve(M_hid, R_hid > TOL)
    axr.plot(fpr_o, tpr_o, lw=2, color="#1f4e79", label=f"delta_obs  AUC={auc_obs:.2f}")
    axr.plot(fpr_h, tpr_h, lw=2, color="#cc5500", label=f"delta_hid  AUC={auc_hid:.2f}")
    axr.plot([0, 1], [0, 1], ls=":", color="#404040", lw=1.0, label="chance")
    axr.set_xlabel("false-positive rate")
    axr.set_ylabel("true-positive rate")
    axr.set_title("Detection ROC for {R > tol}")
    axr.legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_c_c_monitor_honesty.pdf"))
    plt.close(fig)


def _fig_d(ladder_obs, op_delta):
    """(d) regret ladder {stale, gated, oracle} across delta_obs, operating point marked."""
    d = np.array([r[0] for r in ladder_obs])
    eu_or = np.array([r[1] for r in ladder_obs])
    eu_st = np.array([r[2] for r in ladder_obs])
    eu_gt = np.array([r[3] for r in ladder_obs])
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    ax.plot(d, eu_or, lw=2, marker="o", ms=3.5, color="#2e8b57", label="oracle (deployment-optimal)")
    ax.plot(d, eu_gt, lw=2, marker="^", ms=3.5, color="#1f4e79", label="gated recovery (M-triggered)")
    ax.plot(d, eu_st, lw=2, marker="s", ms=3.5, color="#cc5500", label="stale correction")
    if op_delta is not None:
        ax.axvline(op_delta, ls="--", color="#404040", lw=1.0,
                   label=f"monitor fires at delta_obs={op_delta:.2f}")
    ax.set_xlabel("observable shift  delta_obs")
    ax.set_ylabel("expected decision utility  EU")
    ax.set_title("(d) Gated recovery claws back the stale-correction regret")
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_c_d_regret_ladder.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    summary = main()
