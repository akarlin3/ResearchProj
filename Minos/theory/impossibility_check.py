"""Theorem 2(i) — machine verification of the hidden-channel impossibility (GATE 2(i)).

The impossibility (best label-free detection AUC = 1/2 for the hidden channel) is a
data-processing / sufficiency argument (``impossibility.md``). Its single load-bearing premise is
machine-checkable, and so is the conclusion it forces. This script verifies both, turning the former
"drafted, not machine-verified" status into a gated check:

  (i.a) PATHWISE INVARIANCE (the premise).  For the pure hidden channel (delta_obs = 0), the observable
        batch O = {mu} -- and therefore every label-free statistic M = f(O) -- is *bit-for-bit*
        identical across the whole delta_hid sweep, for every seed.  This is the exact-invariance
        P_O^stale = P_O^fresh on which the proof rests, checked across the full sweep (not one point).

  (i.b) FRESH-vs-STALE DETECTION AUC = 1/2 (the conclusion).  Because M(fresh, seed) = M(stale, seed)
        exactly, the monitor scores of fresh and stale batches are the same multiset, so every detector
        built from M has TPR = FPR at every threshold and the fresh-vs-stale AUC is exactly 1/2.  We
        compute it (Mann-Whitney U, ties at average ranks) and bootstrap over seeds: the CI is the
        degenerate [0.5, 0.5].  This is the finite-sample face of the population AUC = 1/2.

  (i.c) REGRET-DETECTION AUC ~ 1/2.  For completeness we reproduce the v3 framing (AUC of M against the
        oracle label {R > tol} on the hidden channel): also at chance, matching RESULTS_C.md (0.500).

Run:  ``.venv-theory/bin/python theory/impossibility_check.py``   (set MINOS_FAST=1 to shrink).
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "minos-core"))

from minos.config import MinosConfig, gaussian_latent_config       # noqa: E402
from minos.correction import fit_loss_calibration, stale_regret     # noqa: E402
from minos.gate import _auc                                         # noqa: E402
from minos.generative import make_population, realise_deploy        # noqa: E402
from minos.monitor import build_reference, detection_auc, monitor   # noqa: E402
from minos.seeding import make_rng                                  # noqa: E402

_FAST = os.environ.get("MINOS_FAST") == "1"
KAP, LAM, RHO = 3.0, 3.0, 0.5
N_CAL = 300_000 if _FAST else 2_000_000        # calibration set for tau_hat / reference
N_DET = 40_000 if _FAST else 120_000           # per-batch detection size
DET_SEEDS = list(range(8001, 8005 if _FAST else 8017))   # 4 (fast) / 16 (full) detection seeds
TOL = 0.02
N_BOOT = 2000


def hr(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def main():
    hr("CP — GATE 2(i): machine verification of the hidden-channel impossibility (Theorem 2(i))")
    print(f"FAST={_FAST}  N_CAL={N_CAL}  N_DET={N_DET}x{len(DET_SEEDS)} seeds  default cell "
          f"kappa={KAP} lambda={LAM} rho={RHO}")

    cfg = gaussian_latent_config(rho=RHO, kappa=KAP, lam=LAM, base=MinosConfig(n_voxels=N_CAL))
    base_cal = make_population(cfg, make_rng(cfg.seed))
    tau_hat = fit_loss_calibration(base_cal, cfg)
    ref = build_reference(base_cal, cfg, tau_hat)
    cfg_det = cfg.replace(n_voxels=N_DET)
    deltas = np.round(np.arange(0.0, 0.2401, 0.03), 3)

    # ------------------------------------------------------------------------------------------
    # (i.a) PATHWISE INVARIANCE — O and M are bit-identical across the whole delta_hid sweep.
    # ------------------------------------------------------------------------------------------
    hr("(i.a) Pathwise invariance: O={mu} and M=f(O) bit-identical across the delta_hid sweep")
    max_dmu = 0.0
    max_dM = 0.0
    M_fresh_by_seed = {}
    M_stale_by_seed = {}     # M at the largest delta_hid (the 'stale' regime), per seed
    for sd in DET_SEEDS:
        base = make_population(cfg_det, make_rng(sd))
        mu0, _ = realise_deploy(base, cfg_det, delta_obs=0.0, delta_hid=0.0)
        M0 = monitor(mu0, cfg_det, ref)
        M_fresh_by_seed[sd] = M0
        for d in deltas[1:]:
            mu_d, _ = realise_deploy(base, cfg_det, delta_obs=0.0, delta_hid=float(d))
            Md = monitor(mu_d, cfg_det, ref)
            max_dmu = max(max_dmu, float(np.max(np.abs(mu_d - mu0))))
            max_dM = max(max_dM, abs(Md - M0))
        # 'stale' = largest shift in the sweep
        mu_s, _ = realise_deploy(base, cfg_det, delta_obs=0.0, delta_hid=float(deltas[-1]))
        M_stale_by_seed[sd] = monitor(mu_s, cfg_det, ref)
    invariant = (max_dmu == 0.0) and (max_dM == 0.0)
    print(f"  max_seed,delta |mu(delta_hid) - mu(0)|  = {max_dmu:.3e}   (must be exactly 0)")
    print(f"  max_seed,delta |M(delta_hid)  - M(0)|   = {max_dM:.3e}   (must be exactly 0)")
    print(f"  => observable law (hence every M=f(O)) is exactly invariant to delta_hid: {invariant}")

    # ------------------------------------------------------------------------------------------
    # (i.b) FRESH-vs-STALE DETECTION AUC = 1/2 (the conclusion), + bootstrap CI over seeds.
    #       Scores = M; labels = {batch is stale}. Because M(fresh,s)==M(stale,s), the two label
    #       groups share the same value multiset -> AUC = 1/2 exactly, every bootstrap resample.
    # ------------------------------------------------------------------------------------------
    hr("(i.b) Fresh-vs-stale detection AUC (any label-free detector) = 1/2, with bootstrap CI")
    seeds = np.array(DET_SEEDS)
    Mf = np.array([M_fresh_by_seed[s] for s in DET_SEEDS])
    Ms = np.array([M_stale_by_seed[s] for s in DET_SEEDS])
    scores = np.concatenate([Mf, Ms])
    labels = np.concatenate([np.zeros(len(Mf), bool), np.ones(len(Ms), bool)])   # stale = positive
    auc_fs = _auc(scores, labels)
    # bootstrap over seeds (resample paired fresh/stale by seed)
    boot = []
    rng = make_rng(99991)
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(seeds), len(seeds))
        sc = np.concatenate([Mf[idx], Ms[idx]])
        lb = np.concatenate([np.zeros(len(idx), bool), np.ones(len(idx), bool)])
        boot.append(_auc(sc, lb))
    lo, hi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))
    print(f"  fresh-vs-stale AUC (Mann-Whitney, ties=avg-rank) = {auc_fs:.6f}")
    print(f"  bootstrap 95% CI over {len(seeds)} seeds            = [{lo:.6f}, {hi:.6f}]")
    print(f"  => every label-free detector has TPR=FPR; ROC is the diagonal; AUC = 1/2 exactly.")

    # ------------------------------------------------------------------------------------------
    # (i.c) REGRET-DETECTION AUC ~ 1/2 (the v3 framing: M vs the oracle label {R>tol}).
    # ------------------------------------------------------------------------------------------
    hr("(i.c) Regret-detection AUC on the hidden channel (M vs {R>tol}) -- v3 framing")
    Ms_all, Rs_all = [], []
    for sd in DET_SEEDS:
        base = make_population(cfg_det, make_rng(sd))
        for d in deltas:
            mu, _ = realise_deploy(base, cfg_det, delta_obs=0.0, delta_hid=float(d))
            Ms_all.append(monitor(mu, cfg_det, ref))
            Rs_all.append(stale_regret(base, cfg_det, tau_hat, delta_obs=0.0, delta_hid=float(d)))
    auc_reg = detection_auc(np.array(Ms_all), np.array(Rs_all), TOL)
    print(f"  hidden-channel regret-detection AUC = {auc_reg:.4f}  (v3 RESULTS_C.md: 0.500, at chance)")

    # ------------------------------------------------------------------------------------------
    # GATE 2(i) verdict.
    # ------------------------------------------------------------------------------------------
    hr("GATE 2(i) — verdict")
    auc_ok = abs(auc_fs - 0.5) < 1e-9 and abs(lo - 0.5) < 1e-9 and abs(hi - 0.5) < 1e-9
    reg_ok = 0.40 <= auc_reg <= 0.60
    ok = invariant and auc_ok and reg_ok
    print(f"  (i.a) pathwise invariance exact         : {invariant}")
    print(f"  (i.b) fresh-vs-stale AUC = 1/2 (CI deg.) : {auc_ok}")
    print(f"  (i.c) regret-detection AUC ~ 1/2         : {reg_ok}  ({auc_reg:.4f})")
    print()
    if ok:
        print("GATE 2(i) PASS: the hidden shift leaves O (hence every M=f(O)) exactly invariant, so the")
        print("  fresh-vs-stale detection AUC is exactly 1/2 -- the impossibility of Theorem 2(i), now")
        print("  machine-verified. Conditional only on the definitional premise that 'hidden' = the")
        print("  component leaving the observable law invariant (true by construction in the Minos split).")
    else:
        print("GATE 2(i) FAIL — the structural premise or its consequence did not hold; investigate.")
    assert ok, "GATE 2(i): hidden-channel impossibility machine-check failed"
    return dict(invariant=invariant, max_dM=max_dM, auc_fs=auc_fs, ci=(lo, hi), auc_reg=auc_reg)


if __name__ == "__main__":
    main()
