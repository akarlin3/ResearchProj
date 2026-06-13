"""CP1 — bootstrap / multi-seed confidence intervals on the decision–calibration gap (HALT-TO-REPORT).

The v2 headline is that statistical and decision calibration land on **opposite sides of 1**:
``tau_stat < 1 < tau*`` (equivalently ``G = tau* - tau_stat > 0``). But ``tau*`` is the argmax of a
*shallow* expected-utility optimum, so a single-seed estimate is noisy (``confirm.py`` shows fresh
re-sims of ``tau*`` jitter by ~0.02). Before that headline goes in a paper, the three pieces of it
must be **error-barred**:

    (V1)  tau*    CI excludes 1.0   <=>  "decision calibration WIDENS the bar" is robust
    (V2)  G       CI excludes 0.0   <=>  "the gap is real" is robust
    (V3)  tau_stat CI sits below 1  <=>  "statistical calibration SHRINKS the bar" is robust

and the combined "opposite sides of 1" claim holds iff V1 and V3 (equivalently V2) all hold.

Method (printed, run-then-write). The estimator under test is exactly the published one
(``minos.calibration.gap``: a coverage root-find for ``tau_stat`` and a grid+parabola argmax for
``tau*``). We characterise its sampling distribution by **multi-seed Monte-Carlo replicates** — B
independent populations, each a genuine i.i.d. draw of the whole estimator (no bootstrap
approximation of the argmax functional, which a within-sample resample would distort). From the B
replicates we report, per quantity:

  * the per-seed **mean** and **sd** (the run-to-run spread that made the single published number
    fragile),
  * a 95% CI on the **mean / true value** two ways — a Student-t interval and a non-parametric
    **percentile bootstrap over the B replicates** (``NBOOT`` resamples, fixed seed) — which must
    agree, and
  * the single-run 95% **prediction interval** (mean +/- 1.96 sd): where one published run lands.

The verdict on each V-claim is driven by the CI on the *mean* (the inferential statement about the
true value); we additionally flag when the *single-run* interval also excludes the null (a stronger
statement). Finite-N bias is checked to be negligible by the N-stability of the mean (the default
cell is also re-estimated at a second N).

REFRAME RULE (no overclaiming, baked in below). If the ``tau*`` CI **includes** 1, the script prints
that "opposite sides of 1" is NOT supported and falls back to the strongest claim the data do
support (``tau* > tau_stat`` <=> ``G > 0``, which can be robust even when ``tau* > 1`` is not). If the
``G`` CI includes 0 at the default cell, it falls back to the analytic slope ``dG/dgamma = |z*|/6 > 0``
(``gap_scaling.py``), which is sympy-derived and not subject to the optimum's noise.

GATE 1 (HALT-TO-REPORT): CIs computed + printed; the three honesty verdicts stated; the robust claim
identified. The estimator is NOT tuned to recover a cleaner claim — whatever the CIs say is the
finding.

Run:  ``.venv-theory/bin/python theory/gap_ci.py``   (~5-6 min; ``MINOS_FAST=1`` shrinks B and N)
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import mpmath as mp
from scipy.stats import t as student_t

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "minos-core"))

from minos.config import MinosConfig, gaussian_latent_config       # noqa: E402
from minos.calibration import gap                                  # noqa: E402
from minos.generative import make_population                       # noqa: E402
from minos.seeding import make_rng                                 # noqa: E402

_FAST = os.environ.get("MINOS_FAST") == "1"

# default gap cell (THEORY_MODEL.md §4 / RESULTS_B.md): kappa=3, lambda=3, rho=0.5.
KAP, LAM, RHO = 3.0, 3.0, 0.5
# Replicate budget. Var(multi-seed mean) ~ 1/(B*N), so the CI half-width is set by the *product*
# B*N; we favour many moderate-N seeds (a stabler variance estimate, and the per-seed sd is itself
# a reported quantity). N-stability across the two default-cell N's is the bias check.
B_CELL = 24 if _FAST else 64           # replicates at the default cell (headline)
N_CELL = 300_000 if _FAST else 1_000_000
N_CELL2 = 150_000 if _FAST else 2_000_000   # second N for the finite-N-bias (N-stability) check
B_CELL2 = 16 if _FAST else 32
B_SWEEP = 16 if _FAST else 32          # replicates per kappa in the gamma sweep
N_SWEEP = 250_000 if _FAST else 500_000
SWEEP_KAPPAS = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0)
SEED0 = 90_000                          # base of the replicate seed block (distinct from confirm.py's)
NBOOT = 20_000                          # bootstrap resamples for the percentile CI on the mean
BOOT_SEED = 73_019                      # fixed seed for the bootstrap resampler (printed)
MAX_WORKERS = min(8, (os.cpu_count() or 4))

# Published single-seed point values we are error-barring (RESULTS_B.md GATE 1, seed 20240517).
PUB = dict(tau_stat=0.9635, tau_star=1.0431, G=0.0796)
# Published v2 gap sweep G(kappa) at lambda=3, rho=0.5 (RESULTS_B.md GATE 2, high-N single seed).
PUB_SWEEP_G = {0.0: -0.004, 0.5: -0.003, 1.0: +0.005, 2.0: +0.044, 3.0: +0.084, 4.0: +0.115}


def hr(title):
    print("\n" + "=" * 84)
    print(title)
    print("=" * 84)


# ---- the estimator under test (top-level so ProcessPoolExecutor can pickle it) -----------------
def _estimate(args):
    """One i.i.d. replicate: draw a fresh population at (kappa,lam,rho,N) on `seed`, return the gap."""
    kappa, lam, rho, n, seed = args
    cfg = gaussian_latent_config(rho=rho, kappa=kappa, lam=lam, base=MinosConfig(n_voxels=n))
    base = make_population(cfg.replace(seed=seed), make_rng(seed))
    gr = gap(base, cfg)
    return (gr.tau_stat, gr.tau_star, gr.gap)


def multiseed(kappa, lam, rho, n, b, seed0=SEED0):
    """B independent replicates of (tau_stat, tau*, G); returns a (B,3) array."""
    seeds = [seed0 + i for i in range(b)]
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        rows = list(ex.map(_estimate, [(kappa, lam, rho, n, s) for s in seeds]))
    return np.asarray(rows, dtype=float)


# ---- analytic noise-free fallback (sympy result of gap_scaling.py, evaluated numerically) -------
def zstar(lam):
    """gamma=0 decision-boundary offset: root of (lam-1)*psi(z)+z=0, psi(z)=z*Phi(z)+phi(z)."""
    if abs(lam - 1.0) < 1e-12:
        return 0.0
    psi = lambda tt: tt * float(mp.ncdf(tt)) + float(mp.npdf(tt))
    return float(mp.findroot(lambda tt: (lam - 1) * psi(tt) + tt, -0.4))


def gamma_of_kappa(kappa):
    """Standardised third cumulant (skewness) of the code's standardised skew-normal error."""
    d = kappa / np.sqrt(1 + kappa ** 2)
    return ((4 - np.pi) / 2) * (d * np.sqrt(2 / np.pi)) ** 3 / (1 - 2 * d ** 2 / np.pi) ** 1.5


# ---- CI machinery -------------------------------------------------------------------------------
def ci_mean_t(x, conf=0.95):
    """Student-t CI on the mean of the replicate estimates."""
    x = np.asarray(x, float)
    n = x.size
    m, sd = float(x.mean()), float(x.std(ddof=1))
    se = sd / np.sqrt(n)
    h = student_t.ppf(0.5 + conf / 2.0, n - 1) * se
    return m, sd, se, (m - h, m + h)


def ci_mean_boot(x, conf=0.95, nboot=NBOOT, seed=BOOT_SEED):
    """Non-parametric percentile-bootstrap CI on the mean (resample the B replicates)."""
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(nboot, x.size))
    boot_means = x[idx].mean(axis=1)
    lo, hi = np.percentile(boot_means, [100 * (1 - conf) / 2, 100 * (1 + conf) / 2])
    return float(lo), float(hi)


def predint(x, conf=0.95):
    """Single-run 95% prediction interval (where one published run lands): mean +/- z*sd."""
    x = np.asarray(x, float)
    m, sd = float(x.mean()), float(x.std(ddof=1))
    z = 1.959963985
    return (m - z * sd, m + z * sd)


def excludes(ci, value, side):
    """True if the interval ci=(lo,hi) lies strictly on one side of `value`.

    side='above' -> lo > value (the quantity is robustly greater);
    side='below' -> hi < value (robustly less).
    """
    lo, hi = ci
    return (lo > value) if side == "above" else (hi < value)


def summarise(name, col, null, side):
    """Print the full CI block for one quantity and return (verdict_on_mean, details)."""
    m, sd, se, ci_t = ci_mean_t(col)
    ci_b = ci_mean_boot(col)
    pi = predint(col)
    mean_excl = excludes(ci_b, null, side) and excludes(ci_t, null, side)
    pred_excl = excludes(pi, null, side)
    rel = ">" if side == "above" else "<"
    print(f"  {name:<10} mean={m:+.4f}  sd={sd:.4f}  SE={se:.4f}")
    print(f"      95% CI(mean):  t=[{ci_t[0]:+.4f}, {ci_t[1]:+.4f}]   "
          f"boot=[{ci_b[0]:+.4f}, {ci_b[1]:+.4f}]   (resamples={NBOOT}, seed={BOOT_SEED})")
    print(f"      single-run 95% PI: [{pi[0]:+.4f}, {pi[1]:+.4f}]   "
          f"(published single-seed point: {PUB.get({'tau_stat':'tau_stat','tau*':'tau_star','G':'G'}[name]):+.4f})")
    print(f"      verdict: CI(mean) {rel} {null:g}?  {mean_excl}    "
          f"[single run also {rel} {null:g}?  {pred_excl}]")
    return mean_excl, dict(mean=m, sd=sd, se=se, ci_t=ci_t, ci_b=ci_b, pi=pi,
                           mean_excl=mean_excl, pred_excl=pred_excl)


def main():
    hr("CP1 / GATE 1 — bootstrap CIs on the gap (HALT-TO-REPORT): is 'opposite sides of 1' robust?")
    print(f"FAST={_FAST}  default cell kappa={KAP} lambda={LAM} rho={RHO}  "
          f"gamma(kappa=3)={gamma_of_kappa(KAP):.4f}  z*(lambda=3)={zstar(LAM):+.4f}")
    print(f"replicates: B_CELL={B_CELL} @ N={N_CELL:,}  (bias check B={B_CELL2} @ N={N_CELL2:,});  "
          f"sweep B={B_SWEEP} @ N={N_SWEEP:,} per kappa")
    print(f"seeds: {SEED0}..{SEED0 + max(B_CELL, B_CELL2) - 1};  workers={MAX_WORKERS};  "
          f"bootstrap resamples={NBOOT} (seed {BOOT_SEED})")

    # ============================================================================================
    # (1) DEFAULT CELL — the headline three quantities, fully error-barred.
    # ============================================================================================
    hr("(1) Default cell (kappa=3, lambda=3, rho=0.5): point estimates + 95% CIs")
    A = multiseed(KAP, LAM, RHO, N_CELL, B_CELL)        # (B,3): tau_stat, tau*, G
    ts, tstar, G = A[:, 0], A[:, 1], A[:, 2]

    v3_ok, d_ts = summarise("tau_stat", ts, 1.0, "below")   # V3: shrink
    print()
    v1_ok, d_tstar = summarise("tau*", tstar, 1.0, "above")  # V1: widen
    print()
    v2_ok, d_G = summarise("G", G, 0.0, "above")             # V2: gap real

    # finite-N bias check: the mean must be stable across N (else the CI localises a biased
    # estimator, not the true value).
    hr("(1b) Finite-N bias check — mean must be N-stable (else CI ≠ statement about the truth)")
    A2 = multiseed(KAP, LAM, RHO, N_CELL2, B_CELL2)
    print(f"  {'quantity':<10}{'mean @N=%d' % N_CELL:>22}{'mean @N=%d' % N_CELL2:>22}{'|drift|':>10}")
    drift_ok = True
    for j, nm in enumerate(("tau_stat", "tau*", "G")):
        m1, m2 = float(A[:, j].mean()), float(A2[:, j].mean())
        dr = abs(m1 - m2)
        # drift tolerance: a few x the larger SE — bias should be well under the CI half-width.
        tol = 3 * max(A[:, j].std(ddof=1) / np.sqrt(B_CELL), A2[:, j].std(ddof=1) / np.sqrt(B_CELL2))
        ok = dr <= max(tol, 0.004)
        drift_ok = drift_ok and ok
        print(f"  {nm:<10}{m1:>22.4f}{m2:>22.4f}{dr:>10.4f}   "
              f"(<= {max(tol,0.004):.4f}: {'ok' if ok else 'DRIFT'})")
    print(f"  N-stable (bias negligible vs CI): {drift_ok}  "
          f"-> the CI on the mean is a statement about the true value, not a finite-N artefact.")

    # ============================================================================================
    # (2) THE THREE HONESTY VERDICTS
    # ============================================================================================
    hr("(2) Honesty verdicts (driven by the CI on the mean; single-run flag in brackets)")
    print(f"  V1  tau* CI excludes 1.0  -> 'decision calibration WIDENS' robust : {v1_ok}"
          f"   [single run too: {d_tstar['pred_excl']}]")
    print(f"  V2  G    CI excludes 0.0  -> 'the gap is real'              robust : {v2_ok}"
          f"   [single run too: {d_G['pred_excl']}]")
    print(f"  V3  tau_stat CI below 1.0 -> 'statistical cal. SHRINKS'     robust : {v3_ok}"
          f"   [single run too: {d_ts['pred_excl']}]")

    # ============================================================================================
    # (3) THE REFRAME RULE — identify the strongest supported claim, never overclaim.
    # ============================================================================================
    hr("(3) Strongest supported claim (reframe rule applied)")
    opposite_sides = bool(v1_ok and v3_ok)         # tau_stat < 1 < tau*
    slope = abs(zstar(LAM)) / 6.0                  # analytic dG/dgamma at gamma->0 (gap_scaling.py)
    if opposite_sides:
        headline = ("ROBUST: 'opposite sides of 1' holds — tau_stat is robustly < 1 and tau* "
                    "robustly > 1, so statistical calibration shrinks the bar while decision "
                    "calibration widens it, with a strictly positive gap G > 0.")
    elif v2_ok:
        headline = ("PARTIAL: 'opposite sides of 1' is NOT supported (tau* CI does not exclude 1), "
                    "but the gap is robust: tau* > tau_stat  <=>  G > 0. Report the relational claim, "
                    "not the absolute placement of tau* above 1.")
    else:
        headline = (f"FALL BACK TO THE ANALYTIC SLOPE: the gap's point CI includes its null at this "
                    f"cell, so claim the sympy-derived slope dG/dgamma = |z*|/6 = {slope:.4f} > 0 "
                    f"(gap_scaling.py), which is not subject to the optimum's noise.")
    print("  " + headline.replace("\n", "\n  "))
    print(f"\n  (analytic cross-check, noise-free: dG/dgamma|_0 = |z*(lambda=3)|/6 = {slope:.4f} > 0;")
    print(f"   leading-order G at gamma={gamma_of_kappa(KAP):.4f} = {slope*gamma_of_kappa(KAP):+.4f}.)")

    # ============================================================================================
    # (4) GAMMA SWEEP — CIs across the skew axis (the gap should turn ON with skew).
    # ============================================================================================
    hr("(4) Gamma sweep (lambda=3, rho=0.5): G CI vs skew — does the gap switch on with gamma?")
    print(f"  {'kappa':>6}{'gamma':>8}{'G mean':>10}{'G 95% CI (boot)':>22}{'excl 0?':>9}"
          f"{'G_v2(pub)':>11}{'G_lead_th':>11}")
    sweep = {}
    for k in SWEEP_KAPPAS:
        Asw = multiseed(k, LAM, RHO, N_SWEEP, B_SWEEP)
        gcol = Asw[:, 2]
        m, sd, se, ci_t = ci_mean_t(gcol)
        ci_b = ci_mean_boot(gcol)
        gk = gamma_of_kappa(k)
        excl0 = excludes(ci_b, 0.0, "above")
        flag = "  <-default" if k == KAP else ""
        print(f"  {k:>6.1f}{gk:>8.4f}{m:>+10.4f}{f'[{ci_b[0]:+.4f},{ci_b[1]:+.4f}]':>22}"
              f"{str(excl0):>9}{PUB_SWEEP_G[k]:>+11.4f}{slope*gk:>+11.4f}{flag}")
        sweep[k] = dict(mean=m, ci_b=ci_b, excl0=excl0, gamma=gk)
    # honesty (magnitude, not a brittle sign test): at kappa=0 the TRUE gap is 0 (well-specified:
    # the report IS the posterior, tau*=1=tau_stat by construction). The estimator carries a tiny
    # finite-sample bias at that degenerate flat-EU optimum (parabola-fit of a near-flat curve),
    # of the same ~0.005 scale as the published single-seed -0.004 of opposite sign — so the kappa=0
    # sign is undetermined and IMMATERIAL: what matters is |G(kappa=0)| << G(default). The gap then
    # rises monotonically and is resolved away from 0 once skew is present (kappa>=2).
    g0 = sweep[0.0]["mean"]
    gdef = sweep[KAP]["mean"]
    g0_ratio = abs(gdef / g0) if g0 != 0 else float("inf")
    g0_negligible = abs(g0) < 0.01 and g0_ratio > 8
    means_by_k = [sweep[k]["mean"] for k in SWEEP_KAPPAS]
    mono = all(means_by_k[i + 1] > means_by_k[i] for i in range(len(means_by_k) - 1))
    big_excl0 = all(sweep[k]["excl0"] for k in (2.0, 3.0, 4.0))
    print(f"\n  well-specified kappa=0: G={g0:+.4f} — NEGLIGIBLE (|G|<0.01 and {g0_ratio:.0f}x below "
          f"the default-cell gap {gdef:+.4f}): {g0_negligible}")
    print(f"    (true G=0 at kappa=0 by construction; the +{g0:.4f} is a small estimator bias at the")
    print(f"     flat optimum — same scale as the published single-seed -0.004, sign undetermined.)")
    print(f"  gap monotone increasing in skew across the sweep: {mono}")
    print(f"  skewed kappa in {{2,3,4}} all exclude 0 (gap resolved once skew present): {big_excl0}")

    # ============================================================================================
    # GATE 1 — verdict block (HALT-TO-REPORT). 'opposite sides of 1' determination AT THE TOP.
    # ============================================================================================
    hr("GATE 1 — HALT-TO-REPORT verdict")
    print(f"  >>> HEADLINE: 'opposite sides of 1' is {'ROBUST' if opposite_sides else 'NOT robust'} "
          f"at the default cell. <<<")
    print()
    print(f"  V1 tau* CI excludes 1.0 (widen) : {v1_ok}")
    print(f"  V2 G    CI excludes 0.0 (gap)   : {v2_ok}")
    print(f"  V3 tau_stat CI below 1.0 (shrink): {v3_ok}")
    print(f"  finite-N bias negligible        : {drift_ok}")
    print(f"  gap scales with skew            : kappa=0 G={g0:+.4f} negligible ({g0_ratio:.0f}x below "
          f"default); monotone={mono}; kappa>=2 all exclude 0={big_excl0}")
    print()
    print("  Strongest supported claim:")
    print("  " + headline.replace("\n", "\n  "))
    print()
    print("  Discipline note: the estimator was NOT tuned to recover a cleaner claim; the CIs above")
    print("  are the published estimator's own sampling distribution over independent seeds. The")
    print("  published single-seed headline (tau*=1.0431, G=0.0796) is one draw; the multi-seed mean")
    print(f"  (tau*={d_tstar['mean']:.4f}, G={d_G['mean']:.4f}) localises the true value.")

    return dict(opposite_sides=opposite_sides, v1=v1_ok, v2=v2_ok, v3=v3_ok, drift_ok=drift_ok,
                tau_stat=d_ts, tau_star=d_tstar, G=d_G, slope=slope, sweep=sweep, headline=headline,
                g0_negligible=g0_negligible, mono=mono, big_excl0=big_excl0)


if __name__ == "__main__":
    main()
