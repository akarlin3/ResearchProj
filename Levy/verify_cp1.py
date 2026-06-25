"""CP1 gate -- the joint (alpha, beta) structural-identifiability object is built AND the
degeneracy characterization has a definite, scoped verdict.

Six independent checks, all green => CP1 gate passes:

  (1) UNIT TESTS    the levy-core pytest suite passes (Mittag-Leffler vs closed forms,
                    forward_joint reduction, joint Fisher PSD/limit, joint MLE recovers truth
                    with two diffusion times, etc.).
  (2) ML KERNEL     E_alpha(-t): exp at alpha=1, erfcx at alpha=1/2 (the two closed forms).
  (3) REDUCTION     forward_joint(alpha=1) == the CP0 stretched-exponential (exact anchor):
                    the joint model contains CP0 as its time-order-1 face.
  (4) DEGENERACY    the joint FIM (alpha,beta) degeneracy is characterised over the
                    physiological grid: a median |rho_alpha_beta| and condition number, with
                    bootstrap CIs on the headline (alpha_hat, beta_hat).
  (5) BOUNDARY      a SECOND diffusion time provably breaks the degeneracy (|rho| and cond drop,
                    cv_alpha recovers) -- the scope limit, reported not hidden.
  (6) VERDICT       a definite, honest answer (degenerate xor separable) with the regime scope.

Run:  <proteus python> verify_cp1.py            # FAST (smoke; default)
      <proteus python> verify_cp1.py --full     # full-N bootstrap
Exit 0 = CP1 gate green; nonzero = a check failed (the script names which).
"""
from __future__ import annotations

import subprocess
import sys
import time

import numpy as np
from scipy import special

import _paths  # noqa: E402

FULL = "--full" in sys.argv


def _hr(title: str) -> None:
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78)


def check_unit_tests() -> None:
    _hr("CP1 check 1/6 -- levy-core unit test suite")
    core = _paths.LEVY_CORE
    t0 = time.time()
    r = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=str(core),
                       capture_output=True, text=True)
    out = r.stdout + r.stderr
    print(f"  pytest -> {'PASS' if r.returncode == 0 else 'FAIL'} (rc={r.returncode}, {time.time()-t0:.1f}s)")
    if r.returncode != 0:
        print("\n".join(out.splitlines()[-25:]))
        raise SystemExit("UNIT TESTS FAILED")
    print("  " + out.strip().splitlines()[-1])
    print("  UNIT TESTS: PASS")


def check_ml_kernel() -> None:
    _hr("CP1 check 2/6 -- Mittag-Leffler kernel vs the two closed forms")
    from levy import mittag_leffler as ml
    t = np.linspace(0.0, 8.0, 30)
    e1 = np.max(np.abs(ml.mlf_neg(1.0, t) - np.exp(-t)))
    e_half = np.max(np.abs(ml.mlf_neg(0.5, t) - special.erfcx(t)))
    print(f"  max|E_1(-t)  - exp(-t)|   = {e1:.2e}")
    print(f"  max|E_.5(-t) - erfcx(t)|  = {e_half:.2e}")
    assert e1 < 1e-9 and e_half < 1e-5
    print("  ML KERNEL: PASS")


def check_reduction() -> None:
    _hr("CP1 check 3/6 -- forward_joint(alpha=1) == CP0 stretched-exponential (exact anchor)")
    from levy import forward
    b = np.array([0.0, 300.0, 1000.0, 2000.0, 3000.0])
    for S0, D, beta in ((1.0, 1.5e-3, 1.7), (2.0, 1.0e-3, 1.4)):
        joint = forward.forward_joint(b, np.array([S0, D, 1.0, beta]))
        cp0 = forward.signal(b, np.array([S0, D, beta / 2.0]))  # heterogeneity exponent == beta/2
        err = np.max(np.abs(joint - cp0))
        print(f"  S0={S0} D={D} beta={beta}: max|joint(a=1) - CP0(a=beta/2)| = {err:.2e}")
        assert err < 1e-9
    print("  REDUCTION: PASS")


def check_degeneracy(rep) -> None:
    _hr("CP1 check 4/6 -- joint (alpha,beta) degeneracy characterised (with bootstrap CIs)")
    print(f"  physiological grid: median |rho_alpha_beta| = {rep.rho_median:.3f} "
          f"(range [{rep.rho_min:.3f}, {rep.rho_max:.3f}]), median FIM cond = {rep.cond_median:.2e}")
    print(f"  headline cell (alpha=0.80,beta=1.70): rho_ab={rep.headline_rho:+.3f}, "
          f"cond={rep.headline_cond:.2e}, cv_alpha={rep.headline_cv_alpha:.2f}, cv_beta={rep.headline_cv_beta:.2f}")
    print(f"  bootstrap alpha_hat 95% CI = [{rep.boot_alpha_ci[0]:.3f}, {rep.boot_alpha_ci[1]:.3f}] (truth 0.80)")
    print(f"  bootstrap beta_hat  95% CI = [{rep.boot_beta_ci[0]:.3f}, {rep.boot_beta_ci[1]:.3f}] (truth 1.70)")
    print(f"  empirical corr(alpha_hat, beta_hat) = {rep.boot_corr:+.3f} (the ridge)")
    print(f"  n_b persistence |rho|: " + ", ".join(f"n_b={nb}:{r:.3f}" for nb, r in zip(rep.nb_list, rep.nb_rho)))
    assert 0.0 <= rep.rho_median <= 1.0
    assert np.isfinite(rep.headline_cond) and rep.headline_cond > 1.0
    assert rep.boot_alpha_ci[0] <= rep.boot_alpha_ci[1]
    assert rep.boot_beta_ci[0] <= rep.boot_beta_ci[1]
    # persistence: the degeneracy must NOT vanish as n_b grows at a single diffusion time
    assert rep.nb_rho.min() >= 0.5, "degeneracy unexpectedly collapses with n_b at single dt"
    print("  DEGENERACY: PASS")


def check_boundary(rep) -> None:
    _hr("CP1 check 5/6 -- a second diffusion time provably breaks the degeneracy")
    print(f"  |rho_ab|: {rep.break_rho_single:.3f} -> {rep.break_rho_two:.3f}")
    print(f"  FIM cond: {rep.break_cond_single:.2e} -> {rep.break_cond_two:.2e}")
    print(f"  cv_alpha: {rep.break_cv_alpha_single:.2f} -> {rep.break_cv_alpha_two:.2f}")
    assert rep.break_rho_two < rep.break_rho_single
    assert rep.break_cond_two < rep.break_cond_single
    assert rep.break_cv_alpha_two < rep.break_cv_alpha_single
    print("  BOUNDARY: PASS")


def check_verdict(rep) -> None:
    _hr("CP1 check 6/6 -- definite, scoped verdict")
    from levy import degeneracy
    print(f"  degenerate = {rep.degenerate}  (pre-registered threshold median|rho|>={degeneracy.DEGEN_RHO})")
    for n in rep.notes:
        print(f"  - {n}")
    # a definite characterization: the verdict matches the threshold logic exactly
    assert rep.degenerate == (rep.rho_median >= degeneracy.DEGEN_RHO)
    print("  VERDICT: PASS (characterisation, scoped to single-diffusion-time clinical acquisition)")


def main() -> int:
    print("CP1 verification -- joint CTRW (alpha, beta) structural identifiability + degeneracy")
    _paths.add_all()
    check_unit_tests()
    check_ml_kernel()
    check_reduction()

    from levy import degeneracy, seeding
    n_boot = 200 if FULL else 60
    t0 = time.time()
    rep = degeneracy.cp1_report(rng=seeding.make_rng(), n_boot=n_boot)
    print(f"\n  (cp1_report computed in {time.time()-t0:.1f}s, n_boot={n_boot})")

    check_degeneracy(rep)
    check_boundary(rep)
    check_verdict(rep)

    _hr("CP1 GATE: PASS")
    print("  Joint (alpha,beta) Fisher/CRLB structural-identifiability object built (net-new);")
    print("  degeneracy characterised over the physiological grid with bootstrap CIs; the")
    print("  constructive two-diffusion-time boundary is reported. See results/RESULTS_CP1.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
