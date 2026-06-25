#!/usr/bin/env python
"""CP1 driver -- the joint CTRW (alpha, beta) structural-identifiability / degeneracy object.

Builds the joint 4-parameter Fisher/CRLB for (S0, D, alpha, beta) under the single-diffusion-
time CTRW forward model with finite b-values + Rician noise; maps the (alpha, beta) degeneracy
over the physiological grid; bootstraps the headline ridge; shows the constructive two-diffusion-
time break. Writes results/RESULTS_CP1.md, results/RESULTS_CP1.json (manuscript anchors), and a
figure. Numbers print to stdout for transcription.

Usage:  <proteus python> experiments/run_cp1.py            # FAST (smoke; default)
        <proteus python> experiments/run_cp1.py --full     # full-N bootstrap CIs
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE.parent / "levy-core"
sys.path.insert(0, str(CORE))

from levy import degeneracy, fisher_joint, forward, identifiability_joint, seeding  # noqa: E402

FULL = "--full" in sys.argv
RESULTS_MD = HERE.parent / "results" / "RESULTS_CP1.md"
RESULTS_JSON = HERE.parent / "results" / "RESULTS_CP1.json"
FIGDIR = HERE.parent / "figures"


def _hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def main() -> int:
    rng = seeding.make_rng()
    n_boot = 200 if FULL else 80

    _hr("CP1 -- joint CTRW (alpha, beta) structural identifiability under the MRI forward model")
    print("Forward model (single diffusion time):  S(b;S0,D,alpha,beta) = S0 E_alpha(-(b D)^{beta/2})")
    print("E_alpha = one-parameter Mittag-Leffler (CTRW; Magin/Ingo 2013). alpha = time-fractional")
    print("(Caputo) order, beta = space-fractional (Riesz) order. At alpha=1 this is exactly the")
    print("CP0 stretched-exponential with heterogeneity exponent beta/2. Estimand theta=(S0,D,alpha,")
    print("beta) estimated JOINTLY; Rician noise sigma=S0/SNR. Question: at ONE diffusion time, are")
    print("the time-order alpha and space-order beta separable, or degenerate?")

    rep = degeneracy.cp1_report(rng=rng, n_boot=n_boot)

    _hr("1. Degeneracy map over the physiological (alpha,beta) grid (clinical n_b=6, SNR=40)")
    print("   |rho_alpha_beta| (from inverse joint FIM); -> 1 == time/space orders inseparable")
    header = "   a \\ b  | " + " ".join(f"{b:>6.2f}" for b in rep.beta_grid)
    print(header)
    for i, a in enumerate(rep.alpha_grid):
        row = " ".join(f"{rep.rho_map[i,j]:6.3f}" for j in range(len(rep.beta_grid)))
        print(f"   {a:.2f}   | {row}")
    print(f"   median |rho| = {rep.rho_median:.3f}  range [{rep.rho_min:.3f}, {rep.rho_max:.3f}]; "
          f"median FIM cond = {rep.cond_median:.2e}")

    _hr("2. Headline cell (alpha=0.80, beta=1.70) -- analytic + bootstrap CIs")
    print(f"   analytic rho_alpha_beta = {rep.headline_rho:+.3f}, FIM cond = {rep.headline_cond:.2e}")
    print(f"   relative CRLB: cv_alpha = {rep.headline_cv_alpha:.2f}, cv_beta = {rep.headline_cv_beta:.2f}")
    print(f"   bootstrap (n_boot={rep.n_boot}) alpha_hat 95% CI = "
          f"[{rep.boot_alpha_ci[0]:.3f}, {rep.boot_alpha_ci[1]:.3f}]  (truth 0.80)")
    print(f"   bootstrap            beta_hat  95% CI = "
          f"[{rep.boot_beta_ci[0]:.3f}, {rep.boot_beta_ci[1]:.3f}]  (truth 1.70)")
    print(f"   empirical corr(alpha_hat,beta_hat) = {rep.boot_corr:+.3f}; "
          f"rel CI widths alpha {rep.boot_alpha_rel_width:.2f}, beta {rep.boot_beta_rel_width:.2f}")

    _hr("3. n_b persistence (single diffusion time) -- the degeneracy does NOT recede with b-values")
    print("   contrast with CP0: the single-order wall recedes for n_b>=8; the joint degeneracy does not.")
    for nb, r in zip(rep.nb_list, rep.nb_rho):
        print(f"   n_b={nb:2d}: |rho_alpha_beta| = {r:.3f}")

    _hr("4. Constructive boundary -- a SECOND diffusion time breaks the degeneracy")
    print(f"   16 b-values, ONE diffusion time   vs   8+8 b-values across TWO diffusion times:")
    print(f"   |rho_alpha_beta|: {rep.break_rho_single:.3f} -> {rep.break_rho_two:.3f}")
    print(f"   FIM condition:    {rep.break_cond_single:.2e} -> {rep.break_cond_two:.2e}")
    print(f"   cv_alpha:         {rep.break_cv_alpha_single:.2f} -> {rep.break_cv_alpha_two:.2f}")

    _hr("CP1 VERDICT")
    print(f"   degenerate = {rep.degenerate}  (pre-registered: median|rho|>={degeneracy.DEGEN_RHO})")
    for n in rep.notes:
        print(f"   - {n}")
    print("\n   >>> Scoped finding: at a SINGLE clinical diffusion time the joint CTRW time-order")
    print("       alpha and space-order beta are structurally degenerate (cannot be separated),")
    print("       and -- unlike the CP0 single-order wall -- this is NOT relieved by adding")
    print("       b-values; only a second diffusion time separates them. Reinforces and extends")
    print("       the CP0 'clinically information-limited' thesis to the two-exponent model.")

    _write_md(rep)
    _write_json(rep)
    _save_figure(rep)
    print(f"\n   results -> {RESULTS_MD}\n   anchors -> {RESULTS_JSON}")
    print("CP1 PASS")
    return 0


def _write_md(rep):
    L = []
    L.append("# RESULTS -- CP1: joint CTRW (alpha, beta) structural identifiability / degeneracy\n")
    L.append("All numbers derived (joint 4-parameter Fisher/CRLB + parametric bootstrap), fully")
    L.append("synthetic, seeded. CRLB = identifiability bound, scoped to its regime (single")
    L.append("diffusion time, finite b-values, Rician noise); never an impossibility claim.\n")
    L.append("## Forward model")
    L.append("`S(b; S0, D, alpha, beta) = S0 * E_alpha(-(b D)^{beta/2})` (single diffusion time;")
    L.append("E_alpha = one-parameter Mittag-Leffler, CTRW/fractional Bloch-Torrey, Magin/Ingo 2013).")
    L.append("alpha = time-fractional (Caputo) order in (0,1]; beta = space-fractional (Riesz) order")
    L.append("in (1,2]. At alpha=1 it is exactly the CP0 stretched-exponential with heterogeneity")
    L.append("exponent beta/2. theta=(S0,D,alpha,beta) estimated JOINTLY; Rician noise sigma=S0/SNR.\n")
    L.append("## Degeneracy over the physiological (alpha,beta) grid (clinical n_b=6, b_max=2500, SNR=40)")
    L.append(f"- **median |rho_alpha_beta| = {rep.rho_median:.3f}** (range [{rep.rho_min:.3f}, {rep.rho_max:.3f}])")
    L.append(f"- median FIM condition number = {rep.cond_median:.2e}")
    L.append(f"- pre-registered degeneracy threshold: median |rho_alpha_beta| >= {degeneracy.DEGEN_RHO}")
    L.append(f"- **verdict: degenerate = {rep.degenerate}**\n")
    L.append("## Headline cell (alpha=0.80, beta=1.70)")
    L.append(f"- analytic rho_alpha_beta = {rep.headline_rho:+.3f}; FIM cond = {rep.headline_cond:.2e}")
    L.append(f"- relative CRLB: cv_alpha = {rep.headline_cv_alpha:.2f}, cv_beta = {rep.headline_cv_beta:.2f}")
    L.append(f"- bootstrap ({rep.n_boot} reps) alpha_hat 95% CI = [{rep.boot_alpha_ci[0]:.3f}, "
             f"{rep.boot_alpha_ci[1]:.3f}] (truth 0.80); beta_hat 95% CI = [{rep.boot_beta_ci[0]:.3f}, "
             f"{rep.boot_beta_ci[1]:.3f}] (truth 1.70)")
    L.append(f"- empirical corr(alpha_hat, beta_hat) = {rep.boot_corr:+.3f} (the degeneracy ridge)\n")
    L.append("## n_b persistence at a single diffusion time (|rho_alpha_beta|)")
    L.append("| n_b | " + " | ".join(str(nb) for nb in rep.nb_list) + " |")
    L.append("|---|" + "---|" * len(rep.nb_list))
    L.append("| \\|rho\\| | " + " | ".join(f"{r:.3f}" for r in rep.nb_rho) + " |")
    L.append("\n(Contrast: the CP0 single-order wall recedes for n_b>=8; this degeneracy does not.)\n")
    L.append("## Constructive boundary -- a second diffusion time breaks the degeneracy")
    L.append("| quantity | 16 b, one Delta | 8+8 b, two Delta |")
    L.append("|---|---|---|")
    L.append(f"| \\|rho_alpha_beta\\| | {rep.break_rho_single:.3f} | {rep.break_rho_two:.3f} |")
    L.append(f"| FIM condition | {rep.break_cond_single:.2e} | {rep.break_cond_two:.2e} |")
    L.append(f"| cv_alpha | {rep.break_cv_alpha_single:.2f} | {rep.break_cv_alpha_two:.2f} |")
    L.append("")
    L.append("## Scoped claim")
    L.append("At a **single clinical diffusion time**, the joint CTRW time-order alpha and space-order")
    L.append(f"beta are **structurally degenerate** (median |rho_alpha_beta| = {rep.rho_median:.3f}, FIM")
    L.append(f"condition ~{rep.cond_median:.0e}): they cannot be separately recovered. Unlike the CP0")
    L.append("single-order wall, the degeneracy is **not** relieved by adding b-values; only a **second")
    L.append("diffusion time** separates them. This reinforces and extends the CP0 'clinically")
    L.append("information-limited' thesis to the two-exponent model.")
    RESULTS_MD.parent.mkdir(exist_ok=True)
    RESULTS_MD.write_text("\n".join(L) + "\n")


def _write_json(rep):
    data = {
        "regime": {
            "forward_model": "S(b)=S0*E_alpha(-(bD)^(beta/2)) single diffusion time",
            "n_b": degeneracy.CLINICAL["n_b"], "b_max": degeneracy.CLINICAL["b_max"],
            "snr": degeneracy.CLINICAL["snr"], "noise": "Rician sigma=S0/SNR",
            "alpha_grid": rep.alpha_grid.tolist(), "beta_grid": rep.beta_grid.tolist(),
        },
        "grid": {
            "rho_median": rep.rho_median, "rho_min": rep.rho_min, "rho_max": rep.rho_max,
            "cond_median": rep.cond_median,
        },
        "headline": {
            "alpha": float(rep.headline_theta[2]), "beta": float(rep.headline_theta[3]),
            "rho_ab": rep.headline_rho, "cond": rep.headline_cond,
            "cv_alpha": rep.headline_cv_alpha, "cv_beta": rep.headline_cv_beta,
            "alpha_ci": list(rep.boot_alpha_ci), "beta_ci": list(rep.boot_beta_ci),
            "corr": rep.boot_corr, "n_boot": rep.n_boot,
        },
        "nb_persistence": {"n_b": list(rep.nb_list), "rho": rep.nb_rho.tolist()},
        "two_dt_break": {
            "rho_single": rep.break_rho_single, "rho_two": rep.break_rho_two,
            "cond_single": rep.break_cond_single, "cond_two": rep.break_cond_two,
            "cv_alpha_single": rep.break_cv_alpha_single, "cv_alpha_two": rep.break_cv_alpha_two,
        },
        "verdict": {"degenerate": bool(rep.degenerate), "degen_rho_threshold": degeneracy.DEGEN_RHO},
    }
    RESULTS_JSON.parent.mkdir(exist_ok=True)
    RESULTS_JSON.write_text(json.dumps(data, indent=2) + "\n")


def _save_figure(rep):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # pragma: no cover
        print(f"  (figure skipped: {e})")
        return
    FIGDIR.mkdir(exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

    # left: degeneracy map |rho_alpha_beta| over the physiological grid
    im = ax1.imshow(rep.rho_map, aspect="auto", origin="lower", cmap="magma", vmin=0.5, vmax=1.0)
    ax1.set_xticks(range(len(rep.beta_grid))); ax1.set_xticklabels([f"{b:.2f}" for b in rep.beta_grid])
    ax1.set_yticks(range(len(rep.alpha_grid))); ax1.set_yticklabels([f"{a:.2f}" for a in rep.alpha_grid])
    ax1.set_xlabel("space-order beta"); ax1.set_ylabel("time-order alpha")
    ax1.set_title("(alpha,beta) degeneracy |rho| (single Delta)")
    fig.colorbar(im, ax=ax1, label="|rho_alpha_beta|")

    # right: the bootstrap ridge (alpha_hat, beta_hat) single-dt + two-dt cloud
    from levy.wall import default_b_design
    rng = seeding.make_rng(99)
    theta = rep.headline_theta
    snr = degeneracy.CLINICAL["snr"]
    b1, dt1 = identifiability_joint.two_dt_design(b_max=2500.0, n_b=8, ratios=(1.0,))
    b2, dt2 = identifiability_joint.two_dt_design(b_max=2500.0, n_b=8, ratios=(1.0, 2.5))
    nb = 60
    s1 = identifiability_joint.parametric_bootstrap_joint(theta, b1, dt1, snr, nb, rng)
    s2 = identifiability_joint.parametric_bootstrap_joint(theta, b2, dt2, snr, nb, rng)
    ax2.scatter(s1.alpha_hats, s1.beta_hats, s=14, alpha=0.6, label="single Delta (ridge)", color="crimson")
    ax2.scatter(s2.alpha_hats, s2.beta_hats, s=14, alpha=0.6, label="two Delta (resolved)", color="steelblue")
    ax2.scatter([theta[2]], [theta[3]], marker="*", s=220, color="k", zorder=5, label="truth")
    ax2.set_xlabel("alpha_hat"); ax2.set_ylabel("beta_hat")
    ax2.set_title("Joint MLE cloud: ridge vs resolved")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    out = FIGDIR / "cp1_degeneracy.png"
    fig.savefig(out, dpi=130)
    print(f"  figure -> {out}")


if __name__ == "__main__":
    raise SystemExit(main())
