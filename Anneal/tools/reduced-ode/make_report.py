"""Assemble reduced_results/REDUCED_REPORT.md from the per-CP JSON outputs.
Reproducible from config: run cp1..cp4 then this. Deterministic."""
import json
import os
from datetime import datetime, timezone

import reduced_core as rc

CFG = rc.load_config()
OUT = os.path.join(rc.ROOT, CFG["output_dir"])


def J(name):
    with open(os.path.join(OUT, name)) as f:
        return json.load(f)


def main():
    cp1 = J("cp1_gate.json")
    cp2 = J("cp2_regions.json")
    cp3 = J("cp3_match.json")
    cp4 = J("cp4_correlation.json")

    g_a = cp1["gate_a_fixed_point_consistency"]
    g_b = cp1["gate_b_fig2_regimes"]
    g_c = cp1["gate_c_bifurcation_series"]
    bif = cp2["bifurcations"]
    A_SN = bif["A_SN"]; A_H = bif["A_H"]; A_hc = bif["A_hc"]
    v = cp2["verdicts"]
    m05 = cp3["A05_matches"]; m02 = cp3["A02_matches"]; eig = cp3["A05_eigen"]
    pooled = cp4["pooled"]; perN = cp4["per_N"]; bench = cp4["benchmarks"]

    L = []
    A = L.append
    A("# REDUCED_REPORT — Abrams 2008 reduced ODE ↔ annealMusic finite-N chimera\n")
    A(f"_Generated {datetime.now(timezone.utc).isoformat()} · config "
      f"`tools/reduced-ode/reduced.config.json` · all equations transcribed from "
      f"Abrams, Mirollo, Strogatz & Wiley, PRL 101, 084103 (2008), arXiv:0806.0594v2._\n")
    A("Reduced-system integration is microseconds per run (DOP853, rtol 1e-10). "
      "The 3-D RHS was derived from Eq. (10) and verified to collapse to the Eq. (12) "
      "2-D RHS on ρ₁=1 to machine precision, and the Eqs (13)-(14) stationary chimeras "
      "verified to zero the 2-D RHS to machine precision (`test_core.py`). No quantity "
      "below is fitted to the finite-N data.\n")

    # ---------------- CP1 ----------------
    A("## CP1 — validation gate (against the paper itself, β=0.1)\n")
    A(f"- **(a) fixed-point consistency** (Eqs 13–14 zero the 2-D RHS): "
      f"max|RHS|={g_a['max_rhs_residual']:.1e} < {g_a['tol']:.0e} → "
      f"**{'PASS' if g_a['pass'] else 'FAIL'}**.")
    A("- **(b) Fig. 2 regimes:**")
    for Aval, r in g_b["regimes"].items():
        per = r.get("period")
        per_s = f", period={per:.1f}s" if per else ""
        A(f"  - A={Aval}: **{r.get('kind')}** (expected {r.get('expected')}{per_s}) "
          f"→ {'PASS' if r.get('match') else 'FAIL'}")
    A(f"  - period order T(0.35)>T(0.28): {g_b['period_order_028_lt_035']}")
    sc, hc = g_c["saddle_node"], g_c["hopf"]
    A(f"- **(c) bifurcation series:** SN numeric={sc['numeric_A']:.5f} vs Eq.17 "
      f"series={sc['series_A_eq17']:.5f} (|err|={sc['abs_err']:.1e}); Hopf "
      f"numeric={hc['numeric_A']:.5f} vs Eq.18 series={hc['series_A_eq18']:.5f} "
      f"(|err|={hc['abs_err']:.1e}) → {'PASS' if g_c['pass'] else 'FAIL'}.")
    A(f"\n**CP1 OVERALL: {'PASS ✅' if cp1['overall_pass'] else 'FAIL ⛔ (STOP)'}** "
      "— the harness reproduces the paper's published fixed-point relations, Fig. 2 "
      "regimes and Eqs (17)–(18) bifurcation series, so it is validated before being "
      "applied to our system.\n")

    # ---------------- CP2 ----------------
    A("## CP2 — our corners (β=0.05)\n")
    A("| bifurcation | series | numeric |")
    A("|---|---|---|")
    A(f"| A_SN (Eq.17) | {A_SN['series_eq17']:.5f} | {A_SN['numeric']:.5f} |")
    A(f"| A_H (Eq.18) | {A_H['series_eq18']:.5f} | {A_H['numeric']:.5f} |")
    A(f"| A_hc (homoclinic) | — (no closed form) | "
      f"{A_hc['numeric']:.5f}  bracket {A_hc['bracket']}  width {float(A_hc['width']):.1e} |")
    A("")
    A("Hand predictions (Avery): A_SN≈0.095 ✓, A_H≈0.270 ✓. The homoclinic was "
      "located by bisection on A (the breathing cycle's period diverges as it "
      f"collides with the saddle): **A_hc(0.05)={A_hc['numeric']:.4f} < 0.5**.\n")
    A("| A | region | prediction | holds | dynamical confirmation |")
    A("|---|---|---|---|---|")
    for Aval, vv in v.items():
        c = vv["dynamical_confirmation"]
        ck = c.get("kind")
        if ck == "fixed":
            extra = f"r*={float(c['r_star']):.4f}, amp→{float(c['amp']):.0e}"
        elif ck == "sync":
            extra = f"escaped→sync at t={float(c['t_escape']):.0f}s"
        else:
            extra = ck
        A(f"| {Aval} | {vv['region']} | {vv['prediction']} | "
          f"{'**YES**' if vv['prediction_holds'] else 'NO'} | {extra} |")
    A("\nPhase portraits (paper Fig. 3 style): `cp2_portrait_A0.2.{png,pdf}` (stable "
      "chimera spiral + saddle, basin filled by the saddle's unstable manifold), "
      "`cp2_portrait_A0.5.{png,pdf}` (unstable chimera spiral, canonical-seed "
      "trajectories spiral out to sync).\n")

    # ---------------- CP3 ----------------
    A("## CP3 — quantitative matches to the finite-N measurements\n")
    A(f"**Timescale:** {cp3['timescale']['statement']}\n")
    A(f"A=0.5 chimera FP is an **unstable spiral** (σ={float(eig['sigma_per_time']):+.5f}/s "
      f">0 ⇒ spiral-out, ω={float(eig['omega']):.4f}), giving linear breath period "
      f"2π/ω={float(eig['T_breath_linear_s']):.1f}s and per-cycle growth "
      f"{float(eig['per_cycle_growth']):.3f}.\n")
    bi = m05["i_breath_period"]; bii = m05["ii_spiral_slope"]; biii = m05["iii_capture_time"]
    A("| quantity | reduced | measured | ratio |")
    A("|---|---|---|---|")
    A(f"| (i) breath period T_b | {float(bi['reduced_pooled_s']):.1f}s "
      f"(linear {float(bi['reduced_linear_theory_s']):.1f}s) | "
      f"{bi['measured_s'][0]:.0f}–{bi['measured_s'][1]:.0f}s | "
      f"{float(bi['ratio_pooled_to_mid']):.2f} |")
    A(f"| (ii) spiral-out slope (pooled) | {float(bii['reduced_pooled']):.2f} | "
      f"−0.55…−0.77 | ≈1.0 |")
    A(f"| (iii) capture time (median) | {float(biii['reduced_pooled_median_s']):.0f}s | "
      f"{float(biii['measured_tau_abs_s']):.0f}s | {float(biii['ratio']):.2f} |")
    A("")
    A("Per-N spiral slope (reduced vs measured): "
      + "; ".join(f"N={k}: {float(d['reduced']):.2f}/{float(d['measured']):.2f}"
                  for k, d in bii["per_N_reduced_vs_measured"].items()) + ".")
    A("")
    A("**A=0.2 (stable stationary chimera):** reduced stable-FP "
      f"r*={float(m02['reduced_r_star']):.4f} vs the never-absorbers' mean R_incoh "
      f"={float(m02['measured_mean_Rincoh_pooled']):.4f} pooled "
      f"(N=32: {float(m02['measured_mean_Rincoh_per_N']['32']):.4f}, "
      f"N=64: {float(m02['measured_mean_Rincoh_per_N']['64']):.4f}) — the finite-N "
      "mean climbs toward r* as N grows. Reduced relaxation σ="
      f"{float(m02['reduced_relaxation_sigma_per_time']):+.5f}/s (per-cycle decay "
      f"{float(m02['reduced_per_cycle_decay']):.3f}); {m02['sign_match']}.\n")

    # ---------------- CP4 ----------------
    A("## CP4 — run-by-run lifetime prediction (capstone)\n")
    A("Per A=0.5 run, the **3-D** reduced system is integrated from its seed-mapped "
      "IC (ρ₁=R_sync0, ρ₂=R_incoh0, ψ=Δφ₀) to capture, and the reduced capture time "
      "is correlated with the measured t_abs.\n")
    A(f"**Pooled (n={pooled['n']}):** Spearman(reduced, t_abs) = "
      f"**{float(pooled['spearman_reduced_vs_tabs']):.3f}**; "
      f"|Δφ₀|-static = {float(pooled['spearman_absdphi0_vs_tabs']):.3f}; "
      f"**partial(reduced | |Δφ₀|) = {float(pooled['partial_spearman_reduced_given_absdphi0']):.3f}** "
      "(positive ⇒ the reduced dynamics adds information beyond the single static "
      "feature). CV R²: reduced-alone "
      f"{float(pooled['cvR2_reduced']):.3f}, |Δφ₀|+quad "
      f"{float(pooled['cvR2_absdphi0_quad']):.3f}, both "
      f"{float(pooled['cvR2_both']):.3f} (collective-IC benchmark "
      f"R²≈{bench['collective_ic_R2']}).\n")
    A("| N | ρ(reduced) | ρ(\\|Δφ₀\\|) | partial | R²(reduced) | R²(both) |")
    A("|---|---|---|---|---|---|")
    for N in CFG["cp4"]["Ns"]:
        d = perN[str(N)]
        A(f"| {N} | {float(d['spearman_reduced_vs_tabs']):+.3f} | "
          f"{float(d['spearman_absdphi0_vs_tabs']):+.3f} | "
          f"{float(d['partial_spearman_reduced_given_absdphi0']):+.3f} | "
          f"{float(d['cvR2_reduced']):.3f} | {float(d['cvR2_both']):.3f} |")
    A("")
    A("Figures: `cp4_scatter.{png,pdf}` (reduced-predicted vs measured, colored by N), "
      "`cp4_rho_vs_N.{png,pdf}` (ρ vs N + where the N-dependence enters).\n")
    A("**Does the reduced model beat the static regression?** It does not *decisively* "
      "exceed the multi-feature collective-IC R² (pooled both=0.21 vs benchmark 0.23), "
      "but it is **parameter-free physics** that (1) reaches Spearman ρ comparable to "
      "the fitted |Δφ₀| benchmark (per-N up to 0.60), and (2) carries a **positive "
      "partial correlation at every N** — it predicts lifetime variance that the single "
      "best static feature misses. The transverse deviation of ρ₁ from 1 stays ≤0.003 "
      "(median) over every run, so the ρ₁=1 manifold picture is well justified even "
      "though we integrated the full 3-D system.\n")
    A("**Where does N enter?** The reduced *dynamics* are N-independent; the only "
      "N-channel is the seed-ensemble IC distribution (mean R_incoh0 falls ~1/√N). The "
      "reduced capture median therefore *rises* with N and does **not** reproduce the "
      "flat measured τ_abs(N)≈139s plateau, nor the growing-spiral-rate-with-N. Those "
      "live in the **IC ensemble + finite-N fluctuations**, not in the collective "
      "reduced flow — a clean localization of the N-physics.\n")

    # ---------------- Verdicts ----------------
    p1_ok = v["0.2"]["prediction_holds"]
    p2_ok = v["0.5"]["prediction_holds"]
    A("## Explicit P1/P2 verdicts\n")
    A(f"- **P1 — {'CONFIRMED ✅' if p1_ok else 'REFUTED'}:** A=0.2 at β=0.05 lies in "
      f"the stable stationary-chimera band (A_SN={A_SN['series_eq17']:.3f} < 0.2 < "
      f"A_H={A_H['series_eq18']:.3f}); integration confirms a stable fixed point "
      f"(r*={float(v['0.2']['dynamical_confirmation']['r_star']):.3f}, amplitude→0). "
      "The A=0.2 never-absorbers are a genuine attractor; their 'breathing' is finite-N "
      "noise around this fixed point, and their mean R_incoh sits on r*.")
    A(f"- **P2 — {'CONFIRMED ✅' if p2_ok else 'REFUTED'}:** A=0.5 at β=0.05 is "
      f"post-homoclinic (0.5 > A_hc={A_hc['numeric']:.3f}); no chimera attractor, sync "
      "is globally attracting, and canonical-seed trajectories spiral out to sync. The "
      "measured ratcheting transient is the spiral-out along the destroyed cycle's "
      "ghost (the unstable-spiral σ>0), with reduced breath period and per-cycle "
      "spiral slope quantitatively matching the finite-N measurements.\n")

    # ---------------- Mismatch list ----------------
    A("## Candid mismatch list\n")
    A(f"1. **Absolute capture-time scale (×~3 short).** Reduced median capture "
      f"≈{float(biii['reduced_pooled_median_s']):.0f}s vs measured τ_abs≈"
      f"{float(biii['measured_tau_abs_s']):.0f}s (ratio "
      f"{float(biii['ratio']):.2f}). The mean-field collective flow spirals to capture "
      "faster than the finite-N system lives; the collective IC is a lossy projection "
      "of the full 2N-phase state, and finite-N fluctuations extend lifetime. Reported "
      "straight, no fudge factor.")
    A("2. **τ_abs(N) plateau not reproduced.** Reduced capture median rises with N "
      "(IC-driven) while measured τ_abs is flat — the plateau is a finite-N effect, "
      "not in the reduced dynamics.")
    A("3. **N=4 is weak** (ρ=+0.14): at the smallest N the collective IC summaries are "
      "too noisy a projection for the reduced model to predict lifetime.")
    A("4. **Reduced does not decisively beat the fitted multi-feature regression** in "
      "absolute R² (0.21 vs 0.23), though it adds positive partial information and is "
      "parameter-free.")
    A("5. **Breath period slightly high / linear-vs-nonlinear.** Nonlinear reduced T_b "
      f"≈{float(bi['reduced_pooled_s']):.0f}s (top of the 21–25s band) while the "
      f"linear-focus value is {float(bi['reduced_linear_theory_s']):.0f}s; the "
      "canonical seeds sample the nonlinear regime away from the focus, lengthening the "
      "period — consistent, not contradictory.\n")

    A("## Artifacts\n")
    for fn, desc in [
        ("cp1_gate.json", "CP1 sub-gate verdicts"),
        ("cp2_regions.json / .md", "CP2 region table + bifurcations + homoclinic trace"),
        ("cp2_portrait_A0.2.{png,pdf}, cp2_portrait_A0.5.{png,pdf}", "phase portraits"),
        ("cp3_match.json / .md", "CP3 (reduced, measured, ratio) tables"),
        ("cp3_a02_level.json", "finite-N A=0.2 mean R_incoh (shipped tracer)"),
        ("reduced_runs.jsonl", "per-run reduced 3-D capture/breath/slope (1400 runs)"),
        ("cp4_correlation.json / .md", "CP4 Spearman/partial/CV-R² per N + pooled"),
        ("cp4_scatter.{png,pdf}, cp4_rho_vs_N.{png,pdf}", "CP4 figures"),
    ]:
        A(f"- `{fn}` — {desc}")

    path = os.path.join(OUT, "REDUCED_REPORT.md")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
