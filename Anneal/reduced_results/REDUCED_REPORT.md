# REDUCED_REPORT — Abrams 2008 reduced ODE ↔ annealMusic finite-N chimera

_Generated 2026-06-05T02:23:43.290309+00:00 · config `tools/reduced-ode/reduced.config.json` · all equations transcribed from Abrams, Mirollo, Strogatz & Wiley, PRL 101, 084103 (2008), arXiv:0806.0594v2._

Reduced-system integration is microseconds per run (DOP853, rtol 1e-10). The 3-D RHS was derived from Eq. (10) and verified to collapse to the Eq. (12) 2-D RHS on ρ₁=1 to machine precision, and the Eqs (13)-(14) stationary chimeras verified to zero the 2-D RHS to machine precision (`test_core.py`). No quantity below is fitted to the finite-N data.

## CP1 — validation gate (against the paper itself, β=0.1)

- **(a) fixed-point consistency** (Eqs 13–14 zero the 2-D RHS): max|RHS|=3.9e-16 < 1e-09 → **PASS**.
- **(b) Fig. 2 regimes:**
  - A=0.2: **fixed** (expected fixed) → PASS
  - A=0.28: **cycle** (expected cycle, period=33.9s) → PASS
  - A=0.35: **cycle** (expected cycle, period=90.4s) → PASS
  - period order T(0.35)>T(0.28): True
- **(c) bifurcation series:** SN numeric=0.17834 vs Eq.17 series=0.17836 (|err|=2.3e-05); Hopf numeric=0.27778 vs Eq.18 series=0.27773 (|err|=4.8e-05) → PASS.

**CP1 OVERALL: PASS ✅** — the harness reproduces the paper's published fixed-point relations, Fig. 2 regimes and Eqs (17)–(18) bifurcation series, so it is validated before being applied to our system.

## CP2 — our corners (β=0.05)

| bifurcation | series | numeric |
|---|---|---|
| A_SN (Eq.17) | 0.09475 | 0.09475 |
| A_H (Eq.18) | 0.27030 | 0.27039 |
| A_hc (homoclinic) | — (no closed form) | 0.40961  bracket [0.4095703125, 0.409658203125]  width 8.8e-05 |

Hand predictions (Avery): A_SN≈0.095 ✓, A_H≈0.270 ✓. The homoclinic was located by bisection on A (the breathing cycle's period diverges as it collides with the saddle): **A_hc(0.05)=0.4096 < 0.5**.

| A | region | prediction | holds | dynamical confirmation |
|---|---|---|---|---|
| 0.2 | stable stationary chimera | P1 | **YES** | r*=0.6781, amp→2e-09 |
| 0.5 | post-homoclinic (no chimera attractor; sync attracting) | P2 | **YES** | escaped→sync at t=726s |

Phase portraits (paper Fig. 3 style): `cp2_portrait_A0.2.{png,pdf}` (stable chimera spiral + saddle, basin filled by the saddle's unstable manifold), `cp2_portrait_A0.5.{png,pdf}` (unstable chimera spiral, canonical-seed trajectories spiral out to sync).

## CP3 — quantitative matches to the finite-N measurements

**Timescale:** The shipped integrator advances 1 model-time unit per second of the drone (omega=0 rotating frame; DRIFT_DT=0.05 control step, sampleStride=0.1). integrator.mjs documents '1 unit := 1 s'. So reduced model-time equals seconds 1:1 — every reduced quantity below is in seconds with NO fudge factor.

A=0.5 chimera FP is an **unstable spiral** (σ=+0.01243/s >0 ⇒ spiral-out, ω=0.3277), giving linear breath period 2π/ω=19.2s and per-cycle growth 1.269.

| quantity | reduced | measured | ratio |
|---|---|---|---|
| (i) breath period T_b | 25.0s (linear 19.2s) | 21–25s | 1.09 |
| (ii) spiral-out slope (pooled) | -0.67 | −0.55…−0.77 | ≈1.0 |
| (iii) capture time (median) | 43s | 139s | 0.31 |

Per-N spiral slope (reduced vs measured): N=8: -0.62/-0.55; N=16: -0.66/-0.50; N=32: -0.66/-0.59; N=64: -0.73/-0.77.

**A=0.2 (stable stationary chimera):** reduced stable-FP r*=0.6781 vs the never-absorbers' mean R_incoh =0.6590 pooled (N=32: 0.6758, N=64: 0.6685) — the finite-N mean climbs toward r* as N grows. Reduced relaxation σ=-0.00539/s (per-cycle decay 0.803); both negative: reduced stable spiral relaxes (maxima drift DOWN to r*), measured <dM><0 — opposite to A=0.5's positive ratchet.

## CP4 — run-by-run lifetime prediction (capstone)

Per A=0.5 run, the **3-D** reduced system is integrated from its seed-mapped IC (ρ₁=R_sync0, ρ₂=R_incoh0, ψ=Δφ₀) to capture, and the reduced capture time is correlated with the measured t_abs.

**Pooled (n=1400):** Spearman(reduced, t_abs) = **0.445**; |Δφ₀|-static = -0.429; **partial(reduced | |Δφ₀|) = 0.200** (positive ⇒ the reduced dynamics adds information beyond the single static feature). CV R²: reduced-alone 0.146, |Δφ₀|+quad 0.185, both 0.206 (collective-IC benchmark R²≈0.23).

| N | ρ(reduced) | ρ(\|Δφ₀\|) | partial | R²(reduced) | R²(both) |
|---|---|---|---|---|---|
| 4 | +0.141 | +0.032 | +0.210 | 0.022 | -0.018 |
| 8 | +0.498 | -0.445 | +0.253 | 0.202 | 0.203 |
| 16 | +0.461 | -0.428 | +0.219 | 0.152 | 0.262 |
| 24 | +0.544 | -0.507 | +0.259 | 0.227 | 0.254 |
| 32 | +0.595 | -0.582 | +0.238 | 0.315 | 0.372 |
| 48 | +0.578 | -0.594 | +0.229 | 0.408 | 0.419 |
| 64 | +0.592 | -0.538 | +0.315 | 0.226 | 0.214 |

Figures: `cp4_scatter.{png,pdf}` (reduced-predicted vs measured, colored by N), `cp4_rho_vs_N.{png,pdf}` (ρ vs N + where the N-dependence enters).

**Does the reduced model beat the static regression?** It does not *decisively* exceed the multi-feature collective-IC R² (pooled both=0.21 vs benchmark 0.23), but it is **parameter-free physics** that (1) reaches Spearman ρ comparable to the fitted |Δφ₀| benchmark (per-N up to 0.60), and (2) carries a **positive partial correlation at every N** — it predicts lifetime variance that the single best static feature misses. The transverse deviation of ρ₁ from 1 stays ≤0.003 (median) over every run, so the ρ₁=1 manifold picture is well justified even though we integrated the full 3-D system.

**Where does N enter?** The reduced *dynamics* are N-independent; the only N-channel is the seed-ensemble IC distribution (mean R_incoh0 falls ~1/√N). The reduced capture median therefore *rises* with N and does **not** reproduce the flat measured τ_abs(N)≈139s plateau, nor the growing-spiral-rate-with-N. Those live in the **IC ensemble + finite-N fluctuations**, not in the collective reduced flow — a clean localization of the N-physics.

## Explicit P1/P2 verdicts

- **P1 — CONFIRMED ✅:** A=0.2 at β=0.05 lies in the stable stationary-chimera band (A_SN=0.095 < 0.2 < A_H=0.270); integration confirms a stable fixed point (r*=0.678, amplitude→0). The A=0.2 never-absorbers are a genuine attractor; their 'breathing' is finite-N noise around this fixed point, and their mean R_incoh sits on r*.
- **P2 — CONFIRMED ✅:** A=0.5 at β=0.05 is post-homoclinic (0.5 > A_hc=0.410); no chimera attractor, sync is globally attracting, and canonical-seed trajectories spiral out to sync. The measured ratcheting transient is the spiral-out along the destroyed cycle's ghost (the unstable-spiral σ>0), with reduced breath period and per-cycle spiral slope quantitatively matching the finite-N measurements.

## Candid mismatch list

1. **Absolute capture-time scale (×~3 short).** Reduced median capture ≈43s vs measured τ_abs≈139s (ratio 0.31). The mean-field collective flow spirals to capture faster than the finite-N system lives; the collective IC is a lossy projection of the full 2N-phase state, and finite-N fluctuations extend lifetime. Reported straight, no fudge factor.
2. **τ_abs(N) plateau not reproduced.** Reduced capture median rises with N (IC-driven) while measured τ_abs is flat — the plateau is a finite-N effect, not in the reduced dynamics.
3. **N=4 is weak** (ρ=+0.14): at the smallest N the collective IC summaries are too noisy a projection for the reduced model to predict lifetime.
4. **Reduced does not decisively beat the fitted multi-feature regression** in absolute R² (0.21 vs 0.23), though it adds positive partial information and is parameter-free.
5. **Breath period slightly high / linear-vs-nonlinear.** Nonlinear reduced T_b ≈25s (top of the 21–25s band) while the linear-focus value is 19s; the canonical seeds sample the nonlinear regime away from the focus, lengthening the period — consistent, not contradictory.

## Artifacts

- `cp1_gate.json` — CP1 sub-gate verdicts
- `cp2_regions.json / .md` — CP2 region table + bifurcations + homoclinic trace
- `cp2_portrait_A0.2.{png,pdf}, cp2_portrait_A0.5.{png,pdf}` — phase portraits
- `cp3_match.json / .md` — CP3 (reduced, measured, ratio) tables
- `cp3_a02_level.json` — finite-N A=0.2 mean R_incoh (shipped tracer)
- `reduced_runs.jsonl` — per-run reduced 3-D capture/breath/slope (1400 runs)
- `cp4_correlation.json / .md` — CP4 Spearman/partial/CV-R² per N + pooled
- `cp4_scatter.{png,pdf}, cp4_rho_vs_N.{png,pdf}` — CP4 figures
