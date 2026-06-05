# Poisson-manifold probe + Weibull aging — MANIFOLD REPORT

Mechanism probe behind the finite-N collapse-time plateau: does distance from the Ott–Antonsen/Poisson submanifold **predict (CP2)**, **precede (CP3)**, and **— under intervention — cause (CP4)** collapse of the δω=0 two-population Sakaguchi–Kuramoto chimera? Plus a **Weibull/aging (CP5)** characterization of the survival shoulder. Built on the merged campaign; shipped integrator reused bit-identical; existing campaign data read-only.

Every physics claim below is tagged **VERIFIED-NUMERICALLY** or **PENDING-READING**.

## Headline synthesis

- **CP1 (precondition): PASS.** The Poisson-manifold distance D is a valid, numerically-verified observable and the manifold is invariant under this coupling.
- **CP2 (predicts?): NO.** D₀ does not predict lifetime — 0 of 7 A=0.5 points show a significant negative D₀–lifetime correlation (|ρ|≤0.08; all p≫0.05); same at A=0.2.
- **CP3 (precedes?): NO.** Event-aligned ⟨D⟩ is flat at the finite-N floor through collapse (paired ΔD≈0, p≈0.6) and _relaxes_ after the merger; D does not lead the order-parameter collapse.
- **CP4 (causes?): NO / at most WEAK.** Injecting manifold distance into the incoherent population leaves τ̂ essentially flat; the family-paired test is null (Δτ≈−0.8 s, p≈0.7).
- **CP5 (aging): YES.** The survival shoulder is real — the Weibull shape k(N) rises from ≈1 (small N) to ≈2 (A=0.5) / ≈5 (A=0.2) at large N (increasing hazard).

**Bottom line (VERIFIED-NUMERICALLY):** the _manifold-escape_ hypothesis — that distance from the Poisson submanifold predicts/precedes/causes collapse — is **not supported** in this δω=0 two-population system over the tested range. Distance from the manifold is neither predictive (CP2), nor precedent (CP3), nor causal (CP4). Separately, the collapse-time distribution is genuinely **aging** (CP5: k>1, rising with N), which is what produces the KM shoulder the campaign flagged. The null mechanism result is itself consistent with the manifold being _invariant but not attracting_ (see PENDING-READING #2): finite-N collapse appears to be a diffusive finite-size effect along the leaves rather than an escape _off_ the leaf — a hypothesis to confirm against the WS/OA literature.

## CP1 — observable + physics validation gate (VERIFIED-NUMERICALLY)

Observable: per population, D = Σ\_{m=2}^{M} |Z_m − (Z_1)^m|² with M=4 (Z_m = ⟨e^{imθ}⟩). On the Poisson submanifold D=0 in the continuum, O(1/N) at finite N.

| Gate           | Test                                                                                             | Verdict                |
| -------------- | ------------------------------------------------------------------------------------------------ | ---------------------- |
| (a) identity   | Z_m=(Z_1)^m for a Möbius-pushed grid (machine precision); D=O(1/N) for random-uniform pre-images | **PASS** (unit-tested) |
| (b) invariance | on-manifold IC evolved 300 s: max D_incoh=0.191, mean=0.048 ≪ O(1); no secular growth            | **PASS** (JS gate)     |
| (c) contrast   | bimodal/clustered population ⇒ large D (≫ on-manifold floor)                                     | **PASS** (unit-tested) |

**CP1 gate: PASS.** The closure identity holds numerically and the sinusoidal mean-field coupling preserves the manifold in this implementation — the precondition for CP2–CP4. Figure: `cp1_invariance.png`.

## CP2 — retrospective correlation: D₀ vs lifetime

Spearman ρ(D₀_incoh, ln lifetime), uncensored runs (ρ<0 ⇒ larger initial distance predicts shorter life):

| A   | N   | uncensored | censored | ρ      | p       |
| --- | --- | ---------- | -------- | ------ | ------- |
| 0.5 | 4   | 200        | 0        | -0.081 | 2.5e-01 |
| 0.5 | 8   | 200        | 0        | -0.041 | 5.6e-01 |
| 0.5 | 16  | 200        | 0        | -0.039 | 5.8e-01 |
| 0.5 | 24  | 200        | 0        | +0.072 | 3.1e-01 |
| 0.5 | 32  | 200        | 0        | -0.049 | 4.9e-01 |
| 0.5 | 48  | 200        | 0        | -0.072 | 3.1e-01 |
| 0.5 | 64  | 200        | 0        | -0.024 | 7.3e-01 |
| 0.2 | 4   | 99         | 1        | +0.069 | 4.9e-01 |
| 0.2 | 8   | 96         | 4        | -0.148 | 1.5e-01 |
| 0.2 | 16  | 100        | 0        | -0.007 | 9.5e-01 |
| 0.2 | 24  | 100        | 0        | -0.138 | 1.7e-01 |
| 0.2 | 32  | 100        | 0        | -0.015 | 8.8e-01 |
| 0.2 | 48  | 100        | 0        | -0.008 | 9.3e-01 |
| 0.2 | 64  | 100        | 0        | -0.010 | 9.2e-01 |

**CP2 headline:** D₀ is a significant negative predictor of lifetime at 0/7 A=0.5 points. ρ(N) trend in `cp2_rho_vs_N.png`; scatter + KM terciles in `cp2_scatter_*.png` / `cp2_km_terciles_*.png`.

## CP3 — escape precedence + lead time

Event-aligned, per-run paired test of ⟨D_incoh⟩ (56 runs): pre-collapse [−25,−10]s = **0.146**, peri [−3,0]s = **0.125**, post-merger = **0.088**. Paired ΔD(peri−pre) median = **-0.004** (Wilcoxon p=5.7e-01).

**Verdict: NO precedence** — D is flat at the finite-N floor through collapse and relaxes after the merger; it does not lead the R-collapse. (Per-run 'lead times' are a breathing-fluctuation artifact; see `cp3_escape.md`.) Figures `cp3_aligned.png`, `cp3_spaghetti.png`.

## CP4 — interventional causal test

Family-paired ε-level intervention (N=16, A=0.5), injected D₀ over [0.183, 0.357]:

| ε      | realized D₀ (med) | τ̂ (s) | KM median (s) |
| ------ | ----------------- | ----- | ------------- |
| 0.00\* | 0.183             | 65.9  | 54.3          |
| 0.10   | 0.214             | 67.5  | 50.8          |
| 0.20   | 0.275             | 66.5  | 48.3          |
| 0.30   | 0.314             | 64.6  | 45.3          |
| 0.40   | 0.352             | 65.3  | 50.1          |
| 0.50   | 0.357             | 65.9  | 46.4          |

Causal readouts: level-trend ρ=**-0.600** (p=2.1e-01); paired Δτ median=**-0.8 s** (Wilcoxon p=7.1e-01); pooled run-level ρ=-0.102 (p=5.3e-02).

See `cp4_intervention.md` for the causal verdict; figure `cp4_tau_vs_d0.png`.

## CP5 — Weibull / aging characterization

Weibull shape k (k≈1 memoryless → k>1 aging) and AIC model choice:

| A   | N   | k [95% CI]        | best model |
| --- | --- | ----------------- | ---------- |
| 0.5 | 4   | 0.82 [0.74, 0.91] | weibull    |
| 0.5 | 8   | 1.26 [1.15, 1.42] | weibull    |
| 0.5 | 16  | 1.49 [1.34, 1.69] | weibull    |
| 0.5 | 24  | 1.61 [1.45, 1.81] | weibull    |
| 0.5 | 32  | 1.52 [1.38, 1.72] | weibull    |
| 0.5 | 48  | 1.85 [1.65, 2.16] | weibull    |
| 0.5 | 64  | 2.19 [1.96, 2.52] | weibull    |
| 0.2 | 4   | 0.64 [0.50, 1.44] | lognorm    |
| 0.2 | 8   | 0.57 [0.47, 0.78] | lognorm    |
| 0.2 | 16  | 2.19 [1.82, 2.73] | weibull    |
| 0.2 | 24  | 2.63 [2.17, 3.30] | weibull    |
| 0.2 | 32  | 3.58 [2.93, 4.57] | weibull    |
| 0.2 | 48  | 5.31 [4.51, 6.61] | weibull    |
| 0.2 | 64  | 4.75 [3.86, 7.74] | weibull    |

Full table `cp5_weibull.md`; k(N) `cp5_kN.png`; fits overlaid `cp5_km_fits_*.png`.

## PENDING-READING — interpretations requiring P&R 2008 / MMS 2009

The numerics above stand on their own. The following **interpretive** statements are written as HYPOTHESES and must be checked against Pikovsky & Rosenblum, PRL 101, 264103 (2008) and Marvel, Mirollo & Strogatz, Chaos 19, 043104 (2009) before appearing as conclusions in the paper:

1. **WS-constant framing.** That D measures drift in the Watanabe–Strogatz constants of motion (vs. a generic moment defect) presumes the WS reduction applies to each finite-N population here — verify the WS variables/constants decomposition (MMS 2009).
2. **Manifold attractivity.** CP1(b) shows the manifold is _invariant_ numerically. Whether it is _attracting_, _neutrally stable_, or _repelling_ (and hence whether finite-N noise should drive escape) is a stability claim — OA/WS theory says the OA manifold is invariant and, for δω=0, the dynamics is _partially integrable_ (P&R 2008), NOT attracting. Confirm before asserting any 'escape' is dynamically driven rather than diffusive.
3. **Partial integrability / foliation.** The picture of motion foliated by conserved quantities with the chimera living on the Poisson leaf is P&R 2008 language — confirm it describes the _two-population_ model with this coupling.
4. **Continuum limit of the mechanism.** Any claim that the CP2/CP4 trend extrapolates to N→∞ (or vanishes there) needs the OA continuum equations, not just the finite-N regression.

— end of report —
