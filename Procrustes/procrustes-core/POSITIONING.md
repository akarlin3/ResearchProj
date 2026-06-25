# Procrustes — positioning & novelty gate

What Procrustes **is not**, and why the surviving wedge is distinct. This file is
the durable record of the novelty gate; the gate tests enforce its claims.

## Verdict trail

- **Original wedge (DEAD).** "Bi-exp-on-non-bi-exp breaks *conditional coverage on
  the D\* tercile axis*." → already done by **Gauge §altmodel + §envelope**, which
  regenerates from the *same* non-bi-exp generators (gamma/log-normal/stretched/
  tri-exp), shows marginal restored (D\*_eff 0.898) and high-D\*_eff tercile-
  conditional coverage fails (0.796/0.800), and **folds it into the same
  identifiability wall** ("a property of the perfusion-estimation problem more
  broadly… survives a generator Eq.(IVIM) does not even contain"), with continuity
  at the bi-exp limit. That fires the original refute condition. **Restated on the
  D\* axis, Procrustes = Gauge.**

- **Repositioned wedge (SURVIVES).** Move *off* the D\* identifiability axis. In
  every Lattice non-bi-exp family only the **perfusion** compartment is made
  non-bi-exp; the tissue term `(1-f)·exp(-b·D)` stays intact, so **D is an exact,
  well-identified ground-truth parameter** — the one Gauge says to *trust*. The
  surviving claim: **bi-exp misspecification breaks the conditional coverage of D
  along the perfusion-departure axis, inside the well-identified D\* regime** —
  a failure mode Gauge's wall cannot explain and that contradicts Gauge's triage
  rule. Refute-first probes (single- and multi-seed) cleared all three refute
  conditions for the heavy-tail (stretched) family.

## Separation from neighbours

| neighbour | what it owns | how Procrustes differs |
|---|---|---|
| **Gauge** (within-model identifiability) | high-**D\*** conditional-coverage wall; bi-exp model assumed *correct*; misspec probe folded into the same wall; "trust D" | Procrustes breaks coverage of the **trusted D**, in the **well-identified D\*** regime, *because* the model is **wrong** — orthogonal axis, opposite parameter |
| **Lei et al. 2018** (CP marginal robustness) | distribution-free **marginal** coverage under any/biased/misspecified base | Procrustes *confirms* this (marginal holds) and is explicitly **not** a marginal claim |
| **Barber–Candès–Ramdas–Tibshirani 2021** | impossibility of distribution-free **conditional** coverage | the theoretical reason the departure axis *can* break; Procrustes is its IVIM-misspecification instance |
| **Wang–Tamir–Bush 2026** (qMRI misspecification, **ASL**) | misspecified-CRB + two-subset variance-consistency test | Procrustes' diagnostic is **CP-coverage-native** (residual structure ranks coverage-failure), and **IVIM/diffusion**, not ASL; different mechanism |
| **IVIM model selection** (AIC/BIC/F-test, mono/bi/tri/stretched) | choosing model *order* by fit criterion | Procrustes uses the observable to **predict where conditional coverage of the trusted parameter fails**, not to pick an order |
| **Casali et al. 2025** (IVIM-UQ) | aleatoric/epistemic UQ, calibration curves/CRPS | no conformal coverage, no non-bi-exp truth, no misspecification diagnostic |

**External cell** (CP-coverage × IVIM × non-bi-exp-truth misspecification diagnostic):
not located in prior art. The binding constraint was internal (Gauge); the
repositioned axis clears it.

## What hardening established (Phases 1–3; GATE B/C/D, see RESULTS.md)

The CP0 probe was hardened under three pre-registered gates. The well-identified-D\*
subset was **re-defined to match Gauge's identifiable region exactly** (true-D\*
terciles via `np.quantile([1/3,2/3])` + `np.digitize`); load-bearing gaps carry
two-level cluster-bootstrap CIs (16-seed headline).

- **GATE B (distinctness) — PASS.** Marginal holds 0.909 [0.900,0.919]; conditional
  gap 0.105 [0.093,0.117]; **well-ID-D\* gap 0.148 [0.131,0.164]**, intensifying to
  0.196 [0.170,0.220] in the strict bottom tercile — while Gauge's OWN high-D\* wall
  shows gap 0.014 [0.001,0.027]. The gradient is **opposite** to Gauge's
  identifiability wall (gap largest where Gauge says "trust D", near zero in Gauge's
  ill-posed regime). Bias 7.7× [5.8,9.5] (high-b aliasing).
- **GATE C (diagnostic scope) — heavy-tail DETECTOR, not universal.** Heavy-tail AUC
  0.684 [0.673,0.694] beats the naive drift monitor 0.550 [0.543,0.556]; pure
  dispersion 0.578 is **below the pre-registered 0.60 floor** (near-hidden, like the
  monitor); tri-exp 0.627 is detectable-but-harmless.
- **GATE D (robustness) — PASS with one honest boundary.** Survives 11/12 conditions
  (SNR 35–100, intensifying with SNR; n 300–2000; b-schemes). **Fails at SNR 25**
  (noise-dominated) — reported, not buried.

→ the wedge is **mechanism-specific (high-b aliasing)**, not generic.

## Pre-registered refute conditions (enforced by gate tests)

- **R1** conditional gap ≈ 0 → no misspecification-driven failure → dead.
- **R2** the gap vanishes inside the well-identified D\* subset → it was Gauge → dead.
- **R3** the bi-exp-limit placebo row is itself broken → artefact, not misspec → dead.
- **boundary** tri-exp must stay null; bi-exp limit must be exact continuity.

## Honest open risk

The **diagnostic** is moderate (AUC ≈ 0.684) and **family-dependent** — strong for
the heavy-tail channel (beats the naive monitor 0.550), near-hidden for pure
dispersion (0.578, below the pre-registered 0.60 floor). Procrustes' contribution
is the *coverage-failure* characterisation + a working IVIM misspecification
diagnostic that beats the hidden-channel baseline for the heavy-tail channel,
**not** a universal detector. The separation also has a stated **SNR boundary**:
it fails at SNR 25 (noise-dominated), holding for SNR ≥ 35. Pushing the diagnostic
(structured-residual GLR, multi-segment consistency) and an in-vivo demonstration
are future work, not claimed here.

## IP / clean-room

Fully synthetic. Ground truth is the Lattice DRO (seed-generated, no data files);
no pancData3 or any clinical data is touched. `[confirm venue]` before submission.
