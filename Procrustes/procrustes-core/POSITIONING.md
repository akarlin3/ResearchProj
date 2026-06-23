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

## What hardening established (8-seed probe)

- **Heavy-tail (stretched) — survives, tight CIs:** marginal 0.894 [0.885,0.903];
  conditional gap 0.126 [0.116,0.136]; **well-ID-D\* gap 0.172 [0.162,0.183]**
  (larger than marginal); bias 7.4×; diagnostic AUC **0.67** (vs Gauge's 0.501).
- **Tri-exp — null:** gap −0.008 [−0.021,0.005] (faster pool decays off high-b).
- **Log-normal — weak break, hidden diagnostic:** gap 0.046; AUC ≈ 0.52.

→ the wedge is **mechanism-specific (high-b aliasing)**, not generic.

## Pre-registered refute conditions (enforced by gate tests)

- **R1** conditional gap ≈ 0 → no misspecification-driven failure → dead.
- **R2** the gap vanishes inside the well-identified D\* subset → it was Gauge → dead.
- **R3** the bi-exp-limit placebo row is itself broken → artefact, not misspec → dead.
- **boundary** tri-exp must stay null; bi-exp limit must be exact continuity.

## Honest open risk

The **diagnostic** is moderate (AUC ≈ 0.67) and **family-dependent** — strong for
the heavy-tail channel, hidden for pure dispersion. Procrustes' contribution is
the *coverage-failure* characterisation + a *first* working IVIM misspecification
diagnostic that beats the hidden-channel baseline, **not** a universal detector.
Pushing the diagnostic (structured-residual GLR, multi-segment consistency)
is future work, not claimed here.

## IP / clean-room

Fully synthetic. Ground truth is the Lattice DRO (seed-generated, no data files);
no pancData3 or any clinical data is touched. `[confirm venue]` before submission.
