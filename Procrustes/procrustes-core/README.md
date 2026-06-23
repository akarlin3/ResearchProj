# procrustes-core

The clean-room core of **Procrustes**: *misspecification-aliasing of a calibrated
error bar*. When a **bi-exponential** IVIM model is fit on **non-bi-exponential**
truth, distribution-free **marginal** coverage survives (Lei 2018) — but the
**conditional** coverage of the *well-identified tissue-diffusion map D* breaks
along the latent misspecification axis, in exactly the regime Gauge's triage rule
calls safe ("trust D"). Ground truth is the **Lattice** digital reference object;
fully synthetic, seed-generated, no data files. `[confirm venue]`

## First-class claim

> Fitting bi-exp on non-bi-exp truth leaves **marginal** coverage of D intact yet
> breaks its **conditional** coverage — a failure that **survives inside the
> well-identified D\* subset** and is therefore mechanistically **distinct from
> Gauge's within-model identifiability wall**. The mechanism is **high-b
> aliasing**: heavy-tailed perfusion departures leak into the high-b tissue slope
> and bias D̂.

**Throughline:** misspecification is invisible from inside the model — another
observable/hidden instance.

## Scope (pre-registered, honest)

| departure family | knob | expectation | why |
|---|---|---|---|
| stretched-exp | `beta` | **breaks** D-coverage | heavy tail aliases into high-b slope |
| log-normal dispersion | `cv` | **weak** break, diagnostic **hidden** | mild high-b leakage; Gauge's hidden channel |
| tri-exponential (faster 3rd pool) | `g` | **null** | a *faster* pool decays away from high-b → D unbiased |

The tri-exp **null** is a feature, not a miss: it shows the wedge is
mechanism-specific, not generic "misspecification breaks everything."

## The CP0 separation (what the gate proves)

For one departure family, observing the **same voxel population** under
increasing departure (base `(D,D*,f)` seed-locked across knob values):

1. **(a) marginal holds** — one departure-blind conformal radius → ~nominal pooled coverage.
2. **(b) conditional fails** — coverage of D degrades monotonically with departure;
   an oracle that knows the departure recalibrates it away (so the axis is *latent*).
3. **(c) distinct from Gauge** — the gap **persists, undiminished, inside the
   well-identified D\* subset** (bottom-2 terciles).

**Refute conditions** (any one ⇒ the wedge is dead): **R1** conditional gap ≈ 0;
**R2** the gap vanishes in the well-identified D\* subset (it was Gauge all along);
**R3** the bi-exp-limit placebo row is itself broken (artefact, not misspecification).
These are encoded as gate tests.

## The diagnostic (the risk-bearing half)

Can a deployer *detect* the misspecification — and thus that D's conditional
coverage is at risk — from observables alone? Procrustes reads the bi-exp fit's
**residual structure** (reduced-χ² / lag-1 autocorrelation / longest same-sign
run). It clears chance for the heavy-tail channel (beating Gauge's AUC≈0.5
naive-transfer monitor) while confirming pure dispersion stays a hidden channel —
distinct from model-order GoF (AIC/BIC) and from the misspecified-CRB variance
test (Wang–Tamir–Bush 2026, ASL).

## Install & run

```bash
pip install -e .            # numpy, scipy
# Lattice DRO must be importable: pip install -e ../../Lattice,
# or `git submodule update --init Lattice`, or set LATTICE_PATH.
pytest                      # gate tests (CP0 separation, boundaries, diagnostic)
python experiments/run_cp0.py        # the headline numbers (-> RESULTS.md)
```

## Sanity gates (must hold)

- **CP0-a** marginal coverage within ±0.03 of nominal under departure-blind calibration.
- **CP0-b** stretched conditional gap > 0.06 (R1 guard).
- **CP0-c** stretched well-identified-D\* gap > 0.03 (R2 guard — survives Gauge).
- **boundary** tri-exp gap ≈ null; bi-exp limit is exact continuity.
- **diagnostic** heavy-tail AUC clears chance and exceeds the dispersion channel.

## Layout

```
procrustes-core/
  procrustes/
    deps.py         locate the Lattice DRO (env / submodule / dev fallback)
    seeding.py      GLOBAL_SEED, make_rng
    config.py       ProcrustesConfig + the three DepartureFamily definitions
    generators.py   Lattice cohorts + bi-exp NLLS fit + residual features
    conformal.py    split-conformal radius / coverage / splits
    separation.py   the CP0 three-part separation (a)/(b)/(c) + bias + R2
    diagnostic.py   observable misspecification diagnostic (AUC, |D-err| rank corr)
  experiments/run_cp0.py   seeded multi-seed driver -> RESULTS.md
  tests/          gate tests encoding the refute conditions
  POSITIONING.md  novelty gate: what this is NOT (vs Gauge, Lei, Barber, ASL-MCRB)
```

## Dependency

Procrustes owns **no** ground-truth generator. The non-bi-exp truth is the
clean-room **Lattice** DRO (seed-generated; no data files). Procrustes contributes
the misspecification-coverage *analysis* and *diagnostic*, not the signals.
