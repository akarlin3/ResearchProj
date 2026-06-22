# Lethe — constrained-validation / honest-limitation results for the IVIM uncertainty program

**Lethe** is the home for *honest-limitation* findings: results where a validation is run
faithfully and the verdict is a sharply-scoped negative. Its first (and currently only)
portion is **Echo**, which arrived here by verdict — a speculative gated build whose CP3
validation routed to Lethe. The package, code, and history below are the **Echo portion**
(the directory was renamed `Echo/ → Lethe/` once the verdict made Echo part of Lethe; the
method keeps the name *Echo*). The manuscript is in [`paper/lethe.tex`](paper/lethe.tex).

---

## The Echo portion — repeatability as a *scale* check on conformal intervals

**Status: speculative gated build — VERDICT RENDERED. CP0 legitimacy PASSED; CP1 scaffold
up; CP2 data gate PASSED (real ACRIN-6698, n=76); CP3 → LETHE (constrained validation).**
On real data the conformal interval is ~4× too narrow to cover real test–retest variation of
D (coverage 0.263 [0.158, 0.355] vs 0.755 target) — a valid honest-limitation verdict, not a
failure. See [`LETHE.md`](LETHE.md), [`paper/lethe.tex`](paper/lethe.tex), and
[`results/RESULTS_VALIDATION.md`](results/RESULTS_VALIDATION.md).

Echo asks one ground-truth-free question about a deployed conformal IVIM interval:

> Is the interval the right **size** to capture a measurement's own irreproducibility?

It answers with **test–retest interval coverage** — does one scan's parameter estimate
fall inside the *other* scan's deployed conformal interval — reported per parameter with a
**BCa bootstrap CI**, on public same-day scan–rescan data (ACRIN-6698, n≈76).

## What makes Echo legitimate (the two hard constraints)

1. **Precision, not accuracy — provably.** Writing each estimate as
   `est = θ_true + bias + ε`, the test–retest discrepancy `Δ = est_B − est_A = ε_B − ε_A`
   **cancels bias exactly**. Echo's statistic is therefore invariant to any systematic
   error common to both scans (demonstrated in `echo_repeat/harness.py` and `tests/`). Echo
   certifies an interval is correctly **sized to measurement noise**; it is **blind to
   accuracy/bias** and makes **no ground-truth coverage claim**. A perfectly
   measurement-scaled 90% interval is *expected* to show ≈76% test–retest coverage, **not
   90%** — the derivable gap between accuracy-coverage and repeat-coverage is exactly why
   the two cannot be conflated.

2. **Distinct from and beyond Gauge — provably.** Gauge's published check (paper §4.2.2;
   the "§3.7" of the protocol) is a Spearman **rank** correlation: does the interval *width*
   widen where scan–rescan scatter is larger (D `r=+0.60`, D\* null)? Echo measures **scale**
   (a coverage rate). These are mathematically independent: rescaling every width by a
   constant leaves Spearman invariant but moves Echo's coverage arbitrarily. Gauge asks
   *"does the band widen where noise is larger?"*; Echo asks *"is the band the right size to
   capture that noise?"* — the question Gauge explicitly declined.

If on real data the coverage signal collapsed to Gauge's rank check, saturated (no scale
content), or under-scaled, Echo routes to **Lethe** (honest-limitation regime) — a valid
verdict. See `VERIFICATION.md` for the locked PASS/FAIL thresholds.

## Reverb — the constructive counterexample (SOLID)

On real data, "precision ≠ coverage" can only be **argued** (no ground truth). **Reverb**
(`echo_repeat/reverb.py`) *shows* it on synthetic ground truth. It is **SOLID** — built only on
**Lattice** (synthetic ground-truth IVIM cohorts, read-only) and **Caliper** (estimator +
conformal ruler, read-only); it depends on no upstream paper.

Reverb draws a Lattice cohort, acquires it **twice from one truth** (test–retest), reduces each to
a whole-tumor **ROI-mean** (Lethe's region level; the √n_vox precision boost), fits bi-exp, and
deploys a bi-exp-calibrated conformal interval. Because the truth is known it measures *both*
repeatability and coverage-of-truth per true-D\* regime. The headline: under realistic
perfusion-model mismatch (dispersed perfusion fit as bi-exp), the perfusion fraction **f at low
D\*** is **excellently repeatable yet badly under-covers the truth** (coverage ≈0.61 [BCa
0.57, 0.64] vs a matched correctly-specified control at ≈0.80, with **identical** repeatability) —
precision blind to a structural bias, visible only because truth is known. A matched bi-exp
control and a family×ROI sensitivity surface show the divergence is the *mismatch*, not IVIM per se.

**Scope (load-bearing):** a synthetic *possibility-and-mechanism* proof — the divergence *can*
occur in IVIM and here is *why*; it does **not** quantify any real-world miscalibration magnitude.
Run `python scripts/run_reverb.py`; the `consistency.py` gate enforces the counterexample.

## Layout

```
echo_repeat/
  statistic.py    the core: test-retest coverage, standardized-residual scale check,
                  analytic reference, Spearman (for contrast), numpy-only BCa bootstrap
  harness.py      synthetic test-retest generator + the CP1 method self-test (abstract scalar)
  reverb.py       the constructive precision-vs-coverage counterexample (SOLID): region-level
                  test-retest on Lattice + Caliper conformal, repeatability vs known-truth coverage
  invivo.py       IVIM forward + segmented plug-in fit + Caliper-conformal deployer +
                  real test-retest signal loader
  provenance.py   download-on-demand provenance manifest writer (mirrors Gauge's posture)
  _paths.py       read-only import chokepoint (Caliper, Lattice SOLID; Gauge/Fashion/Minos PROVISIONAL)
scripts/
  run_harness.py     CP1 method self-test (SOLID) -> results/RESULTS_HARNESS.*
  run_reverb.py      Reverb constructive counterexample (SOLID) -> results/RESULTS_REVERB.*
  fetch_invivo.py    CP2 download-on-demand fetch (reuses Gauge's data-handling template)
  run_validation.py  CP3 real-data validation, locked gate -> results/RESULTS_VALIDATION.*
paper/            CP4 manuscript (ebgaramond + microtype), built PASS-only
tests/            unit tests for the statistic, harness, and fit
ASSUMPTIONS.md    pinned Fashion/Minos/Gauge inputs + SOLID/PROVISIONAL split
PROMOTION.md      promotion (PASS) / Lethe-fold / Reverb-fold paths
VERIFICATION.md   the CP gates and locked thresholds
reproduce.sh      one-command re-validation (CP1 -> CP2 -> CP3 -> CP4)
```

## Reproduce

```bash
pip install -e .            # numpy-only core
bash reproduce.sh          # CP1 self-test always; CP2/CP3 run iff data present; CP4 PASS-only
```

## Data & IP

Echo's own repo and history are **synthetic + open only**. Real repeatability data
(ACRIN-6698, CC-BY-4.0, DOI 10.7937/tcia.kk02-6d95) is **download-on-demand**: no pixel
data is committed, only a provenance manifest. Echo redistributes nothing and mirrors
Gauge's in-vivo data posture exactly.

## License

MIT (see `LICENSE`). Echo imports Caliper (MIT) read-only and depends on Gauge/Fashion/Minos
by read-only import; those dependencies are **PROVISIONAL** (in review) — see `ASSUMPTIONS.md`.
