# reduced-ode — Abrams 2008 reduced ODE ↔ annealMusic finite-N chimera

Maps the two-population Sakaguchi–Kuramoto chimera campaign onto the **reduced
order-parameter ODEs** of Abrams, Mirollo, Strogatz & Wiley, PRL **101**, 084103
(2008) (arXiv:0806.0594v2), and tests whether the paper's bifurcation structure
explains our finite-N measurements. New code only; all prior `*_results/` dirs are
read-only; deterministic; driven by the committed `reduced.config.json`.

## What it does

| CP | question | output |
|----|----------|--------|
| **CP1** | does the harness reproduce the *paper's own* results? (gate) | `cp1_gate.json` |
| **CP2** | which dynamical region are our corners A∈{0.2,0.5} at β=0.05 in? | `cp2_regions.{json,md}`, `cp2_portrait_A0.{2,5}.{png,pdf}` |
| **CP3** | do reduced breath period / spiral slope / capture time match finite-N? | `cp3_match.{json,md}`, `cp3_a02_level.json` |
| **CP4** | run-by-run: does the reduced model predict measured lifetime? | `reduced_runs.jsonl`, `cp4_correlation.{json,md}`, `cp4_scatter.*`, `cp4_rho_vs_N.*` |
| — | final synthesis | `reduced_results/REDUCED_REPORT.md` |

## The model

Conventions match the shipped integrator exactly (verified in campaign CP0):
`μ=(1+A)/2`, `ν=(1−A)/2`, `α=π/2−β`, `ω=0`.

- **2-D system** (invariant manifold ρ₁≡1, paper **Eq. 12**), state `(r=ρ₂, ψ=φ₁−φ₂)`.
- **3-D system** (full rotational reduction of **Eq. 10**), state `(ρ₁, ρ₂, ψ)`,
  derived here and verified to collapse to the 2-D RHS on ρ₁=1 to machine precision.
- Stationary chimeras: **Eqs 13–14**. Saddle-node series **Eq. 17**, Hopf series
  **Eq. 18**. Homoclinic located numerically (no closed form in the paper).

`reduced_core.py` carries all the math; `test_core.py` proves the 3-D→2-D collapse
and the Eqs (13–14) fixed-point consistency to ~1e-16.

## Run

```
pip install numpy scipy matplotlib      # if not present
python3 tools/reduced-ode/run_all.py     # full pipeline (CP1 gate stops on failure)
```

or step by step: `test_core.py`, `cp1_validate.py`, `cp2_classify.py`,
`node trace_a02_level.mjs > reduced_results/cp3_a02_level.json`, `compute_runs.py`,
`cp3_match.py`, `cp4_predict.py`, `make_report.py`.

## Inputs (read-only)

- `transient_results/cp2_features.jsonl` — per-run seed-mapped ICs (R_incoh0,
  R_sync0, Δφ₀) and measured t_abs for every A=0.5 campaign run (transient-tests CP2).
- `absorption_results/absorption_campaign.jsonl` — A=0.2 persistent seeds, re-traced
  by `trace_a02_level.mjs` via the shipped-identical tracer for the finite-N R_incoh.

## Headline results (β=0.05)

- **CP1 gate: PASS** against the paper (fixed points, Fig. 2 regimes, Eqs 17–18).
- A_SN=0.0948, A_H=0.2703, **A_hc=0.4096**. **P1 confirmed** (A=0.2 = stable chimera),
  **P2 confirmed** (A=0.5 = post-homoclinic, sync globally attracting).
- Reduced breath T_b≈25s (meas 21–25s), spiral slope −0.67 (meas −0.55…−0.77),
  capture ≈43s vs τ_abs 139s (×3 short — reported straight).
- CP4: reduced capture time predicts measured t_abs at pooled Spearman ρ=0.45
  (per-N to 0.60), with positive partial correlation beyond |Δφ₀| at every N. The
  N-dependence enters through the seed-ensemble ICs, not the reduced dynamics.
