# Poisson-manifold probe + Weibull aging

Mechanism probe behind the finite-N chimera collapse-time **plateau** measured by the
merged campaign (`tools/chimera-campaign/`). With δω = 0 the two-population
Sakaguchi–Kuramoto system is in Watanabe–Strogatz territory: the dynamics is foliated
by constants of motion with an **invariant (not attracting)** Poisson/Ott–Antonsen
submanifold. This probe asks whether **distance from that submanifold predicts (CP2),
precedes (CP3), and — under intervention — causes (CP4) collapse**, and separately
characterizes the **non-exponential (aging) survival shape (CP5)**.

Everything here is **new** code that imports the shipped integrator primitives from
`../chimera-campaign/integrator.mjs` (verified bit-identical to `src/audio/chimera.ts`).
The shipped voice/supervisor are untouched and the campaign JSONL is read-only input.
Dynamics only — no audio path.

## Observable (CP1)

Per population, circular moments `Z_m = (1/N) Σ_j exp(i·m·θ_j)` and the Poisson-manifold
distance

```
D = Σ_{m=2}^{M} | Z_m − (Z_1)^m |²        (M = 4 by default)
```

On the submanifold `D = 0` in the continuum; at finite N it sits at an O(1/N) sampling
floor. The closure identity `Z_m = (Z_1)^m` and the Möbius push are from MMS 2009;
`moments.mjs` cites the provenance and `moments.test.ts` verifies it numerically.

## Layout

| File                   | Role                                                            |
| ---------------------- | --------------------------------------------------------------- |
| `moments.mjs`          | circular moments, D observable, Möbius push (core CP1)          |
| `moments.test.ts`      | CP1 gate (a) identity + O(1/N), (c) contrast — vitest           |
| `cp1_invariance.mjs`   | CP1 gate (b) invariance — dynamical, printed verdict + JSONL    |
| `trace.mjs`            | time-resolved integration with D(t); custom-IC runner (CP3/CP4) |
| `cp2_d0.mjs`           | regenerate D₀ for every campaign run from its seed              |
| `cp3_traces.mjs`       | quantile-stratified escape traces until collapse                |
| `cp4_intervention.mjs` | paired interventional seed families (causal test)               |
| `analysis.py`          | CP1 plot + CP2/CP3/CP4/CP5 stats/figures + `MANIFOLD_REPORT.md` |
| `manifold.config.json` | committed config (determinism)                                  |

`analysis.py` reuses the campaign's survival utilities (`kaplan_meier`,
`exp_mle_censored`, `km_median`) from `../chimera-campaign/analysis.py`.

## Run

```bash
# CP1 gate
npx vitest run tools/manifold-probe/moments.test.ts   # (a) + (c)
node tools/manifold-probe/cp1_invariance.mjs          # (b) — must PASS before CP2

# Data generation (deterministic, fast)
node tools/manifold-probe/cp2_d0.mjs                  # D₀ per campaign run
node tools/manifold-probe/cp3_traces.mjs             # escape traces
node tools/manifold-probe/cp4_intervention.mjs       # interventional families

# Analysis + report (needs numpy, scipy, matplotlib)
python3 tools/manifold-probe/analysis.py             # all CPs + MANIFOLD_REPORT.md
python3 tools/manifold-probe/analysis.py cp5 report  # a subset
```

Outputs land in `manifold_results/`. Every figure is reproducible from the logged
seeds + `manifold.config.json`.

## Headline result

The **manifold-escape hypothesis is not supported** over the tested range: D₀ does not
predict lifetime (CP2, all p≫0.05), D does not rise before collapse (CP3, paired ΔD≈0;
it _relaxes_ after the merger), and injecting D₀ does not causally shorten lifetime
(CP4, paired test null). Separately, the survival is genuinely **aging** (CP5: Weibull
shape k(N) rises from ≈1 to ≈2 at A=0.5 and ≈5 at A=0.2), which is the KM shoulder the
campaign flagged. See `manifold_results/MANIFOLD_REPORT.md`; interpretive
(WS/attractivity/integrability) claims are quarantined in its **PENDING-READING** list
for verification against Pikovsky–Rosenblum 2008 and Marvel–Mirollo–Strogatz 2009.
