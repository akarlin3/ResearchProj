# `theory/` — Plumbline: analytic hardening of the Minos findings

Two proved statements that turn the Minos-Core simulation findings (v1–v3) into theory, plus the
machine-checked derivations and the theory-vs-sim confirmation. Built on top of the v3 model; v1–v3
are untouched.

| file | what | gate |
|---|---|---|
| `THEORY_MODEL.md` | the analytic model, locked against the simulation code; provenance of every v2/v3 number | **GATE 0** (feasibility, PASS) |
| `gap_scaling.py` | **Theorem 1** — the gap scaling law `G = (1/6)·\|z*(λ)\|·γ`, all coefficients sympy-derived | **GATE 1** (PASS) |
| `gap_ci.py` | **CP1** — bootstrap / multi-seed CIs on `τ_stat, τ*, G`; error-bars the "opposite sides of 1" headline | **GATE 1** (halt-to-report; ROBUST) |
| `detectability.py` | **Theorem 2(ii)** — achievable bound `R_obs ≤ L·δ`, `L` identified from `U`, checked vs v3 | **GATE 2** (PASS) |
| `impossibility.md` | **Theorem 2(i)** — hidden channel undetectable (AUC = ½); pinned definition + partial-leak proposition; **DRAFTED, human proof-review required** | (not machine-verified) |
| `confirm.py` | **CP3** — theory vs v2/v3 for every compared quantity; halts if it does not reproduce | **GATE 3** (PASS, halt-able) |
| `figures/make_figures.py` | **CP4** — the two publication figures (vector PDF, deterministic) | — |
| `consistency_check.py` | **CP4** — re-derives the note's constants and checks every claim traces to a gate | **GATE 4** (PASS) |
| `plumbline.{md,tex}` | the note: intro → Thm 1 + proof + CP1 CIs → Thm 2 → confirmation tables → figures → discussion → positioning | **GATE 4** |

## Reproduce

```bash
python3 -m venv .venv-theory && .venv-theory/bin/pip install numpy scipy sympy mpmath matplotlib pytest
.venv-theory/bin/python theory/gap_scaling.py             # GATE 1 — symbolic, ~5 s
.venv-theory/bin/python theory/gap_ci.py                  # CP1 — bootstrap CIs, ~9 min (MINOS_FAST=1 to shrink)
.venv-theory/bin/python theory/detectability.py           # GATE 2 — bound + L, ~30 s
.venv-theory/bin/python theory/confirm.py                 # GATE 3 — theory vs v2/v3, a few min (MINOS_FAST=1)
.venv-theory/bin/python theory/figures/make_figures.py    # CP4 — figures (vector PDF)
.venv-theory/bin/python theory/consistency_check.py       # GATE 4 — final consistency, ~5 s
tectonic theory/plumbline.tex                             # render the note PDF
```

Scripts import the v3 model from `../minos-core/minos` (added to `sys.path`). Deterministic (seeded);
nothing is fit to the simulation.

## Results in one line

- **Theorem 1:** the decision–calibration gap is a leading-order scaling law in posterior skew,
  `G(γ) = (1/6)|z*(λ)|·γ + O(γ^{4/3})`, with `z*(λ)` the decision-boundary offset. `τ_stat` is
  first-order skew-insensitive (`a=0`); the whole first-order gap is decision-side. The clean line is
  `≈57%` of the gap at the operating point `γ=0.667` (the rest is a second-order coverage shrink).
  Reproduces the v2 sweep to `≤ 0.0065` and the default cell `τ*` to `0.003`.
- **CP1 (error bars):** `τ*` is a shallow optimum, so the published single seed is noisy; multi-seed
  CIs (`B=64`) show `τ_stat = 0.9639 [0.9636,0.9641] < 1 < τ* = 1.0514 [1.0493,1.0536]` and
  `G = 0.0876 [0.0855,0.0897] > 0` — the "opposite sides of 1" headline is **ROBUST** (not a lucky draw).
- **Theorem 2:** observable staleness is bounded-detectable (`R_obs ≤ k_under·W₁`), hidden staleness
  is undetectable by any label-free monitor (AUC = ½) — and the partial-leak proposition extends the
  observable/hidden split to *any* shift (the principled case for labeled repeatability spot-checks).
  Matches v3 (AUC obs `1.00`, hidden `0.50`; `corr(M,R)=+0.96`).
