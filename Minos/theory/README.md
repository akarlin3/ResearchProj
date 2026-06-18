# `theory/` — Plumbline: analytic hardening of the Minos findings

Two proved statements that turn the Minos-Core simulation findings (v1–v3) into theory, plus the
machine-checked derivations and the theory-vs-sim confirmation. Built on top of the v3 model; v1–v3
are untouched.

| file | what | gate |
|---|---|---|
| `THEORY_MODEL.md` | the analytic model, locked against the simulation code; provenance of every v2/v3 number | **GATE 0** (feasibility, PASS) |
| `gap_scaling.py` | **Theorem 1** — the gap scaling law `G = (1/6)·\|z*(λ)\|·γ`, all coefficients sympy-derived | **GATE 1** (PASS) |
| `detectability.py` | **Theorem 2(ii)** — achievable bound `R_obs ≤ L·δ`, `L` identified from `U`, checked vs v3 | **GATE 2** (PASS) |
| `impossibility.md` | **Theorem 2(i)** — hidden channel undetectable (AUC = ½); full data-processing proof + finite-sample statement | **proved** (see `impossibility_check.py`) |
| `impossibility_check.py` | **Theorem 2(i)** machine check — exact pathwise invariance of `O`/`M` across the sweep; fresh-vs-stale AUC = ½ (bootstrap CI `[0.5,0.5]`) | **GATE 2(i)** (PASS) |
| `confirm.py` | **CP3** — theory vs v2/v3 for every compared quantity; halts if it does not reproduce | **GATE 3** (PASS, halt-able) |
| `plumbline.md` | the note: intro → Thm 1 + proof → Thm 2 → confirmation table → discussion → positioning | **GATE 4** |

## Reproduce

```bash
python3 -m venv .venv-theory && .venv-theory/bin/pip install numpy scipy sympy matplotlib
.venv-theory/bin/python theory/gap_scaling.py         # GATE 1 — symbolic, ~5 s
.venv-theory/bin/python theory/detectability.py       # GATE 2 — bound + L, ~20 s
.venv-theory/bin/python theory/impossibility_check.py # GATE 2(i) — invariance + AUC=1/2, ~30 s (MINOS_FAST=1)
.venv-theory/bin/python theory/confirm.py             # GATE 3 — theory vs v2/v3, ~30 s (MINOS_FAST=1 to shrink)
```

Scripts import the v3 model from `../minos-core/minos` (added to `sys.path`). Deterministic (seeded);
nothing is fit to the simulation.

## Results in one line

- **Theorem 1:** the decision–calibration gap is a leading-order scaling law in posterior skew,
  `G(γ) = (1/6)|z*(λ)|·γ + O(γ^{4/3})`, with `z*(λ)` the decision-boundary offset. `τ_stat` is
  first-order skew-insensitive (`a=0`); the whole first-order gap is decision-side. Reproduces the v2
  sweep to `≤ 0.0065` and the default cell `τ*` to `0.003`.
- **Theorem 2:** observable staleness is bounded-detectable (`R_obs ≤ k_under·W₁`), hidden staleness
  is undetectable by any label-free monitor (AUC = ½) — the principled case for labeled repeatability
  spot-checks. Matches v3 (AUC obs `1.00`, hidden `0.50`; `corr(M,R)=+0.96`).
