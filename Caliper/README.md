# Caliper

**An IVIM uncertainty-quantification calibration toolkit.**

Caliper is a small, reviewer-oriented Python package for *measuring and
correcting* the calibration of uncertainty estimates in intravoxel incoherent
motion (IVIM) diffusion MRI. It ships three composable pieces:

1. **`caliper.metrics`** — a **model-agnostic calibration ruler** (numpy-only):
   coverage, quantile ECE, sharpness, pinball/interval score, and
   group-conditional coverage. It knows nothing about IVIM or any particular
   estimator — feed it true values and predicted quantiles.
2. **`caliper.estimator_maf`** — a conditional **masked-autoregressive-flow**
   posterior over `(D, f, D*)` given the multi-b signal decay (optional torch
   extra).
3. **`caliper.conformal`** — a **split-conformal / CQR** wrapper that
   coverage-corrects *any* estimator exposing `predict_quantiles`.

All data is **synthetic and PHI-free**, generated in-repo with fixed seeds
(`caliper.forward`). There are no clinical-data dependencies.

> Every number in this README is produced by `examples/demo.py` in this repo
> (fixed seeds). Re-run it to reproduce them exactly.

---

## Install

```bash
# core ruler + forward model + conformal wrapper (numpy only)
pip install -e .

# add the MAF estimator (pulls in torch)
pip install -e ".[estimator]"
```

Python 3.10–3.12. The core (`metrics`, `forward`, `conformal`) is numpy-only;
torch is required only for `estimator_maf`.

---

## Quickstart

```bash
python examples/demo.py
```

This runs the full pipeline — synthetic IVIM → MAF posterior → split-conformal →
calibration scorecard — under a realistic **deployment shift** (the flow is
trained for high-SNR fitting but evaluated at lower SNR, where model-based IVIM
UQ is known to be over-confident).

### The model-agnostic ruler API

```python
import numpy as np
from caliper import metrics as M

# y_true:   (n, n_params)
# q_pred:   (n, n_params, n_levels)   <- from any estimator's predict_quantiles
# q_levels: (n_levels,)               ascending in (0, 1)
scores = M.score_quantiles(y_true, q_pred, q_levels, alpha=0.10,
                           param_names=["D", "f", "Dstar"],
                           conditioning=y_true)      # tercile-conditional probe
print(M.format_scorecard(scores))
# each ParamScore exposes: coverage, coverage_gap, ece, sharpness,
# mean_pinball, mean_interval_score, conditional (per-tercile coverage)
```

### Conformal wrapper over any estimator

```python
from caliper.conformal import SplitConformalQuantile

cq = SplitConformalQuantile(q_levels).calibrate(q_cal, y_cal)
q_corrected = cq.apply(q_test)   # coverage-corrected quantiles, same shape
```

---

## Results (from `examples/demo.py`)

Nominal central coverage **0.900** (90% intervals, α = 0.10). Held-out synthetic
test set, deployment shift train-SNR 60 → test-SNR 25.

### Raw MAF is over-confident (the known model-based UQ result)

| param | coverage | gap | ECE | sharpness |
|-------|---------:|------:|------:|----------:|
| D     | 0.528 | −0.372 | 0.139 | 0.281 |
| f     | 0.550 | −0.350 | 0.133 | 0.072 |
| D\*   | 0.555 | −0.345 | 0.124 | 26.86 |

The raw posterior intervals are far too tight: ~53–56% empirical coverage
against a 90% target. This is expected and is reported honestly — **not** tuned.

### Split-conformal restores **marginal** coverage

| param | raw coverage | conformal coverage | raw \|gap\| | conformal \|gap\| |
|-------|-------------:|-------------------:|------------:|------------------:|
| D     | 0.528 | 0.890 | 0.372 | **0.010** |
| f     | 0.550 | 0.876 | 0.350 | **0.024** |
| D\*   | 0.555 | 0.903 | 0.345 | **0.003** |

Marginal coverage is restored to within ≤0.024 of nominal for every parameter.

### Honest caveat: **conditional** coverage is *not* restored

Conformal applies a single marginal offset, so it cannot fix coverage that
varies across the parameter range. Post-conformal conditional coverage by
true-D\* tercile:

```
   Dstar  g0(low)=0.972  g1(mid)=0.929  g2(high)=0.810
```

The high-D\* tercile still under-covers (0.810 vs 0.900) while low-D\*
over-covers (0.972). This is the **irreducible identifiability limit** of IVIM
`D*` — pseudo-diffusion is weakly constrained by the signal — and it is a
property of the data, not a bug in the wrapper. Caliper's job is to *measure*
this faithfully, which it does.

---

## What's in the box

```
caliper/
  metrics.py        # numpy-only calibration ruler (the canonical core)
  forward.py        # bi-exponential IVIM model + synthetic cohorts
  estimator_maf.py  # conditional MAF posterior over (D, f, D*)  [torch]
  conformal.py      # split-conformal / CQR coverage correction
examples/demo.py    # one-command end-to-end pipeline (fixed seeds)
tests/              # pytest: metrics, forward, conformal (numpy) + estimator (torch)
```

Run the tests:

```bash
pip install -e ".[dev]"
pytest -q          # 30 tests (estimator tests auto-skip without torch)
```

## License

MIT — see [LICENSE](LICENSE).

## Roadmap

See [ROADMAP.md](ROADMAP.md). Value-of-information, decision-gap, and
deployment validity-monitor functionality are **deliberately deferred** and not
implemented here.
