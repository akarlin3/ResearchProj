# Minos-Core

The formal, **data-independent** core of *Minos*: a framework that prices the
**decision-value of a calibrated per-voxel error bar** — the value of the *uncertainty
itself* — on a treat / spare / escalate action, with a folded-in **trust-gate** that detects
when the reported uncertainty is untrustworthy under deployment shift.

The two first-class quantities are:

- **Value of Calibration** `VoC(tau)` — how much decision utility is lost when the reported
  error bar is mis-*scaled* (over- or under-confident) at a fixed point estimate.
- **Value of the Trust-Gate** `VoTG(delta)` — how much decision utility is *recovered* by
  detecting that the error bar has become *untrustworthy* under shift and acting conservatively.

These deliberately are **not** decision-curve net benefit and **not** population EVPI/EVPPI.
The object being priced is the error bar, not the point estimate. See `POSITIONING.md`.

This build is **100% synthetic** — a toy decision model with a scalar latent severity
`theta`. A real IVIM parameter map + Fashion posterior can later replace the synthetic source
through one marked seam (`minos/generative.py`) without touching the decision / VoI / gate core.

## v2 — the decision–calibration gap (headline)

`v1` proved a calibrated error bar carries positive decision value, but under a **symmetric prior
the value-optimum coincides with `tau=1` by construction**. `v2` asks what happens under a
**misspecified (skewed) posterior** and finds that the two ways to calibrate an error bar
**diverge**:

- **`tau_stat`** — the scale with nominal interval **coverage** (statistical calibration), by a
  coverage root-find.
- **`tau*`** — the scale that maximises **expected decision utility** (decision calibration),
  `argmax_tau E[U]`.
- **decision–calibration gap** `G = tau* − tau_stat` — the quantity none of the neighbours
  computes.

Under a right-skewed posterior with asymmetric cost, statistical calibration says *shrink* the bar
(`tau_stat<1`) while the decision says *widen* it (`tau*>1`): a deliberately under-confident bar
beats the statistically-calibrated one. `G` vanishes in the well-specified / symmetric corner
(recovering v1 at `tau=1`) and grows with posterior skew `kappa`, cost asymmetry `lambda`, and
threshold proximity. The break-even trust-gate shift `delta_be` (when turning the gate on pays off)
is promoted to a first-class number. New code: `minos/calibration.py`
(`tau_stat`, `tau_star`, `gap`, `break_even_shift`); driver `experiments/run_b.py`; see
`DESIGN_B.md`, `RESULTS_B.md`. The v1 symmetric model is retained verbatim (`BASELINE_V1`,
`experiments/run_all.py`) as the degenerate `G≈0` baseline.

```bash
python experiments/run_b.py      # prints the 4 v2 gate blocks, writes figures/fig_gap_*.pdf
```

## v3 — a deployment-validity monitor for a stale loss-calibration correction (headline)

A prior-art sweep found that the *prescriptive* move v2 implies — "tune the error bar so the
**decision** is optimal" — is an existing field: **loss-calibrated Bayesian inference** and
**decision calibration**. So v3 does **not** introduce a calibration method. It treats the
decision-calibrated scale as a **cited baseline** `tau_hat_cal` and contributes the thing that line
does not address: those methods all assume the **calibration set represents deployment**.

The novel object is a **deployment-validity monitor** `M(D_dep)` — a **label-free** statistic that
detects when the learned loss-calibration correction has gone **stale** under distribution shift —
plus a gated-recovery policy that acts on it.

- **`tau_hat_cal`** — the cited loss-calibration baseline, `argmax_tau E[U]` on a labeled
  calibration set (= v2's decision-calibrated scale).
- **stale-correction regret `R`** — utility the stale correction strands vs the deployment-optimal
  scale (oracle, simulation-only).
- **validity monitor `M`** — a utility-weighted divergence between the deployment and calibration
  reported-posterior shape, **regret-targeted** (weighted by the cost stakes near the decision
  threshold), computed from **unlabeled** deployment data. Behind a one-function swappable interface.
- **gated recovery** — where `M > m*`, override the decision-fragile voxels to the conservative arm.

**The honesty constraint is the heart of the build.** A label-free monitor sees only observable
reported summaries. The shift generator exposes two knobs that induce the **same** regret two ways:
`delta_obs` biases the reported point (observable → detectable) and `delta_hid` biases the truth with
the observable summaries held fixed (hidden → **undetectable by construction**). v3 measures detection
separately and **reports the hidden-staleness fraction as a result**: the monitor catches
observable-driven staleness (AUC ≫ 0.5) and is **at chance** on the hidden fraction (AUC ≈ 0.5),
which is what motivates periodic **labeled repeatability spot-checks**.

Two boundaries (see `POSITIONING.md`): (1) **not a calibration method** — that line is cited and
benchmarked; (2) **not a generic OOD detector** — `M` is regret-targeted (utility stakes), not
density-targeted. New code: `minos/correction.py` (`fit_loss_calibration`, `oracle_deploy_scale`,
`stale_regret`), `minos/monitor.py` (`monitor`, `calibrate_threshold`, `gated_recovery_actions`),
`generative.realise_deploy`; driver `experiments/run_c.py`; see `DESIGN_C.md`, `RESULTS_C.md`. v1/v2
are retained verbatim.

```bash
python experiments/run_c.py      # prints the 4 v3 gate blocks, writes figures/fig_c_*.pdf
```

## The math (condensed; full derivation in `DESIGN.md`)

**Actions & utility.** `A = {spare, treat, escalate}`, thresholds `t1 < t2`, under-treatment
slope `k_under` > over-treatment slope `k_over` (asymmetric: under-treating costs more):

```
U(spare,    theta) = - k_under · relu(theta - t1)
U(treat,    theta) = - k_over  · relu(t1 - theta) - k_under · relu(theta - t2)
U(escalate, theta) = - k_over  · relu(t2 - theta)
```

Each action is the unique maximiser on its own severity region, and the best action always
attains `U = 0`, so `max_a U(a, theta) ≡ 0` — the oracle utility is identically zero.

**Generative + measurement.** `theta ~ p(theta)` (mixture, symmetric about `(t1+t2)/2`);
estimate `mu = theta + eta`, `eta ~ N(b, sigma_true^2)`. Reported posterior
`q = N(mu, (tau·s)^2)` where **`tau` is the calibration knob** (`1` = calibrated) and `s` the
intrinsic spread. A shift `delta` inflates `sigma_true`, biases `b` downward, and moves an
observable acquisition feature `w` (the gate's input).

**Bayes step (closed form).** With `EPP(m, sigma) = m·Φ(m/sigma) + sigma·φ(m/sigma)`
(`EPP(m,0)=relu(m)`):

```
EU(spare|q)    = - k_under · EPP(mu - t1, sigma)
EU(treat|q)    = - k_over  · EPP(t1 - mu, sigma) - k_under · EPP(mu - t2, sigma)
EU(escalate|q) = - k_over  · EPP(t2 - mu, sigma)
a*(q)          = argmax_a EU(a | q)
```

**Policies** (scored by the true `U(a, theta_true)`): `point` = `a*(N(mu,0))` (ignore the bar);
`posterior` = `a*(N(mu,(tau s)^2))`; `gated` = `posterior` overridden to `escalate` where the
gate fires; `oracle` = `a*(N(theta_true,0))`.

**Quantities.**
```
EVPI-analog                 = EU(oracle) - EU(posterior)                       (= posterior regret)
value of using the error bar= EU(posterior | tau=1) - EU(point)
VoC(tau)                    = EU(posterior | tau=1) - EU(posterior | tau)      (headline)
VoTG(delta)                 = EU(gated | delta) - EU(posterior | delta)        (headline, under shift)
```

**Trust-gate.** Signal `g(w) = (w - m_w)/s_w` (one-sided OOD / density-ratio proxy),
threshold `g*` at training quantile `q_gate`; where `g > g*`, override to `escalate`.
Detection is scored by `AUC(g, shift-mask)`.

All expectations are seeded Monte Carlo with **common random numbers** across the `(tau, delta)`
sweep, so VoC/VoTG differences are low-variance.

## Install & run

```bash
pip install -e .            # numpy scipy matplotlib (pytest for tests)
pytest                      # 33 tests; the checkpoint gates are assertions
python experiments/run_all.py   # prints the 4 gate blocks, writes figures/*.pdf
```

`run_all.py` reproduces all four gates from the clean seed and writes the four vector-PDF
figures. The exact numbers it prints are transcribed in `RESULTS.md`.

## Sanity gates (must hold)
1. `EU_oracle ≥ EU_posterior ≥ EU_point`; EVPI-analog → 0 as the posterior → point mass.
2. `VoC(tau=1) = 0`, minimal; `VoC(tau) > 0` for `tau ≠ 1`; value of the error bar `> 0`.
3. `VoTG(delta=0) ≈ 0`; `VoTG(delta_test) > 0` with gated regret < posterior regret; detection
   `AUC > 0.5 + margin`.

## Layout
```
minos/
  seeding.py      one global seed -> explicit Generators (no bare np.random)
  config.py       frozen MinosConfig (all parameters)
  utility.py      U(a,theta), Action, EPP, EU(a|q)
  decision.py     bayes_action(q)
  generative.py   prior mixture + measurement + CRN  (# IVIM seam, deferred)
  voi.py          policy EU, EVPI-analog, value-of-error-bar, VoC
  gate.py         gate signal/threshold, gated policy, VoTG, detection AUC
  diagnostics.py  central-interval coverage, ECE
tests/            one file per module; gates encoded as assertions
experiments/run_all.py   seeded driver -> figures + printed numbers
DESIGN.md  RESULTS.md  POSITIONING.md
```

## IVIM seam (deferred)
`minos/generative.py` isolates the synthetic latent source and measurement behind a marked
region (`# IVIM seam — Fashion integration point (deferred)`). Replacing it with a real IVIM
parameter map + Fashion posterior leaves `decision.py`, `voi.py`, `gate.py` untouched.
