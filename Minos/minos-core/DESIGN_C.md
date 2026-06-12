# DESIGN_C — Minos-Core v3: a deployment-validity monitor for a stale loss-calibration correction under shift

Status: **design only — no implementation in this file.** Every quantity below has an estimator
and a degenerate limit; the limits become the CHECKPOINT gates. This document is the CP0
deliverable and locks the design for the v3 contribution. Numbers quoted as *expected* were
confirmed in a throwaway numerical probe on the v2 modules (seed `20240517`, `N=1e6`); the
committed numbers come only from `experiments/run_c.py`.

---

## 0. What changed from v2, and why — read this first

A prior-art sweep found that the **prescriptive** claim "tune the reported error bar so the
*decision* is optimal" is an existing field: **loss-calibrated Bayesian inference**
(Lacoste-Julien, Huszár & Ghahramani 2011; Cobb 2018; Kušmierczyk et al. 2019; post-hoc form,
Vadera et al. 2021) and **decision calibration** (Zhao, Kim, Sahoo, Ma & Ermon 2021). Our v2
decision-calibrated scale `τ*` *is* a one-dimensional instance of that line. So v3 does **not**
introduce a calibration method — that would reinvent existing work.

**v3 contributes the thing that line does not address.** Loss-calibration and decision-calibration
all assume the **calibration set represents deployment**. The novel object here is a
**deployment-validity monitor**: a *label-free* statistic that detects when a learned
loss-calibration correction has gone **stale** under distribution shift — together with the v2 gap
as the motivating result. Two hard boundaries (also in `POSITIONING.md`):

1. **Not a new calibration method.** Loss-calibration / decision-calibration own that; we *cite and
   benchmark against* it (`τ̂*_cal` is the cited baseline).
2. **Not a generic OOD detector.** That field is saturated. The monitor is **regret-targeted**: it
   estimates the staleness of *this correction's decisions*, tied to the utility, **not** to input
   density.

---

## 1. The honesty constraint (non-negotiable; threads every checkpoint)

A label-free monitor can only see **observable** deployment statistics — the reported posterior
summaries `μ` (and the fixed reported spread `s`) and their distribution. A shift that changes the
**truth↔report** relationship *without moving any observable summary* is **undetectable without
labels**. A monitor that appears to catch it is leaking labels or overfitting.

Therefore the deployment shift generator exposes **two orthogonal knobs**, and the cleanest
construction induces the **same** report-understates-truth discrepancy two ways:

- **`delta_obs`** — moves the **observable** channel: biases the reported point **down**,
  `μ = report_center − β·s·δ_obs` (truth `θ` held). The reported-`μ` distribution translates, so the
  staleness is **detectable in principle**.
- **`delta_hid`** — moves the **hidden** channel: biases the **truth up**,
  `θ = report_center + β·s·δ_hid + s·u` (reported `μ = report_center` held). The truth↔report gap
  changes identically to `delta_obs`, but **no observable summary moves**, so the staleness is
  **undetectable in principle**.

At matched `δ`, the two produce the **same** deployment decision problem (same `θ−μ` law → same
`τ*_oracle` → same regret `R`), distinguished *only* by whether the observable `μ`-channel moved.
This is the sharpest possible honesty demonstration: **identical regret, one detectable and one
not.** The build measures detection **separately** for each knob and **reports the hidden-staleness
fraction as a result, not a failure.** Honest claim: *the monitor catches observable-driven
staleness; the hidden fraction requires labeled spot-checks, motivating periodic repeatability
validation.*

*Probe confirmation.* `delta_hid` moves `R` from `0` to `+0.34` while `M` stays pinned at its
`0.0043` noise floor (byte-identical reported `μ`); `delta_obs` moves both `R` (→`+0.18`) and `M`
(→`+0.58`).

---

## 2. Setting and reuse / new map (the CP0 audit)

**Setting.** `D_cal` = a labeled synthetic calibration set (ground-truth `θ` known; the
Fashion-seam stand-in). `D_dep` = an *unlabeled* deployment set, shifted by `(delta_obs, delta_hid)`.
Decision: treat / spare / escalate on scalar severity `θ` (v1 utility `U`, asymmetric `λ=k_under/k_over`).
Reported posterior `q = N(μ, (τ·s)²)`; `τ` is the scalar correction. The calibration cell is the v2
default misspecification — `κ=3, λ=3, ρ=0.5`, posterior-centric, `s=0.5` — so the decision is
`τ`-sensitive (skewed, asymmetric, near the single active threshold `t2`).

| Component | source | v3 status |
|---|---|---|
| Action set, `U(a,θ)`, `EPP`, `EU(a\|q)` | `minos/utility.py` | **reused unchanged** |
| Bayes action `a*(q)` | `minos/decision.py` | **reused unchanged** |
| Per-voxel `realised_utility`, policy EU | `minos/voi.py` | **reused** (primitives for the deploy EU) |
| `τ*` optimiser (grid + parabola-fit refine) | `minos/calibration.py` `tau_star` | **reused** — *is* the loss-calibration baseline `τ̂*_cal` (zero-shift, on `D_cal`) |
| v2 gap `G = τ*−τ_stat`, `τ_stat` | `minos/calibration.py` | **reused** for the motivation figure (a) |
| Trust-gate `g(w)`, `VoTG`, `δ_be` | `minos/gate.py`, `calibration.py` | **reused / retained** (v2 result kept; `_auc` reused) |
| Skew error `u(κ)`, posterior-centric truth, `report_center` | `minos/generative.py` | **reused** (`_unit_skew_error`, `make_population`) |
| Seeding, CRN | `minos/seeding.py` | **reused unchanged** |
| Frozen config | `minos/config.py` | **unchanged** (v3 constants live in the new modules; no v1/v2 default disturbed) |
| `delta_obs`/`delta_hid` deployment realiser | — | **NEW** `generative.realise_deploy` (v1/v2 `realise` untouched) |
| Loss-calibration baseline `τ̂*_cal` (cited) | — | **NEW** `correction.fit_loss_calibration` (thin, over reused `tau_star`) |
| Deploy EU, oracle deploy scale `τ*_oracle`, stale regret `R` | — | **NEW** `correction.py` |
| **Validity monitor `M(D_dep)`** | — | **NEW** `monitor.py` (the contribution; one-function swappable interface) |
| Threshold `m*`, gated-recovery policy, gated regret, detection | — | **NEW** `monitor.py` |

**v1/v2 retained, not replaced.** All v1/v2 modules, tests, drivers, and figures stay; v3 adds new
modules and a new driver only. `BASELINE_V1` and the v2 gap modules are untouched.

---

## 3. The estimands (defined precisely)

### 3.1 Loss-calibration baseline `τ̂*_cal` — cited, not novel
```
τ̂*_cal = argmax_τ  Û(τ; D_cal) = argmax_τ  (1/N) Σ_i U( a*(N(μ_i,(τ s)²)), θ_i )      [labels θ_i used]
```
This is the scalar post-hoc loss-calibration instance (Lacoste-Julien 2011 / Vadera 2021, reduced to
the 1-D scale correction of our toy) — **identical** to v2's `tau_star` on `D_cal` at zero shift.
The module comment-cites the source. It needs the labeled `D_cal` (it reads `θ`).

### 3.2 Deploy EU, oracle deploy scale, stale-correction regret
```
deploy_EU(τ; δ) = (1/N) Σ_i U( a*(N(μ_i(δ),(τ s)²)), θ_i(δ) )         (μ_i(δ), θ_i(δ) from realise_deploy)
τ*_oracle(δ)   = argmax_τ deploy_EU(τ; δ)                              (needs deployment labels — sim only)
R(δ)           = deploy_EU(τ*_oracle(δ); δ) − deploy_EU(τ̂*_cal; δ)  ≥ 0
```
`R` is the utility the **stale** correction leaves on the table versus the deployment-optimal scale.
`R ≥ 0` by definition of argmax. The oracle is known only in simulation (it reads deployment `θ`),
used solely to *validate* the monitor — never inside the monitor or the policy.

### 3.3 Validity monitor `M(D_dep)` — the contribution (DEFAULT)
`M` is a statistic of the **unlabeled** deployment reported points `{μ_i}` plus the **known** utility.
Its signature is `monitor(mu_dep, cfg, ref)` — it takes the reported `μ` array **only**; it is
**structurally** label-free (no `θ`, no `base`). `ref` is a frozen calibration reference built once
at calibration time.

DEFAULT estimand — **utility-weighted divergence** of the decision-relevant reported coordinate:
```
z_i  = (μ_i − t2)/s                                  # signed reported distance to the active threshold (observable)
ω(z) = k_under · ϕ(z)                                # decision-stakes kernel: cost magnitude × near-threshold weight
M(D_dep) = Σ_bins ω(z_b) · | p_dep(z_b) − p_cal(z_b) | · Δz        # utility-weighted L1 (TV-like) divergence
```
`p_cal`, `p_dep` are histogram densities of `{z_i}` on a fixed grid (`p_cal` frozen in `ref`). The
weight `ω` is the **regret-targeting** ingredient: the scale `τ` can only change the treat/escalate
action for reported points within `O(1)` reported-sd of `t2`, and the utility at stake there scales
with the under-treatment cost `k_under`. Far from the threshold (`|z|≫1`), `ω→0`: a distributional
change there *cannot* change which scale is optimal and so **must not** raise the monitor. A plain
OOD/density score would use uniform or density-ratio weights and fire on those irrelevant far-field
changes — that is precisely the distinction from generic OOD (Boundary 2).

**Disjoint code path (asserted in tests).** `M` is a pure function of `mu_dep` and `cfg`; it never
reads `θ_true`. Test: scrambling `base.theta` leaves `M` bit-identical, and the signature exposes no
`θ`.

**Swappable interface (the one design fork).** `monitor(mu_dep, cfg, ref, kind=...)` dispatches over
`{"utility_divergence"` (DEFAULT)`, "action_divergence"}`, with a documented third (conformal).
Alternative **action_divergence**: the cost-scaled total-variation between the *induced action*
distribution `[spare, treat, escalate]` on deployment vs calibration — a coarser, decision-grounded
statistic (3 action bins vs the full reported-coordinate histogram), implemented as a real second
monitor to demonstrate the swap. (A self-consistency residual of the *reported* objective
`Ê(τ)=mean_i EU(a*(q_i)|q_i)` was prototyped and **rejected**: `Ê(τ)` is monotone-decreasing in `τ`
— a decision grading itself under its own posterior always prefers a narrower bar — so its argmax is
degenerate at the boundary and the residual carries no staleness signal. Documented as a dead end,
not shipped.) A **conformal holdout discrepancy** is noted as future. All kinds share the
`(mu_dep, cfg, ref)` signature so swapping touches nothing else.

### 3.4 Threshold `m*` and gated-recovery policy
```
m*  = (1−α) quantile of M under the zero-shift null            (bootstrap deployment-size batches ~ D_cal across seeds; α=0.05)
gated-recovery actions (label-free policy):
    if  M(D_dep) ≤ m*:   a_i = a*(N(μ_i,(τ̂*_cal s)²))                          # trust the correction
    if  M(D_dep) >  m*:   a_i = ESCALATE  where |z_i| < z_guard, else a*(…)     # conservative override near the threshold
```
The override direction is the **conservative arm** dictated by the *known* cost asymmetry
(`k_under > k_over` ⇒ under-treatment is costlier ⇒ escalate), **not** inferred from labels — `M` is a
magnitude detector. The policy reads only `μ_dep` and the scalar `M`; `θ` enters solely to *score*
the resulting utility (regret), never to choose an action (asserted).
```
R_gated(δ) = deploy_EU(τ*_oracle(δ); δ) − EU(gated-recovery actions; δ)
```

---

## 4. Sweep axes

| knob | symbol | meaning | grid | null value |
|---|---|---|---|---|
| observable shift | `delta_obs` | downward bias on reported `μ` | `{0, 0.03, …, 0.30}` | `0` |
| hidden shift | `delta_hid` | upward bias on truth `θ` (`μ` fixed) | `{0, 0.03, …, 0.30}` | `0` |

Detection ROC pools `seeds × δ-grid` batches per knob (independent deployment seeds give `M` its
noise-floor spread and `R` its variation); `corr(M, R)` is computed over the CRN δ-sweep. `tol` for
`{R > tol}` is set as a small documented multiple of the zero-shift regret floor (balanced ROC).

---

## 5. Sanity limits → GATES

**GATE 0 (this document).** Estimands defined (§3); `M`'s no-label code path specified
(`monitor(mu_dep,…)`, scramble-`θ` test); the `delta_obs`/`delta_hid` split specified (§1, the
same-regret-two-ways construction); the swappable monitor interface specified (§3.3); boundaries
written (§0). Sanity limits argued below, **including the at-chance hidden case**.

**GATE 1 — baseline + stale-regret.**
- *Consistency:* `τ̂*_cal → τ*_oracle` as `N_cal → ∞` at zero shift (print bias/variance across
  `N ∈ {1e5,…,4e6}`); the value reproduces v2 (`τ*≈1.04` at `4e6`; v2 `G`, `τ_stat`, `VoC` reproduce).
- `R(0,0) ≈ 0` (at the regret floor); `R` **increases** with `delta_obs` and **separately** with
  `delta_hid`. *Argument:* at zero shift the stale scale is deployment-optimal; each knob makes the
  report understate the truth, moving `τ*_oracle` up and stranding the stale scale.

**GATE 2 — validity monitor `M`.**
- **No-label path:** `monitor` output is invariant to scrambling `θ` (disjoint path).
- *Under `delta_obs`:* `corr(M, R) > thr` (expect Spearman ≈ 1) **and** detection AUC for
  `{R > tol} > 0.5 + margin`.
- *Under `delta_hid`:* detection AUC `≈ 0.5` (**printed — the honest limitation, not a red gate**),
  because `M` is a function of `{μ}` only and `delta_hid` does not move `{μ}`.
- *Zero shift:* false-alarm rate `≤ target` at `m*`.

**GATE 3 — gated recovery + monitor calibration.**
- `m*` controls the zero-shift false-alarm rate.
- *Under `delta_obs`:* gated regret `<` stale regret (recovery); print the ladder
  `{stale, gated, oracle}` across the `delta_obs` sweep with the operating point marked.
- *Zero shift:* gated `≈` stale (no false-alarm harm — the gate doesn't fire).
- *Under `delta_hid`:* gated `≈` stale (**the gate does not spuriously help** — it cannot see the
  hidden shift; gated `<` stale here would signal leakage).

**GATE 4 — figures + reproduction.** Four vector PDFs exist; every `RESULTS_C` number traces to
`run_c.py` stdout; clean-seed rerun identical; tests pass (including the no-label assertion on `M`).

---

## 6. Module map & figures

```
minos-core/
  minos/
    generative.py    + realise_deploy(base, cfg, *, delta_obs, delta_hid) -> (mu, theta)
                       (obs biases mu down; hid biases theta up; CRN; posterior-centric; v1/v2 realise UNCHANGED)
    correction.py    NEW: fit_loss_calibration (cited baseline = tau_star on D_cal),
                       deploy_expected_utility, oracle_deploy_scale, stale_regret R
    monitor.py       NEW: MonitorRef (frozen p_cal + weights + tau_hat_cal); monitor(mu_dep,cfg,ref,kind=)
                       [DEFAULT utility_divergence + action_divergence alt]; calibrate_threshold (m*);
                       gated_recovery_actions; gated_regret; detection helpers (reuse gate._auc)
    utility/decision/voi/gate/diagnostics/calibration/config/seeding.py    UNCHANGED
  tests/
    test_deploy_split.py   NEW: mu-invariance under delta_hid; obs moves mu; hid moves theta; CRN; posterior-centric req
    test_correction.py     NEW: GATE 1 (consistency tau_hat_cal->oracle; reproduces v2; R(0)~0; R up each knob)
    test_monitor.py        NEW: GATE 2/3 (no-label assertion; corr/AUC under obs; AUC~0.5 under hid;
                           false-alarm<=target; gated<stale obs; gated~stale zero & hid)
    test_*.py (v1/v2)      UNCHANGED, still pass
  experiments/
    run_all.py / run_b.py  UNCHANGED
    run_c.py               NEW: CP1 baseline+R; CP2 monitor; CP3 gated+m*; CP4 figures
  DESIGN_C.md  RESULTS_C.md  POSITIONING.md(updated)
  figures/  fig_c_*.pdf  (v1 fig_a..d, v2 fig_gap_* retained)
```

Planned v3 figures (vector PDF, light):
- **(a) motivation:** v2 gap `G` vs `κ` (the misspecification that motivates a loss-calibration correction).
- **(b)** stale-correction regret `R` vs shift, both knobs overlaid.
- **(c) the honesty figure:** monitor `M` vs `R` scatter + detection ROC, **observable vs hidden overlaid**.
- **(d)** regret ladder `{stale, gated, oracle}` across `delta_obs`, monitor operating point marked.

---

## 7. Discipline

Run-then-write; no fabricated numbers; deterministic seeded RNG (`GLOBAL_SEED`, CRN across each
sweep); each checkpoint ends with a GATE (assertions + printed numbers); stop at any red gate; a
single PR. The hidden-shift at-chance AUC is a **documented result**, not a red gate. v1/v2 retained
verbatim.
