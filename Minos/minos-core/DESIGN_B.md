# DESIGN_B — Minos-Core v2: the decision–calibration gap

Status: **design only — no implementation in this file.** Every quantity below has an
estimator and a degenerate limit; the limits become the CHECKPOINT gates. This document
locks the design for Headline B and is the CP0 deliverable.

---

## 0. What changed from v1, and why

v1 (PR #2, merged) proved a calibrated per-voxel error bar carries positive decision value:
`VoC(τ) = EU(posterior | τ=1) − EU(posterior | τ) ≥ 0`, minimised at `τ=1`. But the v1 prior
is **symmetric about the decision midpoint** `(t1+t2)/2`, and the reported error
`θ − μ ~ N(0, s²)` is **exactly Gaussian**. Under those two facts the value-optimal scale and
the statistically-calibrated scale *coincide at `τ=1` by construction*. A reviewer rightly asks
what is surprising: nothing yet.

**v2 headline.** Statistical calibration and *decision* calibration **diverge** under realistic
(misspecified / skewed) posteriors, by a **quantifiable gap**. Define two scales by two
independent criteria:

- **`τ_stat`** — *statistical* calibration: the scale at which the reported posterior attains
  nominal central-interval coverage of the truth. Computed by a **coverage root-find**. Uses
  only the error distribution and the reported intervals; **the utility is never read**.
- **`τ*`** — *decision* calibration: `τ* = argmax_τ E[U_posterior(τ)]`. Computed by a **dense
  grid + local refine** of the expected-utility curve. Uses only the utility/decision core;
  **coverage is never read**.
- **Decision–calibration gap** `G = τ* − τ_stat` (we also report the ratio `τ*/τ_stat`).

`G` is the object **none of the neighbours compute**: Vickers DCA/net-benefit scores the
*marker's point prediction*; ISPOR VoI scores *resolving a parameter*; ARCliDS is a *policy*.
None measures the *distance between statistical and decision calibration of an error bar*.
The well-specified case (`G ≈ 0`, v1) is the **degenerate baseline inside the result**, not the
result.

---

## 1. Reuse / new map (the CP0 audit)

| Component | v1 module | v2 status |
|---|---|---|
| Action set, `U(a, θ)`, `EPP`, `EU(a\|q)` | `minos/utility.py` | **reused unchanged** |
| Bayes action `a*(q)` | `minos/decision.py` | **reused unchanged** |
| Policy EU, EVPI-analog, `VoC(τ)`, `posterior_eu_curve` | `minos/voi.py` | **reused unchanged** (operate on the reported Gaussian + true `θ`; agnostic to how `θ−μ` is distributed) |
| Trust-gate `g`, threshold, gated policy, `VoTG(δ)`, detection AUC | `minos/gate.py` | **reused unchanged** (CP3) |
| Central-interval coverage, ECE | `minos/diagnostics.py` | **reused unchanged** (now measures the *misspecified* coverage automatically, because `realise` supplies skewed `μ`) |
| Seeding | `minos/seeding.py` | **reused unchanged** |
| `τ` calibration knob | `voi`/`decision` (`σ_rep = τ·s`) | **reused unchanged** |
| Frozen config | `minos/config.py` | **extended** (new fields, all defaulting to v1 behaviour; one assertion relaxed — see §3) |
| Generative + measurement | `minos/generative.py` | **extended**: error draw becomes a standardised skew-normal (κ knob) that **reduces to v1's `N(0,s²)` bit-for-bit at κ=0**; latent gains an optional single-Gaussian mode for the `ρ` axis |
| `τ_stat` solver (coverage root-find) | — | **new** `minos/calibration.py` |
| `τ*` optimiser (grid + refine) | — | **new** `minos/calibration.py` |
| gap `G`, ratio | — | **new** `minos/calibration.py` |
| break-even shift `δ_be` (root-find on `VoTG`) | — | **new** `minos/calibration.py` (promotes the v1 footnote) |
| `(κ, λ, ρ)` sweep, gap surface | — | **new** `experiments/run_b.py` |
| v1 symmetric setup | `MinosConfig()` default | **retained** as the named baseline `BASELINE_V1` (mixture latent, κ=0) |

**v1 is retained, not replaced.** `MinosConfig()` with its defaults still describes the v1
symmetric model; `experiments/run_all.py` and all v1 tests keep passing unchanged. v2 adds new
fields whose defaults are exactly v1, plus new modules.

---

## 2. The misspecification family (DEFAULT — flag to override)

The decision core only ever sees the **reported** posterior `q = N(μ, (τ·s)²)` and the **true**
severity `θ`. v2's headline is a *misspecified per-voxel posterior*: the **true posterior is
skewed**, the **report is its Gaussian moment-match**. The cleanest realisation makes the report
centre primary (DESIGN_B's *posterior-centric* model), so the report carries the **correct mean and
sd** and only the **shape** is wrong — isolating the gap as a pure skew effect.

**Posterior-centric generative.** Sample report centres `μ` (the "reported point", spanning the
decision region); draw the truth `θ = μ + s·u`, where `u` is a standardised skew error (mean 0, sd
1, shape `κ`). Then:

- The **true posterior** `θ | μ` is skewed with **mean `μ`, sd `s`** — its Gaussian moment-match is
  exactly `N(μ, s²)`, so **`τ=1` is the moment-matched report by construction**; `τ` scales the
  reported sd (`σ_rep = τ·s`).
- `θ − μ = s·u` is **right-skewed** for `κ>0`: a heavy **under-treatment** tail (truth occasionally
  far more severe than the report). This is the clinically dangerous failure of a bounded severity
  whose Gaussian report cannot represent the upper tail.
- **`κ=0` ⟹ `u = z_eta`, `θ − μ ~ N(0, s²)` exactly ⟹ well-specified.** Crucially, because
  `E[θ|μ]=μ` (no prior shrinkage), at `κ=0` the report **equals** the true posterior, so `τ*=1`
  *rigorously for any report-centre distribution and any utility asymmetry* — the gap at `κ=0` is
  exactly 0, with no prior-curvature confound. (Contrast the v1 *forward* model `μ=θ+η`: with an
  informative, off-centre prior the report `N(μ,s²)` is **not** the Bayes posterior at `τ=1` —
  prior shrinkage drifts `τ*` off 1 even at `κ=0`. That drift, the v1 §6.2 footnote, is a
  *different* effect from misspecification, and the posterior-centric model removes it.)

**Why skew alone, decision-relevant.** A symmetric central credible interval sees only `|θ−μ|`
(the *even* part of the error law), so coverage is **blind to which tail is heavy**. The asymmetric
decision sees the **signed** costly tail. That divergence between what *coverage* measures and what
the *decision* pays for is precisely why `τ_stat` and `τ*` separate.

**CRN-preserving sampler (exact v1 limit).** Two standard-normal base draws `z_eta`, `z_skew`
(drawn once, reused across the whole sweep) and `δ(κ)=κ/√(1+κ²)`:

```
raw = δ·|z_skew| + √(1−δ²)·z_eta           # skew-normal(shape +κ): heavy RIGHT tail for κ>0
u   = (raw − δ·√(2/π)) / √(1 − 2δ²/π)       # standardise: mean 0, sd 1
θ   = μ + s · u                            # posterior-centric truth (θ−μ right-skewed)
```

At `κ=0`: `δ=0 ⟹ u = z_eta`, so `θ = μ + s·z_eta` and, in the v1 forward model, `μ = θ + s·z_eta`
— **identical to v1's `η = s·z_eta`**. `z_skew` is drawn *after* v1's `(θ, z_eta, z_w)`, so v1's
variates are byte-for-byte unchanged and every v1 number reproduces. Flags: `family ∈ {"skewnorm"
(default), "gaussian"}`, `posterior_centric ∈ {True (gap), False (v1 forward)}`. Per the override
note, swapping the family (heavy-tailed / bias-driven) touches **only this generator** — the gap
machinery (`τ_stat`, `τ*`, `G`, sweep, `δ_be`) is mechanism-agnostic.

**Shift reused for the trust-gate.** The v1 shift transform (`σ_true=s(1+αδ)`, bias `−βsδ`, feature
`w~N(δ,1)`) is reused verbatim by CP3; the gap experiments (CP1–CP2) run at `δ=0`.

---

## 3. The two calibration scales — defined independently

The independence is structural: the two estimators import **disjoint** parts of the package.

### 3.1 `τ_stat` — statistical calibration (coverage root-find; never reads utility)

For the reported `q = N(μ, (τs)²)`, the central `L`-interval is `[μ − z_L τ s, μ + z_L τ s]`,
`z_L = Φ⁻¹(½ + L/2)`. Empirical coverage `C(τ; L) = mean( |θ − μ| ≤ z_L τ s )` is exactly
`minos.diagnostics.central_interval_coverage`. `C(·; L)` is monotone increasing in `τ`, with
`C→0`/`1` at the ends, so

```
τ_stat(L) = brentq root of  C(τ; L) − L = 0      (bracket τ ∈ [0.2, 5])
```

is unique. **Reference level** `L_ref = 0.90` (a documented modelling choice; the full
`L ∈ {0.5,…,0.99}` profile and the ECE-minimising scale are reported in RESULTS_B). Code path:
`calibration.tau_stat` → `diagnostics.central_interval_coverage` only. **No `utility`/`decision`/
`voi` import.**

- *κ=0 limit:* `C(τ;L)=P(|N(0,1)|≤z_L τ)`, `C(1;L)=L` at **every** `L` ⟹ `τ_stat=1` (level-free).
- *κ>0:* the `|θ−μ|` law is non-Gaussian, so `C(1;L)≠L`. At central levels (`L≤0.95`) the
  right-skew makes the symmetric interval **over-cover** ⟹ `τ_stat < 1` ("calibration says
  *shrink*"). Only at the extreme tail (`L=0.99`) does the heavy upper tail force under-coverage ⟹
  `τ_stat > 1`. `τ_stat` is **independent of the utility** (it never reads `k_under`/`k_over`).

### 3.2 `τ*` — decision calibration (grid + parabola-fit refine; never reads coverage)

```
EU(τ)  = expected_utility("posterior", base, cfg, tau=τ)        # voi.py; CRN ⟹ smooth in τ
τ*     = argmax_τ EU(τ),  τ ∈ [0.5, 2.5]:
         dense grid (step 0.05) → fit a parabola to EU(τ) in a window around the grid peak →
         vertex (clamped to the window)
```

CRN makes `EU(τ)` a smooth deterministic function for a given sample, but near a **flat** optimum
(the well-specified `κ=0` case) a raw argmax chases finite-sample bumps; the **parabola-fit vertex**
estimates the optimum from the curve's shape and is robust to those bumps. **Flat-EU guard:** if the
EU span over the grid is below `1e-4` the decision is τ-insensitive (e.g. symmetric utility `λ=1`,
where the escalate boundary is τ-independent) — there is no decision-relevant optimum, so by
convention `τ*=1`. Code path: `calibration.tau_star` → `voi.expected_utility` → `decision.bayes_
action` + `utility`. **No `diagnostics`/coverage import.** Equivalently `τ* = argmin_τ VoC(τ)`.

### 3.3 The gap (two components)

```
G      = τ* − τ_stat  =  (τ* − 1)  −  (τ_stat − 1)
ratio  = τ* / τ_stat
```

`G` decomposes into a **decision-side inflation** `(τ*−1) ≥ 0` (needs `κ>0` *and* `λ>1`) minus a
**coverage-side offset** `(τ_stat−1) ≤ 0` at central levels. Under skew the two pull **opposite
ways** — statistical calibration says *shrink* (`τ_stat<1`), decision calibration says *hold or
widen* (`τ*≥1`) — so `G>0`: a *deliberately underconfident* error bar beats the
statistically-calibrated one. `G>0` = "underconfident-favored", the documented sign.

---

## 4. Sweep axes

| axis | symbol | meaning | how varied | well-specified value |
|---|---|---|---|---|
| posterior skew | `κ` | skew-normal shape of `θ−μ` | `κ ∈ {0, 0.5, 1, 2, 3, 4}` | `κ=0` |
| utility asymmetry | `λ` | under:over cost ratio `k_under/k_over` | `k_under ∈ {1,2,3,4,6}`, `k_over=1` | `λ=1` (symmetric) |
| threshold proximity | `ρ` | distance from `E[μ]` to the active threshold `t2`, in units of `s` | report centres `μ ~ N(t2 − ρs, θ_std²)`, `θ_std=0.5`; swept over the resolved regime `ρ ∈ {0,0.5,1,1.5,2}` | large `ρ` (deep in a region) |

For the gap configs the report centres are a **single Gaussian** (`latent_mode="gaussian"`,
`posterior_centric=True`); the spare threshold is pushed far below `t2` so the decision is a clean
**treat vs escalate** around the single active threshold `t2`, leaving room for the `ρ` sweep. The
v1 **mixture** latent is retained for the v1-reproduction baseline.

---

## 5. Sanity limits → GATES

**GATE 0 (this document).** Definitions independent (disjoint imports, §3); reuse/new map written
(§1); limits argued (below).

**GATE 1 — gap estimator.**
- *Well-specified `κ=0`* (posterior-centric gaussian config): `|G| < 0.02`, `τ*≈1`, `τ_stat=1`
  (level-free). Separately, the v1 **mixture** config at `κ=0` reproduces v1: `VoC` argmin at
  `τ=1`. *Argument:* `κ=0 ⟹ θ−μ~N(0,s²)` ⟹ report = true posterior at `τ=1` ⟹ `τ_stat=1` (coverage
  nominal at every level) and `τ*=1` (Bayes-optimal scale), §2.
- *Default misspecified* (`κ=3`, `λ=3`, `ρ=0.5`): `G > 0` with the documented sign `τ* > τ_stat`
  (underconfident-favored). *Argument:* the right-skewed `θ−μ` puts costly under-treatment mass in
  the upper tail the symmetric interval cannot see; coverage over-covers at `L_ref` (`τ_stat<1`)
  while the asymmetric decision hedges up (`τ*>1`).
- *Independence:* `τ_stat` is invariant to the utility (`k_under`/`k_over`); `τ*` takes no coverage
  level — disjoint code paths, asserted in the tests and named in the gate printout.

**GATE 2 — gap map.**
- *Well-specified slice `κ=0`*: `G ≈ 0` across all `λ, ρ` — at `κ=0` both scales are 1 regardless of
  `λ` or `ρ` (no skew ⟹ nothing for the decision to hedge against and coverage is nominal).
- *Monotonicity in `κ`:* for fixed `λ>1, ρ`, `G` increases with `κ` (more skew ⟹ heavier unseen
  upper tail ⟹ `τ_stat` falls and `τ*` rises).
- *`λ` dependence:* for fixed `κ>0`, `G` increases with `λ`; at `λ=1` the decision-side inflation
  vanishes (flat EU, `τ*=1`) leaving the pure **coverage offset** `G = 1 − τ_stat > 0`.
- *`ρ` dependence:* `G` is **largest at `ρ=0`** (report centres straddling the threshold, where the
  skewed tail crosses the boundary) with a **clear downward trend** as `ρ` grows (the threshold
  interaction weakens, so `τ*` relaxes toward 1). The gate checks `argmax_ρ G = ρ=0`, a negative
  least-squares slope, and a net decrease `> 0.01` — robust to the MC wiggle that a strict
  adjacent-pair test would trip on. `ρ` is swept over the **resolved regime** `{0,0.5,1,1.5,2}`:
  beyond `ρ≈2` the population leaves the threshold, the EU curve is degenerately flat, and `τ*` is
  unidentified — that far field is excluded by design, not papered over (§3.2 note).
- Corner table of `G` at the grid extremes printed.

**GATE 3 — break-even shift `δ_be`.** Promote the v1 footnote (`VoTG` dips slightly negative for
small shifts, crosses zero, then grows) to a number: `δ_be = brentq root of VoTG(δ)=0`. Assert
`VoTG(δ < δ_be) < 0`, `VoTG(δ > δ_be) > 0`, `δ_be` finite; print `VoTG` just below/above. Reuses
v1's `gate.votg` verbatim on the baseline config.

**GATE 4 — figures + reproduction.** Four vector PDFs exist; every RESULTS_B number traces to
`run_b.py` stdout; the sweep reproduces from the clean seed.

---

## 6. Module map & figures

```
minos-core/
  minos/
    config.py          + kappa, family, latent_mode, posterior_centric, theta_mean, theta_std
                       (defaults = v1); gaussian_latent_config / DEFAULT_MISSPEC / BASELINE_V1;
                       assertion relaxed: k_under >= k_over (was strictly >)  to allow λ=1
    generative.py      + z_skew base draw; standardised skew-normal error (κ, right tail);
                       gaussian latent; posterior-centric direction (κ=0, mixture == v1 exactly)
    calibration.py     NEW: tau_stat (coverage root-find), tau_star (grid + parabola-fit refine
                       + flat-EU guard), gap, break_even_shift   (disjoint import sets — see §3)
    utility/decision/voi/gate/diagnostics/seeding.py    UNCHANGED
  tests/
    test_misspec.py        NEW: κ=0 == v1; skew sign; standardisation; gaussian latent
    test_calibration.py    NEW: GATE 1 (κ=0 ⟹ G→0, reproduces v1; misspec ⟹ G>0, right sign),
                           independence of paths, GATE 3 (δ_be brackets a sign change)
    test_*.py (v1)         UNCHANGED, still pass
  experiments/
    run_all.py    v1 driver (unchanged)
    run_b.py      NEW: CP1 gap on default+baseline; CP2 (κ,λ,ρ) sweep; CP3 δ_be; CP4 figures
  DESIGN_B.md  RESULTS_B.md  POSITIONING.md(updated)
  figures/  fig_b_*.pdf  (v1 fig_a..d retained)
```

Planned v2 figures (vector PDF, light):
- **(a)** `VoC(τ)` well-specified vs misspecified, with `τ_stat` and `τ*` marked — *the gap made
  visible* (well-specified: min at `τ=1`, `VoC≥0`; misspecified: dips below 0, min at `τ*>τ_stat`).
- **(b) headline:** gap `G` vs `κ` at several `λ`.
- **(c)** `G` surface over `(κ, ρ)`.
- **(d)** `VoTG(δ)` with `δ_be` marked.

---

## 7. Discipline

Run-then-write; no fabricated numbers; deterministic seeded RNG (`GLOBAL_SEED`, CRN across the
sweep); each checkpoint ends with a GATE (assertions + printed numbers); stop at any red gate; a
single PR. The v1 symmetric baseline is retained as `BASELINE_V1` and as the unchanged v1 driver.
