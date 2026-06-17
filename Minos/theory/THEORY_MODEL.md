# THEORY_MODEL.md — Checkpoint 0: the analytic model, locked against the code

**Purpose.** Plumbline proves two analytic statements *about the Minos-Core simulation*. For
"theory predicts sim" to mean anything, the analytic model must be the *actual* model in the
repository — not a convenient idealisation. This file reads the code, writes the model down in the
notation the theorems use, records the v1–v3 numbers the theorems must reproduce (with provenance),
states the assumptions, and runs **GATE 0**: a halt-able feasibility check that the code's model
admits a leading-order expansion.

All numeric claims in §4 and §5 were printed in-session by the CP0 audit script (the block run
before this file was written) or are quoted verbatim from the repo's own RESULTS files with the
commit recorded. No number here was hand-invented.

Provenance of the working tree: branch `theory-plumbline`, based at commit `72406ed` (Minos-Core
v3). v2 numbers are from `minos-core/RESULTS_B.md` (commit `a94e25f`, PR #3); v3 numbers from
`minos-core/RESULTS_C.md` (commit `72406ed`); v1 numbers from `minos-core/RESULTS.md`.

---

## 1. The code, read (not assumed)

### 1.1 Utility `U(a, θ)` — `minos-core/minos/utility.py:43`

Three actions `a ∈ {SPARE, TREAT, ESCALATE}`, two thresholds `t1 < t2`, two cost slopes
`k_over, k_under > 0`. `relu(x) = max(0, x)`. The utility is piecewise-linear and `≤ 0`:

```
U(SPARE,    θ) = − k_under · relu(θ − t1)
U(TREAT,    θ) = − k_over · relu(t1 − θ) − k_under · relu(θ − t2)
U(ESCALATE, θ) = − k_over · relu(t2 − θ)
```

It penalises *under-treatment* (truth above the action's safe zone) with slope `k_under` and
*over-treatment* (truth below) with slope `k_over`. The asymmetry is `λ := k_under / k_over`; the
config enforces `k_under ≥ k_over > 0` (`config.py:70`), so `λ ≥ 1`. The **under:over cost ratio**
`lam` in the gap sweep is realised as `k_under = lam, k_over = 1` (`config.py:114`), so `λ = lam`.

**Closed-form expected utility** under a *reported* Gaussian posterior `q = N(μ, σ²)` with
`σ = τ·s` (`utility.py:57`) uses the expected-positive-part
`EPP(m, σ) = E[max(0, Y)], Y ~ N(m, σ²) = m·Φ(m/σ) + σ·φ(m/σ)` (`utility.py:26`):

```
EU(SPARE   | q) = − k_under · EPP(μ − t1, σ)
EU(TREAT   | q) = − k_over · EPP(t1 − μ, σ) − k_under · EPP(μ − t2, σ)
EU(ESCALATE| q) = − k_over · EPP(t2 − μ, σ)
```

The decision rule is `a*(q) = argmax_a EU(a | q)` (`decision.py:13`).

### 1.2 Generative / skew model — `minos-core/minos/generative.py`

Gap experiments run **posterior-centric** (`config.py:50`, `gaussian_latent_config` sets
`posterior_centric=True`, `generative.py:104`): the latent draw is the *reported centre* `μ`, and
the truth is

```
θ = μ + s · u ,
```

where `u` is a **zero-mean, unit-variance, skew-normal-shaped** error (`_unit_skew_error`,
`generative.py:61`). With `d(κ) = κ / √(1+κ²)`, `raw = d·|z_skew| + √(1−d²)·z_eta` is a
skew-normal(shape `+κ`) variate, standardised to mean 0, sd 1. At `κ = 0` it returns `z_eta`
unchanged, so `u ~ N(0,1)` exactly (the well-specified v1 limit). For `κ > 0`, `u` has a heavy
**right** tail — the heavy *under-treatment* tail (truth occasionally far more severe than the
symmetric report can represent).

Report centres in the gap sweep are `μ ~ N(θ_mean, θ_std²)` with `θ_mean = t2 − ρ·s`
(`config.py:111`), so `ρ` is the **threshold proximity**: distance from `E[μ]` to the active upper
threshold `t2`, in units of `s`.

### 1.3 `τ_stat` — statistical calibration — `minos-core/minos/calibration.py:52`

The scale at which the reported central-`L` interval of `N(μ, (τs)²)` has nominal coverage `L` of
the truth. Coverage (`diagnostics.py:17`) at `δ=0`, posterior-centric, is

```
C(τ) = P( |θ − μ| ≤ z_L · τ · s ) = P( |u| ≤ z_L · τ ),   z_L = Φ⁻¹(½ + L/2),
```

because `θ − μ = s·u`. `τ_stat` is the root of `C(τ) − L = 0` (Brent on `[0.2, 5]`). **It never
reads the utility** (independent code path; the audit confirms `τ_stat` is invariant to `k_under×3`).
Reference level `L_ref = 0.90` (`calibration.py:34`), so `z_L = Φ⁻¹(0.95) ≈ 1.6449`.

### 1.4 `τ*` — decision calibration — `minos-core/minos/calibration.py:81`

`τ* = argmax_τ E[ U( a*(N(μ, (τs)²)), θ ) ]` over the population, by a dense grid + parabola-fit
refine. **It never reads coverage** (independent code path; takes no `L`). Where the decision is
`τ`-insensitive (flat EU, e.g. `λ = 1`, or report centres far from any threshold) it returns `τ*=1`
(`calibration.py:92`, the `FLAT_EU_TOL` / non-concave guard). The gap is `G = τ* − τ_stat`
(`calibration.py:133`).

### 1.5 v3 deployment split + monitor — `generative.py:172`, `correction.py`, `monitor.py`

`realise_deploy` (`generative.py:172`) splits the deployment shift into two orthogonal channels at
the **same** gain `β·s`:

```
OBSERVABLE  δ_obs:  μ = μ₀ − β·s·δ_obs   (reported point moves DOWN; truth θ put)
HIDDEN      δ_hid:  θ = μ₀ + β·s·δ_hid + s·u   (truth moves UP;  reported point μ put)
```

At matched `δ` the two impose the *same* truth↔report discrepancy — the same deployment decision
problem, hence the same stale-correction regret — differing only in whether the **observable**
`μ`-channel moved. `β = 5, s = 0.5 ⇒ β·s = 2.5` (`config.py:55`).

The loss-calibration **baseline** `τ̂_cal = fit_loss_calibration` (`correction.py:46`) is exactly
v2's `τ*` on a labeled calibration set (a cited method: Lacoste-Julien 2011 / Vadera 2021 / Zhao
2021). The **stale-correction regret** is

```
R(δ) = max_τ deploy_EU(τ; δ) − deploy_EU(τ̂_cal; δ) ≥ 0     (correction.py:149)
```

The **label-free monitor** `M` (`monitor.py:144`, default `utility_divergence`) is a function of the
*reported points* `{μ}` only:

```
M = Σ_b ω(z_b) · |p_dep(z_b) − p_cal(z_b)| · dz,   z = (μ − t2)/s,   ω(z) = k_under · φ(z)
```

`ω` is a utility-stakes kernel localised to within ~1 reported sd of `t2` and scaled by the
under-treatment cost — `M` is **regret-targeted, not a density-ratio OOD score** (`monitor.py:78`).
`M` **never reads `θ`** (asserted in code: scrambling `θ` leaves `M` bit-identical).

---

## 2. Code → analytic model (the notation the theorems use)

**Posterior as Gaussian × skew correction.** The truth's law given the report is `θ = μ + s·u`,
so the *true* per-report posterior of `θ` is a location-`μ`, scale-`s` skew-normal. To leading
order in its standardised third cumulant `γ`, its density is the **Gram-Charlier / Edgeworth**
expansion

```
p(θ | μ) = (1/s) · φ((θ−μ)/s) · [ 1 + (γ/6)·He₃((θ−μ)/s) ] + O(γ^{4/3}),   He₃(x) = x³ − 3x.
```

The *reported* posterior is the Gaussian `N(μ, (τs)²)` — at `τ=1` the moment-match (mean and, at
`κ=0`, variance) of the true law. The skew `γ` is the misspecification.

**Skew knob → cumulant.** The skew-normal standardised skewness gives the map (audit §1)

```
γ(κ) = [(4−π)/2] · (d√(2/π))³ / (1 − 2d²/π)^{3/2},   d = κ/√(1+κ²).
```

Verified in-session: `γ(formula)` matches `scipy.stats.skewnorm`'s skewness to 4 dp at every `κ`,
and the *actual code path* `_unit_skew_error(κ=3)` produces a sample with mean `−0.0010`, variance
`0.9996`, skewness `+0.6679` against the target `γ(3) = +0.6670`. The map is monotone:
`κ: 0,0.5,1,2,3,4 → γ: 0, 0.024, 0.137, 0.454, 0.667, 0.784`.

**Utility / thresholds in the gap config.** `t1 = t2 − 20` (`config.py:97`, `GAP_T1_OFFSET`), so only
the **treat/escalate** boundary at `t2` is active; SPARE and the `t1` boundary are inert (audit §2:
`P(SPARE) = 0` exactly at every `τ`). Writing `m = μ − t2` and `x = θ − t2 = m + s·u`, near `t2`:

```
U(TREAT,    θ) ≈ − k_under · (x)_+        (under-treat cost when θ > t2)
U(ESCALATE, θ) = − k_over  · (−x)_+       (over-treat cost when θ < t2)
```

**`τ_stat` in this notation.** `C(τ) = P(|u| ≤ z_L·τ)`; `τ_stat` solves `C(τ_stat) = L`. At `γ=0`,
`u ~ N(0,1) ⇒ τ_stat = 1`.

**`τ*` in this notation.** For a fixed `τ`, the reported rule `a*(μ; τ)` is a **pure threshold on
`μ`**: ESCALATE iff `μ > μ*(τ)`, where `μ*(τ) = t2 + τ·s·z*` and `z* = z*(λ)` solves the
`σ`-independent equation `(λ−1)·ψ(z) + z = 0` with `ψ(z) = z·Φ(z) + φ(z)` (derived from the
EU-difference; `ψ` is the standardised `EPP`). Audit §4: `z*(λ) = 0, −0.276, −0.436, −0.549` for
`λ = 1,2,3,4` — monotone, and `z* < 0` (the boundary sits *below* `t2`: the asymmetric cost makes
escalation fire even somewhat below threshold). Because the family of rules `{ESCALATE iff μ > μ*}`
is swept exactly by `τ`, **choosing `τ` is choosing the threshold `μ*`**, and `τ* = μ*_opt/(s·z*)`
where `μ*_opt` is the realized-utility-optimal threshold.

---

## 3. Assumptions (and where the repo's default regime sits)

1. **Small-skew / leading order.** Expansions are to first order in `γ`. The first neglected terms
   are `O(γ^{4/3})` (skew-normal excess kurtosis ∝ `δ⁴` ∝ `γ^{4/3}`) and `O(γ²)`. The default cell
   sits at `γ(κ=3) = 0.667` — *not* small — so leading order is expected to capture **sign and the
   small-`γ` slope** of the gap, with second-order error at the default cell. CP3 documents the
   regime of validity by sweeping `γ` upward.
2. **Single dominant active threshold.** Only `t2` is active in the gap config (audit §2). ✓
3. **Bounded / Lipschitz utility.** `U` is piecewise-linear, hence Lipschitz in `θ` with constant
   `L_U = max(k_under, k_over) = k_under` (audit §3: `= 3` at the default cell). ✓
4. **`τ*` identifiability.** `τ*` is well-defined only where the EU curve is non-flat — report
   centres within the resolved proximity band. The sweep keeps `ρ ∈ [0, 2]` (`run_b.py:48`); beyond
   that `τ*` defaults to 1 by code. This is a *regime* condition, not a leading-order coefficient
   (see the `ρ` note in §6).

---

## 4. The numbers the theorems must reproduce (from the repo's RESULTS)

**v2 — `RESULTS_B.md` (commit `a94e25f`; seed `20240517`, `RUN_N = 4e6`).** Default misspecified
cell `κ=3, λ=3, ρ=0.5`:

| quantity | value |
|---|---|
| `τ_stat` (nominal 90% coverage) | **0.9635** |
| `τ*` (`argmax_τ E[U]`) | **1.0431** |
| `G = τ* − τ_stat` | **+0.0796** |
| ratio `τ*/τ_stat` | **1.0826** |
| well-specified `κ=0`: `τ_stat, τ*, G` | `0.9998, 0.9959, −0.0039` |
| v1 reproduction `VoC(0.5), VoC(2.0)` | `+0.001169, +0.004361` (= v1 RESULTS.md exactly) |

`G(κ, λ)` at `ρ=0.5` and `G(κ, ρ)` at `λ=3` (the slope checks for CP3) — verbatim from RESULTS_B:

```
G(κ,λ), ρ=0.5:            l=1     l=2     l=3     l=4
                   k=0  +0.001  -0.009  -0.004  +0.004
                   k=1  +0.001  -0.010  +0.005  +0.015
                   k=2  +0.015  +0.016  +0.044  +0.057
                   k=3  +0.036  +0.050  +0.084  +0.100
                   k=4  +0.055  +0.078  +0.115  +0.133

G(κ,ρ), λ=3:             r=0.0   r=0.5   r=1.0   r=1.5   r=2.0
                   k=0  +0.017  -0.004  -0.004  -0.008  +0.004
                   k=3  +0.091  +0.084  +0.082  +0.083  +0.070
                   k=4  +0.120  +0.115  +0.113  +0.114  +0.099
```

Monotone increasing in `κ` and in `λ`; **nearly flat in `ρ` over `[0, 1.5]`** then falling at
`ρ=2` (slope `−0.0084`). `τ_stat(L)` profile @ default: `0.5:1.002, 0.68:0.980, 0.8:0.966,
0.9:0.964, 0.95:0.983, 0.99:1.081` (bulk shrinks, extreme tail widens).

**v3 — `RESULTS_C.md` (commit `72406ed`; seed `20240517`, `CAL_N=2e6, DEP_N=1e6`).** Calibration
cell = the v2 default:

| quantity | value |
|---|---|
| loss-calibration baseline `τ̂_cal` | **1.0480** |
| v2 gap on the cell `τ_stat, τ*, G` | `0.9639, 1.0480, +0.0841` |
| `R(0,0)` (zero-shift regret) | `0.000000` |
| `corr(M, R)` under `δ_obs` | **+0.965** |
| detection AUC `{R>tol}`, `δ_obs` | **1.000** |
| detection AUC `{R>tol}`, `δ_hid` | **0.500** — at chance |
| threshold `m*` (`α=0.05`), zero-shift false-alarm | `0.00876`, `0.062` |

`R` ladders (regret vs each shift knob, comparable by construction):
```
R vs δ_obs: 0.00:0.0000 0.03:0.0039 0.06:0.0147 0.09:0.0320 0.12:0.0561 0.15:0.0845 0.18:0.1164 0.21:0.1500 0.24:0.1838
R vs δ_hid: 0.00:0.0000 0.03:0.0034 0.06:0.0147 0.09:0.0333 0.12:0.0601 0.15:0.0944 0.18:0.1366 0.21:0.1865 0.24:0.2443
```

---

## 5. GATE 0 — feasibility (HALT-ABLE)

| GATE-0 requirement | finding | status |
|---|---|---|
| skew knob maps to a cumulant `γ` | `γ(κ)` = scipy skew-normal skewness to 4 dp; code path matches | **PASS** |
| single dominant active threshold in default regime | `P(SPARE)=0` exactly; only `t2` active | **PASS** |
| bounded / Lipschitz utility | piecewise-linear, `Lip_θ(U) = k_under` | **PASS** |
| posterior = Gaussian × skew correction | `θ = μ + s·u`, Gram-Charlier in `γ` | **PASS** |
| leading-order expansion exists | `τ_stat`, `τ*` both expandable in `γ` about 1 (CP1) | **PASS** |

**GATE 0: PASS.** The code's model can be put in a form admitting a leading-order expansion in the
posterior skew `γ`. Plumbline proceeds with the locked model above.

---

## 6. Two honest refinements the audit forces (carried into the theorems)

These are places where the *actual code* sharpens (not contradicts) the mega-prompt's anticipated
form `G = (a + b·Λ)·γ`. They are stated here so CP1/CP3 test them rather than assume them.

- **`τ_stat` is first-order skew-insensitive: `a = 0`.** `τ_stat` is fixed by a **symmetric**
  (central, two-sided) coverage of the truth about its mean `μ`. The leading Gram-Charlier
  correction to a symmetric probability is the `He₃` term, which is *odd* and integrates to zero
  over `[−c, c]`. Hence `τ_stat = 1 + o(γ)` — the `O(γ)` coefficient `a` vanishes by symmetry, and
  the entire first-order gap comes from the **decision** side. (RESULTS_B's `τ_stat = 0.964` at
  `γ=0.667` is then a second-order effect; the `τ_stat(L)` profile's bulk-vs-tail behaviour is the
  fingerprint.) CP1 derives this; CP3 confirms `(τ_stat−1)/γ → 0`.
- **`Λ` carries utility asymmetry; threshold proximity is a *regime* condition, not a leading
  coefficient.** The realized-utility-optimal boundary obeys a **pointwise** first-order condition
  `k_under·B(μ*) = k_over·A(μ*)` whose population density `f(μ*)` cancels — so to leading order
  `τ*` depends on `λ` (through `z*(λ)`) and `γ`, but **not** on `ρ`. RESULTS_B confirms this
  directly: `G` is flat to ±0.002 over `ρ ∈ [0, 1.5]` at `κ=3`, falling only at `ρ=2` where `τ*`
  loses identifiability. So `Λ ≡` utility-asymmetry factor (`λ`, `z*(λ)`); `ρ` sets *whether* `τ*`
  is resolved and contributes only at second order. CP1 derives the `λ`-dependence; CP3 confirms
  `ρ`-flatness in the resolved band.
