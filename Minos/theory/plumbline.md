# Plumbline — analytic hardening of the Minos decision–calibration findings

*A theory note for the Minos paper (theory section, not standalone). Every number traces to a script
in this directory printed **this session**; every analytic claim to a sympy derivation in
`gap_scaling.py` / `detectability.py`. The model is locked against the simulation code in
`THEORY_MODEL.md`. The gap headline is **error-barred** in `gap_ci.py` (CP1). The impossibility
(Theorem 2, part i) is **drafted and flagged for human proof-review** — see `impossibility.md` — and
is **not** machine-verified.*

---

## 1. Why analytic results

Minos-Core establishes two simulation findings. **(v2)** Statistical and decision calibration of a
reported error bar are different objectives that diverge under a misspecified (skewed) posterior, by
a measurable *decision–calibration gap* `G = τ* − τ_stat`. **(v3)** A loss-calibration correction fit
on a calibration set goes stale under deployment shift; a *label-free validity monitor* recovers the
*observable* fraction of that staleness but is at chance on a *hidden* fraction that induces the same
regret.

Simulations show *that* these effects occur at the parameters tested. They do not, by themselves,
say *why*, *how the effect scales*, *whether the single-seed headline survives error bars*, or
*whether the monitor's blind spot is an artefact of the estimator or a law*. Plumbline supplies the
missing statements:

- **Theorem 1** turns the gap into a **scaling law**: to leading order in the posterior skewness `γ`,
  `G` is linear in `γ` with a coefficient that is a pure decision-boundary quantity. The gap is a
  real scaling phenomenon, not a numerical accident of one cell.
- **CP1** then **error-bars the headline**: `τ*` is a *shallow* optimum, so the published single-seed
  numbers are noisy. Multi-seed CIs (`gap_ci.py`) show the v2 headline — `τ_stat < 1 < τ*` (the two
  scales on **opposite sides of 1**) — is **robust**, not a lucky draw.
- **Theorem 2** turns the monitor's blind spot into a **hard floor**: the observable component of
  staleness is bounded-detectable (`R_obs ≤ L·δ`), and the hidden component is provably undetectable
  by any label-free monitor (best AUC = ½), with the observable/hidden split now defined cleanly for
  *general* shifts (not just the by-construction Minos channels). Label-free monitoring has a
  principled limit, which is the argument for labeled repeatability spot-checks (Echo).

---

## 2. The model (locked; full detail in `THEORY_MODEL.md`)

Posterior-centric gap regime: report centres `μ ~ N(t2 − ρs, θ_std²)`; the truth is `θ = μ + s·u`
with `u` a zero-mean, unit-variance **skew-normal** error of shape `κ`, whose standardised third
cumulant is

```
γ(κ) = [(4−π)/2] · (d√(2/π))³ / (1 − 2d²/π)^{3/2},   d = κ/√(1+κ²).
```

To leading order the true per-report posterior of `θ` is the Gram-Charlier expansion
`p(θ|μ) = (1/s)φ((θ−μ)/s)[1 + (γ/6)He₃((θ−μ)/s)]`, `He₃(x)=x³−3x`. The reported posterior is the
Gaussian `N(μ,(τs)²)`. The piecewise-linear utility has thresholds `t1<t2` and slopes `k_over`
(over-treatment) and `k_under` (under-treatment); the asymmetry is `λ = k_under/k_over ≥ 1`. In the
gap regime `t1` is pushed far below `t2`, so the **only active threshold is `t2`** (the treat/escalate
boundary — verified: `P(SPARE)=0` exactly). `τ_stat` is the coverage scale; `τ*` the decision scale;
`G = τ* − τ_stat`.

**Assumptions (and where the default cell sits).** Expansions are leading-order in `γ`; the first
neglected terms are `O(γ^{4/3})` (kurtosis) and `O(γ²)`. One dominant active threshold. Bounded,
Lipschitz utility (`Lip_θ U = k_under`). `τ*` is identified only where the EU curve is non-flat
(report centres within the resolved proximity band `ρ ∈ [0,2]`). The default cell sits at
`γ(κ=3)=0.667` — **moderate, not small** — so leading order nails the decision-side scale and the
sign/slope of the gap, with a named second-order correction at the operating point (§3.1, §5).

---

## 3. Theorem 1 — the gap scaling law

> **Theorem 1 (leading order in `γ`).** In the locked model, to leading order in the posterior
> skewness `γ`,
>
> ```
> τ_stat(γ) = 1 − a·γ + O(γ²),        a = 0
> τ*(γ)     = 1 + b·Λ(λ)·γ + O(γ^{4/3}),   b = 1/6,  Λ(λ) = −z*(λ) = |z*(λ)|
> G(γ)      = τ* − τ_stat = (1/6)·|z*(λ)|·γ + O(γ^{4/3}),
> ```
>
> where `z*(λ) < 0` is the decision-boundary offset (in reported sd) solving
> `(λ−1)·ψ(z) + z = 0`, with `ψ(z) = z·Φ(z) + φ(z)` the standardised expected-positive-part. The
> gap is positive (the asymmetric cost *widens* the bar, `τ* > 1 = τ_stat`), vanishes as `γ → 0`,
> and is monotone increasing in `γ` and in the cost asymmetry `λ`.

**Proof (derived symbolically in `gap_scaling.py`).**

*`τ_stat`: the coefficient `a = 0`.* Coverage is the **symmetric** probability
`C(τ) = P(|u| ≤ z_L·τ)`. Its `O(γ)` Gram-Charlier correction is
`(1/6)∫_{−c}^{c} φ(u)He₃(u)du`; `He₃` is **odd** and the interval is symmetric about the mean, so the
integral is identically zero (sympy returns `0`). Hence the coverage condition `C(τ_stat)=L` gives
`τ_stat = 1 + o(γ)`: the calibration scale is *first-order skew-insensitive*, and the **entire
first-order gap is decision-side**.

*`τ*`: the slope `b·Λ`.* For fixed `τ` the reported rule is a pure threshold on `μ`:
ESCALATE iff `μ > μ*(τ) = t2 + τ·s·z*(λ)`. So choosing `τ` is choosing the threshold; the
realized-utility-optimal threshold `z_opt(γ)` solves the first-order condition
`h(z;γ) = (λ−1)·g(z;γ) + z = 0`, `g(z;γ) = E_u[(z+u)_+] = ψ(z) + (γ/6)·gc(z) + O(γ²)`. Sympy
evaluates the skew correction in closed form:

```
gc(z) = ∫_{−z}^{∞} (z+u)·He₃(u)·φ(u) du = −z·φ(z).
```

The implicit-function theorem gives
`dz_opt/dγ|₀ = −[(λ−1)gc(z*)/6] / [(λ−1)Φ(z*)+1]`, and since `τ* = z_opt/z*`,

```
dτ*/dγ|₀ = (λ−1)·φ(z*) / (6·[(λ−1)Φ(z*)+1]) = b·Λ(λ),   b = 1/6.
```

Substituting the boundary condition `(λ−1)ψ(z*) = −z*` collapses `Λ` to a pure threshold quantity
(sympy-verified):

```
Λ(λ) = (λ−1)φ(z*) / [(λ−1)Φ(z*)+1] = −z*(λ).
```

So **`G(γ) = (1/6)·|z*(λ)|·γ`**: the gap's slope in skewness is the magnitude of the decision
boundary offset, over six. `Λ` carries the **utility asymmetry** (via `λ` and `z*(λ)`); the threshold
proximity `ρ` does **not** enter the leading coefficient (the optimality condition is pointwise, so
the population density cancels) — it sets only whether `τ*` is *identified* (and a second-order
magnitude). At `λ=1`, `z*=0` ⇒ `Λ=0` ⇒ `τ*=1`: a symmetric cost has no gap, matching the code's
flat-EU guard. ∎

**GATE 1 (printed by `gap_scaling.py` this session).** `gc(z)=−zφ(z)` (sympy); `a=0` (the odd integral
is `0`); `b=1/6`, `Λ(λ)=−z*(λ)` (residual `0`); `G→0` as `γ→0` (linear, no constant); sign `+`; slope
`d(τ*)/dγ` strictly increasing in `λ`: `[0, 0.046, 0.073, 0.092]` for `λ=[1,2,3,4]`
(`z* = [0, −0.276, −0.436, −0.549]`).

### 3.1 Scope: this is a small-`γ` law; at the operating point it is ~half the gap

The clean line `G = (1/6)|z*|γ` is **leading order**. Two numbers fix its scope (computed this
session; see §5 and Fig 1):

- **Regime of validity.** The *pure* leading-order `G_lead = (1/6)|z*|γ` tracks the full theory
  `G_full = τ*_lead − τ_stat^exact` to within `0.012` only up to **`γ ≈ 0.14`** (`κ=1`). Beyond that
  the second-order coverage shrink in `τ_stat` is no longer negligible and must be added.
- **At the operating point `γ = 0.667`,** the leading decision term is `G_lead = 0.0485`, which is
  only **≈57 %** of the full analytic gap `G_full = 0.0847` (and **≈55 %** of the CP1 error-barred gap
  `G = 0.0876`). The remaining **≈43 %** is the coverage shrink `1 − τ_stat^exact = 0.0362` — a
  *separately validated* second-order effect (`τ_stat^exact = 0.9638`, reproducing the sim
  `0.9640±0.0006` to `0.0002`), **not** carried by the `a=0` leading line.

So the clean theorem explains the gap's **sign, slope, and existence as a scaling law**, but the
leading line must **not** be quoted as the operating-point magnitude. `G_lead` and `G_full` are shown
side by side in §5 and Fig 1; the headline magnitude is the error-barred `G` of §3.2 / CP1.

### 3.2 CP1 — the headline is error-barred and robust

`τ*` is the argmax of a *shallow* expected-utility optimum, so a single-seed estimate is noisy. The
v2 headline rests on `τ* > 1` **and** `G > 0`; both are error-barred in `gap_ci.py` by multi-seed
Monte-Carlo replicates (the published estimator re-run on `B = 64` independent populations at
`N = 10⁶`, with a percentile bootstrap of `20 000` resamples, seed `73019`). At the default cell
(`κ=3, λ=3, ρ=0.5`):

| quantity | multi-seed mean | 95 % CI (mean) | single-run 95 % PI | published 1-seed | robust? |
|---|---|---|---|---|---|
| `τ_stat` | `0.9639` | `[0.9636, 0.9641]` | `[0.9619, 0.9658]` | `0.9635` | **< 1 (shrinks): yes** |
| `τ*` | `1.0514` | `[1.0493, 1.0536]` | `[1.0346, 1.0683]` | `1.0431` | **> 1 (widens): yes** |
| `G` | `0.0876` | `[0.0855, 0.0897]` | `[0.0711, 0.1041]` | `0.0796` | **> 0 (gap): yes** |

> **CP1 verdict (GATE 1, halt-to-report).** *"Opposite sides of 1" is **ROBUST.*** All three honesty
> checks pass at the level of the CI on the true value **and** at the level of a single run (the
> prediction interval also excludes the null). The finite-N bias is negligible (the mean is
> N-stable to `0.0002` between `N=10⁶` and `N=2·10⁶`), so the CI is a statement about the true value,
> not a finite-sample artefact. The true `τ* ≈ 1.051` in fact sits **above** the published single
> seed `1.0431`, so the published headline was a representative — if slightly low — draw, not a fluke.
> The estimator was **not** tuned to recover a cleaner claim.

Across the skew axis (`gap_ci.py` sweep, `B = 32` at `N = 5·10⁵`; Fig 1), the gap **switches on with
`γ`**: at `κ=0` it is `G = +0.0046` — negligible (true `G=0` by construction; the `+0.005` is a small
estimator bias at the flat optimum, `≈20×` below the operating-point gap and of the same tiny scale
as the published single-seed `−0.004`, sign undetermined) — rising monotonically to `G = 0.0892` at
`κ=3`, with the CI excluding `0` for every `κ ≥ 1`. The analytic slope `dG/dγ = |z*|/6 = 0.0727 > 0`
is sympy-derived and not subject to the optimum's noise — the noise-free fallback for the gap's sign.

---

## 4. Theorem 2 — label-free detectability

### 4.1 Achievable bound (observable component; gated, `detectability.py`)

> **Theorem 2(ii).** Let `P_O^cal, P_O^dep` be the calibration and deployment laws of the observable
> reported point `μ`. The stale-correction regret of the observable channel obeys
>
> ```
> R_obs ≤ L · δ(P_O^dep, P_O^cal),
> ```
>
> with two consistent forms: **(W₁)** `δ = W₁` and `L = L_U = k_under` (the utility's Lipschitz
> constant in `θ`); **(TV)** `δ = TV` and `L = 2·osc(g_τ)` (twice the oscillation of the per-report
> expected utility `g_τ(μ) = E[U(a*(μ;τ),θ)|μ]`).

**Proof sketch (`detectability.py`).** Optimality chain: `τ̂` is optimal at calibration, so
`R(P_dep) ≤ 2·sup_τ |EU(τ;P_dep) − EU(τ;P_cal)| =: 2ε`. Bounding `ε`: writing
`EU(τ;P) = ∫ g_τ dP(μ)`, `ε ≤ osc(g_τ)·TV(P_O)` (TV form); and since `U` is `L_U`-Lipschitz in `θ`
while the observable shift only displaces the decision boundary by `Δ` in report-space,
`ε ≤ L_U·W₁(P_O)` (W₁ form). For the observable shift `P_O^dep` is `P_O^cal` translated by
`Δ = β·s·δ_obs`, so `W₁(P_O^dep,P_O^cal) = Δ` exactly. **`L` printed this session:** `L_U = k_under =
3.0`; `L_TV = 2·osc(g_τ) = 0.682`. Both bounds hold across the v3 observable sweep
(`δ_obs ∈ [0,0.24]`); the realized regret is sub-linear (≈ quadratic near `δ=0`), so the linear bound
is a conservative envelope (observed `R`-vs-`W₁` slope `0.32 ≤ L_U = 3`; Fig 2).

### 4.2 Impossibility (hidden component; **DRAFTED — REQUIRES HUMAN PROOF-REVIEW**)

> **Theorem 2(i) — drafted, not machine-verified.** For any *purely hidden* shift
> (`P_O^dep = P_O^cal`; the Minos `δ_obs=0` channel realises this), every label-free monitor
> `M = f({μ})` has identical law under *fresh* and *stale* deployment. Hence every detector built from
> the observables has TPR = FPR at every operating point: its ROC is the diagonal and its detection
> AUC is exactly **½**.

**Pinned definitions and the partial-leak proposition (`impossibility.md`, hardened CP2).** The note
fixes the formalisation the proof rests on, so a reviewer evaluates a definition rather than a buried
assumption:

- **The observable σ-algebra** `O = σ({μ_i}, fixed known constants)`; a **label-free monitor** is any
  `O`-measurable `M = f({μ})` (it may not read the truths `{θ_i}`).
- **"Hidden" is a *definition*, not a derived property:** the hidden component of a deployment shift
  is its `O`-**invariant** component — the part that leaves `P_O` (hence every label-free statistic's
  law) unchanged. The Minos `δ_hid` channel is then simply the canonical purely-hidden shift
  (`P_O^dep = P_O^cal` *by construction*).
- **Partial-leak proposition (general shifts).** Factor any shift `cal → dep` through the intermediate
  law `P_mid = P_O^dep ⊗ K_cal` (deployment's observable marginal carried with calibration's `θ|μ`
  kernel): step **A** (`cal→mid`) moves only `P_O` — the **observable component**, on which the
  monitor's discriminating power lives and whose regret obeys part (ii)'s envelope **verbatim**
  (because `g_τ` stays the calibration `g_τ`); step **B** (`mid→dep`) moves only `θ|μ` at fixed
  `P_O` — the **hidden component**, which leaves `M`'s law invariant, so part (i) gives AUC = ½ for it.
  *Net:* a label-free monitor is **boundedly good on exactly the part of any shift that touches `O`,
  and provably powerless on exactly the part that does not.* This answers the "real shifts leak"
  objection: leakage *is* component A (and is what makes a shift detectable, priced by part ii); the
  irreducible `O`-invariant remainder is component B, undetectable. Both bounds are cross-referenced
  to parts (i)/(ii).

**Argument (data-processing/sufficiency).** The hidden shift perturbs only `θ`; the observable law
`P_O^stale = P_O^fresh =: P_O` is exactly invariant (the simulation realises this bit-for-bit:
`M(δ_hid=0) = M(δ_hid=0.20) = 0.007248`, verified this session). Any `M = f(O)` inherits the
invariance (`Law(M|stale) = Law(M|fresh)`, equivalently `I(M;δ_hid) = 0`); by Neyman–Pearson the
likelihood ratio is `≡ 1`, so no test has power and AUC = ½. **This is a sufficiency proof, not
algebra; it carries the human-review flag and is consistent with — not verified by — the v3 result
AUC = 0.500.**

---

## 5. Numerical confirmation (theory vs v2/v3) — `confirm.py`, GATE 3 (HALT-ABLE)

All sim figures below were printed by `confirm.py` (and the error-barred gap by `gap_ci.py`) **this
session**; v2/v3 references are from `RESULTS_B.md` / `RESULTS_C.md`. **No theory constant was tuned
to the simulation.**

### Theorem 1 — gap (default cell `κ=3, λ=3, ρ=0.5`, `γ=0.667`): `G_lead` vs `G_full`

| quantity | theory | sim (this session) | v2 RESULTS_B | CP1 error-barred (95 % CI) |
|---|---|---|---|---|
| `τ*` (decision scale) | `1.0485` (leading order) | `1.0516 ± 0.0043` | `1.0431` | `1.0514` `[1.0493, 1.0536]` |
| `τ_stat` (leading order, `a=0`) | `1.0000` | `0.9640 ± 0.0006` | `0.9635` | `0.9639` `[0.9636, 0.9641]` |
| `τ_stat` (exact coverage model) | `0.9638` | `0.9640 ± 0.0006` | `0.9635` | — |
| **`G_lead`** (leading-order, decision side) | **`0.0485`** | `0.0876 ± 0.0041` | `0.0796` | — |
| **`G_full`** (`τ*_lead − τ_stat^exact`) | **`0.0847`** | `0.0876 ± 0.0041` | `0.0796` | `0.0876` `[0.0855, 0.0897]` |

The **decision-side scaling law nails `τ*`** to `0.0031` at this moderate `γ` (theory `1.0485`; sim
`1.0516`; v2 `1.0431`; v3 `τ̂_cal = 1.0480`). The leading-order `τ_stat = 1` (`a=0`) is correct to
first order; the observed `−0.036` shrink is a higher-order coverage effect, reproduced to `0.0002`
by the skew-normal coverage model (`0.9638` vs sim `0.9640`). **`G_lead` (`0.0485`) is only ≈57 % of
`G_full` (`0.0847`)** at the operating point — the clean line carries the gap's *slope*, not its
operating magnitude; `G_full` (and the CP1 error-barred `0.0876`) carry the magnitude. Theory `G_full`
sits within the `τ*`-optimum jitter of both simulation estimates (sim `0.0876`, v2 `0.0796`;
`|G_full − sim| = 0.0029`).

### Theorem 1 — gap vs skew (theory traces the published v2 sweep, `λ=3, ρ=0.5`; Fig 1)

| `κ` | `γ` | `G` v2 (pub) | `G_lead` | `G_full` | `\|full−v2\|` | `G` CP1 (95 % CI) |
|---|---|---|---|---|---|---|
| 0.0 | 0.000 | `−0.004` | `+0.000` | `+0.000` | `0.004` | `+0.005 [+0.001,+0.008]` |
| 1.0 | 0.137 | `+0.005` | `+0.010` | `+0.012` | `0.007` | `+0.015 [+0.011,+0.018]` |
| 2.0 | 0.454 | `+0.044` | `+0.033` | `+0.048` | `0.004` | `+0.050 [+0.047,+0.054]` |
| 3.0 | 0.667 | `+0.084` | `+0.049` | `+0.085` | `0.001` | `+0.089 [+0.085,+0.093]` |
| 4.0 | 0.784 | `+0.115` | `+0.057` | `+0.113` | `0.002` | `+0.120 [+0.116,+0.124]` |

Max deviation of `G_full` from the v2 sweep: **`0.0065 < 0.012`**. The leading-order slope at `γ→0` is
`|z*(3)|/6 = 0.0727`. **Regime of validity (§3.1):** the *pure* leading-order `G_lead` tracks `G_full`
within `0.012` up to `γ ≈ 0.14`; beyond, the second-order coverage shrink in `τ_stat` is added — but
the **decision-side scale `τ*` and the gap sign / slope are captured throughout**, and the CP1 CI
excludes `0` for every `κ ≥ 1`.

### Theorem 2 — monitor / bound (calibration cell = v2 default)

| quantity | theory / bound | sim (this session) | v3 RESULTS_C |
|---|---|---|---|
| loss-calibration baseline `τ̂_cal` | `≈ 1.048` (= `τ*`) | `1.0543` | `1.0480` |
| `R_obs ≤ L_U·W₁`, `L_U = k_under` | holds, `L_U = 3.0` | holds across `δ_obs∈[0,0.24]` | — |
| observed `R`-vs-`W₁` slope | `≤ L_U = 3` | `0.32` (`detectability.py`) | — |
| observable `corr(M, R)` | `> 0` (tracks) | `+0.964` | `+0.965` |
| observable detection AUC `{R>tol}` | high | `1.000` | `1.000` |
| hidden detection AUC `{R>tol}` | **½** (Thm 2(i)) | `0.500` | `0.500` |

The observable channel tracks regret and is detectable; the hidden channel sits at the theoretical
floor (AUC = ½) while inducing comparable regret (Fig 2), matching the impossibility.

**GATE 3: PASS.** Theory reproduces v2/v3 within the stated regime/tolerance, with no fitted constant
(default-cell, sweep, and monitor sub-gates all PASS).

---

## 6. Discussion

**The gap is a scaling law, not a numerical accident.** Theorem 1 shows `G` is, to leading order,
`(1/6)|z*(λ)|·γ` — linear in the posterior skewness with a coefficient that is purely the
decision-boundary offset. Two consequences. First, the gap **must** appear whenever the posterior is
skewed and the cost is asymmetric (`λ>1`); the v2 finding is structural, not particular to one cell —
and CP1 shows the `τ_stat < 1 < τ*` headline survives error bars (it is not a single-seed artefact of
a shallow optimum). Second, the coefficient is **interpretable**: it is the under:over cost asymmetry
expressed through the threshold offset `z*(λ)`, independent of where the population sits relative to
the threshold (proximity sets only identifiability). The surprising v2 observation — that statistical
calibration *shrinks* the bar while the decision *widens* it — is explained: the shrink is a
second-order, *symmetric-coverage* effect (`a=0` to first order), while the widening is the
first-order, *decision-side* response to skew. At the operating point these contribute ≈43 % and
≈57 % of the gap respectively (§3.1).

**Label-free monitoring has a hard floor.** Theorem 2 splits staleness into an observable component,
bounded-detectable with a utility-Lipschitz constant (`R_obs ≤ k_under·W₁`), and a hidden component
that is **undetectable in principle** by any function of the observables (AUC = ½) — and the CP2
partial-leak proposition shows this split is not a property of the engineered Minos channels but of
*any* shift, partitioned by the observable σ-algebra `O`. The monitor is therefore not "imperfect";
it is *optimal up to a floor that no label-free statistic can cross*. The only way to see the hidden
channel is to read a label — which is exactly the role of periodic **labeled repeatability
spot-checks** (Echo). Theorem 2 is the principled justification for pairing a cheap label-free monitor
with a sparse labeled check, rather than trusting either alone: the monitor prices everything that
leaks into observables, and the spot-check is the *only* instrument that can reach the residue that
does not.

**Scope.** Both results are leading-order / idealised: Theorem 1 to first order in `γ` (with the
default cell at `γ=0.667` requiring the named second-order coverage term for the full gap magnitude —
the clean line is ≈57 % of the gap there); Theorem 2(i) assumes the hidden component leaves observables
exactly invariant (now the *definition* of hidden, generalised to a component of any shift). The
impossibility proof is drafted and **awaits human review**.

---

## 7. Positioning

The loss-calibration baseline Minos cites — tuning a posterior so the *decision* is optimal — is
**loss-calibrated Bayesian inference** (Lacoste-Julien, Huszár & Ghahramani, *Approximate inference
for the loss-calibrated Bayesian*, AISTATS 2011) and its post-hoc form (Vadera et al., *Post-hoc loss
calibration*, UAI 2021), and **decision calibration** (Zhao, Kim, Sahoo, Ma & Ermon, NeurIPS 2021).
Plumbline does not propose a new calibration method; it (i) proves *how far* the decision-optimal
scale departs from the coverage-optimal one (the gap scaling law — a quantity these lines do not
compute, since they each target a single object), and (ii) proves the *detectability limit* of
monitoring such a correction for staleness. The decision-curve / net-benefit literature
(Vickers & Elkin, *Decision curve analysis*, 2006) prices a single threshold's net benefit but not
the statistical-vs-decision *gap*; the OOD / distribution-shift-under-uncertainty literature
(Ovadia et al., *Can you trust your model's uncertainty?*, NeurIPS 2019) studies calibration under
shift empirically but not the label-free detectability floor. Theorem 2(i) is, to our knowledge, the
first statement that the hidden component of correction-staleness is information-theoretically
invisible to any label-free monitor — the formal case for labeled repeatability spot-checks.

---

### Figures

- **`figures/fig1_scaling_law.pdf`** — `G` vs `γ`: full-theory curve, leading-order line, v2 empirical
  points, and the CP1 error-barred gap (95 % CI), with the `γ ≲ 0.14` leading-order regime shaded and
  the operating point `γ=0.667` marked.
- **`figures/fig2_monitor_bound.pdf`** — `R_obs` vs `W₁` under the `L_U·W₁` envelope; observable vs
  hidden channel with an inset on the regret scale, annotating observable AUC `1.00` vs hidden `0.50`.

### Reproduce

```
.venv-theory/bin/python theory/gap_scaling.py        # Theorem 1 — symbolic derivation + GATE 1
.venv-theory/bin/python theory/gap_ci.py             # CP1 — bootstrap CIs on the gap (HALT-TO-REPORT)
.venv-theory/bin/python theory/detectability.py      # Theorem 2(ii) — bound, L, GATE 2
.venv-theory/bin/python theory/confirm.py            # CP3 — theory vs v2/v3, GATE 3 (HALT-able)
.venv-theory/bin/python theory/figures/make_figures.py   # CP4 — publication figures (vector PDF)
```

`THEORY_MODEL.md` — the locked model + GATE 0. `impossibility.md` — Theorem 2(i), drafted, **human
proof-review required**.
