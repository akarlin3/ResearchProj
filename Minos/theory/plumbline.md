# Plumbline — analytic hardening of the Minos decision–calibration findings

*A theory note for the Minos paper (theory section, not a standalone). Every number traces to a
script in this directory; every analytic claim to a sympy derivation in `gap_scaling.py` /
`detectability.py`. The model is locked against the simulation code in `THEORY_MODEL.md`. The
impossibility (Theorem 2, part i) is **drafted and flagged for human proof-review** — see
`impossibility.md` — and is **not** machine-verified.*

---

## 1. Why analytic results

Minos-Core establishes two simulation findings. **(v2)** Statistical and decision calibration of a
reported error bar are different objectives that diverge under a misspecified (skewed) posterior, by
a measurable *decision–calibration gap* `G = τ* − τ_stat`. **(v3)** A loss-calibration correction fit
on a calibration set goes stale under deployment shift; a *label-free validity monitor* recovers the
*observable* fraction of that staleness but is at chance on a *hidden* fraction that induces the same
regret.

Simulations show *that* these effects occur at the parameters tested. They do not, by themselves,
say *why*, *how the effect scales*, or *whether the monitor's blind spot is an artefact of the
estimator or a law*. Plumbline supplies the two missing statements:

- **Theorem 1** turns the gap into a **scaling law**: `G` is, to leading order in the posterior
  skewness `γ`, linear in `γ` with a coefficient that is a pure decision-boundary quantity. The gap
  is a real scaling phenomenon, not a numerical accident of one cell.
- **Theorem 2** turns the monitor's blind spot into a **hard floor**: the observable component of
  staleness is bounded-detectable (`R_obs ≤ L·δ`), and the hidden component is provably undetectable
  by any label-free monitor (best AUC = ½). Label-free monitoring has a principled limit, which is
  the argument for labeled repeatability spot-checks (Echo).

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
`γ(κ=3)=0.667` — **moderate, not small** — so leading order is expected to nail the decision-side
scale and the sign/slope of the gap, with a second-order correction at the default cell (§5).

---

## 3. Theorem 1 — the gap scaling law

> **Theorem 1.** In the locked model, to leading order in the posterior skewness `γ`,
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

**GATE 1 (printed by `gap_scaling.py`).** `gc(z)=−zφ(z)` (sympy); `a=0` (the odd integral is `0`);
`b=1/6`, `Λ(λ)=−z*(λ)` (residual `0`); `G→0` as `γ→0` (linear, no constant); sign `+`; slope
`d(τ*)/dγ` strictly increasing in `λ`: `[0, 0.046, 0.073, 0.092]` for `λ=[1,2,3,4]`
(`z* = [0, −0.276, −0.436, −0.549]`).

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
`Δ = β·s·δ_obs`, so `W₁(P_O^dep,P_O^cal) = Δ` exactly. **`L` printed:** `L_U = k_under = 3.0`;
`L_TV = 2·osc(g_τ) = 0.682`. Both bounds hold across the v3 observable sweep (`δ_obs ∈ [0,0.24]`);
the realized regret is sub-linear (≈ quadratic near `δ=0`), so the linear bound is a conservative
envelope (observed `R`-vs-`W₁` slope `0.32 ≤ L_U = 3`).

### 4.2 Impossibility (hidden component; **DRAFTED — REQUIRES HUMAN PROOF-REVIEW**)

> **Theorem 2(i) — drafted, not machine-verified.** For the hidden channel (`δ_obs=0`), every
> label-free monitor `M = f({μ})` has identical law under *fresh* (`δ_hid=0`) and *stale*
> (`δ_hid>0`) deployment. Hence every detector built from the observables has TPR = FPR at every
> operating point: its ROC is the diagonal and its detection AUC is exactly **½**.

**Argument (data-processing/sufficiency; full text + reviewer notes in `impossibility.md`).** The
hidden shift perturbs only `θ`; the observable law `P_O^stale = P_O^fresh =: P_O` is exactly
invariant (the simulation realises this bit-for-bit: `M(δ_hid=0) = M(δ_hid=0.20) = 0.007248`). Any
`M = f(O)` inherits the invariance (`Law(M|stale) = Law(M|fresh)`, equivalently `I(M;δ_hid) = 0`); by
Neyman–Pearson the likelihood ratio is `≡ 1`, so no test has power and AUC = ½. **This is a
sufficiency proof, not algebra; it carries the human-review flag and is consistent with — not
verified by — the v3 result AUC = 0.500.**

---

## 5. Numerical confirmation (theory vs v2/v3) — `confirm.py`, GATE 3 (HALT-ABLE)

All sim figures below were printed by `confirm.py` this session; v2/v3 references are from
`RESULTS_B.md` / `RESULTS_C.md`. **No theory constant was tuned to the simulation.**

### Theorem 1 — gap (default cell `κ=3, λ=3, ρ=0.5`, `γ=0.667`)

| quantity | theory | sim (this session, 5 seeds @ 2e6) | v2 RESULTS_B |
|---|---|---|---|
| `τ*` (decision scale) | `1.0485` (leading order) | `1.0516 ± 0.0043` | `1.0431` |
| `τ_stat` (leading order) | `1.0000` (`a=0`) | `0.9640 ± 0.0006` | `0.9635` |
| `τ_stat` (exact coverage model) | `0.9638` | `0.9640 ± 0.0006` | `0.9635` |
| `G` (leading order, decision side) | `0.0485` | `0.0876 ± 0.0041` | `0.0796` |
| `G` (`τ*_lead − τ_stat_exact`) | `0.0847` | `0.0876 ± 0.0041` | `0.0796` |

The **decision-side scaling law nails `τ*`** to ~0.3 % at this moderate `γ` (theory `1.0485`;
this session `1.0516`, |Δ|=0.0031; v2 `1.0431`; v3 `τ̂_cal = 1.0480`). The leading-order `τ_stat = 1`
(the `a=0` result) is correct to first order; the observed `−0.036` shrink is a higher-order coverage
effect, reproduced to `0.0002` by the skew-normal coverage model (`0.9638` vs sim `0.9640`).
Combining the two analytic pieces reconstructs the full gap (`0.0847` theory; sim `0.0876`, v2
`0.0796` — theory sits between the two simulation estimates, well within the `τ*`-optimum jitter).

### Theorem 1 — gap vs skew (theory traces the published v2 sweep, `λ=3, ρ=0.5`)

| `κ` | `γ` | `G` (v2 published) | `G` leading-order | `G` full theory | `\|full−v2\|` |
|---|---|---|---|---|---|
| 0.0 | 0.000 | `−0.004` | `+0.000` | `+0.000` | `0.004` |
| 1.0 | 0.137 | `+0.005` | `+0.010` | `+0.012` | `0.007` |
| 2.0 | 0.454 | `+0.044` | `+0.033` | `+0.048` | `0.004` |
| 3.0 | 0.667 | `+0.084` | `+0.049` | `+0.085` | `0.001` |
| 4.0 | 0.784 | `+0.115` | `+0.057` | `+0.113` | `0.002` |

Max deviation across the sweep: **`0.0065 < 0.012`** (tolerance). The leading-order slope at
`γ→0` is `|z*(3)|/6 = 0.0727`. **Regime of validity:** the *pure* leading-order `G` tracks the full
theory within `0.012` up to `γ ≈ 0.14`; beyond that the (separately validated, second-order)
coverage shrink in `τ_stat` must be added — but the **decision-side scale `τ*` and the gap sign /
slope are captured throughout** the tested range.

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
floor (AUC = ½), matching the impossibility.

**GATE 3: PASS.** Theory reproduces v2/v3 within the stated regime/tolerance, with no fitted
constant.

---

## 6. Discussion

**The gap is a scaling law, not a numerical accident.** Theorem 1 shows `G` is, to leading order,
`(1/6)|z*(λ)|·γ` — linear in the posterior skewness with a coefficient that is purely the
decision-boundary offset. Two consequences. First, the gap **must** appear whenever the posterior is
skewed and the cost is asymmetric (`λ>1`); the v2 finding is structural, not particular to one cell.
Second, the coefficient is **interpretable**: it is the under:over cost asymmetry expressed through
the threshold offset `z*(λ)`, and it is independent of where the population sits relative to the
threshold (proximity sets only identifiability). The surprising v2 observation — that statistical
calibration *shrinks* the bar while the decision *widens* it — is explained: the shrink is a
second-order, *symmetric-coverage* effect (`a=0` to first order), while the widening is the
first-order, *decision-side* response to skew.

**Label-free monitoring has a hard floor.** Theorem 2 splits staleness into an observable component,
bounded-detectable with a utility-Lipschitz constant (`R_obs ≤ k_under·W₁`), and a hidden component
that is **undetectable in principle** by any function of the observables (AUC = ½). The monitor is
therefore not "imperfect"; it is *optimal up to a floor that no label-free statistic can cross*. The
only way to see the hidden channel is to read a label — which is exactly the role of periodic
**labeled repeatability spot-checks** (Echo). Theorem 2 is the principled justification for pairing a
cheap label-free monitor with a sparse labeled check, rather than trusting either alone.

**Scope.** Both results are leading-order / idealised: Theorem 1 to first order in `γ` (with the
default cell at `γ=0.667` requiring the named second-order coverage term for the full gap magnitude);
Theorem 2(i) assumes the hidden channel leaves observables exactly invariant (the *definition* of
hidden in the model). The impossibility proof is drafted and **awaits human review**.

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

### Reproduce

```
.venv-theory/bin/python theory/gap_scaling.py     # Theorem 1 — symbolic derivation + GATE 1
.venv-theory/bin/python theory/detectability.py   # Theorem 2(ii) — bound, L, GATE 2
.venv-theory/bin/python theory/confirm.py         # CP3 — theory vs v2/v3, GATE 3 (HALT-able)
```

`THEORY_MODEL.md` — the locked model + GATE 0. `impossibility.md` — Theorem 2(i), drafted, **human
proof-review required**.
