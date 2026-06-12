# POSITIONING — what Minos is *not*

Stub. No results claims here — those live in `RESULTS.md`. This file fixes the
must-distinguish neighbours so the novelty (pricing the **error bar** via VoC and VoTG) is
not silently reframed as one of them.

Minos prices **the per-voxel error bar itself** — its *width* (Value of Calibration) and its
*trustworthiness under shift* (Value of the Trust-Gate) — on a treat / spare / escalate
decision. The object being priced is the uncertainty, not the point estimate, and not a
population parameter.

**v2 headline — the decision–calibration gap.** Minos computes two scales for a reported error
bar by **independent** criteria: `τ_stat`, the scale with nominal interval *coverage* (statistical
calibration), and `τ*`, the scale that maximises *expected decision utility* (decision
calibration). Under a misspecified (skewed) posterior these **diverge**, by a quantity
`G = τ* − τ_stat` — the **decision–calibration gap**. None of the neighbours below computes `G`;
each prices a *single* object and has no notion of a *distance between* statistical and decision
calibration of the uncertainty.

| Neighbour | What it prices | How Minos differs |
|---|---|---|
| **Vickers decision-curve analysis / net benefit** | The decision value of a **marker's point prediction** across threshold probabilities; "is acting on this classifier better than treat-all / treat-none?" | Minos holds the point estimate fixed and prices the **calibration and trustworthiness of the error bar around it** (VoC, VoTG). Net benefit has no notion of the uncertainty being mis-scaled or untrustworthy, and **no notion of `τ_stat` vs `τ*`** — the decision–calibration gap `G` is identically invisible to it. |
| **ISPOR VoI canon (EVPI / EVPPI / EVSI)** | The population value of **resolving or reducing parameter uncertainty** by collecting more data — what a perfect/partial/sample study is worth. | Minos prices a **per-voxel, already-reported** error bar at decision time: VoC is the loss from the reported bar being mis-*scaled*; `G` measures how far *decision*-optimal scaling departs from *statistically*-calibrated scaling under misspecification. No data-collection or parameter-learning step is being valued, and VoI never contrasts coverage-calibration with decision-calibration. (Minos carries an EVPI-*analog* as a degenerate-limit sanity anchor only.) |
| **ARCliDS adaptive-RT CDSS** | An end-to-end **clinical decision-support / RL policy** that recommends adaptive radiotherapy actions from patient state. | Minos is not a policy-learning system or a CDSS; it is a **formal accounting of the decision-value of an uncertainty estimate**, including **when statistical calibration is the wrong objective for the decision** (`G ≠ 0`). It asks "what is this error bar worth, when should it be distrusted, and does calibrating it for coverage help or hurt the decision?", independent of how any policy or estimator was trained. |

**One-line guard.** If an output of Minos is ever described as "net benefit," "decision-curve
analysis," or "EVPI/EVPPI of a parameter," it has been mis-framed: the priced object is the
error bar (VoC / VoTG) and the **gap between its statistical and decision calibration** (`G`),
not the point estimate and not a population parameter.

---

## v3 — a deployment-validity monitor for a stale loss-calibration correction

**The correction is a cited baseline, not our contribution.** The prescriptive move "tune the
reported error bar so the *decision* is optimal" is an existing field — **loss-calibrated Bayesian
inference** and **decision calibration**. Our v2 decision-calibrated scale `τ*` is a 1-D instance of
it, and v3 names it as such: the baseline `τ̂*_cal = argmax_τ Û(τ; D_cal)` (`correction.py`).

- Lacoste-Julien, Huszár & Ghahramani (2011), *Approximate inference for the loss-calibrated
  Bayesian*, AISTATS — the loss-calibrated objective.
- Cobb, Roberts & Gal (2018), *Loss-calibrated approximate inference in Bayesian neural networks*.
- Kušmierczyk, Sakaya & Klami (2019), *Variational Bayesian decision-making for continuous utilities*,
  NeurIPS.
- Vadera, Ghosh, Ng & Marlin (2021), *Post-hoc loss-calibration*, UAI — the **post-hoc** scalar form
  we reduce to.
- Zhao, Kim, Sahoo, Ma & Ermon (2021), *Calibrating predictions to decisions: a novel approach to
  multi-class calibration* (decision calibration), NeurIPS.

**v3 contributes the thing that line does not address: a *deployment-validity monitor*.** All of the
above assume the **calibration set represents deployment**. v3 supplies a **label-free** statistic
`M(D_dep)` that detects when the learned correction has gone **stale under distribution shift**, a
**gated-recovery** policy that acts on it, and the explicit characterisation of *what label-free
monitoring can and cannot catch* (the observable/hidden split).

**Two hard boundaries (must not be silently crossed):**

1. **Not a new calibration method.** Loss-calibration / decision calibration own that; Minos *cites
   and benchmarks against* it. `M` never proposes a better scale — it scores the *staleness of the
   cited correction's decisions*. If an output of v3 is described as "a calibration method," it is
   mis-framed.
2. **Not a generic OOD detector.** That field is saturated. `M` is **regret-targeted**, not
   density-targeted: the deployment↔calibration divergence is weighted by the **utility stakes near
   the decision threshold** (`ω(z)=k_under·ϕ(z)`), so distributional change *far from* the threshold —
   which cannot move the decision-optimal scale — is deliberately ignored. A density-ratio / typicality
   OOD score weights by input density and fires on exactly those irrelevant changes. If an output of
   v3 is described as "an OOD/novelty detector," it is mis-framed.

**The honesty boundary (stated as a result, not hidden).** A label-free monitor sees only observable
reported summaries. v3's shift generator exposes two knobs that induce the **same** regret two ways:
`delta_obs` moves the observable reported point (detectable; AUC ≫ 0.5) and `delta_hid` changes the
truth↔report relationship with the observable summaries held fixed (**undetectable; AUC ≈ 0.5 by
construction**). v3 reports the hidden-staleness fraction as a documented limitation — the honest
claim is *"the monitor catches observable-driven staleness; the hidden fraction requires labeled
repeatability spot-checks,"* which is what motivates periodic labeled validation rather than pretending
label-free monitoring is complete.

| Neighbour (v3) | What it is | How Minos differs |
|---|---|---|
| **Loss-calibrated / decision calibration** (refs above) | Methods that *learn* a decision-optimal posterior/scale on a calibration set | Minos does **not** learn a calibration; it **monitors an already-fitted one for staleness under shift** and recovers utility when it is stale. The correction is the baseline, `M` is the new object. |
| **Generic OOD / covariate-shift detection** (density-ratio, typicality, conformal novelty) | Flags inputs atypical of training, by **input density** | `M` is **regret-targeted**: weighted by utility stakes at the decision boundary, deliberately blind to far-field density change that cannot alter the optimal decision. It answers "is *this correction's decision* stale?", not "is this input unusual?". |
