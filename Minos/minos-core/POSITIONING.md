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
