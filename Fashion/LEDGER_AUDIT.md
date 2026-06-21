# LEDGER_AUDIT.md — railing-first redraft vs the Gnomon claims ledger

Confirms each disposition in `Gnomon/handoff/CLAIMS_LEDGER.md` was honored in the
railing-first redraft for NMR in Biomedicine. Every KEEP/REFRAME carries run-evidence
+ CIs + reworded text; every DROP/OUT-OF-SCOPE is absent from the manuscript body.
Numbers verified against a fresh `build_handoff.py` run (seed 20260620) and frozen in
`NUMBERS_FROZEN.txt`.

| Ledger item | Disposition | Honored? | Where in manuscript | Run-evidence (frozen) |
|---|---|---|---|---|
| **K1** NLLS D\* railing, real data | KEEP (lead) | ✅ | Abstract; §3.1; Fig 1; Table S2 | 54.2% [52.0,56.4] (1e-3); 56.2% (1e-2); 58.7% (SNR>25); prior 54.7% in CI |
| **K2** quantile restores *marginal* D\* coverage | KEEP | ✅ | §3.3; Fig 3A | quantile D\* 0.90 [0.89,0.91]; D 0.95; f 0.96 |
| **K3** flow beats railed NLLS | KEEP | ✅ | §3.3; Fig 3B | cov 0.98 vs 0.76; ECE 0.069 vs 0.121; sharp 0.112 vs 0.181; gaps' CIs exclude 0 |
| **R1** Gaussian under-coverage → conditional | REFRAME | ✅ | §3.2; Fig 2; Table 1; Table S3 | high-D\* honest 0.63 [0.60,0.67] Laplace / 0.81 [0.78,0.84] MCMC; pooled 0.80/0.90 |
| **R2** quantile residual high-D\* wall | REFRAME | ✅ | §3.3; Table 1/S3 | quantile high-D\* 0.81 [0.78,0.84] |
| **R3** below-floor → SD convention | REFRAME | ✅ | §2.4 (Methods); §3.2 | floored Laplace pooled 0.68 / high 0.41; honest recommended |
| **D1** marginal 0.30/0.67 as headline | DROP | ✅ | Appears ONLY in §3.2 as *retired*, naming the floored convention | not a headline anywhere; honest marginal 0.80/0.90 stated |
| OOD self-consistency gate (AUC 0.99/0.59) | OUT OF SCOPE | ✅ absent | §4 out-of-scope note only | not run by Gnomon; not claimed |
| Timing (418 ms / 0.956 s / 635 s) | OUT OF SCOPE | ✅ absent | — | not run by Gnomon; not claimed |
| In-vivo brain held-out-b (0.03) | OUT OF SCOPE | ✅ absent | §4 out-of-scope note only | not run by Gnomon; not claimed |

## Spine decision
The ledger marks K1 and K2/K3 as primary-claim candidates and fixes no order. This
redraft chooses the **boundary-railing-first** spine: K1 leads (real-data
identifiability signature), the conditional table (R1/R3) is the calibration
consequence, and K2/K3 are the resolution — exactly the boundary-railing-first option in
`RETOOL_HANDOFF.md` §5.

## Novelty axis (CP4, binding)
Stated up front in §1 and both cover letters: *first systematic real-data quantification
of D\* boundary-railing across cohorts + honest conditional-coverage characterization +
calibration resolution* — explicitly NOT the bare observation that NLLS fails. Deltas vs
(a) known D\* instability and (b) Casali 2026 are one sentence each in §1 and §4.

## Scope-driven consequence (beyond the ledger)
Because the redraft enforces "every number traces to a Gnomon run," the old
deep-ensemble vs input-perturbation comparison and the CRLB-efficiency table (Table 2) —
which Gnomon did not reproduce and which are not in the keep-set — were also removed from
the body, alongside the explicitly out-of-scope items. Sextant supplies the only
non-Gnomon numbers (cross-cohort railing generalization), clearly labelled as such.

## Verdict
PARTIAL (Gnomon CP3) is represented honestly: the diagnostic (K1) and mechanism (K2/K3)
reproduce; the two marginal coverage targets (T3a/T3b) diverge for the documented
SD-convention reason and are reframed conditionally — not hidden.
