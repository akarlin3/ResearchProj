# Procrustes CP0 — results

Every number below was printed by `python experiments/run_cp0.py 8`
(n=1000, SNR=50, α=0.10 → 90% intervals, 8 seeds = `GLOBAL_SEED + 1009·i`,
target = the well-identified tissue diffusion **D**, well-ID = bottom-2 D\*
terciles). Re-running from the clean seeds reproduces them. Ground truth is the
Lattice DRO (seed-generated, no data files). Brackets are 95% across-seed CIs.

Reference: Gauge's naive-transfer misspecification monitor sits at AUC ≈ 0.501
(a hidden channel for dispersion).

## Headline (the CP0 separation)

| family | expect | marginal | conditional GAP | well-ID-D\* GAP | bias ratio | diag AUC |
|---|---|---|---|---|---|---|
| **stretched-exp (β)** | BREAK | 0.909 [0.898,0.920] | **0.106 [0.092,0.120]** | **0.152 [0.131,0.173]** | 7.7× | **0.684 [0.673,0.694]** |
| log-normal disp. (cv) | WEAK | 0.909 [0.903,0.916] | 0.041 [0.032,0.050] | 0.054 [0.042,0.066] | 3.1× | 0.578 [0.574,0.582] |
| tri-exp 3rd pool (g) | NULL | 0.908 [0.901,0.915] | −0.000 [−0.010,0.010] | −0.009 [−0.023,0.005] | 0.8× | 0.627 [0.617,0.637] |

Reading:
- **(a) marginal holds** — all three families keep pooled coverage ≈ 0.909 under
  the departure-blind conformal radius (Lei 2018; *not* the contribution).
- **(b) conditional fails / (c) distinct from Gauge** — only when the departure
  aliases into the high-b tissue slope (stretched, heavy tail) does D's
  conditional coverage break, and it breaks **harder inside the well-identified
  D\* subset** (gap 0.152 > marginal gap 0.106) — exactly where Gauge's high-D\*
  wall does not operate and Gauge says "trust D".
- **boundary (null)** — the tri-exp *faster* third pool decays away from high-b,
  leaves D unbiased (ratio 0.8×), and produces **no** conditional failure
  (gap ≈ 0). The wedge is mechanism-specific, not generic misspecification.

## Conditional coverage by departure — stretched-exp (the break)

| β | conditional cov | well-ID-D\* cov |
|---|---|---|
| 1.0 (placebo) | 0.958 [0.954,0.963] | 0.947 [0.941,0.954] |
| 0.9 | 0.954 [0.947,0.960] | 0.938 [0.929,0.948] |
| 0.8 | 0.933 [0.922,0.943] | 0.904 [0.891,0.917] |
| 0.7 | 0.897 [0.882,0.911] | 0.852 [0.834,0.871] |
| 0.6 | 0.860 [0.843,0.876] | 0.803 [0.780,0.826] |
| 0.5 | 0.852 [0.835,0.870] | 0.795 [0.771,0.819] |

Monotone, and the placebo (β=1) row holds at/above nominal — the failure is
departure-driven, not a binning artefact (refute R3 clears).

## Tri-exp null — coverage flat across g

| g | conditional cov | well-ID-D\* cov |
|---|---|---|
| 0.0 | 0.905 | 0.873 |
| 0.1 | 0.913 | 0.885 |
| 0.2 | 0.911 | 0.882 |
| 0.3 | 0.907 | 0.877 |
| 0.4 | 0.905 | 0.882 |

Note: tri-exp's diagnostic AUC (0.627) > chance even though coverage does **not**
fail — the third pool is *detectably* non-bi-exp but *harmless* to D. The
diagnostic's job is to flag the families that do break coverage; for those
(stretched) detection and failure coincide.

## Refute conditions — status

| condition | meaning | status |
|---|---|---|
| **R1** | conditional gap ≈ 0 (no misspec-driven failure) | **cleared** (stretched gap 0.106, CI excludes 0) |
| **R2** | gap vanishes in well-ID-D\* subset (= Gauge) | **cleared** (well-ID gap 0.152 > marginal) |
| **R3** | placebo row itself broken | **cleared** (β=1 cond 0.958 ≥ nominal) |
| boundary | tri-exp null; bi-exp limit exact | **holds** (gap ≈ 0; continuity residual 0) |

The repositioned Procrustes wedge survives the gate. Open risk (POSITIONING.md):
the diagnostic is moderate (AUC ≈ 0.68) and family-dependent — strong for the
heavy-tail channel, near-hidden for pure dispersion.
