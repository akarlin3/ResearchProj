# RESULTS_B — Minos-Core v2: the decision–calibration gap

Every number below was printed by `python experiments/run_b.py` in this session
(seed `20240517`; headline cells `RUN_N = 4_000_000`, sweep grids `SWEEP_N = 2_000_000`).
Re-running from the clean seed reproduces them. No number here was hand-entered from
anywhere but that run's stdout.

**Config.** Gap experiment: report centres `mu ~ N(t2 − rho*s, 0.5^2)`, `posterior_centric`
(truth `theta = mu + s*u`, `u` a standardised skew-normal of shape `kappa`); single active
threshold `t2 = 2.0` (spare threshold pushed to `t1 = −18.0`); `s = 0.5`; under:over cost
`lambda = k_under` with `k_over = 1`. `tau_stat` reference level `L_ref = 0.90`.

---

## The headline — statistical and decision calibration diverge

Under a misspecified (skewed) posterior the two calibration scales of the reported error bar
**point in opposite directions**, by a quantifiable gap. At the default misspecified cell
(`kappa=3, lambda=3, rho=0.5`):

| scale | value | what it says |
|---|---|---|
| `tau_stat` (nominal 90% coverage) | **0.9635** | statistical calibration: *shrink* the bar |
| `tau*` (`argmax_tau E[U]`) | **1.0431** | decision calibration: *widen* the bar |
| **gap `G = tau* − tau_stat`** | **+0.0796** | a deliberately under-confident bar beats the calibrated one |
| ratio `tau*/tau_stat` | **1.0826** | |

The statistically-calibrated error bar (which over-covers at the 90% level under right-skew, so
honest coverage-matching shrinks it to `0.96`) is **the wrong width for the decision**: the
asymmetric cost wants it *wider* (`1.04`). `G` is the distance between the two — the quantity none
of the neighbours (DCA/net-benefit, ISPOR VoI, ARCliDS) computes (see `POSITIONING.md`).

The well-specified case is the **degenerate baseline inside the result**: at `kappa=0` both scales
are `1.0` and `G≈0` (below).

## GATE 1 — gap estimator (definitions, well-specified limit, misspecified sign)

| quantity | value |
|---|---|
| well-specified `kappa=0`: `tau_stat`, `tau*`, `G` | `0.9998`, `0.9959`, `−0.0039` |
| v1 mixture reproduction: VoC argmin `tau*` | `0.990` |
| v1 VoC at `n=1e6` (`VoC(0.5)`, `VoC(2.0)`) | `+0.001169`, `+0.004361`  (= v1 RESULTS.md exactly) |
| default misspecified `kappa=3,lambda=3,rho=0.5`: `tau_stat`, `tau*`, `G`, ratio | `0.9635`, `1.0431`, `+0.0796`, `1.0826` |
| independence: `tau_stat` unchanged when utility ×3 | `True` |
| `tau_stat(L)` profile @ default misspec | `0.5:1.002  0.68:0.980  0.8:0.966  0.9:0.964  0.95:0.983  0.99:1.081` |

At `kappa=0` the report equals the true posterior, so `tau_stat=1` (coverage nominal at every
level) and `tau*=1`, giving `G≈0` — and the v1 symmetric baseline reproduces exactly (`VoC(0.5)`,
`VoC(2.0)` match v1 RESULTS.md to six decimals). Under misspecification `G>0` with the documented
sign `tau* > tau_stat`. The two scales come from **disjoint code paths**: `tau_stat` is invariant
to a 3× change in the utility (it never reads the cost), and `tau*` takes no coverage level.

**Robustness of the gap to the coverage criterion.** The `tau_stat(L)` profile shows `tau_stat<1`
("shrink") for every central level `L≤0.95`, so `G>0` is robust across the standard coverage range;
only at the extreme `L=0.99` tail (where the symmetric interval finally *sees* the heavy upper tail
the decision pays for) does `tau_stat` rise to `1.081 ≈ tau*`. The disagreement between coverage and
decision is precisely a *bulk-vs-tail* disagreement.

## GATE 2 — the gap map

`G(kappa, lambda)` at `rho=0.5` (rows `kappa`, cols `lambda`):

```
        l=1.0   l=2.0   l=3.0   l=4.0
 k=0.0  +0.001  -0.009  -0.004  +0.004
 k=0.5  -0.000  -0.012  -0.003  +0.004
 k=1.0  +0.001  -0.010  +0.005  +0.015
 k=2.0  +0.015  +0.016  +0.044  +0.057
 k=3.0  +0.036  +0.050  +0.084  +0.100
 k=4.0  +0.055  +0.078  +0.115  +0.133
```

`G(kappa, rho)` at `lambda=3` (rows `kappa`, cols `rho`):

```
        r=0.0   r=0.5   r=1.0   r=1.5   r=2.0
 k=0.0  +0.017  -0.004  -0.004  -0.008  +0.004
 k=1.0  +0.022  +0.005  +0.006  +0.004  +0.006
 k=2.0  +0.054  +0.044  +0.042  +0.043  +0.034
 k=3.0  +0.091  +0.084  +0.082  +0.083  +0.070
 k=4.0  +0.120  +0.115  +0.113  +0.114  +0.099
```

| sanity check | result |
|---|---|
| `kappa=0` slice ≈ 0 (max \|G\|): `lambda`-row / `rho`-row | `0.009` / `0.017` |
| monotone increasing in `kappa` (`k≥1`, `lambda=3`) | `[0.005, 0.044, 0.084, 0.115]` |
| monotone increasing in `lambda` (`kappa=3`) | `[0.036, 0.050, 0.084, 0.100]` |
| decreasing in `rho`, largest at `rho=0` (`kappa=3`) | `[0.091, 0.084, 0.082, 0.083, 0.070]`, slope `−0.0084` |
| corner table (`kappa,lambda`) | `(0,1)=+0.001 (0,4)=+0.004 (4,1)=+0.055 (4,4)=+0.133` |
| corner table (`kappa,rho`) | `(0,0)=+0.017 (4,0)=+0.120 (4,2)=+0.099` |

The `kappa=0` row is flat-zero across both `lambda` and `rho` (the well-specified baseline). `G`
grows monotonically with the skew `kappa` and with the cost asymmetry `lambda`; at `lambda=1` the
decision-side inflation vanishes (the escalate boundary is `tau`-independent, so `tau*=1`) leaving
the pure coverage offset `G = 1 − tau_stat`. `G` is largest at the threshold (`rho=0`) and falls as
the population moves into the region (`rho` is held in the resolved regime `[0,2]`; beyond that the
EU curve is degenerately flat and `tau*` is unidentified — excluded by design).

## GATE 3 — break-even shift for the trust-gate

The v1 footnote (`VoTG` dips slightly negative for small shifts, then crosses zero and grows) is now
a first-class number on the v1 baseline config:

| quantity | value |
|---|---|
| `VoTG(delta=0)` | `−0.006739` |
| break-even shift `delta_be` | `0.9682` |
| `VoTG(delta_be − 0.1)` | `−0.008858`  (below break-even: gate is a net cost) |
| `VoTG(delta_be + 0.2)` | `+0.030548`  (above break-even: gate pays off) |

The trust-gate is a **net cost** for shifts below `delta_be ≈ 0.97` (its fixed false-positive cost
outweighs the thin rescue) and **pays off** above it. This is the "when is the gate worth turning
on" answer, quantified.

## Figures (vector PDF, `figures/`)

- **`fig_gap_a_voc.pdf`** — (a) `VoC(tau)` well-specified vs misspecified, with `tau_stat` (coverage,
  shrink) and `tau*` (decision, widen) marked and the gap `G` annotated — the misspecified curve dips
  *below zero* at `tau*`, i.e. the under-confident bar beats the moment-matched one.
- **`fig_gap_b_G_vs_kappa.pdf`** — (b) **headline:** gap `G` vs skew `kappa` at `lambda ∈ {1,2,3,4}`;
  near zero at `kappa=0`, fanning upward, steeper with cost asymmetry.
- **`fig_gap_c_surface.pdf`** — (c) `G` surface over `(kappa, rho)`: largest near the threshold and at
  strong skew, with the `G=0` contour bounding the well-specified corner.
- **`fig_gap_d_votg_breakeven.pdf`** — (d) `VoTG(delta)` with the break-even shift `delta_be=0.968`
  marked; the shaded region is where the gate is a net cost.

The v1 figures (`fig_a..d`) are retained unchanged.

## What v2 demonstrates

Statistical calibration and decision calibration of a reported error bar are **different objectives**
that **diverge under realistic (skewed / misspecified) posteriors**, by a quantifiable
**decision–calibration gap** `G = tau* − tau_stat`. Calibrating the bar for nominal coverage can be
the *wrong* thing for the decision: under right-skew the coverage-honest scale *shrinks* the bar
(`tau_stat≈0.96`) while the asymmetric decision wants it *wider* (`tau*≈1.04`), and a deliberately
under-confident bar earns more utility than the statistically-calibrated one. `G` vanishes in the
well-specified / symmetric corner (recovering v1 at `tau=1`) and grows with posterior skew, with
cost asymmetry, and with threshold proximity. None of the neighbours computes this distance — they
each price a single object, not the gap between statistical and decision calibration (see
`POSITIONING.md`). The trust-gate result is sharpened in passing: its break-even shift
`delta_be ≈ 0.97` says exactly when distrusting an over-confident posterior is worth its
false-positive cost.
