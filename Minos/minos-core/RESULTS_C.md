# RESULTS_C — Minos-Core v3: a deployment-validity monitor for a stale loss-calibration correction

Every number below was printed by `python experiments/run_c.py` in this session (seed `20240517`;
`CAL_N = 2_000_000`, `DEP_N = 1_000_000`, detection ROC `DETECT_N = 150_000` × 16 seeds, `m*` null
80 batches). Re-running from the clean seed reproduces them (GATE 4). No number here was hand-entered
from anywhere but that run's stdout.

**Config.** Calibration cell = the v2 default misspecification: report centres `mu ~ N(t2 − rho*s,
0.5^2)`, posterior-centric (truth `theta = mu + s*u`, `u` a standardised skew-normal of shape
`kappa=3`), single active threshold `t2 = 2.0`, `s = 0.5`, under:over cost `lambda = 3`,
proximity `rho = 0.5`. Deployment shift split (gain `beta*s = 2.5`): `delta_obs` biases the reported
point `mu` **down** (observable); `delta_hid` biases the truth `theta` **up** (hidden). Detection
regret tolerance `tol = 0.02`; `m*` false-alarm target `alpha = 0.05`.

---

## The narrative in one line

The v2 gap motivates a **loss-calibration correction** (a cited baseline); that correction is fit on
a calibration set and **goes stale under deployment shift**; a **label-free validity monitor** flags
the *observable* fraction of that staleness; a **gated-recovery** policy claws back the lost utility;
and the *hidden* fraction is quantified and shown to be **undetectable in principle** — which is the
honest case for periodic labeled repeatability spot-checks.

---

## GATE 1 — loss-calibration baseline (cited) + stale-correction regret

The baseline is the scalar post-hoc loss-calibration scale on the labeled calibration set,
**identical to v2's decision-calibrated scale** (`POSITIONING.md` cites Lacoste-Julien 2011 / Vadera
2021 / Zhao 2021):

| quantity | value |
|---|---|
| loss-calibration baseline `tau_hat_cal` | **1.0480** |
| consistency: `std(tau_hat)` across 5 seeds at `N = 1e5 / 5e5 / 2e6` | `0.0234 / 0.0026 / 0.0046` |
| `tau*_oracle@dep` (`N=2e6`, independent draw) vs `mean(tau_hat)@2e6` | `1.0472` vs `1.0527` |
| v2 gap on the cell: `tau_stat`, `tau*`, `G` | `0.9639`, `1.0480`, `+0.0841` |
| `R(0,0)` (zero-shift regret) | `0.000000` |

The estimator is **consistent**: its across-seed variance shrinks ~5× from `N=1e5` to `2e6` and the
large-`N` fit agrees with the deployment oracle to `0.006` (well inside the `0.05`-grid resolution).
The motivating v2 gap reproduces on this cell (`G = +0.084`, `tau* > tau_stat`, matching the v2
headline sign and magnitude). The stale regret is zero at zero shift and grows monotonically with
**each** shift knob:

```
R vs delta_obs:  0.00:0.0000  0.03:0.0039  0.06:0.0147  0.09:0.0320  0.12:0.0561  0.15:0.0845  0.18:0.1164  0.21:0.1500  0.24:0.1838
R vs delta_hid:  0.00:0.0000  0.03:0.0034  0.06:0.0147  0.09:0.0333  0.12:0.0601  0.15:0.0944  0.18:0.1366  0.21:0.1865  0.24:0.2443
```

The two knobs induce **comparable** regret (by construction they impose the same truth↔report
discrepancy) — the only difference is which channel moves. That difference is everything for a
label-free monitor.

## GATE 2 — the validity monitor `M` (label-free): what it catches and what it cannot

| quantity | value |
|---|---|
| no-label code path: `M(scrambled theta) == M(real theta)` | **True** (`M = 0.234732`) |
| under `delta_obs`: `corr(M, R)` | **+0.965** |
| under `delta_obs`: detection AUC for `{R > tol}` | **1.000** (positives 96/144) |
| under `delta_hid`: detection AUC for `{R > tol}` | **0.500** — *at chance* (positives 96/144) |
| threshold `m*` (`alpha = 0.05`) | `0.00876` |
| zero-shift false-alarm rate at `m*` | `0.062` |

The monitor is **structurally label-free** — scrambling `theta` leaves `M` bit-identical, because `M`
is a function of the reported points `{mu}` only. Under the **observable** shift it tracks regret
almost perfectly (`corr = +0.97`) and separates stale from fresh batches with **AUC = 1.00**. Under
the **hidden** shift — *the same regret*, induced by moving the truth instead of the report — the
monitor is at **chance (AUC = 0.50)**, because `delta_hid` does not move `{mu}` and `M` cannot see
it. **This is the honest, documented limitation, not a failure.** It is exactly the regime where
label-free monitoring is blind and labeled repeatability spot-checks are required. The threshold `m*`
holds the zero-shift false-alarm rate at `0.062`, on target for `alpha = 0.05`.

## GATE 3 — gated recovery + the regret ladder

The gated-recovery policy escalates the decision-fragile (near-threshold) voxels **only when `M`
fires**. The regret ladder `{stale, gated, oracle}` across the observable sweep:

```
delta_obs   R_stale    R_gated    fires
  0.00      -0.0000    -0.0000    False
  0.03      +0.0039    +0.0226    True
  0.06      +0.0147    +0.0090    True
  0.09      +0.0320    +0.0015    True
  0.12      +0.0561    +0.0008    True
  0.15      +0.0845    +0.0074    True
  0.18      +0.1164    +0.0206    True
  0.21      +0.1500    +0.0404    True
  0.24      +0.1838    +0.0665    True
```

The monitor's operating point (first `delta_obs` that fires) is `0.03`. Where the shift is
substantial the gate recovers most of the stranded utility (e.g. at `delta_obs = 0.12` it cuts the
regret from `0.056` to `0.001`; at the strongest shift `0.24` from `0.184` to `0.067`, ~64%
recovered). **Honest caveat:** at the *smallest* detected shift (`delta_obs = 0.03`) the conservative
escalation slightly over-corrects (`R_gated = 0.023 > R_stale = 0.004`) — the gate trades mild
over-conservatism at a marginal detection for large recovery everywhere the shift actually bites.

Crucially the gate does no harm where it should be silent:

| scenario | `R_stale` | `R_gated` | gate fires |
|---|---|---|---|
| zero shift | `-0.00002` | `-0.00002` | False |
| `delta_hid = 0.06` | `+0.01467` | `+0.01467` | False |
| `delta_hid = 0.24` | `+0.24432` | `+0.24432` | False |

At zero shift the gate is silent (no false-alarm harm). Under **every** hidden shift the gate is
silent too: `R_gated == R_stale` exactly, for all `delta_hid`. The gate **cannot** recover the
hidden-driven regret — and it does **not** spuriously appear to, which is the leakage check: if gated
beat stale under `delta_hid`, the monitor would be leaking labels. It does not.

## Figures (vector PDF, `figures/`)

- **`fig_c_a_motivation_gap.pdf`** — (a) the v2 decision-calibration gap `G` vs skew `kappa` at
  several cost asymmetries: the misspecification that motivates a loss-calibration correction.
- **`fig_c_b_regret_vs_shift.pdf`** — (b) stale-correction regret `R` vs shift, observable and hidden
  knobs overlaid (comparable regret, two channels).
- **`fig_c_c_monitor_honesty.pdf`** — (c) **the honesty figure:** `M` vs `R` scatter and the detection
  ROC, observable (AUC 1.00) vs hidden (AUC 0.50) overlaid.
- **`fig_c_d_regret_ladder.pdf`** — (d) the regret ladder `{stale, gated, oracle}` across `delta_obs`
  with the monitor's operating point marked.

The v1 (`fig_a..d`) and v2 (`fig_gap_*`) figures are retained unchanged.

## What v3 demonstrates

The v2 decision–calibration gap motivates a **loss-calibration correction** — but that correction is
an **existing method** (loss-calibrated Bayesian inference / decision calibration), so v3 does not
re-invent it; it treats it as a **cited baseline** and asks the question that line leaves open: *what
happens to the correction when deployment drifts away from the calibration set?* It **goes stale**,
stranding utility (`R` grows with shift). The contribution is a **label-free validity monitor** that
predicts that staleness from unlabeled deployment data: it tracks regret (`corr = +0.97`) and detects
it (`AUC = 1.00`) for the fraction of shift that moves **observable** reported summaries, and a
**gated-recovery** policy turns that detection into recovered utility. The monitor is **not** a
calibration method (it proposes no new scale) and **not** a generic OOD detector (it is utility-/
regret-targeted, blind by design to threshold-irrelevant input change) — see `POSITIONING.md`.
Finally, v3 is explicit about the boundary of label-free monitoring: a shift that changes the
truth↔report relationship **without moving any observable summary** induces the *same* regret yet is
**undetectable in principle** (`AUC = 0.50`). That hidden fraction is reported, not hidden, and is the
quantitative case for periodic **labeled repeatability spot-checks** alongside the monitor.
