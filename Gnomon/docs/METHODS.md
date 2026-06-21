# Methods (the complete write-up Fashion was flagged for lacking)

> This is the methods section Gnomon exists to produce — **complete from line one**,
> covering every item Fashion's review flagged. Sections marked `[CP2]`/`[CP3]` are
> filled by the rebuild run; the *specifications* are fixed now so nothing is
> reverse-engineered after seeing results. Numbers, when added, trace to
> `results/reproduction.json`.

## 1. Datasets (every dataset named with a stable ID)

### 1.1 Synthetic substrate — Lattice DRO
- Source: `lattice.make_cohort` (sibling, MIT, synthetic-only), imported read-only.
- Ground truth `(D, Dstar, f)` drawn from Lattice's documented physiological priors;
  seed = `manifest.MASTER_SEED` (20260620), fully reproducible from the seed.
- **Headline 9-cell design:** 3 ground-truth conditions × SNR ∈ {10, 20, 40} × 200
  noise realizations (Rician), matching Fashion's headline set (README.md:58-59).

### 1.2 Open real data — OSIPI abdomen
- **OSIPI TF2.4**, Zenodo record **14605039** (`OSIPI_TF24_data_phantoms.zip`),
  CC-licensed; the in-vivo **abdomen** acquisition (`Data/abdomen.*`).
- Fetched on demand into a gitignored `download/` with a provenance manifest
  (record id, SHA-256 `2a53054d…b3e`). Never redistributed in-tree. `gnomon/osipi.py`.
- **Acquisition (read from the archive):** 144×144×21 volume, **104 measurements**
  over 12 unique b-values {0,10,20,30,40,50,75,100,150,250,400,600} s/mm² (Philips 3T).
- **ROI selection (documented):** the archive's `Data/mask_abdomen_homogeneous.nii.gz`
  homogeneous-tissue mask = **1932 voxels**; each voxel b0-normalized (÷ mean b=0).
  Fashion's "1618 high-SNR ROI voxels" is this ROI under an *unstated* SNR cut — a
  completeness gap. Gnomon reports the D\* railing rate on the **full homogeneous ROI**
  and on a b0-SNR>25 subset (b0-SNR = mean/std over the 15 b=0 repeats), so the
  load-bearing number (the rate) carries its selection sensitivity rather than being
  pinned to an unreproducible voxel count.

## 2. Acquisition (b-value schemes, explicit)

- **clinical-sparse (8 b):** `(0, 10, 20, 40, 80, 200, 400, 800)` s/mm² — Gnomon's
  documented clean-room choice for the *synthetic* design (Fashion's prose did not
  list its 8-b values).
- **dense (16 b):** `(0,10,20,30,50,75,100,150,200,300,400,500,600,700,800,1000)`
  s/mm² — quoted verbatim from Fashion (REVIEWER_RESPONSE.md:126-127).
- **real OSIPI abdomen:** the native 104-point / 12-unique-b scheme above (used as-is
  for the T1 railing fit).

## 3. Forward model

IVIM bi-exponential `S(b)/S0 = (1-f)·exp(-b·D) + f·exp(-b·(D+Dstar))` (Le Bihan
1988), reimplemented from physics (`gnomon/forward.py`). Estimators fit in **scaled**
params `(S0, D3, f, Ds3)` with `D = D3·1e-3`, `Dstar = Ds3·1e-3` so every fit variable
is O(1). Analytic Jacobian (verified vs finite differences). Self-consistency gates:
clean-signal NLLS round-trip recovers truth; continuity `f=0 ⇒ mono-exponential`.

## 4. Estimators & fitting (full detail)

### 4.1 NLLS + boundary-railing
- Box-constrained four-parameter fit (`scipy.optimize.least_squares`, trust-region
  reflective, analytic Jacobian, `max_nfev=400`). **Box** (scaled): `S0∈[0.5,1.5]`,
  `D3∈[0.2,3.0]`, `f∈[0,0.5]`, `Ds3∈[3,150]` (D/Dstar bounds follow Fashion's stated
  NPE prior range). **Init** `(1.0, 1.0, 0.1, 20.0)`.
- Covariance: `σ²·pinv(JᵀJ)` (known σ=1/SNR for the synthetic Laplace; residual-based
  for the real-data fit), each SD capped at the box span so a railed/unidentified D\*
  yields a finite-but-pathological interval.
- **Railing:** parameter railed iff `|x̂ − bound|/(upper−lower) < rail_tol`;
  `rail_tol = 1e-3` primary, `1e-2` sensitivity. D\* railing **rate** is target T1.

### 4.2 Laplace posterior
- Gaussian at the MAP (NLLS fit) with the CRLB covariance (known σ); symmetric SD
  interval → T3a.

### 4.3 MCMC posterior
- Clean-room per-voxel random-walk Metropolis, **vectorized across all voxels**, over a
  Gaussian likelihood (known σ=1/SNR) with a uniform-box prior (the NLLS box).
- Proposal: isotropic Gaussian RW, per-voxel per-param step = `0.6 ×` the CRLB SD;
  init at the NLLS MAP. **burn 1500, keep 2000, thin 2**, seed `MASTER_SEED+3`;
  acceptance ~0.25 (reported per run). From one chain: **SD interval** (T3b) and
  **2.5/97.5 quantile interval** (T3c).

### 4.4 MAF amortized posterior (NPE)
- Conditional autoregressive flow (5 affine AR layers, per-dim MLP conditioners,
  hidden 64, dim order reversed between layers; standard-normal base) over scaled
  `(D3,f,Ds3)` standardized by the prior; context = standardized signal. **Sim budget
  80 000** `(θ,x)` pairs, θ~uniform-box, SNR~U(10,40); Adam lr 1e-3, **40 epochs**,
  batch 512, seed `MASTER_SEED+5`. Inference: 1500 posterior draws/voxel → quantiles.
  Optional `[flow]` extra (torch). No sbi/nflows — written from scratch.

## 5. CRLB assumption + railed-voxel SD convention (the reviewer-flagged item)

Two distinct, load-bearing modelling choices live here. Both are stated; the second is
the one that reconstructs Fashion's marginal headline, so it is documented in full.

**(a) Gaussian/CRLB approximation.** The Laplace/asymptotic covariance
`σ²·pinv(JᵀJ)` is a Gaussian approximation to the posterior, **weakest where D\* is most
skewed** (low SNR / high D\*) — exactly where the under-coverage sits. The MCMC results
do **not** rely on it (the sampler targets the exact posterior under the box prior); it
is used only for the Laplace SD interval and as the MCMC proposal scale.

**(b) Railed-voxel SD convention — the headline-driving choice.** A voxel whose D\* is
railed/unidentified carries **no local Fisher information** about D\* (the Jacobian
column collapses; `JᵀJ` is singular). Two defensible conventions then exist, and they
give *opposite* coverage:

- **Honest CRLB (Gnomon default):** the pseudo-inverse assigns such a voxel a **wide**
  D\* interval (SD capped at the box span). This *over-covers* railed voxels — the
  statistically correct admission that D\* is undetermined there.
- **Overconfident floor (Fashion-implied):** a railed voxel is instead given a small
  **floored** D\* SD (`RAILED_SD_FLOOR = 0.003 mm²/s`), so its interval is narrow and
  *under-covers*. This is "overconfident by design."

The convention is decisive for the marginal number: flooring drops Laplace D\* coverage
from pooled 0.80 → **0.68** (high-D\* tercile 0.63 → 0.41), and the floored pooled value
(0.68) reconstructs Fashion's MCMC-SD claim (0.67). See `gnomon/reframe.py` /
`handoff/conditional_coverage.json` (both conventions, per-tercile, bootstrap CIs) and
`scripts/divergence_diagnostic.py`.

**Recommended:** report under the **honest CRLB** and state the convention explicitly.
The honest convention does not erase the finding — a real, reproducible under-coverage
survives it in the high-D\* tercile (Laplace 0.63 [0.60, 0.67]); it only removes the
manufactured marginal severity. The floored convention should appear, if at all, only
as a labelled illustration of how the original 0.30/0.67 arose.

## 6. Calibration ruler (re-derived independently)

Coverage, ECE (mean |empirical − nominal| over a level grid), sharpness (mean
interval width), pinball loss, interval score — from published definitions
(`gnomon/metrics.py`), numpy only, no Caliper import. Plus a D\*-tercile
conditional-coverage probe (the basis of the reframe, §the conditional table).

## 7. Uncertainty on the numbers

Bootstrap CIs (percentile/BCa, `gnomon/bootstrap.py`, seeded from
`manifest.BOOTSTRAP`) on every load-bearing number; directional gaps (T4) pass only if
the CI excludes 0; per-tercile conditional cells carry their own CIs.

## 8. Reproduction & verdict

`gnomon/reproduce.py` runs the rebuild, compares to the frozen manifest, and writes
`results/reproduction.json` + the verdict. `gnomon/reframe.py` writes the conditional-
coverage table (`handoff/conditional_coverage.json`). One command regenerates the whole
hand-off: `bash Gnomon/reproduce.sh` (gates + verdict) and
`python -m gnomon.reframe` / `python scripts/build_handoff.py` (reframe + hand-off).

## 9. Claims scope

Claims are held to exactly what the rebuild supports; the disposition of each Fashion
claim (KEEP / REFRAME / DROP) is in [`../handoff/CLAIMS_LEDGER.md`](../handoff/CLAIMS_LEDGER.md).
Where Fashion read as overextending — the *marginal* 0.30/0.67 severity, and
single-subject in-vivo generalization — Gnomon scopes the statement to its evidence
(conditional per-D\* coverage; ROI-level railing rate with selection sensitivity) and
records the limitation rather than restating the dramatic figure.

## 10. Completeness checklist (the items Fashion was flagged for) — all present

| Huang-flagged item | Where addressed |
|---|---|
| Dataset IDs | §1.1 (Lattice seed/prior), §1.2 (OSIPI Zenodo 14605039, sha256) |
| Acquisition / b-schemes | §2 (clinical-sparse 8-b, dense 16-b, real 104-pt) |
| Full training (NPE/MAF) | §4.4 (architecture, 80k sims, prior, epochs, seed) |
| Full fitting (NLLS, MCMC) | §4.1 (box/init/solver/railing), §4.3 (sampler spec) |
| CRLB assumption | §5(a) |
| **Railed-voxel SD convention** | §5(b) — both conventions, recommendation justified |
| Calibration ruler definitions | §6 |
| Uncertainty / CIs | §7 |
| Claims scope (no overclaim) | §9 + claims ledger |
