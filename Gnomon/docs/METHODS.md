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
  (record id, file SHA-256, retrieval date). Never redistributed in-tree.
- **High-SNR ROI selection** (the 1618-voxel set): `[CP2]` document the exact b0-SNR
  threshold and ROI mask used; report railing-rate sensitivity to this choice.

## 2. Acquisition (b-value schemes, explicit)

- **clinical-sparse (8 b):** `(0, 10, 20, 40, 80, 200, 400, 800)` s/mm² — Gnomon's
  documented clean-room choice (Fashion's prose did not list its 8-b values).
- **dense (16 b):** `(0,10,20,30,50,75,100,150,200,300,400,500,600,700,800,1000)`
  s/mm² — quoted verbatim from Fashion (REVIEWER_RESPONSE.md:126-127).

## 3. Forward model `[CP2]`

IVIM bi-exponential `S(b)/S0 = (1-f)·exp(-b·D) + f·exp(-b·(D+Dstar))` (Le Bihan
1988), reimplemented from physics. Analytic Jacobian wrt (S0, D, f, Dstar) used for
NLLS and the asymptotic covariance. Self-consistency gate: clean-signal round-trip
recovers truth (residual reported).

## 4. Estimators & fitting (full detail)

### 4.1 NLLS + boundary-railing `[CP2]`
- Box-constrained four-parameter fit (`scipy least_squares`, trust-region
  reflective). Bounds, initialization, and `max_nfev`: documented here when set.
- **Railing definition:** parameter railed iff `|x̂ − bound|/(upper−lower) < rail_tol`;
  `rail_tol = 1e-3` primary, `1e-2` sensitivity. D\* railing **rate** is target T1.

### 4.2 Laplace posterior `[CP2]`
- Gaussian at the MAP with curvature (CRLB-style) covariance; symmetric SD interval.

### 4.3 MCMC posterior `[CP2]`
- Clean-room per-voxel sampler over the Rician/Gaussian likelihood + manifest priors.
- Proposal, chain length, burn-in, thinning, seed, convergence diagnostic:
  documented here when set. From one chain: **SD interval** (T3b) and **2.5/97.5
  quantile interval** (T3c).

### 4.4 MAF amortized posterior (NPE) `[CP2]`
- Architecture/transforms, simulation budget, prior, optimizer, epochs /
  early-stopping, seed: documented here when set. Optional `[flow]` extra (torch).

## 5. CRLB assumption (the reviewer-flagged item) `[CP2]`

The Laplace/asymptotic covariance is a **Gaussian** approximation to the posterior.
It is **weakest exactly where D\* is most skewed** (low SNR), which is where the
central under-coverage claim sits. Stated as a limitation; the direction of the bias
(conservative — it *understates* how overconfident the Gaussian interval is) is noted
so no claim rests on the approximation being tight.

## 6. Calibration ruler (re-derived independently) `[CP2]`

Coverage, ECE (mean |empirical − nominal| over a level grid), sharpness (mean
interval width), pinball loss, interval score — from published definitions, numpy
only, no Caliper import. Plus a D\*-tercile conditional-coverage probe.

## 7. Uncertainty on the numbers `[CP2]`

Bootstrap CIs (percentile/BCa, seeded from `manifest.BOOTSTRAP`) on every
load-bearing number; directional gaps (T4) pass only if the CI excludes 0.

## 8. Reproduction & verdict `[CP3]`

`gnomon/reproduce.py` runs the rebuild, compares to the frozen manifest, and writes
`results/reproduction.json` + the verdict (REPRODUCES / DOES NOT REPRODUCE). One
command: `bash reproduce.sh`.

## 9. Claims scope

Claims are held to what the rebuild supports. Where Fashion was read as overextending
(e.g. single-subject in-vivo generalization), Gnomon scopes the corresponding
statement to its evidence and records the limitation here.
