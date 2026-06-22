# Reverb -- constructive counterexample: precision without coverage (SOLID)

Seed 20260622. Lattice (synthetic truth) + Caliper (estimator/conformal), read-only. ROI region_size=200, single-voxel SNR=40.0, n_cal=n_eval=2000, BCa bootstrap n_boot=2000. Params in Caliper space (D, f, D*).

**Scope:** synthetic possibility-and-mechanism proof -- the divergence CAN occur in IVIM and here is why; it does NOT quantify real-world miscalibration magnitude.

## Headline: f @ D*-lo (perfusion fraction, low-D* regime)
| arm | truth fit as bi-exp | wCV [BCa] | ICC | repeatable | cov(true) [BCa] | broken |
|---|---|---|---|---|---|---|
| matched control | biexp | 0.029 [0.026,0.031] | 0.996 | True | 0.796 [0.765,0.828] | False |
| **mismatch** | **stretched** | 0.032 [0.029,0.035] | 0.995 | True | **0.606 [0.568,0.642]** | **True** |

Precision is identical (wCV 0.029 vs 0.032; ICC ~1.00), yet model mismatch drops true-coverage by 0.190 (0.796 -> 0.606). Marginal f-coverage looks fine in both (0.895 / 0.840) -- the break is hidden until truth is known. **Repeatability is blind to the bias.**

## Verdict
- counterexample_found: **True**
- control_tracks (matched model not broken): True
- precision_blind (wCV ~equal while coverage diverges): True

## Per-regime detail (mismatch arm): repeatability vs conditional true-coverage
| param | regime | wCV | ICC | repeatable | cov(true) [BCa] | cov(mondrian) | counterex |
|---|---|---|---|---|---|---|---|
| D | D*-lo | 0.01 | +1.00 | YES | 0.83 [0.79,0.85] | 0.84 | - |
| D | D*-mid | 0.01 | +1.00 | YES | 0.92 [0.89,0.94] | 0.92 | - |
| D | D*-hi | 0.01 | +1.00 | YES | 0.93 [0.91,0.95] | 0.90 | - |
| f | D*-lo | 0.03 | +1.00 | YES | 0.61 [0.57,0.64] | 0.70 | **YES** |
| f | D*-mid | 0.03 | +1.00 | YES | 0.95 [0.94,0.97] | 0.94 | - |
| f | D*-hi | 0.04 | +1.00 | YES | 0.96 [0.94,0.97] | 0.89 | - |
| Dstar | D*-lo | 0.07 | +0.84 | YES | 1.00 [1.00,1.00] | 0.33 | - |
| Dstar | D*-mid | 0.12 | +0.40 | no | 0.99 [0.98,1.00] | 0.40 | - |
| Dstar | D*-hi | 0.24 | +0.14 | no | 0.46 [0.42,0.50] | 0.97 | - |

## Sensitivity surface: cov(true) of f @ D*-lo by truth-family x ROI size (SNR 40.0)
| truth family | nvox=50 | nvox=100 | nvox=200 | nvox=400 |
|---|---|---|---|---|
| biexp | 0.85 | 0.83 | 0.83 | 0.79 |
| dispersion_gamma | 0.71 | 0.66 | 0.67 | 0.65* |
| dispersion_lognormal | 0.74 | 0.69 | 0.71 | 0.68 |
| stretched | 0.68 | 0.62* | 0.64* | 0.61* |
| triexp | 0.90 | 0.86 | 0.87 | 0.84 |

`*` = repeatable-but-broken (counterexample) cell. The matched bi-exp control never breaks; dispersed-perfusion families break increasingly with ROI size -- precision rises, accuracy does not.