# CP2 — initial manifold distance D₀ vs collapse lifetime

Per (A, N): Spearman rank correlation between the seed's incoherent-population D₀ and ln(lifetime), on the **uncensored** runs only (censored runs are kept in the stratified KM but excluded from the rank correlation). ρ<0 ⇒ larger initial manifold distance predicts shorter life (the manifold-escape hypothesis).

| A   | N   | runs | uncensored | censored | Spearman ρ | p-value  | D₀ median |
| --- | --- | ---- | ---------- | -------- | ---------- | -------- | --------- |
| 0.5 | 4   | 200  | 200        | 0        | -0.081     | 2.51e-01 | 0.630     |
| 0.5 | 8   | 200  | 200        | 0        | -0.041     | 5.61e-01 | 0.354     |
| 0.5 | 16  | 200  | 200        | 0        | -0.039     | 5.82e-01 | 0.169     |
| 0.5 | 24  | 200  | 200        | 0        | +0.072     | 3.09e-01 | 0.123     |
| 0.5 | 32  | 200  | 200        | 0        | -0.049     | 4.88e-01 | 0.087     |
| 0.5 | 48  | 200  | 200        | 0        | -0.072     | 3.09e-01 | 0.060     |
| 0.5 | 64  | 200  | 200        | 0        | -0.024     | 7.35e-01 | 0.040     |
| 0.2 | 4   | 100  | 99         | 1        | +0.069     | 4.95e-01 | 0.733     |
| 0.2 | 8   | 100  | 96         | 4        | -0.148     | 1.49e-01 | 0.323     |
| 0.2 | 16  | 100  | 100        | 0        | -0.007     | 9.45e-01 | 0.153     |
| 0.2 | 24  | 100  | 100        | 0        | -0.138     | 1.72e-01 | 0.119     |
| 0.2 | 32  | 100  | 100        | 0        | -0.015     | 8.80e-01 | 0.099     |
| 0.2 | 48  | 100  | 100        | 0        | -0.008     | 9.33e-01 | 0.057     |
| 0.2 | 64  | 100  | 100        | 0        | -0.010     | 9.20e-01 | 0.040     |
