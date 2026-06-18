# CP4 — run-by-run lifetime prediction

Reduced **3-D** system integrated per run from (ρ₁=R_sync0, ρ₂=R_incoh0, ψ=Δφ₀). Capture = first sustained min(ρ₁,ρ₂)>θ (no-recovery). Correlated with measured t_abs.

**Pooled (n=1400):** Spearman(reduced, t_abs)=**0.445**, |Δφ₀| static=-0.429, partial(reduced | |Δφ₀|)=0.200.
CV R²: reduced=0.146, |Δφ₀|+quad=0.185, both=0.206  (collective-IC benchmark R²≈0.23).

| N | n | ρ(reduced) | ρ(\|Δφ₀\|) | partial(red\|Δφ₀) | R²(red) | R²(\|Δφ₀\|q) | R²(both) | cap med (s) | rho1 dev med |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 200 | +0.141 | +0.032 | +0.210 | 0.022 | -0.037 | -0.018 | 12 | 0.003 |
| 8 | 200 | +0.498 | -0.445 | +0.253 | 0.202 | 0.162 | 0.203 | 42 | 0.003 |
| 16 | 200 | +0.461 | -0.428 | +0.219 | 0.152 | 0.254 | 0.262 | 42 | 0.003 |
| 24 | 200 | +0.544 | -0.507 | +0.259 | 0.227 | 0.234 | 0.254 | 43 | 0.003 |
| 32 | 200 | +0.595 | -0.582 | +0.238 | 0.315 | 0.348 | 0.372 | 43 | 0.003 |
| 48 | 200 | +0.578 | -0.594 | +0.229 | 0.408 | 0.339 | 0.419 | 64 | 0.003 |
| 64 | 200 | +0.592 | -0.538 | +0.315 | 0.226 | 0.205 | 0.214 | 64 | 0.003 |

## Transverse stability of ρ₁=1
Per-run max deviation of ρ₁ from 1 stays small (median column above); the ρ₁=1 manifold's local transverse rate is attracting for a majority of the trajectory at every N (transverse_attracting_frac mean: N=4:0.45, N=8:0.52, N=16:0.58, N=24:0.61, N=32:0.63, N=48:0.64, N=64:0.64).

## Where the N-dependence enters
The reduced **dynamics** (A,β) are N-independent; the only N-channel is the **seed-ensemble IC** distribution. Mean R_incoh0 and mean |Δφ₀| of the canonical seed shift with N (see cp4_rho_vs_N right panel), and the reduced capture median tracks that shift — it does NOT reproduce the flat measured τ_abs(N) plateau. So the finite-N τ-plateau and the growing spiral-rate-with-N both live in the **IC ensemble + finite-N fluctuations**, not in the reduced collective dynamics.
