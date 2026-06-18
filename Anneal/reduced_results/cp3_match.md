# CP3 — quantitative matches to finite-N

**Timescale:** The shipped integrator advances 1 model-time unit per second of the drone (omega=0 rotating frame; DRIFT_DT=0.05 control step, sampleStride=0.1). integrator.mjs documents '1 unit := 1 s'. So reduced model-time equals seconds 1:1 — every reduced quantity below is in seconds with NO fudge factor.

## A=0.5 (post-homoclinic spiral-out)

Chimera FP (unstable spiral): r*=0.334, λ=σ±iω, σ=+0.01243/s, ω=0.3277 ⇒ linear T_breath=2π/ω=19.18s, per-cycle growth exp(σT)=1.269.

| quantity | reduced | measured | ratio |
|---|---|---|---|
| (i) breath period T_b | 25.0s (lin 19.2s) | 21–25s | 1.09 |
| (ii) spiral slope (pooled) | -0.67 | −0.55…−0.77 | 1.02 |
| (iii) capture time (median) | 43s | 139s | 0.31 |

### spiral slope, reduced vs measured, per N
| N | reduced | measured |
|---|---|---|
| 8 | -0.62 | -0.55 |
| 16 | -0.66 | -0.50 |
| 32 | -0.66 | -0.59 |
| 64 | -0.73 | -0.77 |

### breath period & capture, per N (reduced)
| N | T_b (s) | capture median (s) |
|---|---|---|
| 4 | 24.3 | 12 |
| 8 | 24.5 | 42 |
| 16 | 25.0 | 42 |
| 24 | 24.7 | 43 |
| 32 | 24.8 | 43 |
| 48 | 25.6 | 64 |
| 64 | 25.2 | 64 |

## A=0.2 (stable stationary chimera)

Reduced stable-chimera FP r*=0.6781; relaxation σ=-0.00539/s (per-cycle decay 0.803).

| quantity | reduced | measured |
|---|---|---|
| stable-FP r* vs never-absorber mean R_incoh | 0.6781 | 0.6590 (pooled); N=32:0.6758, N=64:0.6685 |
| breath-max drift ⟨ΔM⟩ sign | σ<0 (relaxes, ΔM<0) | −0.0011/cycle (ΔM<0) |

Finite-N mean R_incoh climbs toward r*=0.678 as N grows (N=32: 0.676) — the never-absorbers hover around the reduced stable fixed point, and the reduced relaxation (σ<0) explains the stationary, slightly-negative ⟨ΔM⟩ (opposite sign to the A=0.5 positive ratchet).
