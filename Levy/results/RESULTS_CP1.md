# RESULTS -- CP1: joint CTRW (alpha, beta) structural identifiability / degeneracy

All numbers derived (joint 4-parameter Fisher/CRLB + parametric bootstrap), fully
synthetic, seeded. CRLB = identifiability bound, scoped to its regime (single
diffusion time, finite b-values, Rician noise); never an impossibility claim.

## Forward model
`S(b; S0, D, alpha, beta) = S0 * E_alpha(-(b D)^{beta/2})` (single diffusion time;
E_alpha = one-parameter Mittag-Leffler, CTRW/fractional Bloch-Torrey, Magin/Ingo 2013).
alpha = time-fractional (Caputo) order in (0,1]; beta = space-fractional (Riesz) order
in (1,2]. At alpha=1 it is exactly the CP0 stretched-exponential with heterogeneity
exponent beta/2. theta=(S0,D,alpha,beta) estimated JOINTLY; Rician noise sigma=S0/SNR.

## Degeneracy over the physiological (alpha,beta) grid (clinical n_b=6, b_max=2500, SNR=40)
- **median |rho_alpha_beta| = 0.984** (range [0.950, 0.998])
- median FIM condition number = 1.69e+08
- pre-registered degeneracy threshold: median |rho_alpha_beta| >= 0.9
- **verdict: degenerate = True**

## Headline cell (alpha=0.80, beta=1.70)
- analytic rho_alpha_beta = -0.979; FIM cond = 1.34e+08
- relative CRLB: cv_alpha = 0.56, cv_beta = 0.49
- bootstrap (200 reps) alpha_hat 95% CI = [0.348, 0.980] (truth 0.80); beta_hat 95% CI = [1.194, 2.000] (truth 1.70)
- empirical corr(alpha_hat, beta_hat) = -0.736 (the degeneracy ridge)

## n_b persistence at a single diffusion time (|rho_alpha_beta|)
| n_b | 4 | 6 | 8 | 16 |
|---|---|---|---|---|
| \|rho\| | 0.994 | 0.979 | 0.967 | 0.943 |

(Contrast: the CP0 single-order wall recedes for n_b>=8; this degeneracy does not.)

## Constructive boundary -- a second diffusion time breaks the degeneracy
| quantity | 16 b, one Delta | 8+8 b, two Delta |
|---|---|---|
| \|rho_alpha_beta\| | 0.943 | 0.182 |
| FIM condition | 4.22e+07 | 2.43e+06 |
| cv_alpha | 0.21 | 0.05 |

## Scoped claim
At a **single clinical diffusion time**, the joint CTRW time-order alpha and space-order
beta are **structurally degenerate** (median |rho_alpha_beta| = 0.984, FIM
condition ~2e+08): they cannot be separately recovered. Unlike the CP0
single-order wall, the degeneracy is **not** relieved by adding b-values; only a **second
diffusion time** separates them. This reinforces and extends the CP0 'clinically
information-limited' thesis to the two-exponent model.
