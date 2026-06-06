# Appendix drafts — robustness slice and noise test

Draft text for two referee-pre-empting appendices, plus a suggested edit to the
§7 open-problem paragraph. All numbers are filled from the committed results
JSONs (`absorption_results/beta010_results.json`,
`absorption_results/beta010_reduced_match.json`,
`reduced_results/beta010_placement.json`, `noise_results/noise_results.json`).
Figures: `paper_figures/figA.{pdf,png}` (Appendix A),
`paper_figures/figB.{pdf,png}` (Appendix B).

---

## Appendix A — Robustness of the phenomenology at β = 0.10 (~200 words)

To show the operating-corner phenomenology is not an artifact of the shipped
phase lag β = 0.05, we repeated the finite-N campaign at β = 0.10, the value at
which the reduced model is validated against the published regimes. Reduced-model
placement gives a saddle-node at A_SN = 0.178, a Hopf at A_H = 0.278, and a
homoclinic at A_hc(0.10) = 0.354 (bisection bracket width 9 × 10⁻⁵). We sampled
the post-homoclinic corner A = 0.5 (parallel to the operating corner) and, to
match its distance beyond the homoclinic, a depth-matched point A = A_hc + 0.090
= 0.444; A = 0.2 sits in the stable-chimera band as the never-absorber control.

Every headline result survives. The absorption lifetime is flat in N
(τ̂_abs(64)/τ̂_abs(8) = 0.93 at A = 0.5; 0.87 at A = 0.444; fitted exponential
rate ≈ −1 × 10⁻³ per oscillator, zero censoring). The censored-Weibull shape
stays above one (k_abs = 1.8–2.8 at A = 0.5, 1.4–1.6 at A = 0.444), so hazard
still rises. The per-cycle ratchet is positive (⟨ΔM⟩ = +0.06 to +0.19, bootstrap
CIs excluding zero), and true absorptions remain breath-phase-locked (pooled
Rayleigh p = 6.5 × 10⁻⁵). The reduced flow reproduces the breath period (≈ 23–25
s, measured 20–25 s) and spiral rate. Two honest caveats: the transient is
shallower at β = 0.10 (≈ 2 breath cycles to absorption vs ≈ 6 at the operating
corner), and A = 0.2 is a weaker control here (28–60 % persistent, sitting near
A_SN), with its persistent incoherent level R ≈ 0.66–0.68 close to the reduced
fixed point r\* = 0.73.

---

## Appendix B — A finite-size noise test of fluctuation re-injection (~250 words)

The reduced flow captures ≈ 3× faster than the finite-N system, and that
prolongation is N-independent (§7). The natural candidate — fluctuation-induced
re-injection — predicts that adding finite-size noise to the collective flow
should reproduce both facts. We tested this directly. We integrated the
three-variable reduced flow at the operating corner (A = 0.5, β = 0.05) with the
Euler–Maruyama scheme, adding Gaussian noise of amplitude σ_noise = c/√N to the
collective variables (ρ₁, ρ₂, ψ), with reflection keeping ρ in [0,1]. At c = 0
the scheme reproduces the deterministic capture times to a median 0.25 % per
initial condition. The physical scale, estimated from the measured high-frequency
fluctuation of R_incoh in the baseline data, is c ≈ 0.05. We swept N ∈ {8,16,32,64}
× c ∈ {0, 0.025, 0.05, 0.1, 0.2} × 200 realizations from the §6.4 seed-mapped
initial conditions.

The result is a clean negative. At the physical scale c ≈ 0.05 the noise barely
perturbs capture (prolongation factor 0.97, N-independent). Larger amplitudes do
prolong — at c = 0.2 the trajectory is repeatedly re-injected, accumulating up to
~240 breath cycles before capture — but the prolongation is strongly
N-_dependent_ (factor 2.4 at N = 8 falling to 0.7 at N = 64; coefficient of
variation 0.40 across N), and never reaches the measured ≈ 3 at any cell where it
is N-independent (Fig. B). The rising-hazard constraint holds only at small c
(k > 1 for c ≤ 0.05, k ≈ 1 for c ≥ 0.1), while breath-phase locking survives
throughout. Thus additive 1/√N noise on the collective flow reproduces the
mechanism (re-injection) but with the wrong N-scaling: because its amplitude is
N-dependent, so is the prolongation it produces. The N-independent ×3 factor
requires a different mechanism.

---

## Suggested edit to the §7 open-problem paragraph

Replace the two sentences beginning "A natural candidate is fluctuation-induced
re-injection …" through "… since the fluctuation amplitude is not." with:

> A natural candidate is fluctuation-induced re-injection — finite-size kicks
> scattering the trajectory back toward the ghost region the deterministic spiral
> is leaving. We tested the simplest version of this directly (Appendix B):
> additive Gaussian noise of amplitude c/√N on the collective variables of the
> reduced flow. It does produce re-injection and prolongation, but the
> prolongation is N-_dependent_ — strong at small N and absent by N = 64 — because
> the noise amplitude is. The observed N-independent factor of three therefore
> rules out additive collective noise as the mechanism and points to one whose
> effective re-injection strength does not scale away as 1/√N (for example
> structured or multiplicative fluctuations tied to the breath cycle).

The following sentence ("We report the puzzle as measured …") can stand, with
"leave its theory open" softened to "leave the responsible mechanism open," since
the additive-noise hypothesis is now closed off and the rising-hazard and
breath-locked-phase constraints (still sharp) are joined by a third: N-independence
of the prolongation.
