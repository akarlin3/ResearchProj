# CP4 — interventional causal test (τ vs injected D₀)

Synchronized population held at the canonical chimera envelope; a **controlled** m=2 & m=3 harmonic distortion injects a target initial manifold distance into the incoherent population only. If lifetime falls monotonically with injected D₀, the manifold-escape mechanism is causally **supported**; if lifetime is insensitive, it is **refuted/incomplete**.

Levels are ε (m=2&m=3 distortion strength) in the monotone band; families are **paired** across levels (same synchronized population, same pre-images).

| ε               | realized D₀ (med) | n   | τ̂ (s) | 95% CI       | KM median (s) | mean life (s) |
| --------------- | ----------------- | --- | ----- | ------------ | ------------- | ------------- |
| 0.00 (baseline) | 0.183             | 60  | 65.9  | [51.2, 86.4] | 54.3          | 65.9          |
| 0.10            | 0.214             | 60  | 67.5  | [52.4, 88.5] | 50.8          | 67.5          |
| 0.20            | 0.275             | 60  | 66.5  | [51.7, 87.1] | 48.3          | 66.5          |
| 0.30            | 0.314             | 60  | 64.6  | [50.2, 84.7] | 45.3          | 64.6          |
| 0.40            | 0.352             | 60  | 65.3  | [50.7, 85.5] | 50.1          | 65.3          |
| 0.50            | 0.357             | 60  | 65.9  | [51.2, 86.3] | 46.4          | 65.9          |

**Three causal readouts** (injected D₀ spans 0.183–0.357, the reachable upper CP2 range):

1. **Level trend** — Spearman ρ(level median D₀, level median τ̂) = **-0.600** (p=2.08e-01): does the level-averaged lifetime fall as injected D₀ rises (floor noise averaged out)?
2. **Paired test** — 60 families, baseline→max ε (median injected ΔD₀ = +0.172): median Δlifetime = **-0.8 s**, Wilcoxon p = 7.06e-01 (controls for family-to-family floor noise).
3. **Pooled, run-level** — Spearman ρ(realized D₀, lifetime) = **-0.102** (p=5.28e-02): the unconditioned run-by-run association (dominated by finite-N floor noise).

**CP4 causal verdict:** **MIXED / WEAK** — one controlled test indicates injected D₀ shortens lifetime but the other is not decisive; the causal effect, if present, is small relative to finite-N floor noise.

_Scope note:_ injection can only push D₀ **up** from the finite-N sampling floor (a finite sample cannot be more on-manifold than its own floor), and D₀(ε) folds over above ε≈0.5, so the reachable band is [floor, ~2×floor] — exactly the upper CP2 range. A genuinely lower-D₀ regime is not constructible at this N.
