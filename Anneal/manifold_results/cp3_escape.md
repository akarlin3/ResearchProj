# CP3 — time-resolved manifold escape

60 quantile-stratified traces (N∈{16,32}, A=0.5). The question: does D_incoh **rise before** the order-parameter collapse (manifold escape leads collapse), or only at/after it?

**Primary test — event-aligned, per-run paired comparison of ⟨D_incoh⟩** (pre-collapse window [−25,−10] s vs peri-collapse [−3,0] s; 56 runs):

- ⟨D⟩ pre = **0.146**, peri = **0.125**, post-collapse = **0.088**.
- Per-run paired change peri−pre: median ΔD = **-0.004**, Wilcoxon p = 5.74e-01, fraction rising = 48%.

**Secondary (per-run lead time, robust baseline median+4·MAD over the first half, sustained ≥2 s):** fired in 16/56 runs, median lead 25.3 s.

**CP3 verdict:** **NO precedence** — ⟨D_incoh⟩ is statistically flat from the pre-collapse window into collapse (it sits at the finite-N floor throughout) and in fact RELAXES after the merger (post < peri: the globally-synchronized end state is _more_ on-manifold). D does not lead the order-parameter collapse. The per-run 'lead times' are an artifact of breathing-driven D fluctuations crossing a low-variance baseline, NOT a genuine pre-collapse escape ramp — which is exactly why the primary test uses a paired pre-vs-peri comparison.

Figures: `cp3_aligned.png` (event-aligned ⟨D⟩ and ⟨R⟩ — note D is flat while R rises through θ at collapse), `cp3_spaghetti.png` (per-run traces).
