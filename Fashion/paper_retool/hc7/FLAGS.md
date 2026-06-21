# FLAGS — friction/hc7

Open items, deliberate omissions, and judgement calls a reviewer should see.

## 1. CP2 (real-data railing-aware exclusion) deliberately OMITTED — not clean

The optional Track-A real-data demonstration (excluding railed voxels from the
real OSIPI D\*-map → "improved coherence / fewer implausible D\*") was **not
done**, by design. Without ground truth it is **partly circular**: railed voxels
sit *at the parameter bounds by definition*, so excluding them mechanically
removes the most extreme D\* values and trivially reduces map dispersion — that is
a property of the definition, not independent evidence that the removed voxels were
genuinely unrecoverable. The truth-referenced CP1 (synthetic substrate) already
provides the non-circular version of exactly this exclusion benefit. Per the
prompt ("FLAG and omit if it is not clean; do not over-interpret"), omitted.

## 2. "diagnostic" retained in the TITLE — judged earned, not rewritten

The title keeps "...assumption-free pseudo-diffusion identifiability
**diagnostic**." CP1 demonstrates railing is a per-voxel flag of untrustworthy
D\* (precision/exclusion benefit), so "diagnostic" is earned in the bounded sense
the Discussion now states. Rewriting a submission-ready title was judged a larger,
non-minimal edit than warranted. **Author call:** if the title's "diagnostic"
should be narrowed further (e.g. "reliability flag"), that is a one-line change.

## 3. Synthetic railing magnitude ≠ real 54.2% (regime, pre-existing)

`B4` railed fractions (21% at SNR 20) are unconditioned-prior rates, below the
curated real ROI's 54.2% — the regime-concentration point already documented for
HC2/CS2. The subsection reports the flag's *operating characteristics on truth*,
not the headline magnitude; the two are kept explicitly distinct.

## 4. Precision degrades at high SNR (stated, not hidden)

At SNR 40 precision is exactly 0.50 — at low noise, fewer voxels are unreliable, so
a railed voxel is a coin-flip. The flag is most useful in the low-SNR /
weak-identifiability regime where it is most needed. The table reports all three
SNRs so this is visible, not averaged away.

## 5. "Unreliable D\*" threshold (τ = 0.5) is a choice

Primary definition: normalised D\* error $>0.5$ (off by >50%). A τ = 1.0
sensitivity is stored in the JSON (`sensitivity_tau1.0`) and shows the same
specific-not-sensitive pattern (precision still exceeds base rate; recall still
low). τ = 0.5 is pre-registered in the run, not tuned to the result.

## 6. NOT submission-ready — author framing review required

Final gate item 5: do **not** mark submission-ready until the author confirms the
"diagnostic" language matches the demonstrated precision/recall and that clinical
utility reads as future work, not as an asserted capability. PR is PROVISIONAL; do
not auto-merge.
