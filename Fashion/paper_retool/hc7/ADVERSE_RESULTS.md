# ADVERSE_RESULTS — friction/hc7 (railing utility)

Per the honesty-first contract, the CP1 demonstration is recorded here verbatim —
including the part that does **not** simply confirm "diagnostic." The decisive
verdict is the first entry. Numbers trace to `Sextant/results/phantom_recovery.json`
key `B4_flag_utility` (seed 20260621); see `NUMBERS_FROZEN.txt`.

Legend: **TEMPERED** = claim survives but must be narrowed.
**CONFIRMATORY (concession required)** = survives but a concession must be stated.

---

## CP1 verdict — railing predicts D\* failure, but is SPECIFIC, NOT SENSITIVE — **TEMPERED (the central HC7 result)**

On synthetic ground truth (N=3000/SNR over the trained prior, Fashion's exact
bounded NLLS; "unreliable D\*" $=$ truth-referenced normalised error $>0.5$,
convention-free):

| SNR | base rate | precision $P(\mathrm{unrel}\,|\,\mathrm{railed})$ | recall $P(\mathrm{railed}\,|\,\mathrm{unrel})$ | lift | excl. railed → pooled-relErr reduction |
|----:|----:|----:|----:|----:|----:|
| 10  | 0.538 | **0.68** [0.65, 0.71] | 0.39 [0.37, 0.42] | 1.27 | −0.081 [0.064, 0.103] (14.5%) |
| 20  | 0.414 | **0.61** [0.57, 0.65] | 0.31 [0.28, 0.33] | 1.47 | −0.052 [0.038, 0.064] (13.6%) |
| 40  | 0.294 | **0.50** [0.46, 0.54] | 0.28 [0.26, 0.32] | 1.70 | −0.026 [0.018, 0.034] (10.1%) |

**Not null.** Railing concentrates unreliable D\* well above base rate (lift
1.27–1.70), and *acting* on the flag — excluding railed voxels — lowers the
retained map's pooled D\* error by 10–15% at every SNR, with bootstrap CIs that
exclude zero. The "diagnostic" word is therefore **earned for flagging
untrustworthy D\***.

**But it must be narrowed.** Recall is only 0.28–0.39: most unreliable voxels do
**not** rail, because D\* is weakly identified *broadly* and a fit can carry a
badly recovered interior D\* without reaching a bound. As a binary classifier the
flag's balanced accuracy is ~0.58 — only modestly above chance, consistent with
the near-chance above-median-error AUC already on record (B3: 0.555/0.515). At
SNR 40 precision falls to exactly 0.50: when noise is low, fewer voxels are
unreliable at all, so a railed voxel becomes a coin-flip. **Implication for the
paper:** a *positive* rail is an actionable, high-precision flag; the *absence* of
a rail is not a clean bill of health. The manuscript claims a flag for
untrustworthy D\*, **not** a complete detector — stated as such in §"The railing
flag is actionable" and the Discussion.

---

## Clinical utility is asserted-as-stakes but NOT demonstrated — **CONFIRMATORY (concession required)**

The motivation (perfusion/D\* is clinically meaningful) is retained, but no run in
hand can show that excluding/down-weighting railed voxels changes a downstream
diagnosis or outcome — that needs patient-outcome data a simulation cannot supply.
**Concession (now explicit in the Discussion):** the demonstrated claim is
*measurement reliability* (railing flags untrustworthy D\*), and *clinical
decision utility is deferred to future work.* The friction (HC7: "diagnostic is
partly aspirational because its actionable downstream is deferred") is resolved by
demonstrating the measurement-reliability downstream and scoping the clinical
downstream out, rather than by leaving "diagnostic" to imply a clinical instrument.

---

## The synthetic railing magnitude (~21% at SNR 20) is not the real 54.2% — **TEMPERED (regime, not a defect; pre-existing)**

`B4_flag_utility` railed fractions (0.309/0.210/0.167 at SNR 10/20/40) are the
*unconditioned-prior* rates and are lower than the curated real ROI's 54.2% — the
same regime-concentration point already documented for HC2/CS2 (the real ROI sits
in the low-f/high-D\* corner). The CP1 subsection is explicit that it reports the
flag's *operating characteristics on truth*, not the headline magnitude; the two
are not conflated. No new claim is made about reproducing 54.2% in simulation.
