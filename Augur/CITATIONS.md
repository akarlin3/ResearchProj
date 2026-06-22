# CITATIONS.md — Augur's verified citation list

> **No phantom citations (the Ouroboros lesson).** Every external claim in `synthesis.md` traces to
> a row below. **Tier A** was verified *this build* against the primary source (with a verbatim
> quote). **Tier B** are framing references inherited from Minos/Plumbline §7 and `minos-core`
> POSITIONING (checked there); they must be independently re-verified at submission time (part of the
> publish-gated finalization — see `SUBMISSION_BLOCK.md`).

## Tier A — verified this build (primary source checked + quoted)

### A1. Sun et al. 2019 — `D*`–`Ktrans` weak but significant
- **Citation:** Sun H, Xu Y, Xu Q, Duan J, Zhang H, Liu T, Li L, Chan Q, Xie S, Wang W.
  *Correlation Between Intravoxel Incoherent Motion and Dynamic Contrast-Enhanced Magnetic Resonance
  Imaging Parameters in Rectal Cancer.* **Academic Radiology**, 2019.
- **DOI:** 10.1016/j.acra.2018.08.012 · **PMID:** 30268719
- **Used for:** synthesis §4(3) — cross-modal `D*`–`Ktrans` correlation.
- **Verified quote:** *"relatively weak correlations between D\* and Ktrans (r = 0.389; p < 0.001)"*;
  also *"Moderate correlations were found between f·D\* and Ktrans (r = 0.533; p < 0.001)."*
- **Checked:** PubMed 30268719, this build (2026-06-22).

### A2. Yang et al. 2019 — `D*`–`Ktrans` non-significant; moderate `D*` reproducibility
- **Citation:** Yang X, Xiao X, Lu B, Chen Y, Wen Z, Yu S. *Perfusion-sensitive parameters of
  intravoxel incoherent motion MRI in rectal cancer: evaluation of reproducibility and correlation
  with dynamic contrast-enhanced MRI.* **Acta Radiologica**, 2019.
- **DOI:** 10.1177/0284185118791201 · **PMID:** 30114928
- **Used for:** synthesis §4(3) — cohort-inconsistency of the cross-modal correlation; `D*`
  reproducibility.
- **Verified quote:** *"There was no significant correlation between ve and f, ve and D\*, ve and
  f·D\*, D\* and Ktrans, and D\* and kep"*; reproducibility *"for parameter D\*, ICC = 0.55
  (0.32–0.72), CV = 20.28 ± 3.23%."*
- **Checked:** PubMed 30114928, this build (2026-06-22).

## Tier B — framing references (inherited from Minos/Plumbline; re-verify at submission)

| key | citation | used for | source-of-record |
|---|---|---|---|
| loss-cal-2011 | Lacoste-Julien, Huszár & Ghahramani, *Approximate inference for the loss-calibrated Bayesian*, **AISTATS 2011** | §2 — loss-calibration baseline | Plumbline §7 |
| loss-cal-2021 | Vadera et al., *Post-hoc loss calibration*, **UAI 2021** | §2 — post-hoc loss-calibration | Plumbline §7 |
| dec-cal-2021 | Zhao, Kim, Sahoo, Ma & Ermon, *Calibrating Predictions to Decisions* (decision calibration), **NeurIPS 2021** | §2 — decision calibration | Plumbline §7 |
| dca-2006 | Vickers & Elkin, *Decision curve analysis*, **Med Decis Making 2006** | §2 — net-benefit contrast | Plumbline §7 |
| ood-2019 | Ovadia et al., *Can you trust your model's uncertainty?*, **NeurIPS 2019** | §3 — calibration under shift | Plumbline §7 |
| casali-2026 | Casali et al., *A Comprehensive Framework for UQ of Voxel-wise Supervised Models in IVIM MRI*, **NMR in Biomedicine 2026** (arXiv:2508.04588) | §1/§4 — documented `D*` overconfidence | Gauge `refs.bib` (`casali2026`) |

## In-repo anchors (cross-references, not external citations) — all PROVISIONAL

| anchor | claim used | location |
|---|---|---|
| Fashion | skew-aware posterior restores marginal coverage; residual high-`D*` conditional gap | `Fashion/paper_retool/`, `Gnomon/handoff/CLAIMS_LEDGER.md` |
| Minos/Plumbline | gap law; `O(γ²)` VoI; label-free floor | `Minos/theory/plumbline.md` (Thm 1–2, Prop. 3) |
| Gauge | `D*` identifiability wall; `D*` test–retest `r=−0.17` [−0.39,0.05] | `Gauge/results/conditional_attack_report.txt`, Gauge paper §4.2.2 |
| Lethe | conformal `D` interval ~4× too narrow for test–retest (0.263 vs 0.755) | `Lethe/LETHE.md`, `Lethe/results/RESULTS_VALIDATION.md` |
