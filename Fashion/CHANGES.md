# CHANGES.md — railing-first redraft for NMR in Biomedicine (`redraft/railing-first-nmrbiomed`)

Re-spines Project Fashion to lead with the real-data boundary-railing phenomenon,
built on the **Gnomon** clean-room reference (numeric source of truth) and integrating
**Sextant** (cross-cohort generalization). Audit-first, run-then-write: no number
appears in prose before it is produced by a run. Ledger dispositions are audited in
`LEDGER_AUDIT.md`; open items in `FLAGS.md`; frozen numbers in `NUMBERS_FROZEN.txt`.

## CP0 — enumerate, freeze, plan
- Imported the Gnomon CP4 clean-reference handoff package into the redraft so the paper
  is self-contained and one-command reproducible (`Gnomon/scripts/build_handoff.py`,
  `gnomon/reframe.py`, `handoff/CLAIMS_LEDGER.md`, `handoff/conditional_coverage.json`,
  `RETOOL_HANDOFF.md`, `tests/test_cp4_reframe.py`).
- Re-ran `build_handoff.py` (proteus env; OSIPI zip sha256-verified) — every keep-set and
  reframe number reproduced to the digit; verdict PARTIAL; 17/17 Gnomon tests green.
- Rewrote `NUMBERS_FROZEN.txt` as the railing-first frozen-targets file (was an obsolete
  PDF number-frequency dump from the ruler-first paper).

## CP1 — re-spine to railing-first (KEEP K1)
- New title and thesis: D\* weak identifiability manifests as a measurable, reproducible
  real-data failure. Lead result (Abstract + §3.1 + Fig 1): NLLS rails D\* in
  **54.2% [52.0,56.4]** of high-SNR OSIPI abdomen voxels.
- Moved the amortized-NPE/CRLB-floor headline, the OOD gate, and the brain held-out-b
  limb to a brief explicit out-of-scope note (§4).

## CP2 — marginal → conditional coverage (REFRAME R1/R3, DROP D1)
- Replaced the retracted marginal 0.30/0.67 everywhere with the per-true-D\*-tercile
  conditional table (Table 1; both SD conventions; honest CRLB recommended — high-D\*
  0.63 Laplace / 0.81 MCMC). SD convention documented in Methods §2.4; the floored
  convention is shown only to explain the manufactured severity (pooled 0.68 / high 0.41).
- 0.30/0.67 survive in the body only as *retired* numbers in the reframe narrative (§3.2).

## CP3 — resolution framing (KEEP K2/K3, REFRAME R2)
- Quantile intervals restore marginal D\* coverage (0.90; Fig 3A); amortized flow beats
  railed NLLS on coverage/ECE/sharpness (0.98/0.069/0.112 vs 0.76/0.121/0.181; all gaps'
  CIs exclude 0; Fig 3B) — presented as the *resolution*, not the headline.
- Kept the residual high-D\* wall (0.81) as the identifiability limit (R2).

## CP4 — novelty positioning (binding lever)
- Stated up front in §1 and in both cover letters: first systematic real-data
  quantification of D\* railing across cohorts + honest conditional-coverage
  characterization + calibration resolution; one-sentence deltas vs known D\* instability
  and vs Casali 2026.

## CP5 — Sextant integration
- §3.1 + Table S1: cross-cohort generalization — full-ROI 47.8% [47.1,48.5] and
  independent TCGA-LIHC liver 43.7% (4-b) / 73.4% (3-b), rail direction and generous-bound
  control. Clearly labelled as the Sextant (non-Gnomon) generalization layer.

## CP6 — NMR Biomedicine format
- Graphical abstract: 69-word text (≤80, 2 sentences) + a 50×60 mm figure
  (`graphical_abstract.pdf`).
- Single ~300-word unstructured abstract (converted from the old structured abstract).
- Data-availability statement (repo one-command reproduction + OSIPI Zenodo 14605039 +
  TCGA-LIHC) and a data-ethics statement for the open in-vivo data; author contributions.
- Two cover letters to EIC John R. Griffiths: phenomenon-led and reproducibility/rigor-led.

## CP7 — references
- Converted to Wiley/APA author-year (inline reference list). Added OSIPI (Zenodo
  14605039), Hero & Fessler (1993), and TCGA-LIHC (TCIA; verified DOI from provenance) +
  the TCIA infrastructure citation. Removed 12 references orphaned by the dropped material;
  curated `refs.bib` to the same 18-entry set. Unverifiable journal DOIs were omitted
  rather than fabricated (see `FLAGS.md`).

## Figures & build
- New figures generated directly from frozen run outputs (`make_railing_figures.py`):
  fig1 railing-by-cohort, fig2 conditional coverage, fig3 resolution, graphical abstract.
- Removed the orphaned ruler-first figures (fig1–4, figS1–6) and old mockups.
- `manuscript.tex`, `supplement.tex`, and both cover letters build clean under tectonic.
