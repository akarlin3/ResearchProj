# FLAGS.md — open items for author decision (railing-first redraft)

What is *done* is in `CHANGES.md`; ledger compliance is in `LEDGER_AUDIT.md`. This file
lists only what needs an author decision or external confirmation before submission.

## needs-author-signoff (binding — do not mark submission-ready until resolved)
- **CP4 novelty framing.** The acceptance axis at NMR in Biomedicine is originality. The
  redraft frames the contribution as the *first systematic real-data quantification of
  D\* boundary-railing across cohorts + honest conditional-coverage characterization +
  calibration resolution* (§1, §4, both cover letters). A prompt cannot self-validate
  this "first" claim against the entire literature — the author must confirm the framing
  and the two deltas (vs known D\* instability; vs Casali 2026) hold.

## needs-author-decision
- **Cover-letter variant.** Two variants are provided (`cover_letter_phenomenon.tex`,
  `cover_letter_reproducibility.tex`). Pick one (or merge); only one is submitted.
- **Citation style.** The prompt specified Wiley/APA author-year, which the manuscript
  uses (inline reference list). NMR in Biomedicine's house style is a numbered (Vancouver)
  scheme; if the journal requires numbered citations at submission, convert the inline
  APA list — content is unchanged, only the in-text marker/ordering. `refs.bib` is kept as
  a coherent parallel record to ease that conversion.
- **Abstract length.** The single abstract is ~300 words of readable prose (raw token
  count 316 including LaTeX macros such as `\Dstar{}`). If the venue enforces a hard
  300-word ceiling on the typeset text, trim one clause from the closing sentence.

## needs-external-confirmation (identifiers — not fabricated)
- **Journal-article DOIs.** DOIs that could not be independently verified were *omitted*
  from the manuscript and `refs.bib` rather than asserted (Le Bihan 1988, Koh 2011,
  Lemke 2011, Federau 2012, While 2017, Orton 2014, Meeus 2017, Cranmer 2020). Confirm and
  reinstate at typesetting against each publisher of record. **Verified/retained DOIs:**
  OSIPI Zenodo 10.5281/zenodo.14605039; TCGA-LIHC 10.7937/K9/TCIA.2016.IMMQW8UQ (from
  Sextant provenance); Casali 10.1002/nbm.70227, Jallais 10.7554/eLife.101069,
  Manzano-Patrón 10.1016/j.media.2025.103580 (from the audited bib); Gurney-Champion 2018
  10.1371/journal.pone.0194590 (PLoS-deterministic); Clark 2013 10.1007/s10278-013-9622-7
  (canonical TCIA).
- **Hero & Fessler (1993).** Exact volume/pages/DOI of the MWSCAS proceedings entry to be
  confirmed at typesetting (flagged in `refs.bib`); the citation itself is real.

## out-of-scope by design (not flagged for this paper — recorded for traceability)
- Amortized-NPE CRLB-floor efficiency audit ("overconfident below the floor"), the OOD
  self-consistency gate (AUC 0.99/0.59), inference timing, and the in-vivo brain
  held-out-b limb are deliberately excluded (ledger OUT-OF-SCOPE / not in keep-set). They
  remain available for a separate methods paper.

## reproducibility note
- The 245 MB OSIPI zip is fetched on demand and git-ignored (sha256-verified
  2a53054d…b3e); it is not committed. The TCGA-LIHC arrays are likewise download-on-demand
  and git-ignored. A fresh clone runs `Gnomon/scripts/build_handoff.py` and
  `Sextant/reproduce.sh` to regenerate every number.
