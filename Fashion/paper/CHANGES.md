# CHANGES — assessor-remediation run

Branch `fix/assessor-remediation`. One commit per checkpoint. Every edit verified against the
number-freeze gate (`paper/check_numbers.sh`) — no reported result value changed.

All manuscript/supplement line references are to the files as of the final commit.

| CP | Issue (assessor) | Fix | File / location |
|----|------------------|-----|-----------------|
| 0 | No editable source existed | Reconstructed canonical LaTeX (manuscript + supplement + refs.bib) verbatim from the latest merged MRM PDF; froze all numbers | `paper/manuscript.tex`, `paper/supplement.tex`, `paper/refs.bib`, `Fashion/NUMBERS_FROZEN.txt` |
| 1 | Data Availability missing (Format FAIL, Medium) | Added a Data Availability Statement pointing to the Zenodo reproducibility archive (`doi:10.5281/zenodo.20649669`); lists code, calibration grid, efficiency map; notes OSIPI offer | `manuscript.tex` — Data availability section (Declarations block) |
| 2 | In-vivo datasets unnamed/uncited (Provenance FAIL, Medium) | Named + cited both in-vivo datasets to their single recoverable source, the OSIPI TF2.4 IVIM-MRI Code Collection data (Gurney-Champion et al., Zenodo 2025, `doi:10.5281/zenodo.14605039`); new reference `[22]`; markers at Methods Part 2, Results (brain + abdomen), Limitations, Fig 4 caption; accession named in Supp S3 and S4 | `manuscript.tex` (ref [22]; 5 in-text markers), `supplement.tex` (S3, S4) |
| 3 | N=500 vs ~2,000 unreconciled (Norm gap / R1 Minor) | Verified in code (`npe/run_f_realdata.py`, `npe/run_g_ood_gating.py`): same brain image, same gray-matter mask, same held-out-b partition, same RNG seed (42), differ only in `--n-voxels`. The 500 is **not** a strict subset of the 2,000 for plausible voxel counts, so stated the accurate relationship: overlapping random draws from the same gray-matter ROI, sized differently because a per-voxel ROC (Supp S3) needs more voxels than a coverage curve (Fig 4B) | `manuscript.tex` Methods Part 2; `supplement.tex` S3 |
| 4 | R2 contribution/positioning (prose only) | Reframed to the transferable diagnostic (information-floor / CRLB auditing as a general check for learned qMRI posteriors; aggregate-passes/pointwise-fails as the portable lesson; NPE as worked example). Added a positioning paragraph with one-sentence deltas vs Casali 2026 [19], µGUIDE/Jallais & Palombo 2024 [20], Manzano-Patrón 2025 [21]. Stated explicit scope early (end of Introduction). Zero new numbers, only existing citations | `manuscript.tex` — Introduction |
| 5 | Prose density / long sentences (R3 Minor; MRM return cause) | Split four of the longest compound sentences (Introduction + Discussion) at em-dash/semicolon joins; meaning and numbers preserved exactly | `manuscript.tex` — Introduction, Discussion (amortized-posteriors, limitations) |
| 6 | Venue retarget → MRI (Elsevier) | Title page venue swap (MRM → Magnetic Resonance Imaging (Elsevier)); Elsevier declarations block (CRediT, Declaration of competing interest, Funding, Data availability, Declaration of generative AI — reusing the existing AI-disclosure text); reference style → Elsevier bracketed numeric ([N] in-text + numbered list); fresh cover letter to EIC John C. Gore | `manuscript.tex` (title, declarations, citations); `paper/cover_letter.tex` (new) |

## What was deliberately *not* done

- **No fabrication.** Every identifier (repro DOI, dataset accession) was recovered from the
  repository; **no `[[TODO]]` placeholder was needed** in the manuscript.
- **No numerical change.** All coverage values, CRLB ratios, percentages, Ns, AUCs, and the
  21 reference entries are byte-for-byte as transcribed from the source PDF.
- **Reference style:** in-text citations are Elsevier bracketed numeric and the list is
  numbered, but the 21 entries' internal text was kept verbatim rather than re-flowed into
  the `elsarticle-num` micro-format (initials-first), to honor the "content unchanged"
  constraint. `elsarticle.cls` was verified to build under `tectonic`; a production
  Elsevier-class typesetting can be generated from the committed `refs.bib` without touching
  entry content. The structured abstract was preserved (which `elsarticle`'s abstract
  environment does not natively support), another reason the working `article` build was kept.

## Build

`cd Fashion/paper && make` → `manuscript.pdf`, `supplement.pdf`, `cover_letter.pdf`
(tectonic; clean, no errors, all cross-references resolve).
