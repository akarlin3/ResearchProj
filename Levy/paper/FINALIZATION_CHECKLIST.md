# Levy manuscript -- finalization checklist (pre-submission)

The manuscript compiles offline and is internally number-consistent. The items below are the
deliberate placeholders to resolve at submission. Nothing here is fabricated in the draft.

## Venue / class
- [ ] Swap the wrapper to Springer **Nonlinear Dynamics** class:
      `\documentclass[smallextended]{svjour3}` + `\journalname{Nonlinear Dynamics}`.
      `svjour3.cls` ships with Springer's author kit; it is not in the offline TeX bundle used
      here, so the draft uses the house `article` style (`% TODO-VENUE` in `levy.tex`). All
      content, math, sectioning, and references port unchanged.

## Title / authors (GATE F)
- [ ] Confirm final title (`% TODO-TITLE`).
- [ ] Confirm author block (`% TODO-AUTHORS`): Annealing Signet Institute; Dept. of Applied
      Physics & Applied Mathematics, Columbia University; Dept. of Computer Science, University
      of Colorado Boulder; `ak5232@columbia.edu`.

## Forward-cited tooling: Ouroboros (`% FORWARD-CITE-ouroboros`)
- [ ] Decide the citation posture (author confirms at GATE F):
      (a) **submit-now, forward-cite** Ouroboros as a companion software release (current draft);
      (b) **hold** until Ouroboros is publicly archived, then cite the archive DOI.
- [ ] If (a): replace the `\bibitem{ouroboros}` note with the repository URL / "in preparation".
- [ ] If (b): insert the Zenodo/DOI once minted.
- Levy's net-new contribution (the Fisher/CRLB/identifiability layer) does NOT depend on
      Ouroboros; the dependency is only the reused Gr\"unwald--Letnikov operators and the
      A(alpha) noise-amplification cross-check (read-only, fully synthetic).

## Reference details (`% TODO-REF`) -- confirm against source, do not fabricate
- [ ] Coeurjolly & Istas (2001) -- exact volume/pages/DOI for the fBm/Hurst CRB.
- [ ] Spilling & Barrick (2022) -- exact journal/volume/pages/DOI (PMID 36054778 is verified).
- [ ] Poot et al. (2010) -- exact volume/pages/DOI for the DKI optimal-design CRLB.
- Verified DOIs already in the draft: Bennett 2003 (10.1002/mrm.10581), Magin et al. 2013
      (10.1016/j.micromeso.2013.02.054), Polders et al. 2011 (10.1002/jmri.22554).

## Pre-submission QA
- [ ] `bash paper/build.sh` -> 0 unresolved references, `levy.pdf` renders, figures present.
- [ ] `FULL=1 bash reproduce.sh` green (CP0 + CP1 + CP2 + consistency).
- [ ] `python paper/verify_citations.py` -> all `\cite`/`\bibitem` resolve; Ouroboros forward-cited.
