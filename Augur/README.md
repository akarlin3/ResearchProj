# Augur — the end-stage synthesis of the IVIM uncertainty program

**Augur** is the *perspective* paper that ties the IVIM uncertainty-quantification program together
into one arc — **trust → value-of-information → action** — across **Fashion** (the ruler),
**Minos** (the decision), **Lethe** (the limits), and **Gauge** (the identifiability wall), anchored
by the **"`D*` cross-modally orphaned"** thread. It makes no new measurement; it argues a synthesis.

> **Status: PROVISIONAL · speculative · SUBMISSION-BLOCKED.**
> Augur is built *before* its anchors have published. Every project-anchor is unpublished and pinned
> in [`ASSUMPTIONS.md`](ASSUMPTIONS.md); every external claim cites a **real, checked** source in
> [`CITATIONS.md`](CITATIONS.md). The paper is **not submittable** until **Fashion + Minos** publish
> (Lethe strongly recommended) — the block is enforced by [`check_anchors.py`](check_anchors.py) and
> documented in [`SUBMISSION_BLOCK.md`](SUBMISSION_BLOCK.md).

## The argument in one paragraph

A deployed quantitative-MRI error bar must first be **trusted** (Fashion's calibration ruler: a
skew-aware posterior restores marginal coverage, with a residual high-`D*` hole). Given trust, Minos
prices its **value of information** — the decision–calibration gap `G`, the second-order VoI law
`V = ½|EU″|G² = O(γ²)` (Plumbline Prop. 3, the *Delphi* result), and the label-free detectability
floor. On **action**, two walls appear: the interval is the wrong *size* for real test–retest (Lethe,
~4×) and the wrong *parameter* `D*` is unidentifiable (Gauge's CRLB wall). `D*` is the thread where
all three terminate: unidentifiable from its own signal, un-scalable to its own repeatability
(`r=−0.17`), and only weakly/inconsistently tied to DCE `Ktrans` (`r=0.389` Sun 2019; null Yang
2019). See [`synthesis.md`](synthesis.md).

## Layout

```
Augur/
  README.md            this file (status + arc)
  synthesis.md         the trust -> VoI -> action argument + the D* orphaned thread (cited)
  ASSUMPTIONS.md       pinned anchors (Fashion/Minos/Gauge/Lethe), SOLID vs PROVISIONAL, re-validation
  CITATIONS.md         verified citation list: Tier A (checked this build) + Tier B (inherited)
  SUBMISSION_BLOCK.md  the hard block + the only path to lift it
  check_anchors.py     submission-block gate (exits non-zero while load-bearing anchors unpublished)
  reproduce.sh         one-command re-validation (gate + tests)
  tests/test_augur.py  asserts block engaged + citations/anchors well-formed
```

## Reproduce / re-validate

```bash
bash reproduce.sh          # runs the submission-block gate + tests; prints PROVISIONAL/BLOCKED status
python3 check_anchors.py   # just the block gate (exit 1 = BLOCKED, the expected state today)
python3 -m pytest tests -q
```

## Provenance & IP

Created in-repo (clean, argument-only history; no data in tree or history), mirroring the
`Caliper/` · `Sextant/` · `Gnomon/` pattern in the monorepo. **Clean IP:** Augur touches no
`pancData3` / MSK / private clinical data — it cites in-repo synthetic/open-data results and published
literature only. Each anchor subrepo's own `README.md` is authoritative for its submission status.
