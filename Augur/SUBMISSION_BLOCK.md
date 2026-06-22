# SUBMISSION_BLOCK.md — Augur is NOT submittable yet

**Hard block.** Augur is the end-stage synthesis of the IVIM-UQ program. It synthesizes results from
projects that are **all unpublished**. Submitting Augur before its load-bearing anchors publish would
(a) cite results that may still change in review, and (b) front-run its own dependencies.

## The rule

Augur **may not be submitted** until both load-bearing anchors are published:

- **Fashion** (the trust/ruler anchor) — in review at *NMR in Biomedicine*.
- **Minos** (the value-of-information/decision anchor) — provisional.

**Lethe** (the wrong-size action limit) and **Gauge** (the identifiability wall) are strongly
recommended to be out as well, since §3 of the synthesis rests on them.

## How the block is enforced

`check_anchors.py` pins each anchor's published state and **exits non-zero while any load-bearing
anchor is unpublished**. `reproduce.sh` runs it; `tests/test_augur.py` asserts the block is engaged.
Today, all three return "BLOCKED" — the intended state.

## Lifting the block (the only path)

Per `ASSUMPTIONS.md §3`:

1. When an anchor publishes, set `published=True` + its DOI in `check_anchors.py`'s `ANCHORS` table
   and update `ASSUMPTIONS.md §1`.
2. Run `bash reproduce.sh`.
3. The block lifts only when Fashion **and** Minos are published. Then **re-verify CITATIONS.md
   Tier B** against primary sources before any submission.

Until every step is done, treat Augur as a working draft, not a submittable paper.
