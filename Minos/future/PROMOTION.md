# PROMOTION.md — how `future/` becomes the real Minos paper

> **Documented, not executed.** This file records exactly how `future/` promotes into the real
> Minos manuscript once Fashion and Gauge publish and the re-run passes. **Do not run any of this
> now** — `future/` stays quarantined until the gate below is green against the *published* papers.

## Precondition gate (must be green before any promotion)

1. Fashion **and** Gauge are published (or their revisions are final).
2. `ASSUMPTIONS.md` rows updated to the published artifacts (DOIs / tags / commits).
3. `FULL=1 bash future/reproduce.sh` is **all-green** against those published versions:
   - CP1 theory gates PASS (unchanged — theory half never depended on the assumption);
   - CP2 applied gap PASS on the *published* Fashion ruler;
   - CP3 applied monitor PASS in the *published* Gauge framing;
   - CP4 manuscript consistency PASS (every number still traces).

If any stage fails, the revision genuinely invalidated a dependent result. Fix the result first;
do not promote a stale number.

## File moves (once the gate is green)

| from (`future/`) | to (real Minos) | note |
|---|---|---|
| `future/paper/minos.tex` | `Minos/paper/minos.tex` | the manuscript becomes the canonical Minos paper |
| `future/paper/consistency.py`, `build.sh` | `Minos/paper/` | build + consistency move with it |
| `future/applied/*.py` | `Minos/applied/` (new) | applied half graduates next to `minos-core/` and `theory/` |
| `future/results/*`, `future/figures/*` | `Minos/paper/{results,figures}/` | promoted artifacts |
| `future/ASSUMPTIONS.md` | fold into `Minos/paper/` data-availability + a short "dependencies" note | the manifest's *content* survives as provenance |
| `future/_paths.py` | adjust import roots to the post-move locations | the **only** code edit promotion requires |

`Minos/theory/` and `Minos/minos-core/` **do not move** — they were imported read-only throughout
and remain the validated theory half.

## Flags that clear

- Remove every **PROVISIONAL** marker in figures, numbers, and prose (they were the visible record
  of the assumption; once Fashion/Gauge are published the assumption is discharged).
- Replace "in review / in submission" language and pinned commit hashes with the published
  citations (DOIs) in `ASSUMPTIONS.md` lineage and the manuscript's references.
- Remove the speculative-build banner from `paper/minos.tex`.

## What stays put

- The theory half (`theory/`, `minos-core/`) — untouched.
- This `PROMOTION.md` and the final `ASSUMPTIONS.md` snapshot — kept in git history as the record
  of how the speculative build was discharged.
