# VERIFICATION.md — Gnomon checkpoint ledger

Each checkpoint records what was *run* and what gate it had to clear. Nothing is
claimed here that was not executed.

## CP1 — scaffold + targets manifest  ✅ (this commit)

Gate: **mirrors the sibling subrepos, clean, self-consistent.**

- [x] Embedded as a top-level subrepo with its **own clean root history**, merged
      into the monorepo via `git merge --allow-unrelated-histories` (mirrors
      Lattice/Datum/Lethe). `git log -- Gnomon/` shows Gnomon's own commits.
- [x] Registered in the root `README.md` (Contents, Projects-at-a-glance, Project
      details, How-they-fit-together, Provenance).
- [x] **Reproduction-targets manifest frozen** (`gnomon/manifest.py`): 5 targets
      (T1, T3a, T3b, T3c, T4) with claimed values, frozen tolerances, and
      Fashion-**prose** provenance. `manifest.validate()` passes.
- [x] **Clean-room boundary enforced:** `_paths.py` allows `lattice` only; `caliper`
      in `FORBIDDEN`; static AST test confirms no `gnomon/` module imports Caliper.
- [x] **Clean IP:** no proprietary/clinical data, no data-like files in tree.
- [x] CP1 tests pass: `python -m pytest Gnomon/tests -q` → (recorded at commit).

## CP2 — clean-room rebuild  ⏳

Gate: **runs, self-consistent.** forward/NLLS/Laplace/MCMC/MAF/metrics implemented
from spec; clean-signal round-trip recovers truth; continuity limits hold; every
design choice documented inline + in `docs/METHODS.md`.

## CP3 — the reproduction gate  ⏳  (HARD HALT either way)

Gate: run the rebuild, compare to the frozen manifest with bootstrap CIs.
→ **REPRODUCES** (proceed to CP4) or **DOES NOT REPRODUCE** (divergence report, stop).

## CP4 — package for the retool  ⏳  (reproduce branch only)

Clean implementation + complete methods + reproduced numbers handed to the Fashion
retool; merge-back documented.
