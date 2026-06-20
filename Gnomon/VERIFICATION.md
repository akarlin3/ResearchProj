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

## CP2 — clean-room rebuild  ✅

Gate: **runs, self-consistent.**

- [x] forward (`forward.py`), NLLS + railing (`nlls.py`), Laplace + per-voxel MCMC
      (`bayes.py`), MAF/NPE flow (`flow.py`), independent ruler (`metrics.py`),
      bootstrap (`bootstrap.py`), OSIPI loader (`osipi.py`), CP3 driver
      (`reproduce.py`) all implemented from spec — no Caliper import.
- [x] **Self-consistency gates pass** (`tests/test_cp2_selfconsistency.py`, 7 cases):
      analytic Jacobian matches finite differences; clean-signal NLLS round-trip
      recovers truth (no railing on clean data); continuity (f=0 → mono-exponential)
      holds; cohort draws ground truth from Lattice; ruler coverage exact; the
      end-to-end driver runs and emits a verdict + JSON.
- [x] Full suite green: `python -m pytest Gnomon/tests -q` → **16/16**.
- [x] Every design choice documented inline + in `docs/METHODS.md` (b-schemes, NLLS
      box/init/railing, MCMC sampler spec, MAF training spec, OSIPI ROI selection).

## CP3 — the reproduction gate  ⏳  (HARD HALT either way)

Gate: run the rebuild, compare to the frozen manifest with bootstrap CIs.
→ **REPRODUCES** (proceed to CP4) or **DOES NOT REPRODUCE** (divergence report, stop).

## CP4 — package for the retool  ⏳  (reproduce branch only)

Clean implementation + complete methods + reproduced numbers handed to the Fashion
retool; merge-back documented.
