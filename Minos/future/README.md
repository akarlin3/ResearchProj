# `future/` — the full Minos paper, built now under the Fashion+Gauge-survive assumption

This folder assembles the **complete** Minos paper *before* its two upstream dependencies have
published. It wires the already-validated **theory half** (`Minos/theory` + `Minos/minos-core`,
imported read-only) to a new **applied/decision half** that consumes **Fashion**'s calibrated IVIM
posteriors and is contextualized by **Gauge**'s coverage / high-D\* identifiability wall / label-free
monitor.

**It is speculative by construction.** Fashion and Gauge are in review; this build assumes they land
as submitted. Read [`ASSUMPTIONS.md`](ASSUMPTIONS.md) first — it pins every Fashion/Gauge input and
defines what is **SOLID regardless** vs **PROVISIONAL** (assumption-dependent). Every figure and
number that depends on the assumption is marked PROVISIONAL where it appears.

## Quarantine

- All new work lives **under `future/`**. The validated theory half and the main manuscript files
  outside `future/` are **not modified** — `future/` depends on them by **read-only import** only
  (`_paths.py`).
- Promotion into the real Minos paper is documented in [`PROMOTION.md`](PROMOTION.md) and is **not
  executed** here.

## Layout

```
future/
  README.md          this file
  ASSUMPTIONS.md     pinned Fashion+Gauge inputs; SOLID vs PROVISIONAL split
  PROMOTION.md       how future/ promotes into the real Minos paper (documented, not run)
  reproduce.sh       ONE command: theory gates -> CP2 -> CP3 -> CP4 consistency
  _paths.py          read-only sys.path wiring to theory, minos-core, Fashion, Gauge
  verify_cp1.py      CP1 gate: theory reproduces via the import; wiring resolves
  applied/           CP2/CP3 applied code (gap_applied.py, monitor_applied.py, data.py)
  results/           seeded numeric printouts every paper number traces to (PROVISIONAL)
  figures/           PROVISIONAL-marked figures
  paper/             CP4 full manuscript (minos.tex) + consistency.py + build.sh
```

## Environment

Use the `proteus` conda env (has numpy / scipy / sympy / matplotlib):

```bash
PROT=/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python
$PROT Minos/future/verify_cp1.py          # CP1 gate (FAST); add --full for full-N GATE 3
bash   Minos/future/reproduce.sh          # one-command re-validation of the whole build
```

The theory half is deterministic (seeded); the applied half pins Fashion's `SEED=0` and Gauge's
`DEFAULT_SEED=20260613` (see `ASSUMPTIONS.md`).

## Status

- **CP1 — scaffold + wire theory + manifest:** done. Gate `verify_cp1.py` PASS (theory
  reproduces via read-only import; wiring resolves; quarantine clean).
- **CP2 — applied decision–calibration gap:** done. Gate `applied/gap_applied.py` PASS.
  Honest result: qualitative gap holds on Fashion-calibrated posteriors; quantitative
  leading-order law does **not** transfer (median |G_std/G_theory| ≈ 17×). [`results/RESULTS_CP2.md`](results/RESULTS_CP2.md)
- **CP3 — applied validity monitor:** done. Gate `applied/monitor_applied.py` PASS.
  Observable-fires (AUC up to 1.00) / hidden-blind (AUC 0.52) reproduced on synthetic IVIM.
  [`results/RESULTS_CP3.md`](results/RESULTS_CP3.md)
- **CP4 — full draft assembly:** done. `paper/minos.tex` compiles (tectonic, ebgaramond +
  microtype) and `paper/consistency.py` PASS (20 numbers traced to seeded results).
