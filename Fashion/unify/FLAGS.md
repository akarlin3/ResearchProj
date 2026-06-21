# FLAGS — contradictions, Gnomon divergences, unknowns (Fashion + Sextant unify pass)

Read-only diagnosis. Items ordered by how hard they block a unification route. #1 and #2 are
**blocking**.

---

## #1 — DECISION (was BLOCKING): which railing-first variant is the submission draft?

**Update after inspecting the remote:** the ambiguity is largely *resolved on `origin/main`*.
`origin/main` (tip `6fb8193`) is 4 commits ahead of local `main` (`a3fa8f8`); PRs #31–#34
already replaced the calibration-first Fashion with a **railing-first** Fashion that absorbs
Sextant + Gnomon and retires `0.30/0.67`. So the canonical Fashion *is* railing-first. Two
sub-issues remain:

**(a) Local working tree is stale.** Local `main` still holds the *calibration-first*
Fashion (MRI/Elsevier, with `0.30/0.67`). It is superseded on the remote but is what a local
`cat Fashion/paper/manuscript.tex` shows. Pull `origin/main` to see the real state.

**(b) Two railing-first variants coexist on `origin/main`** — the author must pick ONE:

| | `Fashion/paper/manuscript.tex` (PR #32, 01:08) | `Fashion/paper_retool/` (PR #34 = HEAD, 01:24) |
|---|---|---|
| Title | "Pseudo-diffusion is weakly identifiable in practice…" | "…pseudo-diffusion identifiability diagnostic in IVIM MRI" |
| Venue line | "For submission to NMR in Biomedicine" | (IVIM-MRI framed; Sextant-style subtitle) |
| Figures | graphical abstract + fig1–3 | **none** (text-only, 87 KB PDF) |
| Package | source only (no built PDF in tree) | **complete**: `manuscript.tex`+`.pdf`, `numbers.tex`, 3 cover letters, `consistency.py`, `HUANG_COVERAGE_MAP.md` |
| `0.30/0.67` | retired → conditional 0.63/0.81 | retired → conditional 0.63 |
| Sextant | absorbed (cites JSON) | absorbed |

Both are railing-first and Gnomon-clean. `paper_retool/` is newer and submission-complete but
figure-less and "IVIM MRI"-titled; `paper/` is the figure-rich "weakly identifiable / NMR
Biomed" framing. **Decision for the author:** consolidate to one (recommend the NMR-Biomed
figure-rich framing as the spine + paper_retool's cover-letter/consistency package). The
built final-draft PDF currently on `origin/main` is **`Fashion/paper_retool/manuscript.pdf`**.
The route recommendation (ABSORB) holds for either variant.

## #2 — BLOCKING / GNOMON DIVERGENCE: `0.30 / 0.67` survive in the MRI Fashion

Gnomon (numeric source of truth) **DIVERGES** from Fashion's marginal D\* coverages:

| Target | Fashion (MRI) claims | Gnomon rebuilt | Tol | Result |
|---|---|---|---|---|
| Laplace D\* coverage (marginal) | **0.30** | **0.80** [0.78,0.82] | ±0.10 | **DIVERGES** (Δ ≈ 0.50) |
| MCMC-SD D\* coverage (marginal) | **0.67** | **0.90** [0.89,0.92] | ±0.10 | **DIVERGES** (Δ ≈ 0.23) |

Source: `Gnomon/results/reproduction.json` (`T3a`,`T3b`, `pass:false`), `Gnomon/VERDICT.md`.
Gnomon traces the gap to two **under-documented** choices: (1) an unstated hard high-D\*/low-
perfusion cohort, (2) a *floored* SD convention for railed voxels ("overconfident by design").
Under the **honest** Cramér–Rao convention the marginal numbers are an undramatic 0.80/0.90;
the `0.30/0.67` severity is recovered **only** by flooring (Gnomon: floored pooled Laplace
0.68, high-tercile 0.41).

- **Railing-first Fashion (both variants on `origin/main`):** already **retired** `0.30/0.67`,
  reports conditional per-tercile (0.63/0.81) under a stated honest/floored convention.
  ✅ reconciled — the divergence is fixed on the remote.
- **MRI Fashion (stale local `main` / deprecated calibration version):** still quotes them as
  headline marginal coverages **without** stating the convention. ❌ This is the prompt's
  "retracted marginal numbers surviving" check — they survive **only here**, in the superseded
  version. Status: **historical / local-only** now that `origin/main` carries the railing-first
  fix. Action: discard or archive the MRI calibration version (or split it into the separate
  NPE-CRLB-audit paper, with `0.30/0.67` recast) so it cannot be mistaken for current.
- **Sextant:** never quotes them (ruler demoted). ✅ clean.

## #3 — CONTRADICTION (status, not value): the 54.7% abdomen result is billed oppositely

Identical number, contradictory evidentiary status:
- **MRI Fashion, Supp. Fig S4:** "*N=1* … qualitative illustration … **not used to support any
  quantitative claim**; N=1 precludes population inference."
- **Sextant:** the same n=1618 / **54.7%** [52.2,57.1] is the **REPLICATES-STRONG primary
  statistical claim**.

The point estimate reconciles to Gnomon (54.2% [52.0,56.4]); only the *claimed status*
conflicts. Any unification must pick one billing. Note also Sextant's railing is computed with
Fashion's **own extracted code** (`fashion_reuse`), so Sextant's 54.7% is **not an independent
confirmation** of Fashion — the independent confirmation is Gnomon's 54.2%.

## #4 — UNKNOWN / UNDER-SPECIFIED: the "1618 high-SNR voxels" selection

The 54.7%/1618 figure (Fashion S4 and Sextant) uses an SNR cut that is **under-specified in
Fashion's prose**. Gnomon documents its own ROI as **n=1932** and flags the 1618 selection as
under-specified (`TARGETS.md` T1 note); the railing-first redraft restates the number on the
full 1932-voxel ROI plus a b0-SNR>25 subset (n=799, 58.7%) "so the load-bearing number carries
its selection sensitivity." **Unknown:** the exact SNR threshold/criterion that yields 1618.
Must be pinned in any unified methods section. Not a contradiction — a reproducibility gap.

## #5 — UNKNOWN: Sextant's own target venue

`Sextant/paper/sextant.tex` states **no target venue** (blank title block). It orbits NMR in
Biomedicine (its entire "Differentiation" section is vs Casali 2026, an NMR Biomed paper), and
the railing-first Fashion redraft — which absorbs it — targets NMR Biomed. **Unknown** whether
Sextant was ever intended as an independent submission or always as a Fashion sub-component.
This compounds #1: if Sextant has no venue and Fashion's redraft already contains it, the
ABSORB route is the natural reading.

## #6 — CONSISTENCY (not a divergence): minor CI/point differences are expected

Recorded so they are not mistaken for contradictions:
- Sextant abdomen point **54.7%** (Fashion code, n=1618) vs Gnomon **54.2%** (independent code,
  n=1932): consistent — each inside the other's CI; different code + ROI by design.
- Sextant CI upper 57.1 vs Gnomon 56.4: different n and bootstrap (5000 vs 2000); both contain
  54.7. Not a conflict.
- MCMC-quantile D\* coverage: MRI Fashion **0.94** vs Gnomon **0.90** [0.887,0.914] — within the
  ±0.05 tolerance; Gnomon scores this **REPRODUCES**. Redraft uses 0.90. Consistent.

## #7 — META: Gnomon's own scope constraint

`Gnomon/VERDICT.md` states Gnomon "**does not become a standalone paper**" — it is the clean
technical core / numeric source of truth for the retool. This supports treating Gnomon as the
shared reproducibility spine under any route, not as a third paper to unify.

---

### Lead-with summary

There is **no numeric contradiction between Fashion and Sextant on a shared value** — every
railing rate reconciles to Gnomon. The two items that *looked* blocking are now largely
resolved on `origin/main`: (a) **which Fashion is canonical** (#1) — the remote already adopted
the railing-first version; what remains is a *choice between two railing-first variants*; and
(b) the **`0.30/0.67` Gnomon divergence** (#2) — fixed on the remote, surviving only in the
stale local `main` / deprecated calibration version. The remaining real issues are
**double-counted novelty** (same primary railing claim in Sextant and the railing-first Fashion,
which ABSORB resolves) and the **status-billing contradiction** (#3) on the otherwise-reconciled
54.7%. **Bottom line:** the recommended route is already ~80% executed on `origin/main`; the
author's open decisions are *which railing-first variant ships* and *formally retiring Sextant
as a standalone paper*.
