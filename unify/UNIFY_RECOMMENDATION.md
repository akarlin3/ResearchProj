# UNIFY_RECOMMENDATION — Fashion + Sextant diagnosis (read-only pass)

**Status:** recommendation only. **No manuscript was edited, merged, or deleted.** The
route decision is the author's; a follow-up prompt executes the chosen route.

**Read this first — the recommended route is already substantially executed on `origin/main`.**
`origin/main` (tip `6fb8193`) is **4 commits ahead** of the locally checked-out `main`
(`a3fa8f8`). Those commits (PRs #31–#34, all labelled *"PROVISIONAL — do not auto-merge"*)
**already replaced the calibration-first Fashion with a railing-first Fashion that absorbs
Sextant and Gnomon and retires the `0.30/0.67` coverages.** So the ABSORB route recommended
below is mostly *done* on the remote; the open work is to pick a single final draft and
formally stand Sextant down. Two wrinkles remain (FLAGS.md #1):

1. **The local working tree is stale.** Local `main` still shows the *calibration-first*
   Fashion (MRI/Elsevier) that quotes the Gnomon-diverging `0.30/0.67`. Anyone reading
   `Fashion/paper/manuscript.tex` locally sees a superseded paper.
2. **Two railing-first variants coexist on `origin/main`.** `Fashion/paper/manuscript.tex`
   (PR #32, *"Pseudo-diffusion is weakly identifiable in practice"*, NMR Biomed, graphical
   abstract + 3 figures) and the newer `Fashion/paper_retool/` (PR #34, HEAD; *"…pseudo-
   diffusion identifiability diagnostic in IVIM MRI"*, a complete submission package with
   cover letters, `numbers.tex`, consistency gate, and a built **`manuscript.pdf`**). Both
   are railing-first and Gnomon-clean; **the author must choose which is the submission
   draft.** The built final-draft PDF currently on `origin/main` is
   `Fashion/paper_retool/manuscript.pdf`.

This report treats the **railing-first Fashion** (either variant) as "Fashion" — matching the
prompt's description — and notes (Scenario B) how the diagnosis would differ if the
deprecated calibration-first MRI version were revived.

Everything below is traced to a run/output (Gnomon `results/reproduction.json`, Sextant
`results/railing_results.json`, `Sextant/paper/numbers.tex`, and the two manuscript
sources). No number is cited that was not traced.

---

## Checkpoint 0 — Inventory & profiles

### Artifacts located

| Artifact | Path | Role |
|---|---|---|
| Fashion manuscript (MRI/Elsevier, calibration-first) | `Fashion/paper/manuscript.tex` (on `main`) | Paper — checked-out version |
| Fashion manuscript (NMR Biomed, railing-first) | `Fashion/paper/manuscript.tex` on branch `redraft/railing-first-nmrbiomed` | Paper — newer, unmerged |
| Fashion supplement / refs | `Fashion/paper/supplement.tex`, `refs.bib` | Supporting |
| Sextant manuscript | `Sextant/paper/sextant.tex` (+ `numbers.tex`, `consistency.py`) | Paper |
| Gnomon clean-room reproduction | `Gnomon/results/reproduction.json`, `VERDICT.md`, `TARGETS.md` | **Numeric source of truth** |
| Sextant railing run | `Sextant/results/railing_results.json` | Shared data |
| Shared open data | OSIPI TF2.4 (Zenodo 14605039), TCGA-LIHC (TCIA) | Shared substrate |

### Profile — Fashion (railing-first redraft, NMR in Biomedicine) — *the prompt's "Fashion"*

> *"Pseudo-diffusion is weakly identifiable in practice: a reproducible boundary-railing
> signature of D\* in open in-vivo IVIM, and the shape-aware intervals that resolve it."*

**Thesis:** weak D\* identifiability leaves a *measurable, reproducible, assumption-free
fingerprint* in real data — boundary-railing of the NLLS D\* estimate — and this induces a
*conditional* (not uniform) coverage failure that shape-aware intervals resolve. **Primary
claims:** (1) NLLS D\* rails at **54.2%** of OSIPI-abdomen homogeneous-ROI voxels
(Gnomon clean-room; prior 54.7% inside CI); (2) it **generalises** to 47.8% (full ROI),
43.7% (liver 4-b), 73.4% (liver 3-b) — *these are Sextant's numbers, cited as such*; (3)
the Gaussian under-coverage is conditional (high-D\* tercile 0.63 Laplace / 0.81 MCMC),
**retiring the marginal `0.30/0.67`**; (4) quantile intervals + a clean-room amortized flow
restore calibration/sharpness. **Data/methods:** Gnomon (numeric truth) + Sextant
(cross-cohort layer); OSIPI abdomen + TCGA-LIHC; clean-room NLLS/MCMC/flow. **Venue:** NMR
in Biomedicine (stated).

### Profile — Fashion (calibration-first, MRI/Elsevier) — *the checked-out version*

> *"Calibration and Efficiency of Uncertainty Estimates in IVIM: Quantile Intervals,
> Cross-Paradigm Comparison, and a Cramér–Rao Audit of Amortized Posteriors."*

**Thesis:** IVIM UQ — *aggregate coverage masks pointwise overconfidence*; an amortized NPE
passes marginal coverage yet is overconfident/biased on D\* below the CRLB floor; quantile
intervals and input-perturbation recover calibration. **Primary claims:** cross-paradigm
coverage table (Laplace D\* **0.30**, MCMC-SD **0.67**, MCMC-quantile **0.94**); NPE D\*
claimed SD 0.08–0.67× CRLB; held-out-b collapse to 0.03 in vivo; bias-aware (van Trees)
floor 78–92%. Railing appears **only** as Supp. Fig S4 ("**N=1** abdominal case,
**54.7%** railed, qualitative illustration, *not used for inference*"). **Data/methods:**
Fashion's own NPE/UQ campaign; OSIPI brain + abdomen + liver cohort. **Venue:** Magnetic
Resonance Imaging (Elsevier) (stated).

### Profile — Sextant

> *"Sextant: boundary-railing of conventional NLLS IVIM fits as an assumption-free
> diagnostic."*

**Thesis:** answer the "calibration ruler is *overextended*" critique by **promoting
boundary-railing to the primary claim** — a fact about the optimiser+data that needs no
ground truth, no noise-model trust — and **demoting** the calibration ruler to a scoped
secondary diagnostic. **Primary claims:** OSIPI abdomen homogeneous ROI **54.7%** railed
(CI [52.2, 57.1], n=1618, "REPLICATES-STRONG"), matching "previously reported 54.7%";
full-ROI 47.8%; TCGA-LIHC liver 43.7% (4-b) / 73.4% (3-b); rail-direction + generous-bound
controls. **Data/methods:** **reuses Fashion's exact NLLS code by AST extraction** (not an
independent re-implementation of railing) + new bootstrap/SNR/full-ROI/liver layers; OSIPI
abdomen + TCGA-LIHC. **Venue: unstated** (orbits NMR Biomed — differentiates from Casali
2026, an NMR Biomed paper). See FLAGS.md #5.

---

## Checkpoint 1 — Overlap & contradiction map

### 1A. Claim-by-claim (railing-first Fashion redraft ↔ Sextant)

| # | Claim | Fashion (railing-first) | Sextant | Classification |
|---|---|---|---|---|
| C1 | D\* rails ~54% on OSIPI abdomen | **Primary** (§3.1, Fig 1): 54.2% (Gnomon), prior 54.7% in CI | **Primary** headline: 54.7% [52.2,57.1] | **SHARED — same headline** |
| C2 | Railing generalises (full ROI 47.8%) | §3.1 / Fig 1 (cites Sextant JSON) | Primary result | **SHARED — same number, same run** |
| C3 | Liver replication 43.7% / 73.4% (TCGA-LIHC) | §3.1 / Fig 1 (cites Sextant) | Primary result | **SHARED — same number, same run** |
| C4 | Rail-direction (upper-dominant) + generous-bound control | §2.3/§3.1 ("Sextant additionally reports…") | Sextant §"When and why it rails" | **SHARED (Sextant-owned)** |
| C5 | "First systematic real-data quantification of D\* railing across cohorts" | Stated as the **binding novelty** (§1, §4) | Stated as the **bulletproof primary claim** (abstract, conclusion) | **CONTRADICTION-of-credit / double novelty** |
| C6 | Conditional coverage (high-D\* tercile 0.63/0.81), retire 0.30/0.67 | **Primary** (§3.2, Table 1) | Absent (ruler demoted, not quantified) | **DISTINCT — Fashion only** |
| C7 | Quantile-interval + amortized-flow resolution | **Primary** (§3.3, Fig 3) | Absent | **DISTINCT — Fashion only** |
| C8 | Calibration ruler is "overextended" → demote it | Implicitly agrees (conditional/convention-explicit reframe) | **Central framing** | **OVERLAPPING stance** |
| C9 | Differentiation from Casali 2026 (human-abdomen vs mouse-brain) | §1, §4 | §"Differentiation from Casali" | **SHARED — near-identical argument** |

**Reading:** the railing-first Fashion redraft **is** Sextant's railing result (C1–C5, C9)
plus a calibration-resolution layer (C6–C7) that Sextant lacks. The shared block is the
*entire primary contribution of Sextant*. Sextant's only material beyond Fashion's redraft
is the rail-direction/generous-bound detail (C4) — which the redraft already references.

### 1B. Cross-paper / cross-Gnomon number reconciliation (run-traced)

| Quantity | Gnomon (truth) | Sextant | Fashion railing-first | Fashion MRI | Verdict |
|---|---|---|---|---|---|
| OSIPI abdomen homog. railing | **54.2%** [52.0,56.4] (`reproduction.json` T1, n=1932) | 54.7% [52.2,57.1] (n=1618, `numbers.tex`) | 54.2% (prior 54.7% in CI) | 54.7% (Supp S4, "N=1, not for inference") | ✅ reconciles (54.7 ∈ CI) |
| Full abdominal ROI railing | (Sextant-sourced) | 47.81% [47.10,48.51] | 47.8% [47.1,48.5] | — | ✅ identical run |
| TCGA-LIHC 4-b railing | (Sextant-sourced) | 43.70% [43.21,44.19] | 43.7% [43.2,44.2] | — | ✅ identical run |
| TCGA-LIHC 3-b railing | (Sextant-sourced) | 73.45% [73.02,73.88] | 73.4% [73.0,73.9] | — | ✅ identical run |
| D\* cov, Laplace **SD** (marginal) | **0.80** [0.78,0.82] | not quoted | retired → 0.80 pooled / **0.63** high tercile | **0.30** | ❌ **MRI Fashion DIVERGES from Gnomon** |
| D\* cov, MCMC **SD** (marginal) | **0.90** [0.89,0.92] | not quoted | retired → 0.90 pooled / **0.81** high tercile | **0.67** | ❌ **MRI Fashion DIVERGES from Gnomon** |
| D\* cov, MCMC **quantile** | **0.90** [0.887,0.914] | not quoted | 0.90 [0.89,0.91] | 0.94 | ✅ within ±0.05 tol (Gnomon REPRODUCES) |
| Flow vs railed-NLLS (cov/ECE/sharp) | 0.979/0.069/0.112 vs 0.763/0.121/0.181 | not in Sextant | 0.98/0.069/0.11 vs 0.76/0.121/0.18 | (NPE-vs-CRLB framing) | ✅ matches `reproduction.json` T4 |

**The prompt's specific checks:**
- **Retracted marginal `0.30/0.67` surviving?** → **YES, in the checked-out MRI Fashion**
  (Abstract, Results, Table 1). **NO** in the railing-first redraft (explicitly retired) and
  **NO** in Sextant. → FLAGS.md #2.
- **Convention mismatch (honest vs floored CRLB)?** → Present and *now documented* in the
  railing-first redraft (§2.4 + Table 1 honest vs floored rows; floored Laplace high-tercile
  0.41). The MRI Fashion uses the floored/"overconfident-by-design" baseline that produced
  `0.30/0.67` **without** stating the convention — the exact item Gnomon flagged. → FLAGS.md #2.
- **Conflicting railing/coverage numbers across papers?** → **No numeric railing conflict** —
  all railing rates reconcile to Gnomon. The coverage conflict is **Fashion(MRI)-vs-Gnomon**,
  not Fashion-vs-Sextant.
- **Same result claimed novel in both (double-counting)?** → **YES** (C5): the cross-cohort
  railing quantification is the stated primary novelty of *both* Sextant and the railing-first
  redraft, on identical numbers. → the decisive factor below.

### 1C. The one genuine cross-paper contradiction (status, not value)

The identical **54.7% abdomen result** is asserted with **opposite evidentiary status**:
- **MRI Fashion (Supp S4):** "*N=1* … qualitative illustration … **not used to support any
  quantitative claim**."
- **Sextant:** the same n=1618 / 54.7% is the **REPLICATES-STRONG primary statistical claim**
  with a bootstrap CI.

Same number, contradictory billing. → FLAGS.md #3.

---

## Checkpoint 2 — Relationship diagnosis

**Classification depends on which Fashion is canonical.**

### Scenario A — Fashion = railing-first redraft (the prompt's description): **(a) SAME-PAPER-REDUNDANT**

Evidence:
- The redraft's header literally names Sextant as its "cross-cohort generalization" source
  and pulls Sextant's JSON verbatim (C2–C4; numbers identical to 4 sig figs).
- The redraft's primary novelty sentence (C5) **is** Sextant's primary novelty sentence.
- Sextant's railing computation is **not independent** of Fashion (it AST-extracts Fashion's
  exact `fit_biexp_nlls`); the only independent confirmation is Gnomon, which the redraft
  already uses as truth. So Sextant adds *no independent railing evidence* the redraft lacks.
- Combined scope: one coherent paper (railing phenomenon → conditional coverage → resolution)
  already exists as the redraft; Sextant is a strict subset of its §3.1 + Fig 1.

Novelty-per-paper if both ship: Sextant ≈ 0 *additional* novel claim over the redraft. This
is textbook **salami-slicing of one result**.

### Scenario B — Fashion = calibration-first MRI version: **(b) COMPLEMENTARY-COMPANION, with caveats**

If the MRI calibration paper is canonical, Fashion's thesis (CRLB-floor audit of amortized
posteriors; aggregate-masks-pointwise) and Sextant's thesis (assumption-free railing) are
*distinct contributions on a shared spine* (D\* weak identifiability, shared OSIPI data).
**But** two frictions remain: (i) the 54.7% railing result is still double-published
(Fashion S4 vs Sextant headline, C1/C3); (ii) Sextant explicitly argues Fashion's calibration
ruler is "overextended," i.e. it *undercuts its own companion's central method* — awkward for
a same-author pair. This is the weaker, higher-risk reading.

### Why not (c) distinct/weakly-related

Rejected: the two share the same data, the same 54.7% number, the same Casali differentiation,
and (in the redraft) the same primary novelty sentence. They are not weakly related under
either scenario.

---

## Checkpoint 3 — Route recommendation (the deliverable)

### Recommended route: **ABSORB Sextant into Fashion (railing-first redraft), with Gnomon as the shared numeric core**

Sextant becomes Fashion's **cross-cohort generalization layer** (§3.1 + Figure 1) and a
reproducibility appendix — **not a standalone submission**. This is *already ~80% executed*
in `redraft/railing-first-nmrbiomed`: the redraft cites Sextant's JSON, reproduces its
numbers, and frames railing-first. The remaining work is to formally fold Sextant's
rail-direction/generous-bound detail and reproducibility harness into Fashion and stand Sextant
down as a paper.

**Runner-up:** **COMPANION PAIR** — viable **only** under Scenario B (Fashion stays the
calibration-first MRI/Elsevier paper; Sextant carries railing for an NMR-Biomed-adjacent
venue), with a strict shared claims-ledger so the 54.7% is *owned by exactly one* and cited by
the other, and Sextant's "overextended" framing softened so it does not attack its companion.

### Rubric (1–5; higher = better for that route)

| Criterion | ABSORB (recommended) | COMPANION PAIR (runner-up) | MERGE-as-equals |
|---|---|---|---|
| Scope coherence | **5** — one phenomenon → consequence → fix | 3 — two theses, shared data | 2 — calibration + railing strained into one |
| Salami-slicing risk (lower=better) | **5** (none — one paper) | 1 — same 54.7%/liver result in both | 4 |
| Dilution / Frankenstein risk (lower=better) | **4** — railing-first spine already coherent | 4 | 1 — CRLB-audit + railing dilute each other |
| Venue strategy | **5** — railing-first fits NMR Biomed; Gnomon/Sextant as reproducibility assets | 3 — two submissions, partial overlap to disclose | 2 — no single venue fits both theses |
| Reviewer perception | **5** — honest, reproducible, single claim | 2 — same-author self-overlap + self-undercut | 2 — "two papers stapled" |
| Fit with Fashion's railing-first NMR-Biomed positioning | **5** — *is* that positioning | 2 — only if Fashion reverts to MRI calibration | 3 |
| **Mean** | **≈4.8** | **≈2.5** | **≈2.3** |

### Decisive factor

**The cross-cohort D\* railing quantification is claimed as the primary novelty of *both*
Sextant and the railing-first Fashion redraft, on the *same numbers from the same run* (the
redraft cites Sextant's `railing_results.json` as its source).** Two papers cannot both
headline one result. Because the redraft already contains Sextant's contribution and adds the
calibration resolution Sextant lacks, **absorbing Sextant into Fashion is the only route that
avoids salami-slicing while losing nothing.** A companion pair is defensible solely if Fashion
reverts to the calibration-first MRI paper — which contradicts the prompt's stated
positioning.

**Willing-to-say-separate check:** they should **not** stay two papers. Under the prompt's
railing-first framing they are one paper; only a venue/positioning reversal (Scenario B) makes
two defensible, and even then with a hard claims-ledger.

---

## Checkpoint 4 — Execution sketch (plan only; no edits performed)

### For the recommended route (ABSORB)

- **What merges:** Sextant §"Results" (full-ROI 47.8%, liver 43.7%/73.4%, rail-direction,
  generous-bound) → Fashion §3.1 + Figure 1 (already drafted in the redraft). Sextant's
  pre-registered replication thresholds + `reproduce.sh` → Fashion's Data/Code Availability and
  a reproducibility appendix. Sextant's Casali differentiation → already merged (Fashion §1/§4).
- **What aligns (notation/methods):** single rail definition (`|x̂−bound|/(upper−lower) <
  rail_tol`, `rail_tol=1e-3` primary / `1e-2` sensitivity) — already shared. One bootstrap
  convention to state side by side (Gnomon n=2000 seed 20260621 for the reframe; Sextant n=5000
  seed 20260613) — disclose both seeds. One ROI/SNR-cut definition (resolve the 1618 vs 1932
  selection; FLAGS.md #4).
- **Shared repro/figures:** Gnomon = numeric source of truth for every figure (already declared
  in the redraft header). Sextant = cross-cohort layer only. Figure 1 = Gnomon abdomen point +
  Sextant cohort points (already drafted).
- **Cross-references:** none external once absorbed — Sextant ceases to be a citable paper and
  becomes an internal module; cite Gnomon + Sextant as code/data artifacts (Zenodo/repo), not as
  a companion manuscript.
- **Venue plan:** NMR in Biomedicine (railing-first redraft already targets it; Casali, the
  natural related work, is an NMR Biomed paper).
- **Authorship:** single author (Avery Karlin) throughout — no multi-author reconciliation
  needed; CRediT already identical across all three.
- **Claims-ledger reconciliation (so nothing is double-claimed or contradicted):**
  1. Retire `0.30/0.67` as headline coverages everywhere (done in redraft; **must also be done
     if the MRI version is ever revived**). Replace with the conditional per-tercile table
     under the *stated* honest/floored convention.
  2. Promote 54.7%/54.2% to a single owned claim (Gnomon point estimate, prior 54.7% noted in
     CI); drop the contradictory "N=1, not for inference" billing.
  3. State the honest-vs-floored SD convention explicitly (Gnomon's flagged item).
  4. The amortized-NPE-vs-CRLB-floor audit (the MRI paper's core) is **explicitly scoped OUT**
     of the railing-first paper (redraft §4 already says so) → it is the natural **separate**
     second paper (see below), which *is* a legitimate companion split because it shares no
     primary number with the railing paper.

### For the runner-up (COMPANION PAIR, Scenario B only)

- Keep Fashion = calibration-first MRI/Elsevier (CRLB/NPE audit); keep Sextant = railing for an
  NMR-Biomed-adjacent venue. **Required before either ships:** (i) fix the Gnomon divergence —
  remove/recast `0.30/0.67` in Fashion (FLAGS.md #2); (ii) a shared claims-ledger assigning the
  54.7% railing result to **exactly one** paper, cited by the other; (iii) soften Sextant's
  "overextended" attack on the calibration ruler so it does not undercut its companion;
  (iv) cross-reference each in the other's intro. Salami-slicing risk stays elevated because the
  railing number still appears in both.

### Note on the genuinely-separable second paper

The cleanest *non-redundant* companion split is **not** Fashion-vs-Sextant but
**railing/conditional-coverage paper (absorb Sextant)** vs **amortized-NPE CRLB-floor-audit
paper** (the MRI version's distinctive core: claimed-vs-achieved SD below the floor, held-out-b
0.03 collapse, OOD gate, bias-aware van-Trees floor). These two share *no* primary number and
are a defensible companion pair — whereas Fashion-railing-first vs Sextant share *every*
primary number.

---

## Halt

Deliverables emitted: this report + `FLAGS.md`. **No manuscript source was edited.** The route
decision is the author's. Note that the recommended ABSORB route is **already ~80% executed on
`origin/main`** (PRs #31–#34): the railing-first Fashion there absorbs Sextant + Gnomon and
retires `0.30/0.67`. The author's remaining open decisions are (1) **which railing-first variant
ships** — `Fashion/paper/` (figure-rich, NMR Biomed) vs `Fashion/paper_retool/` (newer,
submission-complete, IVIM-MRI titled; FLAGS.md #1) — and (2) **formally retiring Sextant** as a
standalone paper. The built final-draft PDF on `origin/main` is
`Fashion/paper_retool/manuscript.pdf`.
