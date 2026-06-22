# ASSUMPTIONS.md — Augur's load-bearing assumption

> **Augur is a speculative synthesis.** It assembles the end-stage perspective of the IVIM-UQ
> program *before any of its anchors have published*, under the explicit assumption that they
> survive to publication as currently built. Every anchor is pinned below; every synthesized
> claim that depends on an anchor is **PROVISIONAL**. The submission block (`SUBMISSION_BLOCK.md`,
> enforced by `check_anchors.py`) stays engaged until the load-bearing anchors publish.

Last audited: 2026-06-22, against the `researchProj` monorepo working tree (GitHub
`akarlin3/ResearchProj`), `main` @ the CP5-retool state (latest merged PR #53).

---

## 0. What is solid regardless, and what is assumption-dependent

| | Depends on an unpublished anchor? | Status |
|---|---|---|
| **Minos theory half** — Plumbline Thm 1–2 + Prop. 3 (`Minos/theory`, `Minos/minos-core`) | No — self-contained, machine-verified | **SOLID** |
| The **arc** as a logical argument (trust→VoI→action) | No — it is an argument, not a measurement | SOLID (as argument) |
| External cross-modal `D*`–`Ktrans` citations (Sun 2019; Yang 2019) | No — published, verified | **SOLID** (see `CITATIONS.md`) |
| **§1 Trust** claims (Fashion calibration behaviour) | **Yes** — Fashion in review | **PROVISIONAL** |
| **§2 VoI** applied claims (Minos applied half) | **Yes** — consumes Fashion + Gauge | **PROVISIONAL** (theory half SOLID) |
| **§3 Action** claims (Lethe scale verdict; Gauge wall/monitor) | **Yes** — Lethe + Gauge unpublished | **PROVISIONAL** |
| **Whole-paper submission** | **Yes** | **BLOCKED** until Fashion + Minos publish |

---

## 1. Anchors (pinned; PROVISIONAL until each publishes)

| anchor | role in the arc | pinned status (no manuscript DOI) | the assumption Augur relies on |
|---|---|---|---|
| **Fashion** | Trust — the ruler | **in review, *NMR in Biomedicine*** (retooled, boundary-railing-first; resubmitted from MRM). Code archive Zenodo `10.5281/zenodo.20649669`. | The skew-aware posterior restores *marginal* coverage with a residual **high-`D*` conditional** gap (honest CRLB) — survives review. |
| **Minos** | VoI — the decision | **PROVISIONAL.** Theory half (Plumbline Thm 1–2 + Prop. 3) machine-verified & SOLID; applied half consumes Fashion + Gauge. No DOI. | The decision–calibration gap, the `O(γ²)` VoI law, and the label-free floor hold (theory: independent of publication). |
| **Gauge** | Action — the identifiability wall | **target *MRM*, manuscript assembled** (internal consistency gate PASS, 34/34 numbers trace). No DOI. | The high-`D*` identifiability wall is a real CRLB limit; `D*` test–retest `r=−0.17` (CI [−0.39,0.05]) is reproducible. |
| **Lethe** | Action — the wrong-size limit | **verdict rendered** (constrained-validation; Echo portion; real ACRIN-6698 n≈76). Re-pointed to NMRB. No DOI. | The conformal `D` interval is ~4× too narrow for real test–retest (coverage 0.263 vs 0.755 target) — survives review. |

**Stale-metadata finding (logged, not fixed here):** `Fashion/CITATION.cff` still reads *"in
submission to MRM"* with the pre-retool title; the authoritative status (root `README.md`,
`Minos/future/ASSUMPTIONS.md`) is the retooled NMR-in-Biomedicine submission. Augur pins the
authoritative status. (Fix belongs in a Fashion PR, not Augur.)

---

## 2. Data source — clean (the IP gate)

Augur runs **no data**. It cites in-repo results (synthetic / open ACRIN-6698, CC-BY-4.0) and
published literature only. **No `pancData3`, no MSK, no private clinical data** is touched, imported,
or referenced. IP gate: **PASS** by construction.

---

## 3. Re-validation contract (one command)

When an anchor publishes (or revises):

1. Update its row above with the published DOI / version.
2. Set the matching `published=True` (with the DOI) in `check_anchors.py`'s `ANCHORS` table.
3. Run `bash reproduce.sh`. It (a) re-checks each anchor's status, (b) re-runs the submission-block
   gate, (c) prints what is now clearable. The paper unblocks **only** when Fashion **and** Minos are
   both published (Lethe strongly recommended); until then `check_anchors.py` exits non-zero and the
   block stands.

No PROVISIONAL flag may be cleared except through this contract.
