# Atlas discovery pilot — bounded end-to-end run

**Run date:** 2026-06-13
**Scope:** a BOUNDED pilot of the inverted Atlas pipeline (retrieve → fetch →
biome-tag → screen at the widened operating point). NOT the systematic sweep.
Reproduce: `PYTHONPATH=src python -m proteus.atlas` then `… -m proteus.atlas_screen`.

## What ran (and the honest discovery mode)

The intended retrieval — Foldseek the α/β-hydrolase fold-class query against the
hosted Atlas search — **could not run: the hosted search is DOWN (HTTP 503)** on
every endpoint (`searchStructure`/`searchSequence` `…/ticket`, poll, result), confirmed
live (see `envlog/atlas-endpoints.md` for the 400-vs-503 proof that it's a backend
outage, not a usage error). This is the API-drift case the design anticipates.

Per the fallback policy, and because the multi-GB local DB is **not** downloaded in
this run, the pilot drew an **UNSEEDED 60-accession slice of the HQ clust30
representative set** (the Foldseek DB's own `highquality_clust30.lookup`; seed 1729)
so the fetch→screen machinery runs end-to-end on **real Atlas structures**. This slice
is **not fold-class targeted** — it samples the whole HQ Atlas, so the S4/S5 survival
below is the *unfiltered Atlas base rate*, a conservative floor (see Economics).

## Funnel counts

| Stage | Count | Rate |
|---|---:|---|
| Atlas structures fetched (pilot slice) | **60** | — (all pLDDT > 0.7; mean 0.882, 0–1 scale) |
| → triad-positive (S4 geometry: Ser-His-Asp + oxyanion hole) | **2** | 3.3% of fetched |
| → catalytic pocket (S5, within 12 Å of Ser OG) | **2** | 100% of triad |
| → **above the widened line** (composite ≥ −1.154) | **1** | 1.7% of fetched; 50% of triad |
| → marine subset of the above-line hits | **0** | biome unresolved for all (non-blocking) |

Anchor: IsPETase + LCC-WT, `percentile` mode; controls separated (margin 1.410).
**Widened operating point = −1.154** (production line = −0.296), the line that keeps
all three held-out divergent positives (PET46, Cut190, TfCut2) at precision 1.0,
divergent-recall 1.0 — `envlog/validation-run.md`'s −1.156 line, reproduced.

## Top candidates (screened at the widened line)

| rank | accession | composite | vs widened (−1.154) | vs production (−0.296) | cat. triad (Ser/His/Asp) | pLDDT | len | biome | marine |
|---:|---|---:|:--:|:--:|---|---:|---:|---|:--:|
| 1 | **MGYP000085474821** | **−0.374** | above | below | Ser286 / His308 / Asp288 | 0.819 | 818 | unknown | no |
| — | MGYP000748942012 | −2.720 | below | below | Ser709 / His736 / Asp499 | 0.959 | 759 | unknown | no |

The single above-line candidate (**MGYP000085474821**, composite −0.374) sits *above
the widened line and just below the production line* — squarely in the divergent-
positive band the validation run measured for TfCut2 (−0.12) and the *Thermomonospora*
polyester hydrolase (−0.14). It is a genuine cleft-PETase-like candidate surfaced from
an **unfiltered** Atlas slice with no homology gate. The second triad-bearing structure
is correctly pushed well below the line (−2.720).

**Marine subset:** none. Biome did not resolve for any accession (the Atlas API exposes
no per-accession biome; MGnify returned no mapping for raw representatives). Tagging was
non-blocking — no hit was dropped for a missing biome. Marine enrichment needs the
per-MGYP biome from the metadata sqlite or a working MGnify mapping (future work).

## Economics for the full sweep

Measured pilot rates (on the **unfiltered** Atlas slice — a *floor*, since a real
fold-class search enriches the input):
- triad (S4) survival: **3.3%** (2/60); above-widened-line: **1.7%** (1/60); above-line
  given triad: **50%** (1/2).
- per-structure screen cost: **S4 ≈ 0.10 s/structure** (geometry only); fpocket (S5)
  ≈ 1–3 s but **only fires on triad-positives** (~3.3% unfiltered). One-time control-
  anchor build ≈ 10–15 s.
- **HQ clust30 size ≈ 42.1 M representatives** (estimated: 1.03 GB `highquality_clust30.lookup`
  ÷ ~24.5 B/accession). Full structures ≈ 1 TB.

Projected full local-DB run (`backend: local_db`):
1. **Narrow first with Foldseek, never fpocket the whole DB.** Screening 42 M structures
   structure-by-structure is infeasible (S4 alone ≈ 42 M × 0.1 s ≈ 48 CPU-days; the
   triad-positive fpocket tail adds weeks). The point of `local_db` is to let the
   **Foldseek fold-class search** narrow 42 M → an α/β-hydrolase shortlist (hours, one
   box), and screen **only the shortlist**.
2. **Candidate yield** ≈ `N_shortlist × triad_survival_foldclass × 0.5` (above-line | triad,
   from the pilot). The pilot's 3.3% triad rate is the **unenriched floor**; a fold-class
   Foldseek shortlist (the query *is* the fold) should carry the triad at a much higher
   rate — that enrichment factor is the **key unmeasured quantity** (the hosted search was
   down). At `hit_cap = 2000`, even the floor 1.7% above-line ⇒ ~34 candidates; an enriched
   shortlist would yield proportionally more.
3. **Compute:** one-time ~1 TB HQ Foldseek DB download; Foldseek search of the 9-structure
   query vs 42 M ≈ hours on a multicore box; screening a 10⁴–10⁵ shortlist ≈ minutes–hours
   (S4 cheap; fpocket on the triad-positive tail). Tractable on the existing local/GCE split.

## Reproducibility record (the API drifts — re-resolve before scaling)

- **Atlas version:** `v2023_02` / HQ subset, search DB token `highquality_clust30`.
- **Query set:** `foldclass_union`, **9 structures**, hash **`c5c836fd098b212a`** —
  IsPETase(6EQE), LCC_WT(4EB0), PET46(8B4U), Cut190(4WFI), TfCut2(4CG1), CalB(1TCA),
  AChE_Tc(1EA5), CRL(1CRL), Est2(1EVQ). Trap LCC_ICCG(6THS) excluded. Broad & unseeded.
- **Foldseek params (intended):** mode `3diaa`, `min_bits ≥ 50`, `hit_cap 2000`.
- **Endpoints:** resolved in `envlog/atlas-endpoints.md` (search `POST
  api.esmatlas.com/searchStructure/ticket` — **503 this run**; fetch `GET
  api.esmatlas.com/fetchPredictedStructure/{MGYP}.pdb` — **live**).
- **Discovery mode this run:** `representative_fallback` (hosted search down) — 60 of
  21 402 accessions from the first 524 288 B of `highquality_clust30.lookup`, seed 1729.
- **Operating point:** widened, threshold −1.154 (fpocket jitter ±0.01 vs the −1.156
  validation line). Anchor IsPETase+LCC, percentile mode.
- **Outputs:** `data/interim/atlas_hits.json`, `data/processed/atlas_candidates.{csv,json}`,
  `structures/atlas_hits/*.pdb` (60).

## Honesty note

This is a **pilot on a capped slice**, and — because the hosted Foldseek search was down
— on an **unfiltered representative sample**, *not* a fold-class-targeted search and *not*
the systematic sweep. Consequences, stated plainly:
- The **fold-class enrichment was not exercised** (search 503). The 3.3% triad / 1.7%
  above-line rates are the *unfiltered Atlas base rate*; the real sweep's rates (and the
  enrichment factor) remain unmeasured.
- **Precision against the true metagenomic negative space is unmeasured.** The curated-
  control precision (1.0) does **not** transfer to a corpus this diverse — the Atlas
  contains vast non-hydrolase fold space and triad-bearing non-PET serine hydrolases that
  the 6-control panel never saw. The single above-line candidate is a *lead*, not a
  validated PET hydrolase.
- **Biome/marine tagging is effectively absent** here (all `unknown`); the marine angle
  needs the per-MGYP metadata source wired before it carries weight.

## Scale decision (teed up)

The fetch→screen path and the widened-line scoring are validated end-to-end on real Atlas
structures. The blocker is the hosted search (down) and the unexercised fold-class
enrichment. **Next:** download the HQ clust30 Foldseek DB, flip `atlas.backend: local_db`,
and re-run — the local Foldseek search both removes the dependency on the flaky hosted
endpoint and finally measures the fold-class enrichment the pilot could not.
