# Atlas sweep — CP4/CP5 report: tiered screen, headline, provenance

**Run date:** 2026-06-14 · **Host:** GCE n2-standard-16 (created + torn down this run).
DB: ESMAtlas highquality_clust30 prebuilt Foldseek DB (36,986,627 structures).
Pipeline: foldseek fold-class search → Foldcomp range-fetch → **untouched** S4 geometry
→ S5 cleft → control-anchored composite at the **widened operating point**.

## Headline — enriched above-line rate vs the 1.7 % random floor

Tier-1 = the **top 300 enriched hits by Foldseek score**, screened at the widened line
(composite ≥ **−1.1587**; calibrated on IsPETase/LCC, recovers PET46/Cut190/TfCut2,
precision 1.0 on the controls).

| funnel stage | n |
|---|---|
| enrichment shortlist (CP3) | 217,833 |
| fetched (tier-1) | 300 |
| screened | 300 |
| triad-positive (S4) | 298 |
| pocket-ok (S5) | 298 |
| **above the widened line (S5 composite)** | **23** |

> **Enriched above-line rate = 23/300 = 7.7 %  vs  1.7 % random floor  →  ≈ 4.5× lift.**
> 95 % CI on 23/300 ≈ [5.2 %, 11.3 %], which excludes the 1.7 % floor. **Enrichment
> lifted the signal.** This is the first real read that the fold-class shortlist
> concentrates PETase-like active-site geometry above the metagenomic base rate.

## Important caveat — tier composition (read before over-claiming)

The 23 above-line hits are **21 from the AChE(1EA5) query and 2 from CRL(1CRL); zero
from the PETase/cutinase queries.** Two compounding reasons:

1. **Global-bits tiering crowds out the PETase anchors.** Ranking the union by absolute
   Foldseek bits fills tier-1 with the highest-scoring fold-matches (AChE- and
   PET46-anchored hits reach bits 1800–2810), while the best IsPETase(6EQE)-anchored
   hits top out ~1511 bits and **don't enter the top 300 at all.**
2. **The exposure discriminator is doing its job, inverted from naïve expectation.** The
   composite is computed on the *hit's own* geometry, not the query's. The very
   highest-bits hits are mostly lid-bearing / buried-Ser (lipase-like) and fall *below*
   the line (composite −2 to −3.5); the 23 that clear it are metagenomic α/β-hydrolases
   with **open, exposed catalytic serines** — PETase-like active-site exposure — that
   happened to fold-match the AChE net.

So 7.7 % is a real lift, but tier-1's leads are **exposed-active-site α/β-hydrolases,
not PETase sequence-homologs.** That is consistent with the *unseeded* design (we never
seed on PETase) — but a per-query / score-normalized tiering would surface the
PETase-anchored branch, which this global-bits tier-1 misses. Teed up below.

## Top candidates (above the line, ranked by composite)

`mean_plddt` on the Atlas V0 0–1 scale; biome `unknown` (no metagenome map — see honesty
note); all marine=False (unresolved). Full table: `envlog/atlas_candidates.csv`.

| rank | accession | composite | foldseek_bits | pLDDT | src query |
|---|---|---|---|---|---|
| 1 | **MGYP000470279205** | **+4.292** | 1808 | 0.947 | 1EA5 |
| 2 | MGYP000081732636 | −0.024 | 1780 | 0.969 | 1EA5 |
| 3 | MGYP001597021366 | −0.364 | 1776 | 0.954 | 1EA5 |
| 4 | MGYP000110442245 | −0.401 | 1794 | 0.932 | 1EA5 |
| 5 | MGYP001444862274 | −0.411 | 1790 | 0.928 | 1CRL |
| 6 | MGYP001307280170 | −0.659 | 1801 | 0.951 | 1EA5 |
| 7 | MGYP000942774382 | −0.765 | 1774 | 0.917 | 1EA5 |
| 8 | MGYP001040063565 | −0.817 | 1866 | 0.952 | 1EA5 |
| 9 | MGYP002291380453 | −0.856 | 1795 | 0.908 | 1EA5 |
| 10 | MGYP003342763746 | −0.908 | 1796 | 0.928 | 1EA5 |

`MGYP000470279205` is the standout: composite **+4.29** (far above the line and above
most positive controls), pLDDT 0.95 — a high-confidence, strongly exposed active site.
23 above-line total (the remaining 13 sit between −0.91 and −1.14).

## Recovery check (positive control on real data)

No curated list of **known PETase MGYP accessions** exists to test recovery directly
(the references.csv recovery entries GuaPA/MG8 have no resolved Atlas accession). As a
proxy: the highest-identity hits to the IsPETase(6EQE) query (fident 0.40–0.65, bits
≤ 1511) are genuine metagenomic PETase-fold homologs, but **none entered tier-1** (their
absolute bits fall below the top-300 global-bits cutoff) — so they were neither fetched
nor screened this run. Their screen verdict is a **direct follow-up** for the per-query
tier. No false reassurance is claimed here.

## Provenance (reproducibility record)

- Atlas: `highquality_clust30` (v0), prebuilt Foldseek DB, 36,986,627 entries.
  `_ss` sha256 `1127e5e3…eca84`; base sha256 `8d903b13…7de92`.
- foldseek `718d42176d2f67d36a60866fedfb881f8d5a7ebf`; foldcomp `1.0.0`; fpocket (bioconda).
- Query set hash `c5c836fd098b212a` (9 control structures, unseeded fold-class union).
- Search: `-s 9.5 -e 0.01 --max-seqs 300000 --alignment-type 2` (no `_ca`).
- Screen: widened operating point, threshold **−1.1587** (production −0.2989);
  `screen_budget` = 300.
- Artifacts → `gs://projproteus-fold/atlas-sweep/2026-06-14/`: `provenance.json`,
  `economics.json`, `funnel.json`, `hits_ranked.tsv.gz` (217,833), `result.m8.gz`
  (1.08 M rows), `atlas_candidates.csv`/`.json`.

## Honesty note (carry it)

Candidates are **prioritized, not verified** — no wet-lab. Precision against the true
metagenomic negative space is still **unmeasured**: the 1.7 % floor is the pilot's
random-Atlas base rate, so the 4.5× lift says the fold-class+geometry shortlist
concentrates PETase-like *geometry*, not that any hit is a PET hydrolase. Biome is
**unresolved** (raw ESM Atlas carries no metagenome map here) — the marine/ocean framing
stays aspirational until cross-referenced to AFESM/MGnify. Individual hits are leads for
follow-up, not hits.

## Decision teed up

The lift is real (7.7 % vs 1.7 %), so deeper work is justified — but **not** simply more
of the same global-bits tier. Recommended next, in order:

1. **Per-query / score-normalized tiering** — screen the top-N of *each* query's hits
   (esp. the IsPETase/LCC/Cut190 branches that tier-1 skipped), so the headline reflects
   the PETase-anchored shortlist, not the AChE-dominated high-bits tail. Cheap: reuse
   `hits_ranked.tsv` (in GCS) + Foldcomp-fetch + screen; runs on the Mac or a tiny VM —
   no re-search, no GCE build.
2. **Characterize `MGYP000470279205`** and the other high-composite leads (manual
   inspection, optional docking via the wired Vina P4).
3. **Biome resolution** — cross-reference surfaced MGYPs to AFESM clusters / MGnify to
   test the marine hypothesis on the actual leads.
