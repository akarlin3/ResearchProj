# Atlas sweep — CP3 fold-class enrichment economics

**Run date:** 2026-06-14 · **Host:** GCE n2-standard-16 · **DB:** ESMAtlas
highquality_clust30 prebuilt Foldseek DB (36,986,627 structures).

The measurement the 503'd pilot could not make: **37 M ESM Atlas structures →
α/β-hydrolase fold-class shortlist**, before any screening.

## Query (unseeded fold class)

`foldseek createdb` of the 9 control **structures** (fold-class union minus the
LCC_ICCG trap), query_set_hash `c5c836fd098b212a`:

- PETase/cutinase anchors: IsPETase(6EQE), LCC_WT(4EB0), PET46(8B4U), Cut190(4WFI), TfCut2(4CG1)
- broad α/β-hydrolase (negatives, included on purpose to span the fold): CalB(1TCA), AChE(1EA5), CRL(1CRL), Est2(1EVQ)

Search: `foldseek search -s 9.5 -e 0.01 --max-seqs 300000 --alignment-type 2`
(no `_ca`; 3Di+AA local alignment). Map target → MGYP via headers; `min_bits ≥ 50`; dedup.

## Economics

| quantity | value |
|---|---|
| Atlas structures (denominator) | 36,986,627 |
| total alignment rows | 1,081,416 |
| **unique accessions (min_bits ≥ 50)** | **217,833** |
| fraction of the Atlas | **0.589 %** |
| **search-space enrichment vs random** | **≈ 170×** (37 M → 217,833) |

`min_bits ≥ 50` removed nothing here — the `e ≤ 0.01` alignment gate already kept only
strong hits (min observed bits **192**).

**Foldseek bit-score distribution (unique hits):**

| min | p25 | median | p75 | p90 | p99 | max |
|---|---|---|---|---|---|---|
| 192 | 293 | 403 | 516 | 628 | 1398 | 2810 |

**Hits per query** (union across all 9; an accession can hit several):

| query | hits | | query | hits |
|---|---|---|---|---|
| 8B4U (PET46) | 179,228 | | 4EB0 (LCC_WT) | 136,056 |
| 4CG1 (TfCut2) | 149,645 | | 1EVQ (Est2) | 125,153 |
| 6EQE (IsPETase) | 139,478 | | 4WFI (Cut190) | 113,832 |
| 1TCA (CalB) | 86,212 | | 1CRL (CRL) | 78,358 |
| | | | 1EA5 (AChE) | 73,454 |

## Read

The fold class is common (0.59 % of the Atlas, ~218 k structures) — exactly why the
unseeded search is a **170× concentration**, not a needle-in-haystack miss. The
headline question (does this shortlist actually *enrich for PETase-like active sites*
vs a random Atlas draw?) is answered by screening — see `envlog/atlas-sweep.md`.

Shortlist `hits_ranked.tsv.gz` (217,833 rows: accession, bits, evalue, fident,
best_query) and full `result.m8.gz` at `gs://projproteus-fold/atlas-sweep/2026-06-14/`.
