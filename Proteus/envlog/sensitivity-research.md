# Sensitivity re-search — absence vs. unreachability of the divergent tail

**Run date:** 2026-06-14  
**Search:** GCE n2-standard-16, prebuilt ESMAtlas `highquality_clust30` Foldseek DB (36,986,627 structures), `gce/atlas_sensitivity_search.sh`.  
**Downstream (this file):** local Mac `proteus` env — partition, identity binning, fetch (`fetchPredictedStructure`), screen, seqid. Reuses `screen`/S4/S5 + the seqid path UNTOUCHED at the pinned line **-1.1587**.

## TL;DR

**UNREACHABILITY / CO-FAILURE.** Max-sensitivity structure search retrieves the sub-25% tail in bulk (9,566,468 neighbours; 619 triad-bearing among the screened best), but they clear the pinned line at **11.3%** — roughly a quarter of the 42.9% random floor — and the genuinely divergent (<20%) ones clear it **not at all (0/22)**. **Structure reaches them; geometry cannot discriminate them** — sequence and structure fail together. Demonstrated, not inferred.

## Checkpoint 1 — retrieval at max sensitivity (the threshold test)

Relaxed search params vs the original sweep: `-s 9.5 --max-seqs 4,000,000 -e 10000` (min_bits≈0) vs `-s 9.5 --max-seqs 300000 -e 0.01 min_bits 50`. **20,056,560** PETASE alignment rows → **9,810,546** unique Atlas targets.

**Unique targets by max structural identity (Foldseek fident) to any PETase query** — the bands Figure 1 found empty:

| identity band | this re-search | original (e≤0.01) | NEW (relaxed-only) |
|---|---:|---:|---:|
| <20% | 9,050,775 | 175,252 | 8,879,402 |
| 20-25% | 515,693 | 23,981 | 488,362 |
| 25-30% | 168,245 | 2,555 | 165,225 |
| 30-40% | 63,673 | 450 | 63,161 |
| 40-60% | 11,716 | 56 | 11,658 |
| >60% | 444 | 12 | 432 |

**9,608,240** targets are NEW vs the original e≤0.01 sweep. Sub-25% retrieval (the previously-empty region): <20% = 9,050,775, 20-25% = 515,693.

## Checkpoint 2 — screen the divergent band at the pinned line

Of **9,566,468** targets below 25% structural identity, **2,100** were selected **stratified by fident** (top-N by structural bits within each band, so the deep tail is screened by its strongest — most structurally-convincing — members), fetched and run through the UNTOUCHED S4 → fpocket(triad) → S5 path, scored at **-1.1587**. Sequence identity recomputed via the seqid path (biotite SW, BLOSUM62, coverage floor 0.50).

Selection strata (fident band → available → screened):

- 0-12% fident: 5,217,272 available → top 700 by bits
- 12-18% fident: 3,464,698 available → top 700 by bits
- 18-25% fident: 884,498 available → top 700 by bits

**Funnel (screened sub-cut set):**

| stage | n |
|---|---:|
| screened | 2,100 |
| fetched | 2,100 |
| triad+ (S4) | 1,710 |
| pocket-ok (S5) | 1,702 |
| above line (-1.1587) | 198 |

Of the 198 above-line, **151** have a credible (≥50% query-coverage) sequence homolog to a PETase query (binned by identity below); the other **47** are coverage-floor rejects (median coverage 0.075) — a short high-identity patch over a few % of the query, **not** a credible homolog and not an interpretable identity, so they are excluded from the identity bands (Figure 1's 0.50 coverage guard).

**Genuine divergent tail (Figure 1's primary biotite-SW metric, credible coverage):**

| band | fetched | triad+ | above line |
|---|---:|---:|---:|
| <25% identity | 724 | 619 | 70 |
| <20% identity | 29 | 22 | 0 |

**Above-line│triad by identity band vs the 42.9% floor (the strengthened Figure 1 tail):**

| seq-id bin | n triad | above | rate | vs floor 42.9% |
|---|---:|---:|---:|---|
| <20% | 22 | 0 | 0.0% | p=0.000, RR 0.00× |
| 20-25% | 597 | 70 | 11.7% | p=0.000, RR 0.27× |
| 25-30% | 661 | 78 | 11.8% | p=0.000, RR 0.28× |
| 30-40% | 45 | 3 | 6.7% | p=0.001, RR 0.16× |
| 40-60% | 1 | 0 | 0.0% | p=1.000, RR 0.00× |
| >60% | 0 | 0 | — | *empty* |

## Checkpoint 3 — verdict

**UNREACHABILITY / CO-FAILURE.** The divergent tail is **reachable in bulk** — max-sensitivity structure search retrieves 9,566,468 PETase structural neighbours below 25% identity, and 619 of the 724 screened (credible-coverage) sub-25% hits bear a Ser-His-Asp triad. But the cleft/exposure discriminator **cannot pick PETase-like geometry out of them**: only **70/619 = 11.3%** clear the −1.1587 line — well below the **42.9%** random floor — and in the genuinely divergent **<20%** band, **0/22** clear it. Structure reaches the tail; geometry fails to discriminate it. **Sequence and structure fail together — demonstrated, not inferred.**

Deepest credible-coverage (≥50%) above-line hits, by sequence identity — how far down the line is crossed (none below 20%):

| accession | nearest PETase | seq-id | cov | fident | composite |
|---|---|---:|---:|---:|---:|
| MGYP003629794474 | 8B4U | 20.0% | 0.888 | 15.8% | -1.106 |
| MGYP003311028628 | 8B4U | 20.38% | 0.948 | 17.9% | -0.597 |
| MGYP003294192687 | 4CG1 | 20.65% | 0.618 | 11.9% | 0.278 |
| MGYP001716174841 | 8B4U | 20.68% | 0.829 | 16.6% | -0.866 |
| MGYP003585819079 | 6EQE | 20.86% | 0.66 | 17.9% | -1.120 |
| MGYP001815227229 | 8B4U | 21.05% | 0.963 | 16.7% | -1.051 |
| MGYP000411743521 | 8B4U | 21.27% | 0.933 | 17.6% | -1.059 |
| MGYP000159150033 | 8B4U | 21.4% | 0.777 | 15.6% | 1.622 |
| MGYP000926693233 | 4CG1 | 21.5% | 0.683 | 11.1% | 2.116 |
| MGYP000923013937 | 4WFI | 21.7% | 0.647 | 11.9% | -0.752 |
| MGYP001144797622 | 8B4U | 21.79% | 0.851 | 17.8% | 4.479 |
| MGYP001248854856 | 6EQE | 21.88% | 0.558 | 16.3% | -0.939 |
| MGYP000739964870 | 4CG1 | 21.95% | 0.573 | 11.6% | -0.802 |
| MGYP001277170917 | 8B4U | 22.22% | 0.922 | 17.5% | -1.119 |
| MGYP001428766973 | 8B4U | 22.29% | 0.591 | 11.5% | -0.641 |
| MGYP000909102703 | 8B4U | 22.31% | 0.855 | 21.3% | -0.659 |
| MGYP000067738852 | 6EQE | 22.39% | 0.698 | 17.7% | -0.978 |
| MGYP001076156845 | 4WFI | 22.39% | 0.69 | 24.6% | -1.140 |
| MGYP002796801425 | 8B4U | 22.43% | 0.944 | 21.2% | -0.654 |
| MGYP002852548985 | 6EQE | 22.43% | 0.619 | 15.6% | -0.580 |

_(The 47 low-coverage above-line hits are excluded here — their identity is over a few % of the query, not a credible homolog; including them would not change the verdict.)_

### Scope guard

**Identity measures sequence-search reach, not PET activity.** A sub-25% structural neighbour clearing the line is a geometry/exposure lead, not a verified PET hydrolase. No wet-lab; leads are prioritized, not validated.

## Reproducibility

- **Search:** `gce/atlas_sensitivity_search.sh` (foldseek `-s 9.5 --max-seqs 4000000 -e 10000 --alignment-type 2`), prebuilt `highquality_clust30` 3Di+AA DB.
- **Screen/seqid:** `scripts/sensitivity_screen.py` — reuses `proteus.screen.screen_model` + `build_control_anchor` and the biotite SW seqid method untouched; line pinned to −1.1587; floor 12/28 from `floor.json`.
- **Queries:** 6EQE (IsPETase), 4EB0 (LCC_WT), 8B4U (PET46), 4WFI (Cut190), 4CG1 (TfCut2).
- **Artifacts:** this report + `data/processed/sensitivity_per_hit.csv` / `sensitivity.json`, pushed to `gs://projproteus-fold/sensitivity/2026-06-14/`.
