# Atlas sweep — CP0–CP2 host, toolchain & DB record (GCE, prebuilt Foldseek DB)

**Run date:** 2026-06-14
**Host:** GCE `proteus-atlas` — `n2-standard-16` (16 vCPU, 62 GB RAM), Ubuntu 22.04,
`us-central1-a`, 100 GB `pd-balanced`. Created + torn down within the run (billing
halted at the end; see CP5 in `envlog/atlas-sweep.md`).

This run replaces the planned 120 GB Foldcomp-build sweep: the **prebuilt ESMAtlas30
Foldseek DB is still hosted by Meta** (the prior `atlas-localdb.md` STOP was based on a
mis-read 403 — see its correction banner). So we downloaded the finished DB directly.

---

## CP0 — host + toolchain (verified on the VM)

| component | result |
|---|---|
| arch / vCPU / RAM | x86_64 / 16 / 62 GB ✓ (≥16 vCPU; disk 95 GB free) |
| foldseek | static AVX2 binary `718d42176d2f67d36a60866fedfb881f8d5a7ebf` ✓ |
| foldcomp | pip `1.0.0` ✓ (Foldcomp-read + `decompress` validated) |
| fpocket | bioconda (Miniforge env) ✓ |
| python | 3.11 + numpy 2.4.6, pyyaml, biotite, biopython ✓ |
| proteus pkg | all modules import; screen/S4/S5/calibrate **reused untouched** ✓ |

**Foldcomp-read confirmation:** `foldseek createdb` reads a Foldcomp DB directly
(auto-detected via `.dbtype` + FCZ magic, since Foldseek 5) — so a build *would* avoid
the ~2 TB PDB expansion. Not needed here (prebuilt DB used), but verified as the
fallback.

## CP1/CP2 — DB acquisition (prebuilt, not built)

**Source (live, HTTP 200):**
`https://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/foldseekdb/`
manifest: `raw.githubusercontent.com/facebookresearch/esm/main/scripts/atlas/v0/highquality_clust30/foldseekdb.txt`

Downloaded the **3Di+AA search DB only** (skipped the 88 GB `_ca` C-alpha set, needed
only for TM-align `--alignment-type 1`):

| file | bytes | role |
|---|---|---|
| `highquality_clust30` (base AA) | 7,406,109,749 | amino-acid DB |
| `highquality_clust30_ss` (3Di) | 7,406,109,749 | structure-token DB (searched) |
| `highquality_clust30_h` (headers) | 924,665,675 | MGYP accessions |
| `.index` / `.lookup` / `.source` / `.dbtype` | ~2.4 GB | indices |
| **total downloaded** | **~20 GB** | (vs 120 GB Foldcomp + ~2 TB decompress avoided) |

- **Entries (`.index` line count): 36,986,627** (~37 M representative structures).
- `_ss` sha256: `1127e5e3f6d50243106897416a728f077cf27e8fc3fb8f0bd2dcda8f9c5eca84`
- base sha256: `8d903b13698fd6bdd3b461ee343d9d00cda1e33a575421dbd24a1946a5c7de92`
- Validated: trivial `foldseek search` returns parseable hits (the full CP3 search ran
  end-to-end). **No `_ca`** → foldseek warns and disables the structure-bit re-sort
  ("does not affect E-values"); and `convertalis` must use a **non-structure**
  `--format-output` (the `prob`/`lddt`/`tmscore` columns require `_ca`).

## CP4 structure access — Foldcomp range-fetch (no bulk download)

The prebuilt Foldseek DB has no full-atom coordinates and the Atlas
`fetchPredictedStructure` API is **403 (down)**. Full structures for screening are
pulled **per-accession from the 120 GB Foldcomp archive without downloading it**:

- The Steinegger worker `foldcomp.steineggerlab.workers.dev/highquality_clust30`
  302-redirects to `steineggerlab.s3.amazonaws.com`, which **honors HTTP Range** (206).
- Download only the Foldcomp `.index` + `.lookup` (~2 GB) → map `MGYP → (offset,len)` →
  `Range`-GET each hit's ~5 KB FCZ block → `foldcomp.decompress` → full-atom PDB
  (validated: a 269-res ESMFold model, 2284 atoms, pLDDT in B-factor on 0–1 scale).
- 300 tier-1 structures fetched this way in **6 s**.

## Reproducibility

Search inputs, full `result.m8`, ranked shortlist, and candidates are staged at
`gs://projproteus-fold/atlas-sweep/2026-06-14/`. A deeper tier needs no re-search:
take more rows from `hits_ranked.tsv.gz`, Foldcomp-fetch, screen.
