# ESM Metagenomic Atlas — resolved endpoints & provenance

> Resolved at CP0 of the Atlas-integration run. **The Atlas API drifts**, so these
> were resolved from the live site (not hardcoded guesses) and probed empirically.
> Re-resolve before any large run. Resolution method recorded per row.

**Resolution date:** 2026-06-13
**Atlas version targeted:** `v2023_02` (a.k.a. v0 high-quality subset); search DB
token `highquality_clust30` (from the site's `REACT_APP_FOLDSEEK_DB`).

## How these were resolved (no guessing)

1. Fetched the SPA bundle `https://esmatlas.com/static/js/main.b2c22278.js` and
   read the API constants + request builders.
2. Empirically probed every endpoint with `curl` and (for the search ticket flow)
   the live browser network panel.

## Endpoints

### 1. Foldseek structure search — `searchStructure` (MMseqs2/Foldseek "ticket" API)
- **Submit:** `POST https://api.esmatlas.com/searchStructure/ticket`
  - `Content-Type: application/x-www-form-urlencoded`
  - body: `q=<URL-ENCODED PDB FILE CONTENTS>&mode=3diaa&email=&database[]=highquality_clust30`
  - the **query is a structure** (PDB text). The server does the AA/struct→3Di
    conversion; **no local Foldseek / ProstT5 needed for the api backend.**
- **Poll status:** `GET https://api.esmatlas.com/searchStructure/ticket/{id}` → `{id,status}`
- **Fetch results:** `GET https://api.esmatlas.com/searchStructure/result/{id}/0`
  → MMseqs2-server result JSON (m8-style alignments: target, bits, eval, …).
- Sequence variant: `POST https://api.esmatlas.com/searchSequence/ticket`,
  body `q=%3E<name>%0A<SEQ>&mode=accept&email=&database[]=highquality_clust30`.
- **CORS** `*`; methods `OPTIONS,POST,GET`; no auth token on `…/ticket`; no
  rate-limit/`Retry-After` headers exposed.
- **OBSERVED STATUS (2026-06-13): DOWN — HTTP 503 "Service Temporarily Unavailable"**
  on `POST …/searchStructure/ticket`, `…/ticket/{id}`, `…/result/{id}/0`, and the
  sequence variant. Proof it's a backend outage and not a usage error: the correct
  segment `ticket` returns **503** (routes through API Gateway to the search
  integration, which is unavailable), whereas a wrong segment such as
  `…/searchStructure/highquality_clust30` returns **400 "POST 'x' is not supported."**
  (Lambda path-check) and bare `…/searchStructure` returns API-Gateway
  `{"message":"Missing Authentication Token"}`. → triggers the **local-DB fallback**.

### 2. Fetch predicted structure by accession — LIVE ✅
- `GET https://api.esmatlas.com/fetchPredictedStructure/{MGYP}.pdb` (or `.cif`)
- Returns an ESMFold **V0** model; per-residue **pLDDT is the B-factor on a 0–1
  scale** (HQ subset mean ≈ 0.77). `content-type: chemical/x-pdb`, CORS `*`.
- A non-existent / non-HQ accession returns **HTTP 403 "Forbidden"** (S3
  AccessDenied style), *not* 404. Verified live, e.g. `MGYP002537940442` → 200,
  175 824 bytes.

### 3. ESMFold fold-a-sequence — LIVE ✅ (not used by the pilot; recorded for completeness)
- `POST https://api.esmatlas.com/foldSequence/v1/pdb/` — body = raw AA sequence →
  returns ESMFold PDB text.

### 4. Bulk high-quality clust30 Foldseek DB (for the local_db sweep) — NOT downloaded this run
- **Foldseek DB file manifest:** `https://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/foldseekdb.txt`
  (a list of part URLs; download with `aria2c --input-file=…`). Base:
  `https://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/foldseekdb/…`
- **Tarballs of PDBs:** `…/esmatlas/v0/highquality_clust30/tarballs/`
- **Accession list (DB representatives):** `highquality_clust30.lookup`
  (served via `https://foldcomp.steineggerlab.workers.dev/highquality_clust30.lookup`).
- **Metadata sqlite (per-MGYP):** `https://dl.fbaipublicfiles.com/esmatlas/v2023_02/metadata-rc2.sqlite.gz`
  (range-readable; gzip not random-seekable — needs full download to query).
- **Foldcomp alternative:** `foldcomp.setup('highquality_clust30')`.
- HQ set ≈ **1 TB** as structures. Per the run spec, **do not download in the pilot.**

## Biome metadata
- The Atlas API exposes **no per-accession biome/source endpoint** (only
  `fetchPredictedStructure`, `fetchSequence`, `fetchConfidencePrediction`,
  `fetchMetrics`). Biome is resolved best-effort via the **MGnify API**
  (`https://www.ebi.ac.uk/metagenomics/api/v1/`), non-blocking: unresolved →
  `biome: unknown`, `marine: false`. (The authoritative per-MGYP biome lives in the
  multi-GB `metadata-rc2.sqlite.gz`, which the pilot does not download.)
