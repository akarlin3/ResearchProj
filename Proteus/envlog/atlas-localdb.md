# Atlas local-DB sweep — Checkpoint 0 RESOURCE AUDIT

> ⚠️ **SUPERSEDED / CORRECTED (2026-06-14).** This audit's central premise — that
> Meta deprecated the ESM Atlas hosting and the prebuilt Foldseek DB is gone — was
> **WRONG.** The `403 AccessDenied` seen below was only an S3 **directory-listing**
> denial; the prebuilt files are alive one level down under a `foldseekdb/`
> subdirectory (`…/highquality_clust30/foldseekdb/<file>` → HTTP 200), with the
> manifest hosted on GitHub raw (`facebookresearch/esm/.../foldseekdb.txt`). The
> **~20 GB prebuilt 3Di+AA Foldseek DB** (skipping the 88 GB `_ca`) downloads
> directly — **no 120 GB Foldcomp build, no ~2 TB decompress, no terabyte disk.**
> The sweep was subsequently **completed on a cheap GCE n2-standard-16** using the
> prebuilt DB. See `envlog/atlas-foldcomp.md` (CP0–CP2), `envlog/atlas-enrichment.md`
> (CP3 economics), and `envlog/atlas-sweep.md` (CP4–CP5 results + headline).
> The disk/host reasoning below is retained only as a record of the (mistaken) STOP.

**Run date:** 2026-06-13
**Host audited:** M4 MacBook Air (Apple Silicon, arm64, 16 GB unified, MPS/CPU).
**Verdict:** 🛑 **STOP — do NOT run the local-DB sweep on this host.** Move it to a
big-disk, many-core box (GCE persistent disk ≥ 2 TB, or the Origin desktop).

This file records the audit; it is intentionally written *instead of* downloading
anything (Checkpoint 0 says enumerate with HEAD requests, then STOP if insufficient).

---

## 1. Code / tooling state

| Item | State | Note |
|---|---|---|
| `atlas.py` `local_db` backend | present, **wired but untested** | lives on branch `worktree-atlas-pilot` (commit `771c309`), **not on `main`**. `LocalDbBackend.search` is a validate-on-first-use stub that raises `BackendNotReady` until the DB + foldseek exist. |
| `atlas_screen.py` (S4→S5 widened-line screen) | present, intact | same pilot branch; reuses `proteus.screen.screen_model` + `calibrate` untouched. |
| `screen` / S4 / S5 | intact on `main` | `src/proteus/{screen,s4_geometry,s5_cleft_filter,calibrate}.py`. |
| **foldseek** | ❌ **NOT on PATH** | `which foldseek` → not found. Hard blocker for any search; installable (brew/bioconda/static) but absent now. |
| `envlog/atlas-endpoints.md` | ❌ does not exist | referenced by `atlas.py` docstring but the pilot never wrote it. |

## 2. The DB source has moved — the prebuilt Foldseek DB is GONE

Checkpoint 0 assumed `dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/foldseekdb.txt`
is a live manifest of Foldseek DB-part URLs. It is not — **Meta deprecated the ESM
Atlas bulk hosting.** HEAD/GET probes (2026-06-13):

| URL | Result |
|---|---|
| `…/highquality_clust30/` (dir) | **403 AccessDenied** (XML) |
| `…/highquality_clust30/foldseekdb.txt` | **403 AccessDenied** |
| `…/highquality_clust30/highquality_clust30.lookup` | **403 AccessDenied** |
| `…/highquality_clust30/foldseekdb.tar.gz` | **403 AccessDenied** |
| foldseek `databases.json` registry (gwdg / mmseqs mirrors) | 302→moved / non-JSON; **ESMAtlas30 no longer resolvable** |
| upstream `foldseek/data/databases.sh` | no ESM-Atlas entry matches |
| `api.esmatlas.com/` hosted search | 403 / 400 (the degraded state that triggered this run); `esmatlas.com/` site itself is 200 |

**Only live mirror:** Steinegger lab Foldcomp worker
(`https://foldcomp.steineggerlab.workers.dev/`) — and it serves the **Foldcomp
*structure* archive, not a Foldseek-searchable 3Di DB**:

| File | Content-Length | What it is |
|---|---|---|
| `highquality_clust30` | **120,442,521,198 B ≈ 120.4 GB** | Foldcomp-compressed structures (NOT a `_ss` 3Di DB) |
| `highquality_clust30.index` | 956,142,258 B ≈ 0.96 GB | per-entry offsets |
| `highquality_clust30.lookup` | 1,030,274,634 B ≈ 1.03 GB | internal id → MGYP |
| `highquality_clust30.dbtype` | 4 B | foldcomp dbtype |
| `highquality_clust30_ss` (3Di sub-DB) | **404 — does not exist** | a Foldcomp DB has no Foldseek 3Di sub-DB |

Entry count (from `.index`/`.lookup` size ÷ ~30 B/line): **≈ 32 million representative
structures.**

## 3. Disk math — why this host STOPs

Free disk on the working volume: **205 GiB (~220 GB).**

A **Foldcomp archive is not directly Foldseek-searchable.** The only route from the
live source to a `foldseek search` is:

1. **Download** the Foldcomp DB: 120.4 + 0.96 + 1.03 GB ≈ **122 GB** — *fits, barely.*
2. **`foldcomp decompress`** → ~32 M PDB structures. Foldcomp expands ~10–20×, so the
   decompressed structure set is **≈ 1.2 – 2.4 TB** — *does NOT fit.*
3. **`foldseek createdb`** over those → a tens-of-GB 3Di+AA DB, **plus** Foldseek `tmp`
   workspace (tens of GB), **plus** the run's query DB / result / fetched hit PDBs.

**Peak requirement ≫ 1 TB vs 205 GiB free → deficit ~5–10×.** The download alone fits;
the mandatory decompress→createdb expansion does not.

> Note: the disk shortfall is *not* the simple "free disk < prebuilt-DB size" the CP0
> brief hypothesized — a prebuilt tens-of-GB 3Di DB *would* fit in 205 GiB. The real
> blockers are (a) that prebuilt DB is **no longer downloadable** (Meta 403; registry
> moved), leaving only the 120 GB Foldcomp archive, and (b) converting that archive
> needs a 1–2 TB intermediate that does not fit.

## 4. Compute blocker (independent of disk)

`foldseek createdb` over ~32 M structures is a **multi-day CPU job** on the M4 Air's
8 cores, even before the first search. The Air is the wrong host on both axes.

## 5. Decision — where this should run

**STOP on the M4 Air.** Run on a big-disk, many-core host:

- **GCE** (project already wired: `compute.gce_burst`, project `projproteus`, zone
  `us-central1-a`, bucket `gs://projproteus-fold`). Provision a **≥ 2 TB persistent
  disk** and a high-vCPU machine. To bound peak disk, **stream** `foldcomp decompress`
  → `foldseek createdb` in shards rather than materializing all 32 M PDBs at once, or
  locate a **prebuilt ESMAtlas30 Foldseek 3Di DB** (Zenodo / a current Steinegger
  mirror) which, at tens of GB, sidesteps the decompress entirely and would even fit a
  modest disk.
- **Origin desktop** — viable if it has ≥ 2 TB free and many cores.

**First action on the new host:** install foldseek (`brew install foldseek` /
bioconda / static binary), confirm `foldseek version`, then re-run this audit's HEAD
probes before pulling 120 GB. Do **not** set `atlas.local_db_path` or flip
`atlas.backend: local_db` on the Air — left unchanged by design.

## 6. Reproducibility — probe commands

```bash
# tooling + disk
which foldseek; df -h .

# source enumeration (HEAD only — no download)
curl -sIL https://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/foldseekdb.txt        # -> 403
for f in highquality_clust30 highquality_clust30.index highquality_clust30.lookup \
         highquality_clust30.dbtype highquality_clust30_ss; do
  curl -s -o /dev/null -w "%{http_code} clen=%header{content-length}  $f\n" \
       -IL "https://foldcomp.steineggerlab.workers.dev/$f"
done
```

Recorded: Atlas version `v2023_02` (config `atlas.version`); DB `highquality_clust30`;
Foldcomp archive 120.4 GB / ~32 M entries; no DB downloaded (HEAD-only audit); foldseek
absent; verdict STOP.
