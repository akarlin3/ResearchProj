# Proteus — env resolution status, failures & deviations

Per Checkpoint 2: "Any tool that won't resolve: record in `envlog/` with the
reason and continue — don't abort the env over one package."

**Target host:** M4 MacBook Air (osx-arm64, MPS/CPU, no CUDA).
**Where this scaffold was generated:** Linux x86_64 cloud container (NOT the Mac).

## The one structural fact that gates everything below

The install + verify step (Checkpoint 2) and the positive-output smoke suite
(Checkpoint 3) **must run on the M4** to be meaningful. They could not be executed
here because this container is Linux x86_64 with no conda/mamba and no Apple
Silicon / MPS. So the lockfiles (`requirements-lock.txt`,
`envlog/conda-env-resolved.yml`) are **placeholders pending generation on the Mac**
— not a real freeze. Fabricating osx-arm64 pins from a Linux host would be wrong
(it would pin the wrong-platform wheels), so we explicitly did not.

## Per-tool resolution plan & expected source (resolve on the Mac, then lock)

| Component | Expected source on osx-arm64 | Status / note |
|---|---|---|
| **PyTorch (MPS)** | conda-forge osx-arm64 (`pytorch`), else pip arm64 wheel | resolve on Mac; assert `torch.backends.mps.is_available()`. CPU/MPS only — **no CUDA**. |
| **transformers / sentencepiece** | pip | ProstT5 tokenizer. Some `transformers` releases mis-route ProstT5's SentencePiece model through an incompatible converter; smoke test falls back to loading `spiece.model` via `sentencepiece` directly. Pin a transformers whose T5 converter handles ProstT5, or rely on the fallback. |
| **Foldseek** | bioconda osx-arm64 → **Homebrew** (`brew install foldseek`) | bioconda arm coverage varies; Homebrew is the reliable arm64 fallback. Record which resolved. |
| **MMseqs2** | bioconda osx-arm64 → **Homebrew** (`brew install mmseqs2`) | same fallback chain as Foldseek. |
| **fpocket** | bioconda → **Homebrew** (`brew install fpocket`) → source build | most likely Homebrew or source on arm64. |
| **AutoDock Vina** | conda-forge osx-arm64 (`vina`) + `meeko` (pip) | Meeko preferred for ligand/receptor prep. |
| **ADFRsuite** | **NOT installed — x86-only** | Scripps distribution has no arm64 build. If required, run under Rosetta 2 (`CONDA_SUBDIR=osx-64`) or prefer Meeko. Flagged, not blocking. |
| **numpy/pandas/biopython/biotite/pyyaml/tqdm/pytest** | conda-forge osx-arm64 | standard, expected clean. |
| **ESMFold / fair-esm / Chai-1 / GNINA / DiffDock** | **n/a locally — Vast-only** | Intentionally absent from the local env. Live in `vast/Dockerfile.fold`; pinned there on first build. |

## What WAS validated here (host-agnostic, pure-Python)

- `s3_fold.py --dry-run` — validates a FASTA and emits a valid job manifest
  (`run_location: vast`, per-seq sha256, fold params from config). ✅ runs & passes.
- Smoke suite **collects and runs**; the S3 dry-run test passes; every
  tool-dependent test SKIPS cleanly (tools absent on this container) rather than
  erroring. ✅
- `controls/fetch_controls.py` reachability: RCSB download of 6EQE returned the
  exact byte count recorded in `controls/MANIFEST.json` (614952 bytes). ✅

## Net env status (on this Linux container): YELLOW — scaffold only

Not GREEN, and cannot be from here: no Apple Silicon to exercise MPS, no osx-arm64
conda solve. GREEN is reachable **only on the M4** after `mamba env create`, a
clean import/run of each LOCAL tool, and a smoke suite where every LOCAL tool shows
positive output. There is **no hard GPU blocker** — folding (S3) is offloaded to
Vast.ai by design (`vast/Dockerfile.fold`), so MPS/CPU is sufficient locally.
