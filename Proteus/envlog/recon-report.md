# Proteus — Checkpoint 0 RECON REPORT

**Generated:** 2026-06-08
**Target host (design):** M4 MacBook Air, Apple Silicon (arm64), MPS, no CUDA.
**Host this scaffold was generated on:** Claude Code on the web — ephemeral
Linux x86_64 cloud container (NOT the target Mac).

> ⚠️ **Read this first.** The repo is scaffolded for an Apple-Silicon M4 Air, but
> the scaffolding session itself ran in a Linux x86_64 container. The probe table
> below therefore reports the *container*, not the Mac. The **recon procedure and
> the GO/WARNING verdict logic are the deliverable**; re-run the probe block on the
> actual M4 (it is one shell paste) to capture the real numbers before installing.
> Nothing in the architecture decision changes: folding is offloaded to Vast.ai, so
> there is no hard GPU blocker on either host.

## Probe block (paste on the M4 to regenerate the real numbers)

```bash
sw_vers                                   # macOS version
uname -m                                  # expect: arm64
sysctl -n machdep.cpu.brand_string        # chip (e.g. "Apple M4")
sysctl -n hw.memsize | awk '{print $1/1073741824" GB unified RAM"}'
vm_stat | awk '/free/{print}'             # free pages -> free RAM
python3 --version
which mamba conda || echo "no conda — install Miniforge (osx-arm64)"
which brew    || echo "no Homebrew — install for arm64 tool fallback"
/usr/bin/pgrep -q oahd && echo "Rosetta 2: present" || echo "Rosetta 2: absent"
python3 -c "import torch; print('mps', torch.backends.mps.is_available())" 2>/dev/null \
    || echo "torch not yet installed (expected pre-CP2)"
df -h .                                   # free disk on working volume
```

## Probes — values observed in THIS scaffolding container (Linux, not the Mac)

| Probe | Target M4 Air (expected) | This container (observed) | Verdict |
|---|---|---|---|
| OS / arch | macOS, `arm64` | Ubuntu/Linux, `x86_64`, kernel 6.18.5 | container ≠ target — re-probe on Mac |
| Chip | Apple M4 | Intel Xeon @ 2.80 GHz (4 vCPU) | re-probe on Mac |
| Unified RAM | 16 GB (M4 Air base) | 16 GB total, ~15.8 GB free | see WARNING below |
| Python | 3.11+ (via Miniforge) | 3.11.15 (system) | OK |
| conda / mamba | install Miniforge (osx-arm64) | both absent | install Miniforge on Mac |
| Homebrew | present (arm64 fallback) | absent | install on Mac (`brew`) |
| Rosetta 2 | present (for x86-only tools) | n/a (already x86 Linux) | install on Mac if ADFR needed |
| torch / MPS | `mps.is_available() == True` | torch not installed | resolve in CP2, verify on Mac |
| CUDA | none (by design) | none | not required — folding offloaded |
| Free disk | check working volume | 31 GB free of 252 GB | OK (corpora live in `data/`, gitignored) |

## Verdict logic (per spec)

This env has **no hard GPU blocker**, because ESMFold / Chai-1 folding is
intentionally offloaded to a Vast.ai burst box (Linux + CUDA, scaffolded in
`vast/`, not installed here). The pipeline narrows **locally** — MMseqs2
dereplication, ProstT5 seq→3Di, Foldseek fold-class triage, fpocket/geometry,
Vina — and ships only the S2 shortlist up for folding. So MPS/CPU is sufficient
for everything that runs on this Mac.

Flag a **WARNING** (not a blocker) only if either holds; otherwise **GO**:

1. **RAM < 16 GB** → ProstT5 / MMseqs2 will be tight.
2. **Neither conda nor Homebrew available** → no install path for the bioconda /
   arm64 tools.

## VERDICT: **GO** (with one WARNING to confirm on the Mac)

- No hard GPU blocker (folding offloaded — by design). ✅ GO premise holds.
- **WARNING — RAM is exactly 16 GB on the M4 Air base config.** That is the spec's
  threshold (`< 16 GB`), so the base M4 Air sits right at the edge: ProstT5
  embedding and large MMseqs2 clusters will be tight. Mitigation already wired into
  `config/proteus.yaml`: modest `s1_tokenize.batch_size`, length caps in `corpus`,
  and dereplication (S0) before any heavy step. Prefer a 24 GB M4 Air if available.
- conda/Homebrew are absent **in this container** but are a one-time install on the
  Mac — not a blocker for the target host. Install **Miniforge (osx-arm64)** for the
  conda solver and **Homebrew (arm64)** as the native fallback for tools lacking an
  osx-arm64 conda build.

**Decision:** GO. Continue through Checkpoints 1–5. Install/verify (CP2) and the
positive-output smoke suite (CP3) must be executed **on the M4** to reach a real
GREEN — they cannot be validated against MPS or osx-arm64 conda from this Linux
container (see `envlog/env-failures.md` for exactly what is and isn't validated here).
