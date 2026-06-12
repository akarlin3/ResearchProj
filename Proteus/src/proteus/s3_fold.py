"""S3 — ESMFold batch fold. RUNS ON GCE (Linux+CUDA), NOT on this Mac.

The pipeline narrows LOCALLY (S0–S2) so only the S2 shortlist is folded. ESMFold
is GPU-heavy and is intentionally offloaded to the GCE burst box
(gce/Dockerfile.fold). On Apple Silicon this stage is **dry-run only**: it
validates the input FASTA and emits a job manifest to ship up to GCE. It MUST
NOT attempt to fold on MPS.

Local usage (dry-run), from the repo root:
    PYTHONPATH=src python -m proteus.s3_fold --dry-run \
        --fasta data/interim/s2_shortlist.fasta \
        --out   data/interim/s3_job_manifest.json

The emitted manifest is the contract handed to the GCE burst runner: it lists
each sequence to fold (id, length, sha256) plus the fold parameters resolved from
config/proteus.yaml (plddt_min, chunk_size, max_recycles, random_seed). The real
fold (loading fair-esm, running esmfold_v1 on CUDA, pLDDT filtering) happens
on GCE — see gce/sync.md.

TODO(P3): implement the GCE-side runner in gce/ (load ESMFold, seed torch from
          random_seed, batch-fold with length chunking, write PDBs + per-model
          pLDDT, drop models below plddt_min).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_CONFIG = os.path.join(REPO, "config", "proteus.yaml")

# Standard one-letter amino-acid alphabet (+ X for unknown). Anything else in a
# sequence is a hard validation error — we don't want to ship junk up to GCE.
_AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


def _load_config(path: str) -> dict:
    """Load config/proteus.yaml. Falls back to embedded defaults if PyYAML or the
    file is unavailable, so a dry-run never hard-fails on a thin host."""
    defaults = {
        "random_seed": 1729,
        "s3_fold": {"plddt_min": 70.0, "chunk_size": 400, "max_recycles": 3,
                    "device": "cuda", "run_location": "gce"},
        "corpus": {"min_length": 80, "max_length": 1000},
    }
    try:
        import yaml  # noqa: PLC0415
        with open(path) as fh:
            cfg = yaml.safe_load(fh) or {}
        # shallow-merge so missing keys fall back to defaults
        for k, v in defaults.items():
            cfg.setdefault(k, v)
        for k, v in defaults["s3_fold"].items():
            cfg["s3_fold"].setdefault(k, v)
        return cfg
    except Exception:  # noqa: BLE001 - config is best-effort for a dry-run
        return defaults


def parse_fasta(path: str):
    """Yield (header_id, sequence) pairs. Minimal, dependency-free FASTA reader."""
    rid, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(seq)
                rid = line[1:].split()[0] if len(line) > 1 else ""
                seq = []
            else:
                seq.append(line.strip())
    if rid is not None:
        yield rid, "".join(seq)


def validate_records(records, min_len: int, max_len: int):
    """Validate parsed FASTA records. Returns (valid_entries, errors)."""
    valid, errors = [], []
    seen = set()
    for rid, seq in records:
        seq = seq.upper()
        if not rid:
            errors.append("record with empty/blank header id")
            continue
        if rid in seen:
            errors.append(f"{rid}: duplicate sequence id")
            continue
        seen.add(rid)
        if not seq:
            errors.append(f"{rid}: empty sequence")
            continue
        bad = sorted(set(seq) - _AA)
        if bad:
            errors.append(f"{rid}: non-amino-acid chars {bad}")
            continue
        n = len(seq)
        too_long = n > max_len  # not fatal: folded in chunks on GCE, just flag it
        valid.append({
            "id": rid,
            "length": n,
            "sha256": hashlib.sha256(seq.encode()).hexdigest(),
            "below_min_length": n < min_len,
            "exceeds_max_length": too_long,
        })
    return valid, errors


def build_manifest(fasta_path: str, cfg: dict, valid) -> dict:
    s3 = cfg["s3_fold"]
    return {
        "schema": "proteus.s3_fold.job_manifest/v1",
        "generated": datetime.now(timezone.utc).isoformat(),
        "stage": "S3_esmfold_batch",
        "run_location": "gce",          # ESMFold runs on GCE, never on this Mac
        "note": "Ship this manifest + the shortlist FASTA to the GCE burst box "
                "(gce/Dockerfile.fold). Do NOT fold on MPS. See gce/sync.md.",
        "input_fasta": os.path.relpath(fasta_path, REPO),
        "random_seed": cfg.get("random_seed"),
        "fold_params": {
            "model": "esmfold_v1",
            "device": s3.get("device", "cuda"),
            "plddt_min": s3.get("plddt_min"),
            "chunk_size": s3.get("chunk_size"),
            "max_recycles": s3.get("max_recycles"),
        },
        "n_sequences": len(valid),
        "total_residues": sum(e["length"] for e in valid),
        "sequences": valid,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="validate FASTA + emit the GCE job manifest. The ONLY mode "
                         "supported locally — folding on this Mac is refused.")
    ap.add_argument("--fasta", default=os.path.join(REPO, "data", "interim", "s2_shortlist.fasta"),
                    help="S2 shortlist FASTA to fold on GCE")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "interim", "s3_job_manifest.json"),
                    help="path to write the job manifest")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    if not args.dry_run:
        print("S3 (ESMFold) does not run on this Apple-Silicon Mac — folding is "
              "offloaded to GCE. Use --dry-run to validate input and emit the "
              "job manifest, then run the fold on GCE (see gce/sync.md).",
              file=sys.stderr)
        return 2

    if not os.path.exists(args.fasta):
        print(f"input FASTA not found: {args.fasta}", file=sys.stderr)
        return 2

    cfg = _load_config(args.config)
    corpus = cfg.get("corpus", {})
    valid, errors = validate_records(
        parse_fasta(args.fasta),
        int(corpus.get("min_length", 0)),
        int(corpus.get("max_length", 10**9)),
    )

    for e in errors:
        print(f"[invalid] {e}", file=sys.stderr)
    if not valid:
        print("no valid sequences to fold — manifest not written", file=sys.stderr)
        return 1

    manifest = build_manifest(args.fasta, cfg, valid)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    print(f"[dry-run] validated {len(valid)} sequence(s), "
          f"{manifest['total_residues']} residues; {len(errors)} rejected.")
    print(f"[dry-run] job manifest -> {os.path.relpath(args.out, os.getcwd())}")
    print("[dry-run] next: ship the manifest + FASTA to GCE (gce/sync.md); fold on CUDA there.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
