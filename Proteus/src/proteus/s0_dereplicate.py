"""S0 — Dereplicate the input corpus with MMseqs2.

This is DEREPLICATION, not a homology gate: we collapse near-identical
sequences (>= s0_dereplicate.min_seq_id) to a representative set so downstream
folding isn't wasted on duplicates. We deliberately do NOT filter by similarity
to known PETases — that would discard the divergent dark-tail candidates we are
hunting for. The guard below logs that no homology gate was applied.

Reads thresholds and the random seed from config/proteus.yaml.

Backend: MMseqs2 ``easy-cluster`` (greedy set-cover). It emits, for a prefix P:
  - P_rep_seq.fasta   the cluster representatives  -> data/interim/s0_representatives.fasta
  - P_cluster.tsv     rep<TAB>member membership     -> data/interim/s0_clusters.tsv
  - P_all_seqs.fasta  all members grouped (unused downstream)

Seeding note: MMseqs2 ``easy-cluster`` exposes NO RNG seed — the clustering is a
deterministic greedy set-cover, and the only stochastic knob is ``--shuffle``
(input-order shuffle of the DB). For reproducibility we set ``--shuffle 0`` and
still read + log the single global ``random_seed`` from config for provenance
(see envlog/env-failures.md, "MMseqs2 has no clustering seed").

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.s0_dereplicate \
        --in  tests/data/mini_corpus.fasta \
        --rep data/interim/s0_representatives.fasta \
        --tsv data/interim/s0_clusters.tsv
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_CONFIG = os.path.join(REPO, "config", "proteus.yaml")

DEFAULT_IN = os.path.join(REPO, "data", "interim", "corpus.fasta")
DEFAULT_REP = os.path.join(REPO, "data", "interim", "s0_representatives.fasta")
DEFAULT_TSV = os.path.join(REPO, "data", "interim", "s0_clusters.tsv")


def _load_config(path: str) -> dict:
    """Load config/proteus.yaml via the shared utils loader, falling back to the
    documented S0 defaults if PyYAML / the file is unavailable on a thin host."""
    defaults = {
        "random_seed": 1729,
        "s0_dereplicate": {"min_seq_id": 0.95, "coverage": 0.90, "cov_mode": 1},
    }
    try:
        sys.path.insert(0, os.path.join(REPO, "src"))
        from proteus.utils import load_config  # noqa: PLC0415
        cfg = load_config(path)
        cfg.setdefault("random_seed", defaults["random_seed"])
        cfg.setdefault("s0_dereplicate", {})
        for k, v in defaults["s0_dereplicate"].items():
            cfg["s0_dereplicate"].setdefault(k, v)
        return cfg
    except Exception:  # noqa: BLE001 - config load is best-effort
        return defaults


def parse_fasta(path: str):
    """Yield (id, sequence) pairs. Minimal dependency-free FASTA reader."""
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


def parse_clusters_tsv(path: str) -> dict[str, list[str]]:
    """Read an MMseqs2 ``_cluster.tsv`` (rep<TAB>member) into {rep: [members...]}."""
    clusters: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep, member = line.split("\t")[:2]
            clusters.setdefault(rep, []).append(member)
    return clusters


def dereplicate(in_fasta: str, rep_fasta: str, clusters_tsv: str, cfg: dict,
                tmp_dir: str | None = None, mmseqs_bin: str = "mmseqs") -> dict:
    """Run MMseqs2 easy-cluster and place the representative FASTA + cluster TSV.

    Returns a summary dict: n_input, n_representatives, n_clusters, params, seed.
    Raises FileNotFoundError if mmseqs is absent and RuntimeError on a failed run.
    """
    if shutil.which(mmseqs_bin) is None:
        raise FileNotFoundError(
            f"'{mmseqs_bin}' not found on PATH — install MMseqs2 "
            "(bioconda/Homebrew; see envlog/env-failures.md)")

    s0 = cfg["s0_dereplicate"]
    min_seq_id = float(s0["min_seq_id"])
    coverage = float(s0["coverage"])
    cov_mode = int(s0["cov_mode"])
    seed = int(cfg["random_seed"])

    # --- GUARD: dereplication, NOT a homology gate. -------------------------- #
    # We collapse near-identical sequences only; we never filter by similarity to
    # known PETases / any reference. Every input sequence is eligible to be (or to
    # be represented by) a representative. Logged so the run record proves it.
    print("[S0][guard] DEREPLICATION ONLY — no homology gate applied: sequences "
          "are clustered against EACH OTHER, never filtered by similarity to "
          "known PETases or any reference set (the divergent dark tail is kept).")
    print(f"[S0] params: min_seq_id={min_seq_id} coverage={coverage} "
          f"cov_mode={cov_mode}")
    print(f"[S0] random_seed={seed} (provenance only — MMseqs2 easy-cluster has no "
          "RNG seed; determinism pinned via --shuffle 0)")

    n_input = sum(1 for _ in parse_fasta(in_fasta))
    if n_input == 0:
        raise RuntimeError(f"input FASTA has no sequences: {in_fasta}")

    cleanup = tmp_dir is None
    work = tmp_dir or tempfile.mkdtemp(prefix="s0_mmseqs_")
    try:
        os.makedirs(work, exist_ok=True)
        prefix = os.path.join(work, "clu")
        mmseqs_tmp = os.path.join(work, "tmp")
        cmd = [
            mmseqs_bin, "easy-cluster", in_fasta, prefix, mmseqs_tmp,
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode),
            "--shuffle", "0",            # deterministic input order (no RNG seed exists)
        ]
        print(f"[S0] running: {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"mmseqs easy-cluster failed (rc={proc.returncode}):\n"
                f"{proc.stderr[-1000:]}")

        src_rep = f"{prefix}_rep_seq.fasta"
        src_tsv = f"{prefix}_cluster.tsv"
        for src in (src_rep, src_tsv):
            if not os.path.exists(src):
                raise RuntimeError(f"expected MMseqs2 output missing: {src}")

        for dest in (rep_fasta, clusters_tsv):
            os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
        shutil.copyfile(src_rep, rep_fasta)
        shutil.copyfile(src_tsv, clusters_tsv)
    finally:
        if cleanup:
            shutil.rmtree(work, ignore_errors=True)

    reps = [rid for rid, _ in parse_fasta(rep_fasta)]
    clusters = parse_clusters_tsv(clusters_tsv)
    summary = {
        "n_input": n_input,
        "n_representatives": len(reps),
        "n_clusters": len(clusters),
        "representatives": reps,
        "params": {"min_seq_id": min_seq_id, "coverage": coverage,
                   "cov_mode": cov_mode},
        "random_seed": seed,
        "homology_gate": False,
        "rep_fasta": rep_fasta,
        "clusters_tsv": clusters_tsv,
    }
    collapsed = n_input - len(reps)
    print(f"[S0] {n_input} input -> {len(reps)} representatives "
          f"({len(clusters)} clusters; {collapsed} sequence(s) collapsed as "
          "near-duplicates)")
    print(f"[S0] representatives -> {os.path.relpath(rep_fasta, os.getcwd())}")
    print(f"[S0] cluster membership -> {os.path.relpath(clusters_tsv, os.getcwd())}")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_fasta", default=DEFAULT_IN,
                    help="input corpus FASTA to dereplicate")
    ap.add_argument("--rep", default=DEFAULT_REP,
                    help="output representative FASTA")
    ap.add_argument("--tsv", default=DEFAULT_TSV,
                    help="output cluster-membership TSV (rep<TAB>member)")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--tmp", default=None,
                    help="working dir for MMseqs2 (default: a fresh temp dir)")
    args = ap.parse_args(argv)

    if not os.path.exists(args.in_fasta):
        print(f"input FASTA not found: {args.in_fasta}", file=sys.stderr)
        return 2

    cfg = _load_config(args.config)
    try:
        dereplicate(args.in_fasta, args.rep, args.tsv, cfg, tmp_dir=args.tmp)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 3
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
