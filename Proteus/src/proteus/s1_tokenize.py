"""S1 — Tokenize sequences to the 3Di structural alphabet.

ProstT5 (Rostlab/ProstT5) translates amino-acid sequence into Foldseek's 3Di
alphabet WITHOUT folding, giving a cheap structural representation for the
fold-class triage in S2.

BACKEND (decided in CP0 — see envlog/env-failures.md "S1 backend"):
  * PREFERRED — Foldseek-native ProstT5.  ``foldseek createdb <fasta> <db>
    --prostt5-model <weights>`` builds the Foldseek query DB directly and runs
    the ProstT5 AA->3Di translation on CPU. On Apple Silicon this sidesteps both
    the MPS question and the arm64 transformers/sentencepiece fragility, and the
    output is already a valid Foldseek query DB that S2 consumes directly. We
    also export an inspectable 3Di FASTA from the DB's ``_ss`` records.
  * FALLBACK — transformers ProstT5.  Used only if Foldseek lacks the
    ``--prostt5-model`` flag (older builds). Loads Rostlab/ProstT5, translates
    AA->3Di in batches (s1_tokenize.batch_size, device from config w/ cpu
    fallback), and writes a 3Di FASTA.

Weights path resolves from (in order): --prostt5-model arg, config
paths.prostt5_weights, env PROTEUS_PROSTT5_MODEL, then an auto-download via
``foldseek databases ProstT5`` cached under config paths.models. The resolved
location is logged.

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.s1_tokenize \
        --in  data/interim/s0_representatives.fasta \
        --out data/interim/s1_3di
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_CONFIG = os.path.join(REPO, "config", "proteus.yaml")

DEFAULT_IN = os.path.join(REPO, "data", "interim", "s0_representatives.fasta")
DEFAULT_OUT = os.path.join(REPO, "data", "interim", "s1_3di")

# Files produced inside the S1 output dir.
QUERYDB_NAME = "querydb"          # Foldseek query DB prefix (consumed by S2)
THREEDI_FASTA = "s1_3di.fasta"    # inspectable 3Di artifact


def _load_config(path: str) -> dict:
    defaults = {
        "random_seed": 1729,
        "device": "cpu",
        "paths": {"prostt5_weights": "", "models": "models"},
        "s1_tokenize": {"batch_size": 16, "device": "auto"},
    }
    try:
        sys.path.insert(0, os.path.join(REPO, "src"))
        from proteus.utils import load_config  # noqa: PLC0415
        cfg = load_config(path)
        for k, v in defaults.items():
            cfg.setdefault(k, v)
        for k, v in defaults["paths"].items():
            cfg.setdefault("paths", {}).setdefault(k, v)
        for k, v in defaults["s1_tokenize"].items():
            cfg.setdefault("s1_tokenize", {}).setdefault(k, v)
        return cfg
    except Exception:  # noqa: BLE001
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


# --------------------------------------------------------------------------- #
# Weights resolution
# --------------------------------------------------------------------------- #
def resolve_prostt5_weights(cfg: dict, explicit: str | None = None,
                            foldseek_bin: str = "foldseek",
                            allow_download: bool = True) -> str | None:
    """Resolve the ProstT5 weights path (dir or file). Order: explicit arg ->
    config paths.prostt5_weights -> env PROTEUS_PROSTT5_MODEL -> auto-download via
    ``foldseek databases ProstT5`` under paths.models. Returns a path or None."""
    candidates = [
        explicit,
        cfg.get("paths", {}).get("prostt5_weights") or None,
        os.environ.get("PROTEUS_PROSTT5_MODEL") or None,
    ]
    for cand in candidates:
        if cand and os.path.exists(cand):
            print(f"[S1] ProstT5 weights: {cand}")
            return cand

    if not allow_download or shutil.which(foldseek_bin) is None:
        return None

    models_dir = cfg.get("paths", {}).get("models", "models")
    if not os.path.isabs(models_dir):
        models_dir = os.path.join(REPO, models_dir)
    dest = os.path.join(models_dir, "prostt5")
    gguf = os.path.join(dest, "prostt5-f16.gguf")
    if os.path.exists(gguf) or (os.path.isdir(dest) and os.listdir(dest)):
        print(f"[S1] ProstT5 weights (cached): {dest}")
        return dest

    os.makedirs(models_dir, exist_ok=True)
    tmp = os.path.join(models_dir, "_prostt5_dl_tmp")
    print(f"[S1] downloading ProstT5 weights via foldseek -> {dest} (first run)")
    proc = subprocess.run([foldseek_bin, "databases", "ProstT5", dest, tmp],
                          capture_output=True, text=True)
    shutil.rmtree(tmp, ignore_errors=True)
    if proc.returncode != 0 or not os.path.exists(gguf):
        print(f"[S1] weights download failed (rc={proc.returncode}): "
              f"{proc.stderr[-400:]}", file=sys.stderr)
        return None
    print(f"[S1] ProstT5 weights cached: {dest}")
    return dest


# --------------------------------------------------------------------------- #
# Backend: Foldseek-native ProstT5 (preferred)
# --------------------------------------------------------------------------- #
def foldseek_supports_prostt5(foldseek_bin: str = "foldseek") -> bool:
    """True if this Foldseek build exposes ``createdb --prostt5-model`` (native
    ProstT5 3Di generation). Older builds lack it -> fall back to transformers."""
    if shutil.which(foldseek_bin) is None:
        return False
    proc = subprocess.run([foldseek_bin, "createdb", "--help"],
                          capture_output=True, text=True)
    return "--prostt5-model" in (proc.stdout + proc.stderr)


def _extract_3di_fasta(querydb: str, out_fasta: str, foldseek_bin: str) -> None:
    """Write an inspectable 3Di FASTA from a Foldseek query DB's ``_ss`` records.

    ``createdb --prostt5-model`` writes the 3Di sequence DB as ``<db>_ss`` but no
    matching header DB, so we link the amino-acid headers (``<db>_h``) onto it
    first, then convert2fasta.
    """
    ss_db = f"{querydb}_ss"
    ss_h = f"{querydb}_ss_h"
    if not os.path.exists(f"{ss_h}.index") and not os.path.exists(ss_h):
        subprocess.run([foldseek_bin, "lndb", f"{querydb}_h", ss_h],
                       capture_output=True, text=True, check=True)
    subprocess.run([foldseek_bin, "convert2fasta", ss_db, out_fasta],
                   capture_output=True, text=True, check=True)


def tokenize_foldseek(rep_fasta: str, out_dir: str, weights: str,
                      foldseek_bin: str = "foldseek", threads: int | None = None) -> dict:
    """Build the Foldseek query DB (with 3Di) from `rep_fasta` and export a 3Di
    FASTA. Returns {querydb, threedi_fasta, records:[{id,length,len_3di}]}."""
    os.makedirs(out_dir, exist_ok=True)
    querydb = os.path.join(out_dir, QUERYDB_NAME)
    out_fasta = os.path.join(out_dir, THREEDI_FASTA)
    tmp = os.path.join(out_dir, "_createdb_tmp")

    # Clean any stale DB so createdb starts from a known state.
    for f in os.listdir(out_dir):
        if f.startswith(QUERYDB_NAME):
            os.remove(os.path.join(out_dir, f))

    cmd = [foldseek_bin, "createdb", rep_fasta, querydb,
           "--prostt5-model", weights]
    if threads:
        cmd += ["--threads", str(threads)]
    print(f"[S1] running: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    shutil.rmtree(tmp, ignore_errors=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"foldseek createdb failed (rc={proc.returncode}):\n{proc.stderr[-1000:]}")
    if not os.path.exists(f"{querydb}_ss") and not os.path.exists(f"{querydb}_ss.index"):
        raise RuntimeError(f"foldseek produced no 3Di (_ss) DB for {querydb}")

    _extract_3di_fasta(querydb, out_fasta, foldseek_bin)
    threedi = dict(parse_fasta(out_fasta))
    records = []
    for rid, seq in parse_fasta(rep_fasta):
        td = threedi.get(rid, "")
        records.append({"id": rid, "length": len(seq), "len_3di": len(td)})
    return {"querydb": querydb, "threedi_fasta": out_fasta, "records": records}


# --------------------------------------------------------------------------- #
# Backend: transformers ProstT5 (fallback — only if Foldseek lacks the flag)
# --------------------------------------------------------------------------- #
def _resolve_device(cfg: dict) -> str:
    """s1_tokenize.device 'auto' -> mps if available else cpu; never cuda locally."""
    want = str(cfg.get("s1_tokenize", {}).get("device", "auto")).lower()
    try:
        import torch  # noqa: PLC0415
        if want in ("auto", "mps") and torch.backends.mps.is_available():
            return "mps"
    except Exception:  # noqa: BLE001
        pass
    return "cpu"


def tokenize_transformers(rep_fasta: str, out_dir: str, cfg: dict,
                          weights: str = "Rostlab/ProstT5") -> dict:
    """Fallback AA->3Di via transformers ProstT5. Writes a 3Di FASTA (no Foldseek
    query DB — S2 must build one from it with `foldseek base:createdb`). Returns
    the same summary shape as the Foldseek-native backend.

    This is the documented ProstT5 translation recipe: prefix each (spaced) amino
    -acid sequence with the ``<AA2fold>`` control token and decode 3Di with the
    full encoder-decoder via ``generate``. ProstT5 emits 3Di as LOWERCASE letters
    (to disambiguate from AA); we upper-case them for Foldseek. Sampling is seeded
    from the single global ``random_seed`` for reproducibility.

    NOTE: unlike the Foldseek-native path, generation does not hard-guarantee a
    1:1 length, which is exactly why CP0 prefers Foldseek-native on arm64 (see
    envlog/env-failures.md). We log any length mismatch rather than silently pad.
    """
    import re  # noqa: PLC0415

    import torch  # noqa: PLC0415
    from transformers import T5ForConditionalGeneration, T5Tokenizer  # noqa: PLC0415

    os.makedirs(out_dir, exist_ok=True)
    out_fasta = os.path.join(out_dir, THREEDI_FASTA)
    device = _resolve_device(cfg)
    batch_size = int(cfg.get("s1_tokenize", {}).get("batch_size", 16))
    seed = int(cfg.get("random_seed", 0))
    torch.manual_seed(seed)
    print(f"[S1][transformers] device={device} batch_size={batch_size} "
          f"weights={weights} seed={seed}")

    tok = T5Tokenizer.from_pretrained(weights, do_lower_case=False)
    model = T5ForConditionalGeneration.from_pretrained(weights).to(device).eval()

    gen_kwargs = {  # ProstT5 model-card defaults for the AA->3Di direction
        "do_sample": True, "num_beams": 3, "top_p": 0.95, "temperature": 1.2,
        "top_k": 6, "repetition_penalty": 1.2, "num_return_sequences": 1,
    }
    records, items = [], list(parse_fasta(rep_fasta))
    threedi_out: dict[str, str] = {}
    for start in range(0, len(items), batch_size):
        chunk = items[start:start + batch_size]
        spaced = [" ".join(re.sub(r"[UZOB]", "X", s)) for _, s in chunk]
        seqs = ["<AA2fold> " + s for s in spaced]
        enc = tok.batch_encode_plus(seqs, add_special_tokens=True, padding="longest",
                                    return_tensors="pt").to(device)
        max_len = max(len(s) for _, s in chunk) + 1
        with torch.no_grad():
            gen = model.generate(enc.input_ids, attention_mask=enc.attention_mask,
                                 max_length=max_len, min_length=1,
                                 early_stopping=True, **gen_kwargs)
        decoded = tok.batch_decode(gen, skip_special_tokens=True)
        for (rid, seq), td in zip(chunk, decoded):
            td = td.replace(" ", "").upper()  # 3Di letters, Foldseek-style upper-case
            if len(td) != len(seq):
                print(f"[S1][transformers][warn] {rid}: 3Di length {len(td)} != "
                      f"sequence length {len(seq)} (generation drift)", file=sys.stderr)
            threedi_out[rid] = td
    with open(out_fasta, "w") as fh:
        for rid, seq in items:
            td = threedi_out.get(rid, "")
            fh.write(f">{rid}\n{td}\n")
            records.append({"id": rid, "length": len(seq), "len_3di": len(td)})
    return {"querydb": None, "threedi_fasta": out_fasta, "records": records}


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def tokenize(rep_fasta: str, out_dir: str, cfg: dict, prostt5_model: str | None = None,
             foldseek_bin: str = "foldseek") -> dict:
    """Tokenize S0 representatives to 3Di using the CP0-chosen backend.

    Returns a summary dict: backend, n_records, n_nonempty, records, querydb,
    threedi_fasta, all_lengths_match.
    """
    if foldseek_supports_prostt5(foldseek_bin):
        weights = resolve_prostt5_weights(cfg, prostt5_model, foldseek_bin)
        if weights is None:
            raise RuntimeError(
                "Foldseek-native backend selected but ProstT5 weights could not be "
                "resolved (set paths.prostt5_weights, PROTEUS_PROSTT5_MODEL, or allow "
                "`foldseek databases ProstT5` to download).")
        backend = "foldseek-native"
        res = tokenize_foldseek(rep_fasta, out_dir, weights, foldseek_bin)
    else:
        print("[S1] Foldseek lacks --prostt5-model — using transformers fallback.")
        backend = "transformers"
        res = tokenize_transformers(rep_fasta, out_dir, cfg)

    recs = res["records"]
    nonempty = [r for r in recs if r["len_3di"] > 0]
    lengths_match = all(r["len_3di"] == r["length"] for r in recs) and bool(recs)
    summary = {
        "backend": backend,
        "n_records": len(recs),
        "n_nonempty": len(nonempty),
        "all_lengths_match": lengths_match,
        "records": recs,
        "querydb": res["querydb"],
        "threedi_fasta": res["threedi_fasta"],
    }
    print(f"[S1] backend={backend}: {len(recs)} representative(s) -> "
          f"{len(nonempty)} non-empty 3Di string(s); "
          f"length-match={lengths_match}")
    if res["querydb"]:
        print(f"[S1] Foldseek query DB (for S2) -> "
              f"{os.path.relpath(res['querydb'], os.getcwd())}")
    print(f"[S1] inspectable 3Di FASTA -> "
          f"{os.path.relpath(res['threedi_fasta'], os.getcwd())}")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_fasta", default=DEFAULT_IN,
                    help="S0 representative FASTA to tokenize")
    ap.add_argument("--out", default=DEFAULT_OUT,
                    help="output dir for the Foldseek query DB + 3Di FASTA")
    ap.add_argument("--prostt5-model", default=None,
                    help="ProstT5 weights path (overrides config/env)")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    if not os.path.exists(args.in_fasta):
        print(f"input FASTA not found: {args.in_fasta}", file=sys.stderr)
        return 2

    cfg = _load_config(args.config)
    try:
        summary = tokenize(args.in_fasta, args.out, cfg, prostt5_model=args.prostt5_model)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0 if summary["n_nonempty"] == summary["n_records"] and summary["n_records"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
