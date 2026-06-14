"""Atlas retrieval front-end — Foldseek the fold-class query AGAINST the ESM
Metagenomic Atlas, fetch the pre-folded hit structures, biome-tag them.

The Atlas is PRE-FOLDED structures, so discovery INVERTS the normal pipeline:
instead of folding dark-proteome candidates (S0/S1/S3), we search the broad
alpha/beta-hydrolase fold class against the Atlas and pull back hit structures to
feed straight into `screen` (S4->S5 at the widened operating point). This module
is the retrieval front-end only; `proteus.atlas_screen` wires the hits into the
validated screen.

PLUGGABLE BACKEND behind one interface (`AtlasBackend.search`):
  * api (default, pilot) — submit the fold-class reference STRUCTURES to the Atlas
    hosted Foldseek "ticket" search (server-side 3Di; no local Foldseek needed),
    poll, parse hits. On a hosted-search outage / cap / truncation it raises
    `HostedSearchUnavailable`, the caller logs it and surfaces the local-DB
    fallback (and, for the pilot, an unseeded representative slice so the
    fetch->screen path still runs end-to-end).
  * local_db (built, validate-on-first-use) — `foldseek search <query_db>
    <atlas.local_db_path>`. Wired but UNTESTED until the multi-GB HQ clust30 DB is
    downloaded and Foldseek is on PATH; building the query DB reuses the S1
    Foldseek-native ProstT5 path.

UNSEEDED PRINCIPLE: the query is the BROAD fold class (lipases, esterases,
cutinases, AChE, PETases — `structures/` per controls/references.csv), never
PETase templates, so divergent non-homologous alpha/beta-hydrolases surface;
S4/S5 discriminate downstream.

Resolved endpoints + provenance (the API drifts): envlog/atlas-endpoints.md.

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.atlas \
        --out data/interim/atlas_hits.json --hits-dir structures/atlas_hits
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from proteus.calibrate import read_structure_controls
from proteus.utils import DEFAULT_CONFIG, REPO, get_seed, load_config

MGYP_RE = re.compile(r"MGYP\d{6,}")

DEFAULT_ENDPOINTS = {
    "base": "https://api.esmatlas.com",
    "search_structure": "https://api.esmatlas.com/searchStructure",
    "fetch_structure": "https://api.esmatlas.com/fetchPredictedStructure",
    "search_mode": "3diaa",
}
# A browser-ish UA: some Atlas/CloudFront paths 403 a bare urllib UA.
_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
       "AppleWebKit/537.36 (KHTML, like Gecko) Proteus/atlas")


class HostedSearchUnavailable(RuntimeError):
    """The hosted Atlas search errored, was unreachable, was rate-limited, or
    returned a truncated/empty result. Caller logs it and falls back to local_db."""


class BackendNotReady(RuntimeError):
    """A backend is wired but its preconditions are unmet (e.g. local_db needs the
    downloaded HQ Foldseek DB + Foldseek on PATH)."""


# --------------------------------------------------------------------------- #
# HTTP (stdlib urllib; no extra deps)
# --------------------------------------------------------------------------- #
def _request(url: str, *, data: bytes | None = None, headers: dict | None = None,
             timeout: float = 60, max_bytes: int | None = None) -> tuple[int, bytes]:
    """One HTTP call. Returns (status, body). Never raises on HTTP status — returns
    the code so callers can branch (503 -> fallback, 403 -> missing accession). Caps
    the body at `max_bytes` so a huge file (e.g. the 1 GB .lookup) is never fully
    pulled. Raises only on transport-level failure (DNS/timeout/connection)."""
    hdrs = {"User-Agent": _UA}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, data=data, headers=hdrs,
                                 method="POST" if data is not None else "GET")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read(max_bytes) if max_bytes else resp.read()
            return resp.status, body
    except urllib.error.HTTPError as exc:  # 4xx/5xx still carry a body + code
        try:
            body = exc.read(max_bytes) if max_bytes else exc.read()
        except Exception:  # noqa: BLE001
            body = b""
        return exc.code, body


def _post_form(url: str, fields: list[tuple[str, str]], timeout: float) -> tuple[int, bytes]:
    data = urllib.parse.urlencode(fields).encode()
    return _request(url, data=data, timeout=timeout,
                    headers={"Content-Type": "application/x-www-form-urlencoded"})


# --------------------------------------------------------------------------- #
# Fold-class union query set (UNSEEDED)
# --------------------------------------------------------------------------- #
def foldclass_union(cfg: dict, struct_dir: str) -> list[dict]:
    """The fold-class union query: the curated alpha/beta-hydrolase reference
    STRUCTURES present in `struct_dir` (controls/references.csv structure rows),
    minus any `atlas.query_exclude_ids`. Spans PETases + cutinases + lipases +
    esterases + AChE — the broad fold class, NOT PETase templates."""
    exclude = set(cfg.get("atlas", {}).get("query_exclude_ids", []) or [])
    members = []
    for c in read_structure_controls():
        if c["id"] in exclude:
            continue
        pdb = os.path.join(struct_dir, f"{c['accession']}.pdb")
        if os.path.exists(pdb):
            members.append({"id": c["id"], "accession": c["accession"],
                            "role": c["role"], "pdb": pdb})
    return members


def query_set_hash(members: list[dict]) -> str:
    """Reproducibility hash of the query set: sha256 over each member's sorted
    (id, accession, file-sha256). Pins exactly which structures were queried."""
    h = hashlib.sha256()
    for m in sorted(members, key=lambda x: x["accession"]):
        h.update(m["id"].encode())
        h.update(m["accession"].encode())
        with open(m["pdb"], "rb") as fh:
            h.update(hashlib.sha256(fh.read()).hexdigest().encode())
    return h.hexdigest()[:16]


# --------------------------------------------------------------------------- #
# Backends
# --------------------------------------------------------------------------- #
def _endpoints(cfg: dict) -> dict:
    ep = dict(DEFAULT_ENDPOINTS)
    ep.update(cfg.get("atlas", {}).get("endpoints", {}) or {})
    return ep


def _parse_search_results(payload: dict) -> list[dict]:
    """Tolerant parser for the MMseqs2/Foldseek-server result JSON. The hosted
    search is currently 503 so the exact live schema can't be observed; this reads
    the documented shape ({results:[{alignments:[...]}]}) and tolerates field-name
    variants (score|bits, eval|evalue, tmScore|prob). Returns per-hit dicts."""
    hits = []
    results = payload.get("results") or payload.get("result") or []
    if isinstance(results, dict):
        results = [results]
    for block in results:
        aligns = block.get("alignments") or block.get("alns") or []
        # foldseek-server sometimes nests one list deeper ([[...]])
        if aligns and isinstance(aligns[0], list):
            aligns = [a for grp in aligns for a in grp]
        for a in aligns:
            if not isinstance(a, dict):
                continue
            target = str(a.get("target") or a.get("targetId") or a.get("t") or "")
            m = MGYP_RE.search(target)
            if not m:
                continue
            bits = a.get("bits", a.get("score"))
            evalue = a.get("eval", a.get("evalue"))
            tm = a.get("tmScore", a.get("prob"))
            hits.append({
                "accession": m.group(0),
                "bits": float(bits) if bits is not None else None,
                "evalue": float(evalue) if evalue is not None else None,
                "tm": float(tm) if tm is not None else None,
            })
    return hits


class ApiBackend:
    """Hosted Atlas Foldseek search (the ticket API). Submits each fold-class
    reference STRUCTURE; the server does the 3Di conversion."""

    name = "api"

    def __init__(self, cfg: dict):
        a = cfg.get("atlas", {})
        self.ep = _endpoints(cfg)
        self.db = self.ep.get("search_db", "highquality_clust30")  # FOLDSEEK_DB token
        self.mode = self.ep.get("search_mode", "3diaa")
        self.timeout = float(a.get("request_timeout_s", 60))
        self.poll_interval = float(a.get("poll_interval_s", 3))
        self.poll_max = float(a.get("poll_max_s", 180))

    def _submit(self, pdb_text: str) -> str:
        url = f"{self.ep['search_structure']}/ticket"
        fields = [("q", pdb_text), ("mode", self.mode), ("email", ""),
                  ("database[]", self.db)]
        status, body = _post_form(url, fields, self.timeout)
        if status != 200:
            raise HostedSearchUnavailable(
                f"submit -> HTTP {status}: {body[:160].decode('utf-8', 'replace')}")
        try:
            return json.loads(body)["id"]
        except Exception as exc:  # noqa: BLE001
            raise HostedSearchUnavailable(f"submit: unparsable ticket: {exc}")

    def _await(self, ticket: str) -> None:
        url = f"{self.ep['search_structure']}/ticket/{ticket}"
        deadline = self.poll_max
        waited = 0.0
        while waited <= deadline:
            status, body = _request(url, timeout=self.timeout)
            if status != 200:
                raise HostedSearchUnavailable(f"poll -> HTTP {status}")
            st = (json.loads(body or b"{}").get("status") or "").upper()
            if st in ("COMPLETE", "COMPLETED", "DONE"):
                return
            if st in ("ERROR", "FAILED", "UNKNOWN"):
                raise HostedSearchUnavailable(f"ticket status {st}")
            time.sleep(self.poll_interval)
            waited += self.poll_interval
        raise HostedSearchUnavailable(f"ticket {ticket} not complete after {deadline}s")

    def _results(self, ticket: str) -> list[dict]:
        url = f"{self.ep['search_structure']}/result/{ticket}/0"
        status, body = _request(url, timeout=self.timeout)
        if status != 200:
            raise HostedSearchUnavailable(f"results -> HTTP {status}")
        try:
            payload = json.loads(body)
        except Exception as exc:  # noqa: BLE001
            raise HostedSearchUnavailable(f"results: unparsable JSON: {exc}")
        return _parse_search_results(payload)

    def search(self, query_members: list[dict], log=print) -> list[dict]:
        """Search every fold-class query structure; UNION the hits (best bits per
        accession), tagging which query member surfaced each. Any hosted error on
        any member raises HostedSearchUnavailable (the whole api search is treated
        as unavailable so the caller falls back cleanly)."""
        best: dict[str, dict] = {}
        for qm in query_members:
            with open(qm["pdb"]) as fh:
                pdb_text = fh.read()
            log(f"[atlas:api] submit query {qm['id']} ({qm['accession']}) "
                f"-> {self.ep['search_structure']}/ticket")
            ticket = self._submit(pdb_text)
            self._await(ticket)
            for h in self._results(ticket):
                h = dict(h, query_id=qm["id"], source="search")
                acc = h["accession"]
                if acc not in best or (h["bits"] or -1) > (best[acc]["bits"] or -1):
                    best[acc] = h
        return list(best.values())


class LocalDbBackend:
    """Local Foldseek search of the fold-class query DB against the downloaded HQ
    clust30 Foldseek DB. WIRED but UNTESTED until the DB exists and Foldseek is on
    PATH. Building the query DB reuses the S1 Foldseek-native ProstT5 path."""

    name = "local_db"

    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.db_path = cfg.get("atlas", {}).get("local_db_path", "") or ""

    def search(self, query_members: list[dict], log=print) -> list[dict]:
        import shutil  # noqa: PLC0415
        if not self.db_path or not (os.path.exists(self.db_path)
                                    or os.path.exists(self.db_path + ".dbtype")):
            raise BackendNotReady(
                "local_db backend: atlas.local_db_path is unset/missing — download "
                "the HQ clust30 Foldseek DB (see envlog/atlas-endpoints.md) first.")
        if shutil.which("foldseek") is None:
            raise BackendNotReady("local_db backend: foldseek not on PATH.")
        # Build the query DB from the fold-class union via the S1 ProstT5 path, then
        # `foldseek search`. Implemented for the sweep; intentionally not exercised
        # in the pilot (no DB downloaded). Reuses proteus.s1_tokenize + s2 helpers.
        raise BackendNotReady(
            "local_db search is wired but validate-on-first-use: run it once the HQ "
            "Foldseek DB is downloaded (it reuses S1 ProstT5 query-DB build + "
            "`foldseek search`; see proteus.s2_foldclass_triage for the search call).")


def make_backend(cfg: dict):
    name = cfg.get("atlas", {}).get("backend", "api")
    if name == "api":
        return ApiBackend(cfg)
    if name == "local_db":
        return LocalDbBackend(cfg)
    raise ValueError(f"unknown atlas.backend {name!r} (expected api|local_db)")


# --------------------------------------------------------------------------- #
# Degraded discovery fallback (pilot only) — unseeded representative slice
# --------------------------------------------------------------------------- #
def representative_sample(cfg: dict, n: int, log=print) -> list[dict]:
    """When the hosted search is down AND no local_db is configured, draw a bounded
    UNSEEDED slice of the HQ clust30 representative set (the Foldseek DB's own
    accession list) so the fetch->screen machinery still runs end-to-end on real
    Atlas structures. This is NOT fold-class targeted — it samples the whole HQ
    Atlas, so the downstream S4/S5 survival is the *unfiltered base rate* (a
    conservative floor / economics anchor), clearly labelled as such everywhere.

    Range-fetches only the first chunk of the ~1 GB .lookup (never the whole file),
    then seed-samples `n` accessions reproducibly."""
    a = cfg.get("atlas", {})
    url = a.get("representative_lookup_url")
    if not url:
        return []
    # Pull a bounded window (≈512 KB ~ 20k ids) — hard byte cap so the 1 GB file is
    # never fully downloaded even if the server ignores Range.
    window = 512 * 1024
    status, body = _request(url, timeout=float(a.get("request_timeout_s", 60)),
                            headers={"Range": f"bytes=0-{window - 1}"},
                            max_bytes=window)
    if status not in (200, 206) or not body:
        log(f"[atlas:fallback] representative lookup -> HTTP {status}; no sample")
        return []
    text = body.decode("utf-8", "replace")
    # drop a possibly-truncated final line
    accs = MGYP_RE.findall(text[:text.rfind("\n") + 1] or text)
    uniq = sorted(dict.fromkeys(accs))
    if not uniq:
        return []
    rng = random.Random(get_seed(cfg))
    n = min(n, len(uniq))
    picked = sorted(rng.sample(uniq, n))
    log(f"[atlas:fallback] HQ clust30 representative slice: sampled {n} of "
        f"{len(uniq)} accessions from first {len(body)} B of the .lookup "
        f"(seed={get_seed(cfg)})")
    return [{"accession": acc, "bits": None, "evalue": None, "tm": None,
             "query_id": None, "source": "representative_fallback"} for acc in picked]


# --------------------------------------------------------------------------- #
# Fetch structures + mean pLDDT + biome
# --------------------------------------------------------------------------- #
def mean_plddt(pdb_path: str) -> tuple[float | None, int]:
    """Mean per-residue pLDDT (= mean CA B-factor) and residue count for an ESMFold
    Atlas model. Atlas V0 models carry pLDDT on a 0–1 scale (HQ subset mean ≈0.77)."""
    vals = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    if not vals:
        return None, 0
    return round(sum(vals) / len(vals), 4), len(vals)


def fetch_structure(accession: str, dest_dir: str, endpoints: dict,
                    timeout: float = 60) -> dict:
    """Fetch one pre-folded Atlas structure by accession into dest_dir. Returns a
    record with the path, byte size, mean pLDDT and residue count, or status=missing
    (HTTP 403 = accession not in the served HQ subset) / error (transport)."""
    os.makedirs(dest_dir, exist_ok=True)
    url = f"{endpoints['fetch_structure']}/{accession}.pdb"
    rec = {"accession": accession, "pdb": None, "bytes": 0,
           "mean_plddt": None, "n_res": 0, "fetch_status": None}
    try:
        status, body = _request(url, timeout=timeout)
    except Exception as exc:  # noqa: BLE001 - transport failure is non-fatal per-hit
        rec["fetch_status"] = f"error:{type(exc).__name__}"
        return rec
    rec["fetch_status"] = status
    if status != 200 or not body.startswith(b"HEADER") and b"ATOM" not in body[:400]:
        return rec  # 403 = not in HQ subset; anything non-PDB -> skip
    path = os.path.join(dest_dir, f"{accession}.pdb")
    with open(path, "wb") as fh:
        fh.write(body)
    rec["pdb"] = path
    rec["bytes"] = len(body)
    rec["mean_plddt"], rec["n_res"] = mean_plddt(path)
    return rec


_MARINE_TERMS = ("marine", "ocean", "sea", "seawater", "coast", "reef", "pelagic",
                 "hydrothermal", "saline water", "estuar")


def biome_lookup(accession: str, mgnify_api: str, timeout: float = 30) -> dict:
    """Best-effort biome for an MGnify accession via the MGnify API. NON-BLOCKING:
    any failure (no mapping, network, rate limit) -> biome 'unknown', marine False.
    The Atlas API exposes no per-accession biome; the authoritative per-MGYP biome
    lives in the multi-GB metadata sqlite the pilot does not download — so unknown
    is the expected, acceptable result for most raw representatives."""
    rec = {"biome": "unknown", "marine": False, "biome_source": "unresolved"}
    if not mgnify_api:
        return rec
    url = f"{mgnify_api.rstrip('/')}/proteins/{accession}"
    try:
        status, body = _request(url, timeout=timeout,
                                headers={"Accept": "application/json"})
        if status != 200 or not body:
            return rec
        data = json.loads(body)
    except Exception:  # noqa: BLE001
        return rec
    # Tolerant: scan the JSON for a biome/lineage string.
    blob = json.dumps(data).lower()
    m = re.search(r'root[:>][^"]+', blob)
    if m:
        rec["biome"] = m.group(0)[:120]
        rec["biome_source"] = "mgnify"
    rec["marine"] = any(t in blob for t in _MARINE_TERMS)
    return rec


# --------------------------------------------------------------------------- #
# Driver: retrieve (CP1) + fetch & biome (CP2)
# --------------------------------------------------------------------------- #
def run_retrieval(cfg: dict, struct_dir: str, out_json: str, hits_dir: str,
                  log=print) -> dict:
    """End-to-end retrieval: build the fold-class query, run the configured backend
    (falling back on a hosted-search outage), apply min_bits + hit_cap, fetch each
    hit structure, record mean pLDDT + biome. Writes data/interim/atlas_hits.json."""
    a = cfg.get("atlas", {})
    members = foldclass_union(cfg, struct_dir)
    if not members:
        raise RuntimeError(f"no fold-class query structures found in {struct_dir} "
                           "(need the control PDBs from controls/fetch_controls.py)")
    qhash = query_set_hash(members)
    log(f"[atlas] query_set={a.get('query_set')} ({len(members)} structures, "
        f"hash={qhash}): {[m['id'] for m in members]}")
    log(f"[atlas] backend={a.get('backend')} version={a.get('version')} "
        f"hit_cap={a.get('hit_cap')} min_bits={a.get('min_bits')} "
        f"operating_point={a.get('operating_point')}")

    backend = make_backend(cfg)
    discovery_mode = backend.name
    fallback_note = None
    try:
        hits = backend.search(members, log=log)
        log(f"[atlas] {backend.name} search returned {len(hits)} raw hit(s)")
    except (HostedSearchUnavailable, BackendNotReady) as exc:
        log(f"[atlas][WARN] {backend.name} search unavailable: {exc}")
        fallback_note = str(exc)
        hits = []
        # Surface the documented fallback: download the local DB and flip backend.
        log("[atlas] FALLBACK: hosted search unavailable -> the systematic path is "
            "`backend: local_db` against the downloaded HQ clust30 Foldseek DB "
            "(envlog/atlas-endpoints.md).")
        if a.get("allow_representative_fallback", False) and backend.name == "api":
            n = int(a.get("pilot_sample_size", 60))
            log(f"[atlas] pilot degraded-discovery: drawing an UNSEEDED {n}-accession "
                "representative slice so fetch->screen runs end-to-end on real Atlas "
                "structures (NOT fold-class targeted — base-rate floor).")
            hits = representative_sample(cfg, n, log=log)
            discovery_mode = "representative_fallback"

    # min_bits (skip None — representative slice has no score) + hit_cap.
    min_bits = float(a.get("min_bits", 0))
    kept = [h for h in hits if h.get("bits") is None or h["bits"] >= min_bits]
    dropped_bits = len(hits) - len(kept)
    cap = int(a.get("hit_cap", 2000))
    capped = len(kept) > cap
    kept = sorted(kept, key=lambda h: (h.get("bits") is not None, h.get("bits") or 0),
                  reverse=True)[:cap]
    log(f"[atlas] thresholds: {len(kept)} hit(s) after min_bits>={min_bits} "
        f"(dropped {dropped_bits}) + hit_cap={cap}{' (TRUNCATED)' if capped else ''}")

    # Fetch structures + pLDDT + biome.
    ep = _endpoints(cfg)
    timeout = float(a.get("request_timeout_s", 60))
    mgnify = a.get("mgnify_api", "")
    fetched, n_ok, n_missing = [], 0, 0
    for i, h in enumerate(kept, 1):
        fr = fetch_structure(h["accession"], hits_dir, ep, timeout=timeout)
        rec = dict(h)
        rec.update({k: fr[k] for k in ("pdb", "bytes", "mean_plddt", "n_res",
                                       "fetch_status")})
        if fr["pdb"]:
            n_ok += 1
            rec.update(biome_lookup(h["accession"], mgnify, timeout=min(timeout, 30)))
        else:
            n_missing += 1
            rec.update({"biome": "unknown", "marine": False,
                        "biome_source": "unresolved"})
        fetched.append(rec)
        if i % 10 == 0 or i == len(kept):
            log(f"[atlas] fetched {i}/{len(kept)} ({n_ok} ok, {n_missing} missing)")

    summary = {
        "atlas_version": a.get("version"),
        "query_set": a.get("query_set"),
        "query_set_hash": qhash,
        "query_members": [{"id": m["id"], "accession": m["accession"],
                           "role": m["role"]} for m in members],
        "backend": backend.name,
        "discovery_mode": discovery_mode,
        "fallback_note": fallback_note,
        "min_bits": min_bits,
        "hit_cap": cap,
        "hit_cap_truncated": capped,
        "operating_point": a.get("operating_point"),
        "endpoints": ep,
        "n_hits_raw": len(hits),
        "n_hits_kept": len(kept),
        "n_fetched": n_ok,
        "n_missing": n_missing,
        "n_marine": sum(1 for r in fetched if r.get("marine")),
        "n_biome_resolved": sum(1 for r in fetched
                                if r.get("biome_source") == "mgnify"),
        "hits": fetched,
    }
    os.makedirs(os.path.dirname(os.path.abspath(out_json)), exist_ok=True)
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    log(f"[atlas] {n_ok} structure(s) -> {os.path.relpath(hits_dir, os.getcwd())}/ ; "
        f"manifest -> {os.path.relpath(out_json, os.getcwd())}")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"),
                    help="dir with the fold-class reference (control) structures")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "interim", "atlas_hits.json"))
    ap.add_argument("--hits-dir", default=os.path.join(REPO, "structures", "atlas_hits"))
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if "atlas" not in cfg:
        print("config has no `atlas` block — add it (see config/proteus.yaml).",
              file=sys.stderr)
        return 2
    try:
        summary = run_retrieval(cfg, args.struct_dir, args.out, args.hits_dir)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(f"[atlas] done: {summary['n_fetched']} fetched, "
          f"{summary['n_marine']} marine-tagged, mode={summary['discovery_mode']}")
    return 0 if summary["n_fetched"] > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
