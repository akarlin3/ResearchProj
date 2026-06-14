"""Wire Atlas hits -> the validated `screen` at the WIDENED operating point.

Checkpoint 3 of the Atlas run. Reads the fetched hit structures (structures/
atlas_hits/, manifest data/interim/atlas_hits.json) and runs each through the
EXACT validated path — S4 geometry gate -> S5 cleft -> control-anchored composite
(`proteus.screen.screen_model`) — scored against the SAME IsPETase/LCC anchor as
calibration, at the **widened operating point** (the -1.156 line that recovers all
three divergent positives; `envlog/validation-run.md`). No new metric logic.

Threshold selection (config `atlas.operating_point`):
  * widened (default)  — `calibrate.recovery_screen(...).widened_operating_point`
                         (min over positives + recovered divergent positives).
  * production         — the lowest-positive-control line (`screen`'s default).

Emits a ranked candidate table carrying accession, composite, above/below line,
biome, marine, pLDDT, and triad residues:
  - data/processed/atlas_candidates.csv
  - data/processed/atlas_candidates.json

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.atlas_screen \
        --hits data/interim/atlas_hits.json --struct-dir structures \
        --out data/processed/atlas_candidates
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

from proteus.calibrate import analyze_controls, recovery_screen, score_analysis
from proteus.screen import screen_model
from proteus.utils import DEFAULT_CONFIG, REPO, load_config


def resolve_operating_point(cfg: dict, struct_dir: str) -> dict:
    """Build the IsPETase/LCC control anchor and the requested operating point.

    Returns {anchor, threshold, line, mode, production_threshold, widened, ...}.
    Reuses calibrate untouched: production line = lowest positive composite; widened
    line = the recovery-screen widened_operating_point that also keeps the held-out
    divergent positives."""
    analysis = analyze_controls(cfg, struct_dir)
    cal = score_analysis(analysis, cfg)
    op = cal.get("operating_point") or {}
    if "threshold" not in op:
        raise RuntimeError("calibration produced no operating point (need a positive "
                           "control with a catalytic pocket) — cannot screen.")
    production = op["threshold"]
    line = cfg.get("atlas", {}).get("operating_point", "widened")

    widened = None
    threshold = production
    if line == "widened":
        recov = recovery_screen(cfg, struct_dir, cal)
        widened = recov.get("widened_operating_point")
        if not widened:
            raise RuntimeError("widened operating point requested but recovery_screen "
                               "proposed none (no recovery structure cleared the "
                               "negatives) — fetch the recovery controls or use "
                               "operating_point: production.")
        threshold = widened["threshold"]
    elif line != "production":
        raise ValueError(f"atlas.operating_point must be widened|production, got {line!r}")

    return {
        "anchor": cal["anchor"],
        "mode": cal["mode"],
        "positive_ids": cal["positive_ids"],
        "line": line,
        "threshold": threshold,
        "production_threshold": production,
        "widened": widened,
        "separated": cal["verdict"].get("separated"),
        "margin": cal["verdict"].get("margin"),
    }


def load_hits(hits_json: str) -> tuple[list[dict], dict]:
    """Load the retrieval manifest -> (per-accession meta records, full summary)."""
    with open(hits_json) as fh:
        summary = json.load(fh)
    recs = [h for h in summary.get("hits", []) if h.get("pdb")]
    return recs, summary


def screen_hits(cfg: dict, hits: list[dict], op: dict, struct_dir: str,
                log=print) -> list[dict]:
    """Run each fetched hit through screen_model at the resolved threshold, merging
    the Atlas metadata (biome, marine, pLDDT) onto the screen verdict."""
    anchor, threshold = op["anchor"], op["threshold"]
    out = []
    for h in hits:
        pdb = h["pdb"]
        if not os.path.isabs(pdb):
            pdb = os.path.join(REPO, pdb)
        if not os.path.exists(pdb):
            log(f"[atlas-screen][skip] {h['accession']}: structure missing ({pdb})")
            continue
        rec = screen_model(pdb, cfg, anchor, threshold,
                           mean_plddt=h.get("mean_plddt"), cand_id=h["accession"])
        rec.update({
            "accession": h["accession"],
            "biome": h.get("biome", "unknown"),
            "marine": bool(h.get("marine", False)),
            "biome_source": h.get("biome_source", "unresolved"),
            "atlas_bits": h.get("bits"),
            "atlas_source": h.get("source"),
            "query_id": h.get("query_id"),
            "n_res": h.get("n_res"),
        })
        out.append(rec)
        comp = "-" if rec["composite"] is None else f"{rec['composite']:.3f}"
        verdict = ("HIT" if rec["petase_like_hit"] else
                   (f"below-line({comp})" if rec["pocket_ok"] else
                    rec.get("stage_failed", "rejected")))
        log(f"[atlas-screen] {rec['accession']}: triad="
            f"{'Y' if rec['triad_found'] else 'N'} pocket="
            f"{'Y' if rec['pocket_ok'] else 'N'} composite={comp} "
            f"marine={'Y' if rec['marine'] else 'N'} -> {verdict}")
    return out


def rank_and_summarize(candidates: list[dict], op: dict, summary_in: dict) -> dict:
    scorable = [c for c in candidates if c["composite"] is not None]
    scorable.sort(key=lambda c: c["composite"], reverse=True)
    for rank, c in enumerate(scorable, 1):
        c["rank"] = rank
    hits = [c for c in scorable if c["petase_like_hit"]]
    return {
        "atlas_version": summary_in.get("atlas_version"),
        "query_set_hash": summary_in.get("query_set_hash"),
        "discovery_mode": summary_in.get("discovery_mode"),
        "operating_point": op["line"],
        "threshold": op["threshold"],
        "production_threshold": op["production_threshold"],
        "widened": op["widened"],
        "anchor_mode": op["mode"],
        "positive_ids": op["positive_ids"],
        "controls_separated": op["separated"],
        "n_screened": len(candidates),
        "n_triad": sum(1 for c in candidates if c["triad_found"]),
        "n_pocket": sum(1 for c in candidates if c["pocket_ok"]),
        "n_above_line": len(hits),
        "n_marine": sum(1 for c in candidates if c["marine"]),
        "n_marine_hits": sum(1 for c in hits if c["marine"]),
        "hit_accessions": [c["accession"] for c in hits],
        "ranking": [c["accession"] for c in scorable],
        "candidates": candidates,
    }


CSV_COLS = ["rank", "accession", "composite", "above_line", "petase_like_hit",
            "biome", "marine", "mean_plddt", "n_res", "triad_found",
            "catalytic_ser", "his", "acid", "atlas_bits", "atlas_source", "query_id"]


def write_outputs(summary: dict, out_prefix: str) -> tuple[str, str]:
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    js = out_prefix + ".json"
    with open(js, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    by_acc = {c["accession"]: c for c in summary["candidates"]}
    order = summary["ranking"] + [c["accession"] for c in summary["candidates"]
                                  if c["accession"] not in summary["ranking"]]
    csv_path = out_prefix + ".csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS)
        w.writeheader()
        for acc in order:
            c = by_acc[acc]
            w.writerow({
                "rank": c.get("rank", ""), "accession": c["accession"],
                "composite": c["composite"], "above_line": c["above_threshold"],
                "petase_like_hit": c["petase_like_hit"], "biome": c.get("biome"),
                "marine": c["marine"], "mean_plddt": c.get("mean_plddt"),
                "n_res": c.get("n_res"), "triad_found": c["triad_found"],
                "catalytic_ser": c["catalytic_ser"], "his": c["his"], "acid": c["acid"],
                "atlas_bits": c.get("atlas_bits"), "atlas_source": c.get("atlas_source"),
                "query_id": c.get("query_id"),
            })
    return js, csv_path


def run(cfg: dict, hits_json: str, struct_dir: str, out_prefix: str,
        log=print) -> dict:
    hits, summary_in = load_hits(hits_json)
    if not hits:
        raise RuntimeError(f"no fetched hit structures in {hits_json} — run "
                           "proteus.atlas first")
    op = resolve_operating_point(cfg, struct_dir)
    log(f"[atlas-screen] anchor={op['positive_ids']} mode={op['mode']} "
        f"line={op['line']} threshold={op['threshold']:.4f} "
        f"(production={op['production_threshold']:.4f}; controls "
        f"separated={op['separated']}, margin={op['margin']})")
    candidates = screen_hits(cfg, hits, op, struct_dir, log=log)
    summary = rank_and_summarize(candidates, op, summary_in)
    js, csv_path = write_outputs(summary, out_prefix)
    log(f"[atlas-screen] {summary['n_screened']} screened: {summary['n_triad']} triad, "
        f"{summary['n_pocket']} pocket, {summary['n_above_line']} above the "
        f"{op['line']} line ({summary['n_marine_hits']} marine).")
    log(f"[atlas-screen] candidates -> {os.path.relpath(csv_path, os.getcwd())} ; "
        f"{os.path.relpath(js, os.getcwd())}")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--hits", default=os.path.join(REPO, "data", "interim", "atlas_hits.json"))
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"))
    ap.add_argument("--out", default=os.path.join(REPO, "data", "processed", "atlas_candidates"))
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if not os.path.exists(args.hits):
        print(f"hits manifest not found: {args.hits} (run proteus.atlas first)",
              file=sys.stderr)
        return 2
    try:
        run(cfg, args.hits, args.struct_dir, args.out)
    except (RuntimeError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
