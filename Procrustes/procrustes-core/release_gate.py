#!/usr/bin/env python3
"""Procrustes RELEASE gate -- the HOLD mechanism (separate from reproduction).

Reproduction and release are different concerns. ``reproduce.sh`` proves the
manuscript is COMPLETE and reproduces green. THIS gate reports whether the finished
paper may be SUBMITTED. The load-bearing novelty -- distinctness from Gauge --
forward-cites Gauge, so the natural hold is keyed on GAUGE:

    GAUGE_PUBLISHED   (distinctness anchor)   -- published=true AND a real DOI

ARMING THE HOLD IS THE AUTHOR'S CALL (GATE G). The author may instead SUBMIT NOW,
forward-citing Gauge as in-review (release_config.json -> author_decision.mode =
"submit_now"). This gate therefore reports one of three states:

    CLEAR        -- Gauge is published (DOI present): submit with the published cite.
    SUBMIT-NOW   -- author_decision.mode == "submit_now": submit forward-citing Gauge.
    HELD         -- default: awaiting Gauge publication.

It refreshes the SUBMISSION_HOLD marker that the final report reads. As in the
sibling gates, published=true with a null/empty DOI is a configuration error
(no fabricated DOIs), never a release.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "release_config.json"
MARKER = HERE / "SUBMISSION_HOLD"
LOAD_BEARING = ("GAUGE",)   # the distinctness anchor


def load_config(path: Path = CONFIG) -> dict:
    if not path.exists():
        print(f"RELEASE-GATE FAIL: missing {path.name}", file=sys.stderr)
        raise SystemExit(2)
    return json.loads(path.read_text())


def anchor_published(a: dict) -> bool:
    return bool(a.get("published")) and bool(a.get("doi"))


def evaluate(cfg: dict):
    anchors = cfg["anchors"]
    for name, a in anchors.items():
        if a.get("published") and not a.get("doi"):
            print(f"RELEASE-GATE FAIL: {name} marked published with no DOI "
                  f"(no fabricated/empty DOIs permitted)", file=sys.stderr)
            raise SystemExit(2)
    unmet = [f"{n}_PUBLISHED" for n in LOAD_BEARING if not anchor_published(anchors[n])]
    mode = cfg.get("author_decision", {}).get("mode", "hold")
    if not unmet:
        return "CLEAR", unmet, mode
    if mode == "submit_now":
        return "SUBMIT-NOW", unmet, mode
    return "HELD", unmet, mode


def write_marker(state, unmet, mode):
    if state == "CLEAR":
        body = ("SUBMISSION_HOLD: CLEAR\n"
                "Gauge is published; submit with the published citation. Re-verify refs.bib "
                "(swap the @unpublished Gauge entry for its published DOI) first.\n")
    elif state == "SUBMIT-NOW":
        body = ("SUBMISSION_HOLD: CLEAR (author override)\n"
                "author_decision.mode == 'submit_now': the author has elected to submit now, "
                "forward-citing Gauge as in-review. The distinctness claim is marked PROVISIONAL "
                "in-manuscript; swap to Gauge's DOI on publication.\n")
    else:
        body = ("SUBMISSION_HOLD: ACTIVE\n"
                "HELD -- awaiting Gauge publication (the distinctness anchor).\n"
                f"Unmet release conditions: {', '.join(unmet)}.\n"
                "The manuscript is COMPLETE and reproduces green (reproduce.sh); only SUBMISSION "
                "is held. The author may override to submit-now via release_config.json.\n")
    MARKER.write_text(body)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    dry_run = "--dry-run" in argv
    cfg = load_config()
    state, unmet, mode = evaluate(cfg)
    if not dry_run:
        write_marker(state, unmet, mode)
    print("Procrustes release gate (HOLD mechanism; arming is the author's call)")
    print("=" * 68)
    for name, a in cfg["anchors"].items():
        tag = "RELEASE-CONDITION" if name in LOAD_BEARING else "recommended    "
        st = f"PUBLISHED ({a['doi']})" if anchor_published(a) else f"UNPUBLISHED -- {a['status']}"
        print(f"  {name:8s} [{tag}] {a['role'][:48]:48s} -> {st}")
    print(f"  author_decision.mode = {mode!r}")
    print("=" * 68)
    if state == "CLEAR":
        print("RELEASE GATE: CLEAR -- Gauge is published; submit with the published cite.")
        return 0
    if state == "SUBMIT-NOW":
        print("RELEASE GATE: SUBMIT-NOW (author override) -- submit forward-citing Gauge.")
        return 0
    print("RELEASE GATE: HELD -- awaiting Gauge publication (the distinctness anchor).")
    print(f"Unmet: {', '.join(unmet)}. Manuscript is COMPLETE and reproduces green; only "
          "SUBMISSION is held. Author may override (author_decision.mode='submit_now').")
    return 1


if __name__ == "__main__":
    sys.exit(main())
