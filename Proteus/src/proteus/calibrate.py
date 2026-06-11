"""S4 -> S5 calibration & separation report on the control set (Checkpoint 5).

Runs the catalytic-geometry gate (S4) then the cleft filter (S5) across every control
structure, z-scores the cleft metrics against the positive controls, and asks the
make-or-break question: do the known PET hydrolases (IsPETase, LCC-WT) score above
every non-PET serine hydrolase that shares the fold and triad?

Emits:
  - envlog/calibration-report.md   (table + separation verdict + operating point + LOO)
  - data/processed/s5_scores.csv   (ranked composite table)
  - data/processed/s4_triads.json  (per-model triad detail)

N is tiny (2 positive scaffolds, 4 negatives). This establishes FACE VALIDITY and a
PROVISIONAL operating point only — see the honest-stats section of the report.

Usage:
    PYTHONPATH=src python -m proteus.calibrate
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil

from proteus.s4_geometry import analyze_model
from proteus.s5_cleft_filter import analyze_cleft, score_controls
from proteus.utils import DEFAULT_CONFIG, REPO, load_config

CONTROLS_CSV = os.path.join(REPO, "controls", "references.csv")
MANIFEST = os.path.join(REPO, "controls", "MANIFEST.json")
TRAP_IDS = {"LCC_ICCG"}  # the S165A inactivated control — expected to have no triad
# Controls that MUST be present to calibrate at all (Checkpoint 0 precondition).
REQUIRED_STRUCTURES = ("6EQE", "4EB0", "6THS", "1TCA")


# --------------------------------------------------------------------------- #
# Checkpoint 0 — precondition audit
# --------------------------------------------------------------------------- #
def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def preconditions(cfg: dict, struct_dir: str) -> dict:
    """Enumerate every precondition for calibration. Returns {ok: bool, checks: [...]}.

    Calibration cannot run without (a) a working env (torch/biotite/biopython + fpocket
    on PATH), (b) the control structures present with sha256 matching MANIFEST.json, and
    (c) a config with the s4_geometry/s5_cleft_filter blocks."""
    checks = []

    def add(name, ok, detail):
        checks.append({"name": name, "ok": bool(ok), "detail": detail})

    # (a) environment
    for mod in ("torch", "biotite", "Bio"):
        try:
            __import__(mod)
            add(f"import {mod}", True, "ok")
        except Exception as exc:  # noqa: BLE001
            add(f"import {mod}", False, f"{type(exc).__name__}: {exc}")
    fp = shutil.which("fpocket")
    add("fpocket on PATH", fp is not None, fp or "not found — install fpocket")

    # (b) controls present + sha256 match
    manifest_sha = {}
    if os.path.exists(MANIFEST):
        try:
            man = json.load(open(MANIFEST))
            for s in man.get("structures", []):
                if s.get("sha256"):
                    manifest_sha[s["accession"].upper()] = s["sha256"]
            add("MANIFEST.json loads", True, f"{len(manifest_sha)} sha256 records")
        except Exception as exc:  # noqa: BLE001
            add("MANIFEST.json loads", False, str(exc))
    else:
        add("MANIFEST.json present", False, "run controls/fetch_controls.py")

    for pid in REQUIRED_STRUCTURES:
        path = os.path.join(struct_dir, f"{pid}.pdb")
        if not os.path.exists(path):
            add(f"control {pid}.pdb", False, "missing — run controls/fetch_controls.py")
            continue
        want = manifest_sha.get(pid)
        if not want:
            add(f"control {pid}.pdb", False, "no sha256 in MANIFEST — re-run fetch_controls.py")
            continue
        got = _sha256(path)
        add(f"control {pid}.pdb sha256", got == want,
            "match" if got == want else f"MISMATCH got {got[:12]} want {want[:12]}")

    # (c) config blocks
    add("config random_seed", "random_seed" in cfg, str(cfg.get("random_seed")))
    add("config s4_geometry block", "s4_geometry" in cfg,
        "present" if "s4_geometry" in cfg else "MISSING")
    add("config s5_cleft_filter block", "s5_cleft_filter" in cfg,
        "present" if "s5_cleft_filter" in cfg else "MISSING")

    return {"ok": all(c["ok"] for c in checks), "checks": checks}


def print_precondition_report(audit: dict) -> None:
    print("=== PRECONDITION REPORT (Checkpoint 0) ===")
    for c in audit["checks"]:
        print(f"  [{'PASS' if c['ok'] else 'FAIL'}] {c['name']}: {c['detail']}")
    print(f"  -> {'GO' if audit['ok'] else 'STOP — preconditions unmet'}")


def read_structure_controls(path: str = CONTROLS_CSV) -> list:
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            row = {k: (v or "").strip() for k, v in row.items() if k is not None}
            if row.get("type") == "structure":
                rows.append({"id": row["id"], "accession": row["accession"].upper(),
                             "role": row["role"]})
    return rows


def classify(control: dict, positive_ids) -> str:
    if control["id"] in positive_ids:
        return "positive"
    if control["id"] in TRAP_IDS:
        return "trap"
    if control["role"] == "negative":
        return "negative"
    return "other"


def run_calibration(cfg: dict, struct_dir: str) -> dict:
    positive_ids = list(cfg["s5_cleft_filter"]["positive_controls"])
    controls = read_structure_controls()

    per_control, s4_detail = {}, {}
    for c in controls:
        cls = classify(c, positive_ids)
        pdb = os.path.join(struct_dir, f"{c['accession']}.pdb")
        rec = {"id": c["id"], "accession": c["accession"], "role": c["role"],
               "class": cls, "present": os.path.exists(pdb)}
        if not rec["present"]:
            per_control[c["id"]] = rec
            continue
        s4 = analyze_model(pdb, cfg)
        s4_detail[c["accession"]] = s4
        rec["triad_found"] = s4["triad_found"]
        rec["catalytic_ser"] = s4["best"]["ser"]["res_id"] if s4["best"] else None
        if s4["triad_found"] and rec["catalytic_ser"] is not None:
            s5 = analyze_cleft(pdb, rec["catalytic_ser"], cfg)
            rec["s5"] = s5
            rec["pocket_ok"] = s5["pocket_id"] is not None
        else:
            rec["pocket_ok"] = False
        per_control[c["id"]] = rec

    # score the controls that reached a catalytic pocket
    scored_ids = [cid for cid, r in per_control.items() if r.get("pocket_ok")]
    metrics_by_id = {cid: per_control[cid]["s5"]["metrics"] for cid in scored_ids}
    scores = score_controls(metrics_by_id, positive_ids, cfg)
    for cid in scored_ids:
        per_control[cid]["composite"] = scores[cid]["composite"]
        per_control[cid]["z"] = scores[cid]["z"]

    # rank (highest composite first)
    ranking = sorted(scored_ids, key=lambda i: per_control[i]["composite"], reverse=True)
    for rank, cid in enumerate(ranking, 1):
        per_control[cid]["rank"] = rank

    # separation verdict
    pos = [cid for cid in scored_ids if per_control[cid]["class"] == "positive"]
    neg = [cid for cid in scored_ids if per_control[cid]["class"] == "negative"]
    verdict = _separation(per_control, pos, neg)

    # operating point at recall = 1.0 on positives
    op = _operating_point(per_control, pos, neg, verdict)

    # leave-one-out face validity
    loo = _leave_one_out(metrics_by_id, per_control, positive_ids, neg, cfg)

    return {
        "anchor": scores["_anchor"],
        "positive_ids": positive_ids,
        "per_control": per_control,
        "ranking": ranking,
        "verdict": verdict,
        "operating_point": op,
        "loo": loo,
        "s4_detail": s4_detail,
        "trap": {cid: r for cid, r in per_control.items() if r.get("class") == "trap"},
    }


def _separation(per_control, pos, neg) -> dict:
    if not pos or not neg:
        return {"separated": None, "reason": "need >=1 positive and >=1 negative scored"}
    min_pos = min(per_control[c]["composite"] for c in pos)
    max_neg = max(per_control[c]["composite"] for c in neg)
    lowest_pos = min(pos, key=lambda c: per_control[c]["composite"])
    highest_neg = max(neg, key=lambda c: per_control[c]["composite"])
    return {
        "separated": bool(min_pos > max_neg),
        "min_positive": min_pos, "max_negative": max_neg,
        "margin": round(min_pos - max_neg, 4),
        "lowest_positive": lowest_pos, "highest_negative": highest_neg,
    }


def _operating_point(per_control, pos, neg, verdict) -> dict:
    if not pos:
        return {}
    thr = min(per_control[c]["composite"] for c in pos)  # recall=1.0 on positives
    pos_kept = [c for c in pos if per_control[c]["composite"] >= thr]
    neg_kept = [c for c in neg if per_control[c]["composite"] >= thr]
    tp, fp = len(pos_kept), len(neg_kept)
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    return {
        "threshold": thr, "recall_positives": 1.0,
        "precision": round(precision, 4),
        "true_positives": tp, "false_positives": fp,
        "negatives_above_line": sorted(neg_kept),
        "negative_positions": {c: round(per_control[c]["composite"] - thr, 4) for c in neg},
    }


def _leave_one_out(metrics_by_id, per_control, positive_ids, neg, cfg) -> list:
    out = []
    pos_scored = [p for p in positive_ids if p in metrics_by_id]
    for dropped in pos_scored:
        remaining = [p for p in positive_ids if p != dropped]
        re_scores = score_controls(metrics_by_id, remaining, cfg)
        dropped_comp = re_scores[dropped]["composite"]
        max_neg = max((re_scores[c]["composite"] for c in neg), default=float("-inf"))
        out.append({
            "dropped": dropped, "anchored_on": remaining,
            "dropped_composite": dropped_comp,
            "max_negative_composite": round(max_neg, 4),
            "still_above_negatives": bool(dropped_comp > max_neg),
        })
    return out


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_csv(result: dict, path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = ["rank", "id", "accession", "class", "triad_found", "catalytic_ser",
            "exposure", "aromatics", "druggability", "depth", "volume",
            "hydrophobicity", "polarity", "raw_og_sasa", "composite"]
    rows = []
    for cid in result["ranking"]:
        r = result["per_control"][cid]
        m = r["s5"]["metrics"]
        rows.append({
            "rank": r.get("rank"), "id": cid, "accession": r["accession"],
            "class": r["class"], "triad_found": r["triad_found"],
            "catalytic_ser": r["catalytic_ser"], "exposure": m["exposure"],
            "aromatics": int(m["aromatics"]), "druggability": m["druggability"],
            "depth": m["depth"], "volume": m["volume"],
            "hydrophobicity": m["hydrophobicity"], "polarity": m["polarity"],
            "raw_og_sasa": round(r["s5"]["raw_og_sasa"], 3), "composite": r["composite"],
        })
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def _fmt(x, nd=3):
    return f"{x:.{nd}f}" if isinstance(x, (int, float)) else str(x)


def write_report(result: dict, path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    pc = result["per_control"]
    v = result["verdict"]
    op = result["operating_point"]
    L = []
    L.append("# Proteus S4+S5 calibration & separation report")
    L.append("")
    L.append("Generated by `python -m proteus.calibrate`. S4 = catalytic-geometry gate "
             "(Ser-His-Asp triad + oxyanion hole); S5 = cleft filter (fpocket base "
             "metrics + active-site exposure + aromatic subsites), scored relative to "
             "the positive controls.")
    L.append("")
    L.append(f"**Positive anchor set:** {', '.join(result['positive_ids'])}  "
             f"(IsPETase=6EQE, LCC-WT=4EB0)")
    L.append("")

    # main table
    L.append("## Per-control results")
    L.append("")
    L.append("| rank | control | acc | class | triad | cat.Ser | exposure | aromatics "
             "| drugg. | depth | volume | composite |")
    L.append("|---:|---|---|---|:--:|---:|---:|---:|---:|---:|---:|---:|")
    rows = result["ranking"] + [cid for cid in pc
                                if cid not in result["ranking"]]
    for cid in rows:
        r = pc[cid]
        if r.get("pocket_ok"):
            m = r["s5"]["metrics"]
            L.append(f"| {r.get('rank','')} | {cid} | {r['accession']} | {r['class']} "
                     f"| {'Y' if r['triad_found'] else 'N'} | {r['catalytic_ser']} "
                     f"| {_fmt(m['exposure'],2)} | {int(m['aromatics'])} "
                     f"| {_fmt(m['druggability'],3)} | {_fmt(m['depth'],3)} "
                     f"| {_fmt(m['volume'],0)} | {_fmt(r['composite'],3)} |")
        else:
            tf = r.get("triad_found")
            tf = ("Y" if tf else "N") if tf is not None else "-"
            note = "no triad (S165A trap)" if r.get("class") == "trap" else \
                   ("no catalytic pocket" if r.get("present") else "structure absent")
            L.append(f"| - | {cid} | {r['accession']} | {r['class']} | {tf} "
                     f"| {r.get('catalytic_ser','-')} | — | — | — | — | — | {note} |")
    L.append("")

    # raw OG SASA diagnostic
    L.append("### Why exposure is peripherality, not raw OG SASA")
    L.append("")
    L.append("Raw catalytic-Ser OG SASA does **not** separate the classes — it is "
             "dominated by whether the crystal caught the lid open or closed:")
    L.append("")
    L.append("| control | class | raw OG SASA (A^2) | exposure = OG->centroid (A) |")
    L.append("|---|---|---:|---:|")
    for cid in rows:
        r = pc[cid]
        if r.get("pocket_ok"):
            L.append(f"| {cid} | {r['class']} | {_fmt(r['s5']['raw_og_sasa'],2)} "
                     f"| {_fmt(r['s5']['metrics']['exposure'],2)} |")
    L.append("")

    # separation verdict
    L.append("## Separation verdict")
    L.append("")
    if v.get("separated") is None:
        L.append(f"Cannot judge separation: {v.get('reason')}")
    elif v["separated"]:
        L.append(f"**PASS — positives separate from all negatives.** The lowest positive "
                 f"({v['lowest_positive']}, composite {_fmt(v['min_positive'],3)}) scores "
                 f"above the highest negative ({v['highest_negative']}, composite "
                 f"{_fmt(v['max_negative'],3)}).")
        L.append("")
        L.append(f"**Margin = {_fmt(v['margin'],3)}** composite units between the lowest "
                 f"positive and the highest negative.")
    else:
        L.append(f"**FAIL — STOP AND RETHINK.** The lowest positive "
                 f"({v['lowest_positive']}, {_fmt(v['min_positive'],3)}) does NOT score "
                 f"above the highest negative ({v['highest_negative']}, "
                 f"{_fmt(v['max_negative'],3)}). Margin = {_fmt(v['margin'],3)}. The cleft "
                 f"metrics as configured cannot tell PET hydrolases from these decoys. "
                 f"Do not proceed downstream until the filter or the control set is fixed.")
    L.append("")

    # operating point
    if op:
        L.append("## Provisional operating point (recall = 1.0 on positives)")
        L.append("")
        L.append(f"Threshold set at the lowest positive composite = "
                 f"**{_fmt(op['threshold'],4)}** (keeps every known PETase by construction).")
        L.append("")
        L.append(f"- Recall on positives: {op['recall_positives']:.2f} "
                 f"({op['true_positives']}/{op['true_positives']})")
        L.append(f"- Precision at this line: **{_fmt(op['precision'],3)}** "
                 f"({op['true_positives']} TP / {op['true_positives'] + op['false_positives']} kept)")
        L.append(f"- Negatives above the line (false positives): "
                 f"{op['negatives_above_line'] or 'none'}")
        L.append("")
        L.append("Each negative's composite relative to the line (negative = correctly "
                 "below the threshold):")
        L.append("")
        L.append("| negative | composite - threshold |")
        L.append("|---|---:|")
        for cid, delta in sorted(op["negative_positions"].items(),
                                 key=lambda kv: kv[1], reverse=True):
            L.append(f"| {cid} | {_fmt(delta,3)} |")
        L.append("")

    # leave-one-out
    L.append("## Leave-one-out face validity")
    L.append("")
    L.append("Drop each positive, re-anchor on the remaining positive(s), and check the "
             "dropped positive still scores above every negative:")
    L.append("")
    L.append("| dropped positive | anchored on | dropped composite | max negative | still above? |")
    L.append("|---|---|---:|---:|:--:|")
    for lo in result["loo"]:
        L.append(f"| {lo['dropped']} | {', '.join(lo['anchored_on']) or '(none)'} "
                 f"| {_fmt(lo['dropped_composite'],3)} | {_fmt(lo['max_negative_composite'],3)} "
                 f"| {'YES' if lo['still_above_negatives'] else 'NO'} |")
    L.append("")

    # trap
    L.append("## 6THS S165A trap (expected null)")
    L.append("")
    for cid, r in result["trap"].items():
        tf = r.get("triad_found")
        L.append(f"- **{cid} ({r['accession']})**: triad_found = {tf} "
                 f"-> {'EXPECTED NULL OK (catalytic Ser is mutated to Ala; no triad to find)' if tf is False else 'UNEXPECTED — a triad was detected in the inactivated mutant'}")
    L.append("")

    # honest stats
    L.append("## Honest statistics — read this before trusting the threshold")
    L.append("")
    L.append("- **N is tiny:** 2 distinct positive scaffolds (IsPETase, LCC) and 4 "
             "negatives (CalB, AChE, *C. rugosa* lipase, Est2). This is a FACE-VALIDITY "
             "check and a PROVISIONAL operating point — **not** a trustworthy threshold.")
    L.append("- With only 2 positives the per-metric standard deviation is unreliable, so "
             "z-score scales are floored by half the overall control spread "
             "(`_robust_scale`). Treat z-magnitudes as indicative, not calibrated.")
    L.append("- The separation here is carried by **exposure (catalytic-Ser "
             "peripherality)**; druggability and volume separate the big-pocket lipases "
             "but miss the compact esterase (Est2). Raw OG SASA and raw aromatic counts "
             "do **not** separate on this panel (aromatic-gorge esterases inflate the "
             "aromatic count), which is why aromatics is down-weighted.")
    L.append("- Metric orientations and weights were set from biology a priori, not tuned "
             "to maximise separation on these 6 structures; still, a 6-point fit cannot "
             "be considered validated.")
    L.append("- **Next step:** fold the divergent positives GuaPA and MG8 on the Vast.ai "
             "burst box (S3) to add real, sequence-divergent PETases to the positive set, "
             "expand the decoy panel (more lipase/esterase/cutinase-adjacent folds), then "
             "re-calibrate. Re-check the aromatic-subsite metric on the *folded* models "
             "(rotamer noise) after Chai-1 refinement.")
    L.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"))
    ap.add_argument("--report", default=os.path.join(REPO, "envlog", "calibration-report.md"))
    ap.add_argument("--csv", default=os.path.join(REPO, "data", "processed", "s5_scores.csv"))
    ap.add_argument("--s4-json", default=os.path.join(REPO, "data", "processed", "s4_triads.json"))
    ap.add_argument("--check-only", action="store_true",
                    help="run only the Checkpoint 0 precondition audit and exit")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)

    # Checkpoint 0 — enumerate preconditions; STOP if any fail.
    audit = preconditions(cfg, args.struct_dir)
    print_precondition_report(audit)
    if not audit["ok"]:
        print("Aborting: fix the FAILED preconditions above, then re-run.")
        return 2
    if args.check_only:
        return 0
    print()

    result = run_calibration(cfg, args.struct_dir)

    os.makedirs(os.path.dirname(os.path.abspath(args.s4_json)), exist_ok=True)
    with open(args.s4_json, "w") as fh:
        json.dump(result["s4_detail"], fh, indent=2)
        fh.write("\n")
    write_csv(result, args.csv)
    write_report(result, args.report)

    v = result["verdict"]
    print("=== S4 -> S5 CALIBRATION ===")
    for cid in result["ranking"]:
        r = result["per_control"][cid]
        print(f"  #{r['rank']} {cid:9s} [{r['class']:8s}] composite={r['composite']:+.3f} "
              f"exposure={r['s5']['metrics']['exposure']:.2f}")
    if v.get("separated") is None:
        print(f"SEPARATION: indeterminate ({v.get('reason')})")
    else:
        print(f"SEPARATION: {'PASS' if v['separated'] else 'FAIL — STOP AND RETHINK'} "
              f"(margin={v['margin']:+.3f}; lowest pos {v['lowest_positive']}="
              f"{v['min_positive']:.3f} vs highest neg {v['highest_negative']}="
              f"{v['max_negative']:.3f})")
    if result["operating_point"]:
        op = result["operating_point"]
        print(f"OPERATING POINT: thr={op['threshold']:.3f} recall=1.0 "
              f"precision={op['precision']:.3f} fp={op['false_positives']}")
    print(f"LOO: " + "; ".join(
        f"{lo['dropped']}->{'ok' if lo['still_above_negatives'] else 'FAIL'}"
        for lo in result["loo"]))
    print(f"report -> {os.path.relpath(args.report, os.getcwd())}")
    print(f"csv    -> {os.path.relpath(args.csv, os.getcwd())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
