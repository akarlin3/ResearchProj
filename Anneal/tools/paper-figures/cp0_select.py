"""
CP0 — light audit + deterministic Fig. 1 run selection.

Locates the inputs the two figures compose from (read-only), lists the three
best Fig. 1 candidates (N=16, A=0.5, 1–2 recovered grazes before absorption,
total span ≲300 s), scores them for visual clarity, picks the clearest, and
asserts the pick equals the seed pinned in figures.config.json (so the whole
pipeline is reproducible from config). Writes paper_figures/cp0_audit.json.

Inputs (read-only):
  absorption_results/phase_traces.jsonl   — per-run R_incoh=min(R1,R2) + labels
  absorption_results/absorption_campaign.jsonl — per-run labels (t_graze, t_abs,
                                                  n_grazes_before_abs, T_b)
Machinery referenced:
  tools/absorption-recampaign/breath.mjs       — breath-peak detector (T_b)
  tools/absorption-recampaign/labeling.mjs     — t_graze / t_abs labeler
  tools/reduced-ode/reduced_core.py            — A_SN_series/A_H_series (Eqs 17/18),
                                                 locate_saddle_node/locate_hopf,
                                                 locate_homoclinic (bisection)
The full per-population R1/R2 series is not stored (phase_traces carries the
headline min); Fig. 1 regenerates R1, R2 deterministically from the seed via
fig1_trace.mjs (shipped-identical RK4).
"""
from __future__ import annotations

import json
import os

import style

ROOT = style.ROOT
CFG = style.load_config()
OUT = os.path.join(ROOT, CFG["output_dir"])
SEL = CFG["fig1"]["selection"]


def clarity(c):
    """Deterministic visual-clarity score: clear recovery dip, graze/absorption
    separation, visible breathing, moderate span, room for breaths before graze."""
    s = 0.0
    s += (0.80 - c["rec_depth"]) * 3.0 if c["rec_depth"] < 0.80 else -5.0
    s += min(c["sep"], 150) / 150 * 2.0
    s += min(c["nbp"], 8) / 8 * 2.0
    s += 1.5 if 100 <= c["span"] <= 240 else 0.0
    s += min(c["pre_graze"], 90) / 90 * 1.5
    s += 0.5 if c["ng"] == 1 else 0.0
    return s


def score_run(d):
    if d.get("abs_censored") or d.get("t_abs") is None:
        return None
    ng = d["n_grazes_before_abs"]
    if ng not in SEL["allowed_n_grazes"]:
        return None
    if d["t_abs"] > SEL["max_span_s"]:
        return None
    R = d["R_incoh"]
    gi, ai = d["grazeIndex"], d["absIndex"]
    if gi < 0 or ai < 0 or ai <= gi:
        return None
    rec_depth = min(R[gi:ai])
    return {
        "seed": d["seed"], "t_graze": d["t_graze"], "t_abs": d["t_abs"],
        "ng": ng, "Tb": d["T_b"], "nbp": d["n_breath_peaks"],
        "sep": round(d["t_abs"] - d["t_graze"], 1),
        "rec_depth": round(float(rec_depth), 3),
        "span": round(d["t_abs"], 1), "pre_graze": round(d["t_graze"], 1),
    }


def main():
    N = CFG["fig1"]["N"]
    src = os.path.join(ROOT, SEL["trace_source"])
    cands = []
    with open(src) as f:
        for line in f:
            d = json.loads(line)
            if d["N"] == N and d["A"] == CFG["fig1"]["A"]:
                sc = score_run(d)
                if sc:
                    cands.append(sc)
    for c in cands:
        c["score"] = round(clarity(c), 2)
    cands.sort(key=lambda c: c["score"], reverse=True)

    print("CP0 — light audit")
    print("  inputs (read-only):")
    print("    phase_traces.jsonl  (R_incoh + labels)")
    print("    absorption_campaign.jsonl  (labels: t_graze, t_abs, n_grazes, T_b)")
    print("  breath-peak detector:  tools/absorption-recampaign/breath.mjs")
    print("  reduced-ODE bifurcation: reduced_core.A_SN_series/A_H_series, "
          "locate_saddle_node/locate_hopf, locate_homoclinic")
    print(f"\n  {len(cands)} candidates (N={N}, A={CFG['fig1']['A']}, "
          f"1–2 grazes, span≤{SEL['max_span_s']:.0f}s). Top 3:")
    top3 = cands[:3]
    for c in top3:
        print(f"    seed={c['seed']}  t_graze={c['t_graze']}  t_abs={c['t_abs']}  "
              f"n_grazes={c['ng']}  T_b={c['Tb']}  sep={c['sep']}s  "
              f"reform_dip={c['rec_depth']}  span={c['span']}s  "
              f"score={clarity(c):.2f}")

    chosen = top3[0]
    pinned = CFG["fig1"]["seed"]
    ok = chosen["seed"] == pinned
    print(f"\n  CHOICE: seed {chosen['seed']} (N={N}, A={CFG['fig1']['A']}) — "
          f"clearest: well-separated graze ({chosen['t_graze']}s) and absorption "
          f"({chosen['t_abs']}s), deep reform dip to {chosen['rec_depth']}, "
          f"{chosen['nbp']} visible breaths (T_b={chosen['Tb']}s).")
    print(f"  pinned in figures.config.json: seed {pinned} → "
          f"{'MATCH' if ok else 'MISMATCH'}")
    assert ok, f"selection {chosen['seed']} != pinned {pinned}"

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "cp0_audit.json"), "w") as f:
        json.dump({"N": N, "A": CFG["fig1"]["A"], "n_candidates": len(cands),
                   "top3": top3, "chosen_seed": chosen["seed"],
                   "pinned_seed": pinned, "match": ok}, f, indent=1)
    print("  wrote", os.path.relpath(os.path.join(OUT, "cp0_audit.json"), ROOT))


if __name__ == "__main__":
    main()
