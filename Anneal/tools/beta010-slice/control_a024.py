#!/usr/bin/env python3
"""
A=0.24 mid-band stable-chimera control (beta=0.10) — deferred item 3.

Reads the control campaign produced by
  node tools/absorption-recampaign/sweep.mjs \
       --config tools/absorption-recampaign/absorption.beta010.A024.config.json
and summarizes the never-absorb (persistent) fraction and the absorbing-minority
statistics per N.

Purpose: A=0.24 sits well inside the stable stationary-chimera band
(A_SN=0.178 < 0.24 < A_H=0.278), further from the saddle-node than the A=0.2
secondary control. The A=0.2 control's persistent fraction swings widely with N
(28% at N=8 -> 60% at N=64), which a referee could read as N-dependent persistence
rather than a genuine attractor. This control checks whether a more central
mid-band point gives an N-stable persistent fraction.

Writes absorption_results/beta010_A024_control.json.
"""
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "absorption_results" / "absorption_campaign_beta010_A024.jsonl"
A2 = ROOT / "absorption_results" / "absorption_campaign_beta010.jsonl"
OUT = ROOT / "absorption_results" / "beta010_A024_control.json"


def per_N(rows, A):
    out = {}
    for N in sorted({r["N"] for r in rows if abs(r["A"] - A) < 1e-9}):
        rs = [r for r in rows if r["N"] == N and abs(r["A"] - A) < 1e-9]
        absd = sorted(r["t_abs"] for r in rs if not r["abs_censored"])
        cens = sum(1 for r in rs if r["abs_censored"])
        out[N] = {
            "n": len(rs),
            "persistent_n": cens,
            "persistent_frac": cens / len(rs),
            "absorbers_n": len(absd),
            "t_abs_median": (st.median(absd) if absd else None),
            "t_abs_min": (absd[0] if absd else None),
            "t_abs_max": (absd[-1] if absd else None),
        }
    return out


def main():
    rows = [json.loads(l) for l in RUN.read_text().splitlines() if l.strip()]
    summary = {
        "campaign": "beta010-A024-v1",
        "beta": 0.10,
        "A": 0.24,
        "band": "A_SN=0.178 < 0.24 < A_H=0.278 (stable stationary-chimera band)",
        "by_N": per_N(rows, 0.24),
    }
    # direct comparison to the marginal A=0.2 secondary control, if present
    if A2.exists():
        a2rows = [json.loads(l) for l in A2.read_text().splitlines() if l.strip()]
        summary["compare_A0.2"] = per_N(a2rows, 0.2)

    OUT.write_text(json.dumps(summary, indent=2))

    def fr(d, N):
        return f"{d[N]['persistent_frac']*100:.0f}% (n={d[N]['n']})" if N in d else "—"

    print(f"A=0.24 mid-band control (beta=0.10): {len(rows)} runs")
    a024 = summary["by_N"]
    for N in sorted(a024):
        s = a024[N]
        tt = (
            f"median t_abs={s['t_abs_median']:.1f}s [{s['t_abs_min']:.1f}, {s['t_abs_max']:.1f}]"
            if s["absorbers_n"]
            else "no absorbers"
        )
        print(
            f"  N={N}: persistent {s['persistent_frac']*100:.0f}% "
            f"({s['persistent_n']}/{s['n']}); absorbers {s['absorbers_n']} ({tt})"
        )
    if "compare_A0.2" in summary:
        a2 = summary["compare_A0.2"]
        print("  vs A=0.20 control:  N=8 " + fr(a2, 8) + "  N=64 " + fr(a2, 64))
        print("     A=0.24        :  N=8 " + fr(a024, 8) + "  N=64 " + fr(a024, 64))
    print(f"\nWrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
