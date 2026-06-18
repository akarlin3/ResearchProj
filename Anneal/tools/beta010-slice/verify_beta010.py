"""CP1 verification — spot-check two beta=0.10 trajectories against their labels.

Picks (from the stored phase traces) one absorbed run with >=1 recovered graze
and one direct absorber, plots min(R1,R2) with the theta/recovery thresholds and
t_graze/t_abs marks, and asserts the labels by eye-equivalent checks:
  - recovered-graze run: n_grazes_before_abs > 0 AND R_incoh dips < 0.80 after the
    first theta-crossing and before t_abs (a real recovery), then a final
    no-recovery crossing at t_abs.
  - direct absorber: a sustained theta-crossing leads to absorption.

Writes absorption_results/beta010_verify_traces.png and prints PASS/FAIL.
Run: python3 tools/beta010-slice/verify_beta010.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
TRACES = ROOT / "absorption_results/phase_traces_beta010.jsonl"
OUT = ROOT / "absorption_results"
THETA, REC = 0.85, 0.80

rows = [json.loads(s) for s in TRACES.read_text().splitlines() if s.strip()]
absorbed = [r for r in rows if not r["abs_censored"] and "R_incoh" in r and r["A"] == 0.5]

# recovered-graze example: highest n_grazes_before_abs; direct absorber: n=0.
recov = max((r for r in absorbed if r["n_grazes_before_abs"] >= 1),
            key=lambda r: r["n_grazes_before_abs"])
direct = next(r for r in absorbed if r["n_grazes_before_abs"] == 0)

failures = []


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        failures.append(name)


fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, r, kind in [(axes[0], recov, "recovered graze"), (axes[1], direct, "direct absorber")]:
    dt = r["sampleDt"]
    y = np.asarray(r["R_incoh"], float)
    t = np.arange(len(y)) * dt
    ax.plot(t, y, lw=0.8)
    ax.axhline(THETA, color="orange", ls="--", lw=0.8, label="θ=0.85")
    ax.axhline(REC, color="green", ls=":", lw=0.8, label="recovery=0.80")
    ax.axvline(r["t_graze"], color="gray", ls="-", lw=0.6)
    ax.axvline(r["t_abs"], color="red", ls="-", lw=0.8)
    ax.set_title(f"{kind}: N={r['N']} seed={r['seed']}\n"
                 f"n_graze={r['n_grazes_before_abs']} t_graze={r['t_graze']:.1f} t_abs={r['t_abs']:.1f}")
    ax.set_xlabel("t (s)"); ax.set_ylabel("min(R1,R2)")
    ax.set_xlim(0, min(r["t_abs"] + 30, t[-1])); ax.legend(fontsize=7)

    print(f"\n{kind} (N={r['N']} seed={r['seed']}):")
    abs_idx = r["absIndex"]
    pre = y[:abs_idx]
    crossed = np.where(pre > THETA)[0]
    if kind == "recovered graze":
        check("n_grazes_before_abs > 0", r["n_grazes_before_abs"] > 0)
        # a real recovery: after the FIRST crossing, R dips below 0.80 before t_abs
        first_cross = crossed[0] if len(crossed) else len(pre)
        dip_after = pre[first_cross:].min() if first_cross < len(pre) else 1.0
        check("recovery dip < 0.80 after first graze, before t_abs", dip_after < REC)
        check("t_abs > t_graze (absorption follows the recovered graze)", r["t_abs"] > r["t_graze"])
    else:
        check("direct absorber: n_grazes_before_abs == 0", r["n_grazes_before_abs"] == 0)
        check("R_incoh reaches θ at/by absorption", y[max(0, abs_idx - 1):abs_idx + 50].max() > THETA)

fig.tight_layout()
fig.savefig(OUT / "beta010_verify_traces.png", dpi=130)
print(f"\nWrote {OUT / 'beta010_verify_traces.png'}")
print("\nVERIFICATION:", "ALL PASS" if not failures else f"FAILURES: {failures}")
sys.exit(1 if failures else 0)
