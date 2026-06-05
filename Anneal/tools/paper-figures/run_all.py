"""
Driver — regenerate Figures 1 & 7 and the report, deterministically, from
figures.config.json.

Steps:
  CP0  cp0_select.py            audit + pin-checked Fig.1 run selection
  F1   fig1_trace.mjs (node)    regenerate the example trace from the seed
       fig1.py                  plot Figure 1 + caption
  F7   fig7_curves.py           SN/Hopf series+numeric, homoclinic bisection (slow)
       fig7.py                  plot Figure 7 + caption
  RPT  make_report.py           FIGURES_REPORT.md

Usage:
  python3 tools/paper-figures/run_all.py [--skip-curves]
    --skip-curves   reuse cached paper_figures/fig7_curves.json (skip the ~4-5
                    min homoclinic trace) when only plots/report changed.
"""
from __future__ import annotations

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "paper_figures")


def run(cmd, **kw):
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=ROOT, **kw)


def main():
    skip_curves = "--skip-curves" in sys.argv
    os.makedirs(OUT, exist_ok=True)

    run([sys.executable, os.path.join(HERE, "cp0_select.py")])

    # Figure 1: deterministic trace (node) → plot
    with open(os.path.join(OUT, "fig1_trace.json"), "w") as f:
        run(["node", os.path.join(HERE, "fig1_trace.mjs")], stdout=f)
    run([sys.executable, os.path.join(HERE, "fig1.py")])

    # Figure 7: curves (slow) → plot
    cache = os.path.join(OUT, "fig7_curves.json")
    if skip_curves and os.path.exists(cache):
        print("\n[skip-curves] reusing cached fig7_curves.json")
    else:
        run([sys.executable, os.path.join(HERE, "fig7_curves.py")])
    run([sys.executable, os.path.join(HERE, "fig7.py")])

    run([sys.executable, os.path.join(HERE, "make_report.py")])
    print("\nDone. Artifacts in paper_figures/.")


if __name__ == "__main__":
    main()
