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
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "paper_figures")

# Paper figure <- source result figure (stem, no extension). Figures 1, 7, A and
# B are generated directly into paper_figures/ by their own scripts; the rest are
# produced by the per-tool analysis pipelines and assembled here (previously a
# manual copy step). All source titles are CP-codename-free (publication style).
FIG_MAP = {
    "fig2": "absorption_results/tau_old_vs_new",
    "fig3": "absorption_results/k_abs_vs_N",
    "fig4": "absorption_results/geometric_p",
    "fig5": "absorption_results/absorption_phase_rose",
    "fig6": "transient_results/transient_decider",
    "fig9": "reduced_results/cp4_scatter",
}


def run(cmd, **kw):
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=ROOT, **kw)


def assemble_figs():
    """Copy the per-tool source figures into paper_figures/figN.{pdf,png}.
    Scripted (not manual) so stale copies cannot drift from their sources."""
    print("\n# assembling paper figures from source pipelines")
    for fig, stem in FIG_MAP.items():
        for ext in ("pdf", "png"):
            src = os.path.join(ROOT, f"{stem}.{ext}")
            dst = os.path.join(OUT, f"{fig}.{ext}")
            if os.path.exists(src):
                shutil.copyfile(src, dst)
                print(f"  {fig}.{ext} <- {stem}.{ext}")
            else:
                print(f"  [warn] missing source {stem}.{ext} (run its analysis pipeline first)")


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

    # Appendix figures: beta=0.10 slice (A) and finite-size noise test (B).
    run([sys.executable, os.path.join(HERE, "figA.py")])
    run([sys.executable, os.path.join(HERE, "figB.py")])

    # Assemble Figures 2-6, 9 from their source pipelines (scripted copy step).
    assemble_figs()

    run([sys.executable, os.path.join(HERE, "make_report.py")])
    print("\nDone. Artifacts in paper_figures/.")


if __name__ == "__main__":
    main()
