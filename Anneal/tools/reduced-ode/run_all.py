"""Reproduce the entire reduced-ODE study from config, in order.

  python3 tools/reduced-ode/run_all.py

Steps (deterministic): core self-tests → CP1 gate (STOPS the pipeline if it
fails) → CP2 corners+portraits → A=0.2 finite-N level (node) → compute reduced
runs → CP3 matches → CP4 capstone → assemble REDUCED_REPORT.md.
"""
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, check=True):
    print(f"\n$ {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=HERE)
    if check and r.returncode != 0:
        print(f"STEP FAILED (exit {r.returncode}) — stopping.")
        sys.exit(r.returncode)
    return r.returncode


def main():
    PY = sys.executable
    run([PY, "test_core.py"])
    # CP1 is the gate: a non-zero exit means STOP and await review.
    if run([PY, "cp1_validate.py"], check=False) != 0:
        print("CP1 GATE FAILED — not proceeding to CP2+. Await human review.")
        sys.exit(1)
    run([PY, "cp2_classify.py"])
    # finite-N A=0.2 level uses the shipped node tracer
    out = os.path.join(HERE, "..", "..", "reduced_results", "cp3_a02_level.json")
    with open(out, "w") as f:
        r = subprocess.run(["node", "trace_a02_level.mjs"], cwd=HERE, stdout=f)
    if r.returncode != 0:
        print("node trace_a02_level.mjs failed — stopping.")
        sys.exit(r.returncode)
    run([PY, "compute_runs.py"])
    run([PY, "cp3_match.py"])
    run([PY, "cp4_predict.py"])
    run([PY, "make_report.py"])
    print("\nDONE — see reduced_results/REDUCED_REPORT.md")


if __name__ == "__main__":
    main()
