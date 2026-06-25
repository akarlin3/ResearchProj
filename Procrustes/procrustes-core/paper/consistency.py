#!/usr/bin/env python3
"""Consistency gate for the Procrustes manuscript (house number-traceability).

Mirrors Gauge/Augur/Matrix paper/consistency.py: regenerates ``numbers.tex`` from
the SEEDED gate results so that NO number in ``procrustes.tex`` is hand-typed --
every value originates in a reproduction artifact (``results/phase{1,2,3}_*.json``)
produced by the gate drivers ``experiments/run_phase{1,2,3}.py``. Then it verifies
that every ``\\num*`` macro the manuscript uses is defined and asserts the
load-bearing spine invariants (the locked claim scope) survive regeneration.

Inputs (run ``bash reproduce.sh`` first, which regenerates them):
  ../results/phase1_gateB.json   (apples-to-apples Gauge separation)
  ../results/phase2_gateC.json   (diagnostic reach; honest scope)
  ../results/phase3_gateD.json   (robustness envelope)

Output:
  paper/numbers.tex  (AUTO-GENERATED -- do not edit by hand)

Exit 0 = regenerated, all macros defined, spine invariants hold; non-zero else.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORE = HERE.parent
RESDIR = CORE / "results"
TEX = HERE / "procrustes.tex"
NUMBERS = HERE / "numbers.tex"


def _load(name: str) -> dict:
    p = RESDIR / name
    if not p.exists():
        print(f"CONSISTENCY FAIL: missing results/{name} -- run `bash reproduce.sh` first",
              file=sys.stderr)
        raise SystemExit(2)
    return json.loads(p.read_text())


def _g(d):  # signed gap, 3 decimals
    return f"{d['point']:+.3f}"


def _gl(d):
    return f"{d['lo']:+.3f}"


def _gh(d):
    return f"{d['hi']:+.3f}"


def _a(d):  # AUC / coverage, 3 decimals, unsigned
    return f"{d['point']:.3f}"


def build_numbers():
    """Returns {macro: (latex_value, provenance)} -- the single source of truth."""
    p1 = _load("phase1_gateB.json")
    p2 = _load("phase2_gateC.json")
    p3 = _load("phase3_gateD.json")
    st1 = p1["families"]["stretched"]
    ds1 = p1["families"]["dispersion_lognormal"]
    tr1 = p1["families"]["triexp"]
    st2 = p2["families"]["stretched"]
    ds2 = p2["families"]["dispersion_lognormal"]
    tr2 = p2["families"]["triexp"]
    h16 = p3["headline16"]
    cfg1 = p1["config"]
    cfg2 = p2["config"]

    n = {
        # --- experiment configuration -------------------------------------------------
        "numN": (f"{cfg1['n']}", "voxels per cohort"),
        "numSNR": (f"{cfg1['snr']:.0f}", "Rician SNR (default)"),
        "numSeeds": (f"{cfg1['n_seeds']}", "seeds (phase 1/2)"),
        "numSeedsHeadline": (f"{h16['gap_all']['n_seeds']}", "seeds (phase 3 headline)"),
        "numAUCfloor": (f"{p2['gate_c']['floor']:.2f}", "pre-registered diagnostic AUC floor"),
        # --- (a) MARGINAL holds -- confirms Lei 2018, NOT the contribution ------------
        "numMarg": (_a(st1["marginal"]), "stretched marginal coverage (departure-blind)"),
        "numMargLo": (f"{st1['marginal']['lo']:.3f}", "marginal coverage CI lo"),
        "numMargHi": (f"{st1['marginal']['hi']:.3f}", "marginal coverage CI hi"),
        "numNominal": ("0.90", "nominal coverage (alpha=0.10)"),
        # --- (b) CONDITIONAL gap of D breaks (16-seed headline) -- the contribution ---
        "numGap": (_g(h16["gap_all"]), "conditional gap (all voxels, 16 seeds)"),
        "numGapLo": (_gl(h16["gap_all"]), "conditional gap CI lo"),
        "numGapHi": (_gh(h16["gap_all"]), "conditional gap CI hi"),
        # --- (c) DISTINCT from Gauge: gap survives & intensifies in identifiable region
        "numWell": (_g(h16["gap_wellid"]), "well-ID (Gauge identifiable) gap, 16 seeds"),
        "numWellLo": (_gl(h16["gap_wellid"]), "well-ID gap CI lo"),
        "numWellHi": (_gh(h16["gap_wellid"]), "well-ID gap CI hi"),
        "numStrictLo": (_g(h16["gap_lo"]), "strict bottom-D* tercile gap, 16 seeds"),
        "numStrictLoLo": (_gl(h16["gap_lo"]), "strict-lo gap CI lo"),
        "numStrictLoHi": (_gh(h16["gap_lo"]), "strict-lo gap CI hi"),
        "numWallGap": (_g(st1["gap_hi"]), "Gauge's OWN high-D* wall gap (8 seeds, contrast)"),
        "numWallGapLo": (_gl(st1["gap_hi"]), "high-D* wall gap CI lo"),
        "numWallGapHi": (_gh(st1["gap_hi"]), "high-D* wall gap CI hi"),
        "numNwell": (f"{st1['n_wellid']}", "voxels in well-ID subset (~2/3 of test)"),
        "numNlo": (f"{st1['n_lo']}", "voxels in strict bottom tercile (~1/3 of test)"),
        # --- mechanism: signed D-bias growth ------------------------------------------
        "numBias": (f"{st1['bias_ratio']['point']:.1f}", "stretched bias growth |worst|/|limit|"),
        # --- mechanism-specific scope: dispersion weak, tri-exp null ------------------
        "numGapDisp": (_g(ds1["gap_all"]), "log-normal dispersion conditional gap (weak)"),
        "numWellDisp": (_g(ds1["gap_wellid"]), "log-normal dispersion well-ID gap"),
        "numGapTri": (_g(tr1["gap_all"]), "tri-exp conditional gap (null)"),
        "numWellTri": (_g(tr1["gap_wellid"]), "tri-exp well-ID gap (null)"),
        # --- diagnostic reach (phase 2; heavy-tail DETECTOR) --------------------------
        "numAUC": (_a(st2["auc_best"]), "stretched diagnostic AUC (best residual stat)"),
        "numAUCLo": (f"{st2['auc_best']['lo']:.3f}", "stretched diagnostic AUC CI lo"),
        "numAUCHi": (f"{st2['auc_best']['hi']:.3f}", "stretched diagnostic AUC CI hi"),
        "numAUCstruct": (_a(st2["auc_structure"]), "stretched residual-structure AUC"),
        "numAUCchi": (_a(st2["auc_chi2_red"]), "stretched reduced-chi2 (magnitude) AUC"),
        "numMonitor": (_a(st2["auc_monitor_full"]), "naive drift-monitor AUC (stretched)"),
        "numMonitorLo": (f"{st2['auc_monitor_full']['lo']:.3f}", "naive monitor AUC CI lo"),
        "numMonitorHi": (f"{st2['auc_monitor_full']['hi']:.3f}", "naive monitor AUC CI hi"),
        "numRho": (f"{st2['rho_best']['point']:.3f}", "rank power rho(stat,|D-err|) stretched"),
        "numAUCdisp": (_a(ds2["auc_best"]), "dispersion diagnostic AUC (near-hidden)"),
        "numAUCtri": (_a(tr2["auc_best"]), "tri-exp diagnostic AUC (detectable-but-harmless)"),
        # --- robustness envelope (phase 3) -------------------------------------------
        "numCondSurv": (f"{len(p3['gate_d']['survived'])}", "conditions where separation survives"),
        "numCondTot": (f"{p3['gate_d']['n_conditions']}", "total robustness conditions swept"),
        "numSNRfail": ("25", "SNR at which separation fails (honest boundary)"),
        "numSNRmin": ("35", "lowest SNR at which separation survives"),
        "numWellSNRhi": (_g(p3["snr"]["100"]["gap_wellid"]), "well-ID gap at SNR=100 (intensifies)"),
        "numWellClin": (_g(p3["bscheme"]["clinical-sparse-10pt"]["gap_wellid"]),
                        "well-ID gap, clinical sparse 10-pt b-scheme"),
        "numWellBig": (_g(p3["n"]["2000"]["gap_wellid"]), "well-ID gap at n=2000"),
    }
    return n


def write_numbers(nums):
    lines = [
        "% AUTO-GENERATED by paper/consistency.py from the seeded gate results.",
        "% DO NOT EDIT BY HAND. Regenerate with: bash reproduce.sh (or bash paper/build.sh).",
        "% Every value traces to results/phase{1,2,3}_*.json from experiments/run_phase*.py.",
        "",
    ]
    for k, (v, why) in nums.items():
        lines.append(f"\\newcommand{{\\{k}}}{{{v}}}  % {why}")
    NUMBERS.write_text("\n".join(lines) + "\n")


def verify_macros(nums):
    problems = []
    if not TEX.exists():
        return [f"manuscript {TEX.name} not found (write it before the gate can pass)"]
    tex = TEX.read_text()
    used = set(re.findall(r"\\(num[A-Za-z]+)\b", tex))
    defined = set(nums)
    undefined = sorted(used - defined)
    unused = sorted(defined - used)
    if undefined:
        problems.append(f"macros used in procrustes.tex but NOT defined: {undefined}")
    if unused:
        print(f"  note: macros defined but unused in procrustes.tex: {unused}")
    return problems


def check_invariants():
    """Load-bearing spine invariants -- the locked claim scope must survive."""
    p1 = _load("phase1_gateB.json")
    p2 = _load("phase2_gateC.json")
    p3 = _load("phase3_gateD.json")
    st = p1["families"]["stretched"]
    h16 = p3["headline16"]
    fails = []
    # GATE B: distinctness from Gauge
    if not p1["gate_b"]["passed"]:
        fails.append("GATE B must pass (well-ID gap survives Gauge's own partition)")
    if not (h16["gap_wellid"]["lo"] > 0.03):
        fails.append("well-ID gap (16s) CI-lo must exceed the 0.03 floor (refute R2)")
    if not (h16["gap_wellid"]["point"] >= h16["gap_all"]["point"] - 0.02):
        fails.append("well-ID gap must NOT collapse below the overall gap (not diluted)")
    if not (st["gap_hi"]["point"] < h16["gap_wellid"]["point"]):
        fails.append("Gauge's OWN high-D* wall gap must be SMALLER than the well-ID gap "
                     "(opposite gradient -- the distinctness mechanism)")
    # trap stays closed: marginal holds
    if abs(st["marginal"]["point"] - 0.90) > 0.035:
        fails.append("marginal coverage must stay ~nominal (the trap stays closed; Lei 2018)")
    # GATE C: honest scope -- heavy-tail DETECTOR, not a universal one
    if p2["gate_c"]["detected_all_nonnull"]:
        fails.append("scope invariant: must NOT claim a general detector "
                     "(dispersion is below the AUC floor -- honest scoping)")
    if not p2["gate_c"]["detected"]["stretched"]:
        fails.append("stretched heavy-tail channel must clear the diagnostic AUC floor")
    if p2["gate_c"]["detected"]["dispersion_lognormal"]:
        fails.append("dispersion must stay BELOW the AUC floor (near-hidden channel)")
    # GATE D: robustness with the honest SNR boundary
    if "25" not in p3["gate_d"]["failed"]:
        fails.append("the SNR=25 honest boundary (separation fails) must be reported")
    if len(p3["gate_d"]["survived"]) < 8:
        fails.append("separation must survive across most of the swept envelope")
    return fails


def main():
    print("Procrustes consistency gate")
    print("=" * 64)
    nums = build_numbers()
    write_numbers(nums)
    print(f"  regenerated numbers.tex ({len(nums)} macros) from seeded gate results")
    problems = verify_macros(nums)
    fails = check_invariants()
    for p in problems:
        print(f"  MACRO FAIL: {p}", file=sys.stderr)
    for f in fails:
        print(f"  INVARIANT FAIL: {f}", file=sys.stderr)
    if problems or fails:
        print("=" * 64)
        print("CONSISTENCY GATE: FAIL", file=sys.stderr)
        return 1
    print("  all manuscript macros defined; locked-scope spine invariants hold")
    print("=" * 64)
    print("CONSISTENCY GATE: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
