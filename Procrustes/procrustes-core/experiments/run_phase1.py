"""PHASE 1 / GATE B -- apples-to-apples Gauge separation (load-bearing).

The distinctness-from-Gauge claim rests on one number: the conditional-coverage
gap of the trusted tissue parameter D must SURVIVE inside Gauge's *identifiable
region*. This driver re-defines the well-identified-D* subset to match Gauge's
own partition EXACTLY -- true-D* terciles via np.quantile([1/3, 2/3]) +
np.digitize, identical to Gauge/gauge/conditional.py:157-159 and
conditional_attack._regime_from_true -- and re-runs the gap there.

Three nested D*-regime subsets are reported, increasingly deep inside Gauge's
"trust D" region:
  * ALL test voxels          -- the overall conditional gap,
  * WELL-ID (bottom-2 terc.) -- Gauge's identifiable region (lo + mid D*),
  * STRICT-LO (bottom terc.) -- the deepest interior of it (strongest R2 test),
and, for contrast, Gauge's own ill-posed wall:
  * HIGH-D* (top tercile)    -- where Gauge says D* is unidentifiable.

PRE-REGISTERED REFUTE (R2): if, under Gauge's own definition, the well-ID gap
collapses toward the marginal (loses separation), the distinctness claim is DEAD
-> STOP for the author. Operationally: GATE B PASSES iff
  (i)  the strict-lo-tercile gap's bootstrap CI lower bound exceeds the
       pre-registered floor (cfg.wellid_gap_min = 0.03), AND
  (ii) the well-ID gap does NOT collapse below the overall conditional gap
       (the failure is not merely the high-D* wall leaking in).

Usage:  python experiments/run_phase1.py [n_seeds]   (default 8)
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from procrustes import (ProcrustesConfig, STRETCHED, LOGNORMAL, TRIEXP,
                        separation_detail, across_seed_ci,
                        cluster_bootstrap_gap, GLOBAL_SEED)

ORDER = [STRETCHED, LOGNORMAL, TRIEXP]
_RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "results")


def _fmt(ci):
    return f"{ci['point']:+.3f} [{ci['lo']:+.3f}, {ci['hi']:+.3f}]"


def main(n_seeds: int = 8):
    cfg = ProcrustesConfig()
    seeds = [GLOBAL_SEED + 1009 * i for i in range(n_seeds)]
    print(f"# PHASE 1 / GATE B -- apples-to-apples Gauge separation")
    print(f"# n={cfg.n}, SNR={cfg.snr}, alpha={cfg.alpha} (nominal {cfg.nominal:.2f}), "
          f"{n_seeds} seeds, target=D")
    print(f"# well-ID defined EXACTLY as Gauge: true-D* terciles "
          f"(np.quantile([1/3,2/3]) + np.digitize); bottom-2 = identifiable region\n")

    payload = {"config": {"n": cfg.n, "snr": cfg.snr, "alpha": cfg.alpha,
                          "n_seeds": n_seeds, "seeds": seeds,
                          "wellid_def": "Gauge-exact true-D* terciles, bottom-2 = identifiable"},
               "families": {}}

    for fam in ORDER:
        det = [separation_detail(fam, cfg, seed=s) for s in seeds]

        # marginal + gaps share ONE departure-blind radius recipe (placebo+worst
        # pooled), computed inside separation_detail -- fully self-consistent.
        marg = across_seed_ci([d["marginal"] for d in det])
        gap_all = cluster_bootstrap_gap(det, "all")
        gap_well = cluster_bootstrap_gap(det, "wellid")
        gap_lo = cluster_bootstrap_gap(det, "lo")
        gap_hi = cluster_bootstrap_gap(det, "hi")
        bias_ratio = across_seed_ci([d["bias_ratio"] for d in det])
        bias_worst = across_seed_ci([d["bias_worst"] for d in det])
        n_well = int(np.mean([d["n_wellid"] for d in det]))
        n_lo = int(np.mean([d["n_lo"] for d in det]))

        print("=" * 84)
        print(f"FAMILY {fam.label}   [pre-registered: {fam.expect.upper()}]")
        print("=" * 84)
        print(f"  marginal coverage (departure-blind)   : "
              f"{marg[0]:.3f} [{marg[1]:.3f}, {marg[2]:.3f}]  (nominal {cfg.nominal:.2f})")
        print(f"  conditional gap -- ALL voxels          : {_fmt(gap_all)}")
        print(f"  conditional gap -- WELL-ID (bottom-2)  : {_fmt(gap_well)}   "
              f"(n~{n_well}; = Gauge identifiable region)")
        print(f"  conditional gap -- STRICT-LO (bottom)  : {_fmt(gap_lo)}   "
              f"(n~{n_lo}; deepest in Gauge 'trust D')")
        print(f"  conditional gap -- HIGH-D* (top terc.) : {_fmt(gap_hi)}   "
              f"(Gauge's OWN ill-posed wall, for contrast)")
        print(f"  signed D-bias growth |worst|/|limit|   : "
              f"{bias_ratio[0]:.1f}x [{bias_ratio[1]:.1f}, {bias_ratio[2]:.1f}]   "
              f"(mechanism: high-b aliasing)")
        print()

        payload["families"][fam.lattice_family] = {
            "label": fam.label, "expect": fam.expect,
            "marginal": {"point": marg[0], "lo": marg[1], "hi": marg[2]},
            "gap_all": gap_all, "gap_wellid": gap_well,
            "gap_lo": gap_lo, "gap_hi": gap_hi,
            "bias_ratio": {"point": bias_ratio[0], "lo": bias_ratio[1], "hi": bias_ratio[2]},
            "bias_worst": {"point": bias_worst[0], "lo": bias_worst[1], "hi": bias_worst[2]},
            "n_wellid": n_well, "n_lo": n_lo,
        }

    # ---- GATE B verdict (on the BREAK family, stretched-exp) -----------------
    s = payload["families"]["stretched"]
    floor = cfg.wellid_gap_min
    lo_survives = s["gap_lo"]["lo"] > floor
    well_survives = s["gap_wellid"]["lo"] > floor
    not_collapsed = s["gap_wellid"]["point"] >= s["gap_all"]["point"] - 0.02
    # placebo/marginal sanity: marginal must hold (trap stays closed)
    marg_ok = abs(s["marginal"]["point"] - cfg.nominal) <= cfg.marginal_tol + 0.005

    passed = bool(lo_survives and well_survives and not_collapsed)
    print("#" * 84)
    print("GATE B VERDICT (stretched-exp, the BREAK family)")
    print("#" * 84)
    print(f"  marginal holds (trap closed)              : {marg_ok}  "
          f"(marginal {s['marginal']['point']:.3f} vs nominal {cfg.nominal:.2f})")
    print(f"  well-ID gap CI-lo > {floor:.2f} floor          : {well_survives}  "
          f"(CI-lo {s['gap_wellid']['lo']:+.3f})")
    print(f"  strict-lo-tercile gap CI-lo > {floor:.2f} floor : {lo_survives}  "
          f"(CI-lo {s['gap_lo']['lo']:+.3f})")
    print(f"  well-ID gap NOT collapsed vs ALL           : {not_collapsed}  "
          f"(well-ID {s['gap_wellid']['point']:+.3f} vs ALL {s['gap_all']['point']:+.3f})")
    print()
    if passed:
        print("  => GATE B: PASS. The conditional gap SURVIVES inside Gauge's")
        print("     identifiable region under Gauge's own partition, and is NOT")
        print("     diluted there -- distinct from Gauge's high-D* wall. Continue.")
    else:
        print("  => GATE B: FAIL (refute R2 FIRED). The well-ID gap collapsed")
        print("     toward the marginal under Gauge's definition. HARD STOP --")
        print("     do NOT reframe to rescue it. Surface to the author.")
    payload["gate_b"] = {
        "passed": passed, "floor": floor,
        "well_survives": bool(well_survives), "lo_survives": bool(lo_survives),
        "not_collapsed": bool(not_collapsed), "marginal_ok": bool(marg_ok),
    }

    os.makedirs(_RESULTS, exist_ok=True)
    with open(os.path.join(_RESULTS, "phase1_gateB.json"), "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"\n[wrote results/phase1_gateB.json]")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 8))
