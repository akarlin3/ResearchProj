"""CP1 — implementation + validation gate against the PAPER's own results.

Three sub-gates, ALL must pass (else STOP and report the discrepancy):
  (a) Eqs (13)-(14) stationary chimeras zero rdot=psidot in the 2-D system.
  (b) Paper Fig. 2 regimes at beta=0.1: A=0.2 stable FP; A=0.28 breathing cycle;
      A=0.35 long-period breathing. Classify and report periods.
  (c) Numerically locate the saddle-node and Hopf at beta=0.1 and match the
      series Eqs (17)-(18) within series-truncation error.

Writes reduced_results/cp1_gate.json. Exits non-zero if any sub-gate FAILS.
"""
import json
import math
import os
import sys

import numpy as np

import reduced_core as rc

CFG = rc.load_config()
rc.set_config(CFG)
OUT = os.path.join(rc.ROOT, CFG["output_dir"])
os.makedirs(OUT, exist_ok=True)
BETA = CFG["beta_paper"]  # 0.1


def gate_a():
    fp = CFG["fixed_point"]
    psi_grid = np.linspace(fp["psi_min"], fp["psi_max"], fp["psi_grid_n"])
    psi, r, A = rc.fixed_point_family(BETA, psi_grid)
    worst = 0.0
    n = 0
    for pj, rj, Aj in zip(psi, r, A):
        if not (0.0 < Aj < 0.8):
            continue
        d = rc.rhs_2d([rj, pj], rc.Params(A=Aj, beta=BETA))
        worst = max(worst, float(np.linalg.norm(d)))
        n += 1
    tol = 1e-9
    return {
        "beta": BETA,
        "n_points": n,
        "A_range": [float(A.min()), float(A.max())] if len(A) else None,
        "max_rhs_residual": worst,
        "tol": tol,
        "pass": worst < tol,
    }


def gate_b():
    cls = CFG["classify"]
    results = {}
    expect = {0.2: "fixed", 0.28: "cycle", 0.35: "cycle"}
    for A in cls["A_list"]:
        fp = rc.get_fixed_point(A, BETA, branch="chimera")
        if fp is None:
            results[str(A)] = {"error": "no fixed point on upper branch"}
            continue
        r_star, psi_star, st = fp
        # start just off the fixed point
        r0 = min(0.999, r_star + cls["r0_offset"])
        psi0 = psi_star + cls["psi0_offset"]
        res = rc.classify_2d(
            rc.Params(A=A, beta=BETA),
            r0,
            psi0,
            cls["t_settle"],
            cls["t_window"],
            const_tol=cls["const_tol"],
        )
        res["A"] = A
        res["r_star_eq13_14"] = r_star
        res["fp_stable"] = st["stable"]
        res["fp_eig"] = [complex(e).__repr__() for e in st["eig"]]
        res["expected"] = expect[A]
        res["match"] = res["kind"] == expect[A]
        results[str(A)] = res
    all_match = all(v.get("match", False) for v in results.values())
    # period ordering: A=0.35 should be longer-period than A=0.28
    p28 = results.get("0.28", {}).get("period")
    p35 = results.get("0.35", {}).get("period")
    period_order_ok = (
        p28 is not None and p35 is not None and p35 > p28
    )
    return {
        "beta": BETA,
        "regimes": results,
        "period_order_028_lt_035": period_order_ok,
        "pass": bool(all_match and period_order_ok),
    }


def gate_c():
    sn = rc.locate_saddle_node(BETA)
    hopf = rc.locate_hopf(BETA)
    sn_series = rc.A_SN_series(BETA, CFG)
    h_series = rc.A_H_series(BETA, CFG)
    # series truncation error ~ next-order term magnitude at beta=0.1.
    # A_SN truncated at beta^5; next term O(beta^6) ~ |coef|*1e-6, coef O(1-10) =>
    # tolerance ~ 5e-3 is generous and honest. A_H truncated at beta^4; next is
    # O(beta^6) (only even powers) ~ 1e-6 => the Hopf series is far tighter, but the
    # NUMERICAL location (finite grid + numeric Jacobian) sets the floor ~ a few e-3.
    tol_sn = 5e-3
    tol_h = 5e-3
    sn_num = sn[0] if sn else None
    h_num = hopf[0] if hopf else None
    sn_err = abs(sn_num - sn_series) if sn_num is not None else None
    h_err = abs(h_num - h_series) if h_num is not None else None
    return {
        "beta": BETA,
        "saddle_node": {
            "numeric_A": sn_num,
            "numeric_r": sn[1] if sn else None,
            "numeric_psi": sn[2] if sn else None,
            "series_A_eq17": sn_series,
            "abs_err": sn_err,
            "tol": tol_sn,
            "pass": sn_err is not None and sn_err < tol_sn,
        },
        "hopf": {
            "numeric_A": h_num,
            "numeric_psi": hopf[1] if hopf else None,
            "series_A_eq18": h_series,
            "abs_err": h_err,
            "tol": tol_h,
            "pass": h_err is not None and h_err < tol_h,
        },
        "pass": bool(
            sn_err is not None
            and sn_err < tol_sn
            and h_err is not None
            and h_err < tol_h
        ),
    }


def main():
    a = gate_a()
    b = gate_b()
    c = gate_c()
    overall = bool(a["pass"] and b["pass"] and c["pass"])
    out = {
        "checkpoint": "CP1",
        "purpose": "validation gate against Abrams et al. 2008 (the paper itself)",
        "beta_paper": BETA,
        "gate_a_fixed_point_consistency": a,
        "gate_b_fig2_regimes": b,
        "gate_c_bifurcation_series": c,
        "overall_pass": overall,
    }
    with open(os.path.join(OUT, "cp1_gate.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)

    print("=== CP1 validation gate (beta=0.1, vs the paper) ===")
    print(f"(a) fixed-point consistency: max|RHS|={a['max_rhs_residual']:.2e}  "
          f"-> {'PASS' if a['pass'] else 'FAIL'}")
    print("(b) Fig.2 regimes:")
    for A, v in b["regimes"].items():
        per = v.get("period")
        per_s = f"{per:.2f}" if per else "—"
        print(f"    A={A}: kind={v.get('kind')} (expected {v.get('expected')}) "
              f"period={per_s}  -> {'PASS' if v.get('match') else 'FAIL'}")
    print(f"    period order T(0.35)>T(0.28): {b['period_order_028_lt_035']}  "
          f"-> {'PASS' if b['pass'] else 'FAIL'}")
    print("(c) bifurcation series:")
    sc = c["saddle_node"]
    hc = c["hopf"]
    fmt = lambda x: f"{x:.5f}" if x is not None else "None"
    fe = lambda x: f"{x:.2e}" if x is not None else "None"
    print(f"    SN : numeric={fmt(sc['numeric_A'])}  series(Eq17)={fmt(sc['series_A_eq17'])}  "
          f"|err|={fe(sc['abs_err'])}  -> {'PASS' if sc['pass'] else 'FAIL'}")
    print(f"    Hopf: numeric={fmt(hc['numeric_A'])}  series(Eq18)={fmt(hc['series_A_eq18'])}  "
          f"|err|={fe(hc['abs_err'])}  -> {'PASS' if hc['pass'] else 'FAIL'}")
    print(f"OVERALL CP1: {'PASS' if overall else 'FAIL — STOP, await review'}")
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()
