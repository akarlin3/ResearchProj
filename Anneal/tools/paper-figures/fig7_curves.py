"""
Figure 7 — compute the (β, A) bifurcation curves from our OWN reduced-ODE
machinery (tools/reduced-ode/reduced_core.py), then cache them so the plot is a
fast redraw.

  - Saddle-node A_SN(β): Eq. (17) series (smooth curve) + numeric fold location
    (locate_saddle_node) at ~10 β as check marks; reports max |series−numeric|.
  - Hopf A_H(β):        Eq. (18) series (smooth curve) + numeric trace(J)=0
    location (locate_hopf, windowed) where the spiral is well-defined.
  - Homoclinic A_hc(β): the existing escape-to-sync bisection over an adaptive
    bracket [A_H+δ, cap] at ~12 β in [0.02, 0.20]. SELF-CHECK: β=0.05 must
    reproduce 0.4096 within tolerance (verdict printed + cached).
  - Takens–Bogdanov point (0.2239, 0.3372) carried through from config.

Reads figures.config.json + reduced.config.json (homoclinic t_max/sync_r/tol).
Writes paper_figures/fig7_curves.json. Prior results read-only. Deterministic.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

import style

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = style.ROOT
sys.path.insert(0, os.path.join(ROOT, "tools", "reduced-ode"))
import reduced_core as rc  # noqa: E402

CFG = style.load_config()
RCFG = rc.load_config()
rc.set_config(RCFG)
OUT = os.path.join(ROOT, CFG["output_dir"])
F7 = CFG["fig7"]


def homoclinic_at(beta, a_lo, a_hi, tol):
    """Escape-to-sync bisection (reduced_core._escaped_to_sync) on A: below A_hc
    a trajectory off the chimera FP settles on a bounded breathing cycle; above
    it escapes to sync (r→1). Returns dict with A_hc, bracket, validity, and the
    bracket endpoints' escape flags (so fragility near TB is visible)."""
    h = RCFG["homoclinic"]
    t_max, sync_r = h["t_max"], h["sync_r"]
    lo_esc, lo_info = rc._escaped_to_sync(a_lo, beta, t_max, sync_r)
    hi_esc, hi_info = rc._escaped_to_sync(a_hi, beta, t_max, sync_r)
    if lo_esc is None or hi_esc is None:
        return {"beta": beta, "ok": False, "reason": "FP missing at bracket"}
    if lo_esc or not hi_esc:
        return {"beta": beta, "ok": False, "reason": "bracket does not straddle",
                "lo_escaped": bool(lo_esc), "hi_escaped": bool(hi_esc)}
    lo, hi = a_lo, a_hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        esc, _ = rc._escaped_to_sync(mid, beta, t_max, sync_r)
        if esc:
            hi = mid
        else:
            lo = mid
    return {"beta": beta, "ok": True, "A_hc": 0.5 * (lo + hi),
            "bracket": [lo, hi], "width": hi - lo}


def main():
    t_start = time.time()
    b0, b1 = F7["beta_curve"]
    nb = F7["beta_curve_n"]
    beta_curve = np.linspace(b0, b1, nb)

    # --- series curves (Eqs 17/18) over the dense grid ---
    A_SN_series = np.array([rc.A_SN_series(b, RCFG) for b in beta_curve])
    A_H_series = np.array([rc.A_H_series(b, RCFG) for b in beta_curve])

    # --- numeric SN / Hopf at the check-mark β values ---
    sn_num, hopf_num = [], []
    for b in F7["beta_numeric"]:
        snn = rc.locate_saddle_node(b)
        if snn:
            sn_num.append({"beta": b, "A": float(snn[0]),
                           "series": float(rc.A_SN_series(b, RCFG))})
        ah_s = rc.A_H_series(b, RCFG)
        # windowed numeric Hopf: avoid the spurious far-branch crossing at low β
        a_lo = max(0.05, rc.A_SN_series(b, RCFG) + 0.02)
        hn = rc.locate_hopf(b, a_lo=a_lo, a_hi=ah_s + 0.06)
        if hn:
            hopf_num.append({"beta": b, "A": float(hn[0]), "series": float(ah_s)})

    sn_dev = max((abs(p["A"] - p["series"]) for p in sn_num), default=None)
    hopf_dev = max((abs(p["A"] - p["series"]) for p in hopf_num), default=None)

    # --- homoclinic bisection over β (adaptive bracket) ---
    hcfg = F7["homoclinic"]
    tol = hcfg["bisect_tol"]
    hc_points, hc_degraded = [], None
    print("tracing homoclinic A_hc(β) by escape-to-sync bisection ...")
    for b in F7["beta_homoclinic"]:
        a_lo = rc.A_H_series(b, RCFG) + hcfg["a_lo_above_hopf"]
        res = homoclinic_at(b, a_lo, hcfg["a_hi_cap"], tol)
        if res.get("ok"):
            hc_points.append({"beta": b, "A_hc": res["A_hc"],
                              "bracket": res["bracket"], "width": res["width"]})
            print(f"  β={b:.3f}  A_hc={res['A_hc']:.4f}  "
                  f"(bracket±{res['width']:.1e})  [{time.time()-t_start:.0f}s]")
        else:
            hc_degraded = {"beta": b, **res}
            print(f"  β={b:.3f}  homoclinic bisection degraded: {res.get('reason')}"
                  f" — stopping curve here")
            break

    # --- β=0.05 self-check ---
    sc = F7["homoclinic_selfcheck"]
    sc_hit = next((p for p in hc_points if abs(p["beta"] - sc["beta"]) < 1e-9), None)
    if sc_hit is None:
        res = homoclinic_at(sc["beta"], rc.A_H_series(sc["beta"], RCFG)
                            + hcfg["a_lo_above_hopf"], hcfg["a_hi_cap"], tol)
        sc_hit = {"beta": sc["beta"], "A_hc": res.get("A_hc"),
                  "bracket": res.get("bracket"), "width": res.get("width")}
    sc_err = abs(sc_hit["A_hc"] - sc["expected_A_hc"])
    sc_pass = sc_err <= sc["abs_tol"]
    verdict = (f"β={sc['beta']}: A_hc={sc_hit['A_hc']:.4f} vs expected "
               f"{sc['expected_A_hc']} (|Δ|={sc_err:.4g}, tol={sc['abs_tol']}) → "
               f"{'PASS' if sc_pass else 'FAIL'}")
    print("SELF-CHECK:", verdict)

    out = {
        "beta_curve": beta_curve.tolist(),
        "A_SN_series": A_SN_series.tolist(),
        "A_H_series": A_H_series.tolist(),
        "sn_numeric": sn_num,
        "hopf_numeric": hopf_num,
        "homoclinic": hc_points,
        "homoclinic_degraded": hc_degraded,
        "deviations": {"A_SN_max": sn_dev, "A_H_max": hopf_dev},
        "selfcheck": {"beta": sc["beta"], "A_hc": sc_hit["A_hc"],
                      "expected": sc["expected_A_hc"], "abs_err": sc_err,
                      "bracket": sc_hit.get("bracket"), "tol": sc["abs_tol"],
                      "pass": bool(sc_pass), "verdict": verdict},
        "tb_point": F7["tb_point"],
        "corners": F7["corners"],
        "compute_seconds": round(time.time() - t_start, 1),
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "fig7_curves.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", os.path.relpath(os.path.join(OUT, "fig7_curves.json"), ROOT),
          f"({out['compute_seconds']}s)")
    print(f"deviations: A_SN max={sn_dev:.2e}  A_H max={hopf_dev:.2e}")


if __name__ == "__main__":
    main()
