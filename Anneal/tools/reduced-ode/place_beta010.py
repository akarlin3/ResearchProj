"""beta=0.10 reduced-model placement for the robustness slice.

Computes, at beta=0.10 (the published-regime value the reduced model is
validated against):
  - A_SN(0.10), A_H(0.10) [Eq.17/Eq.18 series + numerical],
  - A_hc(0.10) [homoclinic, bisection], with bracket width,
  - region + fixed-point classification (sigma, omega) for the selected
    A_post / A_stable points.

Reuses tools/reduced-ode/reduced_core.py exclusively (no reimplementation).
Writes reduced_results/beta010_placement.json and .md.

Run: python3 tools/reduced-ode/place_beta010.py
"""
import json
import os

import numpy as np

import reduced_core as rc

CFG = rc.load_config()
rc.set_config(CFG)
OUT = os.path.join(rc.ROOT, CFG["output_dir"])
os.makedirs(OUT, exist_ok=True)

BETA = 0.10
# Candidate corner points on the existing campaign grid.
A_POST = 0.5     # post-homoclinic (parallels the operating corner)
A_STABLE = 0.2   # stable stationary-chimera band (never-absorber control)


def region_of(A, A_SN, A_H, A_hc):
    if A < A_SN:
        return "no-chimera (below saddle-node)"
    if A < A_H:
        return "stable stationary chimera"
    if A_hc is not None and A < A_hc:
        return "breathing chimera (limit cycle)"
    return "post-homoclinic (no chimera attractor; sync attracting)"


def classify_point(A, beta, A_SN, A_H, A_hc):
    """Region + (sigma, omega) of the chimera fixed point (when one exists)."""
    reg = region_of(A, A_SN, A_H, A_hc)
    info = {"A": A, "region": reg}
    try:
        fp = rc.get_fixed_point(A, beta, branch="chimera")
        if fp is None:
            raise ValueError("fixed-point family does not reach this A")
        r_star, psi_star, st = fp
        ev = st["eig"]
        # leading eigenvalue (largest real part)
        lead = ev[int(np.argmax(np.real(ev)))]
        info.update({
            "r_star": float(r_star),
            "psi_star": float(psi_star),
            "sigma": float(np.real(lead)),       # growth rate (real part)
            "omega": float(abs(np.imag(lead))),  # spiral frequency (|imag|)
            "stable_fp": bool(st["stable"]),
            "eig": [complex(e).__repr__() for e in ev],
        })
    except Exception as e:  # post-homoclinic: chimera FP may be irrelevant
        info["fp_note"] = f"no usable chimera fixed point: {e}"
    # dynamical confirmation
    if A < A_H:
        try:
            fp = rc.get_fixed_point(A, beta, branch="chimera")
            cls = rc.classify_2d(
                rc.Params(A=A, beta=beta),
                min(0.999, fp[0] + 0.02), fp[1] + 0.05,
                CFG["classify"]["t_settle"], CFG["classify"]["t_window"],
                const_tol=CFG["classify"]["const_tol"],
            )
            info["dynamical_confirmation"] = cls
        except Exception as e:
            info["dynamical_confirmation"] = {"error": str(e)}
    else:
        esc, einfo = rc._escaped_to_sync(
            A, beta, CFG["homoclinic"]["t_max"], CFG["homoclinic"]["sync_r"])
        info["dynamical_confirmation"] = {"kind": "sync" if esc else "bounded", **(einfo or {})}
    return info


def main():
    A_SN_s = rc.A_SN_series(BETA, CFG)
    A_H_s = rc.A_H_series(BETA, CFG)
    sn_num = rc.locate_saddle_node(BETA)
    h_num = rc.locate_hopf(BETA)

    # Homoclinic bisection. Widen the bracket if the default (tuned for 0.05)
    # does not straddle A_hc(0.10).
    cfg = json.loads(json.dumps(CFG))  # deep copy
    print("locating homoclinic at beta=0.10 (bisection)...")
    hc = rc.locate_homoclinic(BETA, cfg)
    if "A_hc" not in hc:
        for lo, hi in [(0.28, 0.50), (0.28, 0.55), (0.30, 0.60)]:
            cfg["homoclinic"]["A_lo"] = lo
            cfg["homoclinic"]["A_hi"] = hi
            hc = rc.locate_homoclinic(BETA, cfg)
            if "A_hc" in hc:
                break
    A_hc = hc.get("A_hc")

    post = classify_point(A_POST, BETA, A_SN_s, A_H_s, A_hc)
    stable = classify_point(A_STABLE, BETA, A_SN_s, A_H_s, A_hc)

    # depth-matched contingency candidate: operating corner (0.05) depth was
    # A_op - A_hc(0.05). At 0.10 the depth-matched A is A_hc(0.10) + 0.09
    # (the reviewer's matched offset; 0.41-0.32 ~= 0.09 at the operating corner).
    A_depth = round(A_hc + 0.09, 3) if A_hc else None

    out = {
        "beta": BETA,
        "A_SN": {"series_eq17": A_SN_s, "numeric": sn_num[0] if sn_num else None},
        "A_H": {"series_eq18": A_H_s, "numeric": h_num[0] if h_num else None},
        "A_hc": {"numeric": A_hc, "bracket": hc.get("A_hc_bracket"), "width": hc.get("width")},
        "selections": {
            "A_post": A_POST, "A_stable": A_STABLE,
            "A_depth_matched_candidate": A_depth,
            "depth_offset_used": 0.09,
        },
        "post": post,
        "stable": stable,
    }
    with open(os.path.join(OUT, "beta010_placement.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)

    lines = [
        "# beta=0.10 reduced-model placement\n",
        f"A_SN(0.10): series(Eq17)={A_SN_s:.5f}  numeric={(sn_num[0] if sn_num else float('nan')):.5f}",
        f"A_H(0.10):  series(Eq18)={A_H_s:.5f}  numeric={(h_num[0] if h_num else float('nan')):.5f}",
        f"A_hc(0.10): numeric={A_hc:.5f}  bracket={hc.get('A_hc_bracket')}  width={hc.get('width'):.2e}"
        if A_hc else f"A_hc(0.10): FAILED — {hc.get('error')}",
        "",
        f"A_post   = {A_POST}: {post['region']}",
        (f"   sigma={post.get('sigma')}, omega={post.get('omega')}"
         if "sigma" in post else f"   {post.get('fp_note','')}"),
        f"   dynamical: {post['dynamical_confirmation'].get('kind', post['dynamical_confirmation'])}",
        f"A_stable = {A_STABLE}: {stable['region']}",
        (f"   sigma={stable.get('sigma')}, omega={stable.get('omega')}, r*={stable.get('r_star')}"
         if "sigma" in stable else f"   {stable.get('fp_note','')}"),
        f"   dynamical: {stable['dynamical_confirmation'].get('kind', stable['dynamical_confirmation'])}",
        "",
        f"depth-matched contingency A = A_hc(0.10)+0.09 = {A_depth}",
    ]
    with open(os.path.join(OUT, "beta010_placement.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
