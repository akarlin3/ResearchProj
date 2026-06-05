"""CP2 — classify OUR corners (beta=0.05).

- A_SN(0.05), A_H(0.05) [series + numerical]; A_hc(0.05) [homoclinic, numerical].
- Verdict table for A in {0.2, 0.5}: which region, against predictions P1/P2.
- Phase portraits of the 2-D system at both corners (paper Fig. 3 style).
Writes reduced_results/cp2_regions.json, cp2_regions.md, and the portrait PNG/PDF.
"""
import json
import os

import numpy as np

import reduced_core as rc
import portraits

CFG = rc.load_config()
rc.set_config(CFG)
OUT = os.path.join(rc.ROOT, CFG["output_dir"])
os.makedirs(OUT, exist_ok=True)
BETA = CFG["beta_ours"]  # 0.05


def load_seed_ics(n=6):
    """Representative canonical-seed ICs (r0=Rincoh0, psi0=dphi0) — the IC map is
    A-independent, so the A=0.5 feature rows give the canonical-seed neighbourhood
    for both corners. Pick a spread of |dphi0|."""
    path = os.path.join(rc.ROOT, CFG["cp4"]["input_features"])
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    # spread over dphi0 at a mid N for a legible portrait
    rows = [r for r in rows if r["N"] == 16]
    rows.sort(key=lambda r: r["absdphi0"])
    idx = np.linspace(0, len(rows) - 1, n).astype(int)
    return [(rows[i]["Rincoh0"], rows[i]["dphi0"]) for i in idx]


def region_of(A, A_SN, A_H, A_hc):
    if A < A_SN:
        return "no-chimera (below saddle-node)"
    if A < A_H:
        return "stable stationary chimera"
    if A_hc is not None and A < A_hc:
        return "breathing chimera (limit cycle)"
    return "post-homoclinic (no chimera attractor; sync attracting)"


def main():
    A_SN_s = rc.A_SN_series(BETA, CFG)
    A_H_s = rc.A_H_series(BETA, CFG)
    sn_num = rc.locate_saddle_node(BETA)
    h_num = rc.locate_hopf(BETA)
    print("locating homoclinic (bisection)...")
    hc = rc.locate_homoclinic(BETA, CFG)
    A_hc = hc.get("A_hc")

    bif = {
        "beta": BETA,
        "A_SN": {"series_eq17": A_SN_s, "numeric": sn_num[0] if sn_num else None},
        "A_H": {"series_eq18": A_H_s, "numeric": h_num[0] if h_num else None},
        "A_hc": {
            "numeric": A_hc,
            "bracket": hc.get("A_hc_bracket"),
            "width": hc.get("width"),
        },
    }

    # verdicts for A in {0.2, 0.5}
    verdict = {}
    for A in CFG["A_ours"]:
        reg = region_of(A, A_SN_s, A_H_s, A_hc)
        # dynamical confirmation
        if A < A_H_s:
            fp = rc.get_fixed_point(A, BETA, branch="chimera")
            cls = rc.classify_2d(
                rc.Params(A=A, beta=BETA),
                min(0.999, fp[0] + 0.02),
                fp[1] + 0.05,
                CFG["classify"]["t_settle"],
                CFG["classify"]["t_window"],
                const_tol=CFG["classify"]["const_tol"],
            )
            confirm = cls
        else:
            esc, info = rc._escaped_to_sync(
                A, BETA, CFG["homoclinic"]["t_max"], CFG["homoclinic"]["sync_r"]
            )
            confirm = {"kind": "sync" if esc else "bounded", **info}
        pred = "P1" if abs(A - 0.2) < 1e-9 else ("P2" if abs(A - 0.5) < 1e-9 else None)
        pred_text = CFG["predictions"].get(pred, "")
        pred_ok = (
            (pred == "P1" and reg == "stable stationary chimera")
            or (pred == "P2" and reg.startswith("post-homoclinic"))
        )
        verdict[str(A)] = {
            "A": A,
            "region": reg,
            "dynamical_confirmation": confirm,
            "prediction": pred,
            "prediction_text": pred_text,
            "prediction_holds": bool(pred_ok),
        }

    # portraits
    seed_ics = load_seed_ics()
    portraits.make_portrait(
        0.2, BETA, seed_ics, os.path.join(OUT, "cp2_portrait_A0.2"), CFG,
        title="2-D phase portrait  A=0.2, β=0.05  (stable stationary chimera)",
        traj_t=500.0,
    )
    portraits.make_portrait(
        0.5, BETA, seed_ics, os.path.join(OUT, "cp2_portrait_A0.5"), CFG,
        title="2-D phase portrait  A=0.5, β=0.05  (post-homoclinic → sync)",
        traj_t=300.0,
    )

    out = {"checkpoint": "CP2", "beta": BETA, "bifurcations": bif,
           "verdicts": verdict, "homoclinic_trace": hc.get("trace")}
    with open(os.path.join(OUT, "cp2_regions.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)

    # markdown table
    lines = [
        "# CP2 — our corners (β=0.05)\n",
        f"A_SN(0.05): series(Eq17)={A_SN_s:.5f}  numeric={sn_num[0]:.5f}",
        f"A_H(0.05): series(Eq18)={A_H_s:.5f}  numeric={h_num[0]:.5f}",
        f"A_hc(0.05): numeric={A_hc:.5f}  bracket={hc.get('A_hc_bracket')}  width={hc.get('width'):.2e}",
        "",
        "| A | region | prediction | holds | dynamical confirmation |",
        "|---|---|---|---|---|",
    ]
    for A, v in verdict.items():
        c = v["dynamical_confirmation"]
        ck = c.get("kind")
        extra = ""
        if ck == "fixed":
            extra = f"r*={c.get('r_star'):.4f}, amp={c.get('amp'):.1e}"
        elif ck == "cycle":
            extra = f"cycle, period={c.get('period')}"
        elif ck == "sync":
            extra = f"escaped→sync, t_esc={c.get('t_escape')}"
        lines.append(
            f"| {A} | {v['region']} | {v['prediction']} | "
            f"{'YES' if v['prediction_holds'] else 'NO'} | {ck} ({extra}) |"
        )
    with open(os.path.join(OUT, "cp2_regions.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print("\nportraits: cp2_portrait_A0.2.{png,pdf}, cp2_portrait_A0.5.{png,pdf}")


if __name__ == "__main__":
    main()
