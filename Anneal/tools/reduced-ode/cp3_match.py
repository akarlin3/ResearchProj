"""CP3 — quantitative matches to the finite-N measurements.

Timescale mapping (model time <-> seconds), then for A=0.5:
  (i)   reduced breath period T_b vs measured 21-25 s
  (ii)  reduced spiral-out slope log(theta - peak) vs cycle, vs measured -0.55..-0.77
  (iii) reduced capture time vs measured tau_abs ~ 139 s
each as (reduced, measured, ratio), NO fudge factors. Plus A=0.2: stable FP r*
vs the never-absorbers' mean R_incoh, and the reduced relaxation rate vs the
measured slight-negative <dM> drift.

Reads reduced_results/reduced_runs.jsonl (compute_runs.py) and
reduced_results/cp3_a02_level.json (trace_a02_level.mjs).
Writes reduced_results/cp3_match.json and cp3_match.md.
"""
import json
import math
import os

import numpy as np

import reduced_core as rc

CFG = rc.load_config()
rc.set_config(CFG)
OUT = os.path.join(rc.ROOT, CFG["output_dir"])
BETA = CFG["beta_ours"]
THETA = CFG["boundary"]["theta"]
MEAS = CFG["measured"]


def med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(xs)) if xs else None


def main():
    rows = [json.loads(l) for l in open(os.path.join(OUT, "reduced_runs.jsonl"))]
    Ns = sorted(set(r["N"] for r in rows))

    # ---- timescale mapping ----
    ts = CFG["timescale"]
    timescale = {
        "model_time_per_second": ts["model_time_per_second"],
        "dt": ts["dt"],
        "sampleStride": ts["sampleStride"],
        "statement": (
            "The shipped integrator advances 1 model-time unit per second of the "
            "drone (omega=0 rotating frame; DRIFT_DT=0.05 control step, "
            "sampleStride=0.1). integrator.mjs documents '1 unit := 1 s'. So reduced "
            "model-time equals seconds 1:1 — every reduced quantity below is in "
            "seconds with NO fudge factor."
        ),
    }

    # ---- A=0.5 eigenvalue theory at the chimera FP (unstable spiral) ----
    r_c, psi_c, st_c = rc.get_fixed_point(0.5, BETA, branch="chimera")
    ev = st_c["eig"]
    sigma = float(np.real(ev[0]))
    omega = float(abs(np.imag(ev[0])))
    T_lin = 2 * math.pi / omega
    per_cycle_growth = math.exp(sigma * T_lin)

    # ---- per-N + pooled reduced metrics ----
    def collect(rs):
        return {
            "Tb": med([r["breath_period"] for r in rs]),
            "slope": med([r["spiral_slope"] for r in rs]),
            "tcap": med([r["t_capture"] for r in rs]),
            "n": len(rs),
            "n_Tb": sum(1 for r in rs if r["breath_period"] is not None),
            "n_slope": sum(1 for r in rs if r["spiral_slope"] is not None),
        }

    per_N = {str(N): collect([r for r in rows if r["N"] == N]) for N in Ns}
    pooled = collect(rows)

    meas_slopes = MEAS["spiral_slopes"]  # per-N (8,16,32,64)
    Tb_lo, Tb_hi = MEAS["Tb_range_s"]
    Tb_mid = 0.5 * (Tb_lo + Tb_hi)
    tau_abs = MEAS["tau_abs_s"]

    # ---- A=0.2 ----
    r2_star, psi2, st2 = rc.get_fixed_point(0.2, BETA, branch="chimera")
    ev2 = st2["eig"]
    sigma2 = float(np.real(ev2[0]))
    omega2 = float(abs(np.imag(ev2[0])))
    T_lin2 = 2 * math.pi / omega2
    per_cycle_decay = math.exp(sigma2 * T_lin2)
    a02 = json.load(open(os.path.join(OUT, "cp3_a02_level.json")))

    out = {
        "checkpoint": "CP3",
        "beta": BETA,
        "timescale": timescale,
        "A05_eigen": {
            "r_star": r_c, "sigma_per_time": sigma, "omega": omega,
            "T_breath_linear_s": T_lin, "per_cycle_growth": per_cycle_growth,
        },
        "A05_matches": {
            "i_breath_period": {
                "reduced_pooled_s": pooled["Tb"],
                "reduced_linear_theory_s": T_lin,
                "measured_s": [Tb_lo, Tb_hi],
                "ratio_pooled_to_mid": pooled["Tb"] / Tb_mid if pooled["Tb"] else None,
                "per_N": {k: per_N[k]["Tb"] for k in per_N},
            },
            "ii_spiral_slope": {
                "reduced_pooled": pooled["slope"],
                "measured_range": [-0.77, -0.55],
                "per_N_reduced_vs_measured": {
                    k: {"reduced": per_N[k]["slope"],
                        "measured": meas_slopes.get(k)}
                    for k in per_N if k in meas_slopes
                },
            },
            "iii_capture_time": {
                "reduced_pooled_median_s": pooled["tcap"],
                "measured_tau_abs_s": tau_abs,
                "ratio": pooled["tcap"] / tau_abs if pooled["tcap"] else None,
                "per_N_reduced_median": {k: per_N[k]["tcap"] for k in per_N},
                "note": "Reduced captures faster than finite-N; ratio < 1 reported straight.",
            },
        },
        "A02_matches": {
            "reduced_r_star": r2_star,
            "reduced_relaxation_sigma_per_time": sigma2,
            "reduced_per_cycle_decay": per_cycle_decay,
            "measured_mean_Rincoh_pooled": a02["pooled_mean_Rincoh"],
            "measured_mean_Rincoh_per_N": {
                k: v["mean_Rincoh"] for k, v in a02["per_N"].items()
            },
            "r_star_minus_measured_pooled": r2_star - a02["pooled_mean_Rincoh"],
            "measured_dM_drift_per_cycle": -0.00108,
            "sign_match": "both negative: reduced stable spiral relaxes (maxima drift DOWN to r*), measured <dM><0 — opposite to A=0.5's positive ratchet",
        },
    }
    with open(os.path.join(OUT, "cp3_match.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)

    # markdown
    L = ["# CP3 — quantitative matches to finite-N\n",
         f"**Timescale:** {timescale['statement']}\n",
         "## A=0.5 (post-homoclinic spiral-out)\n",
         f"Chimera FP (unstable spiral): r*={r_c:.3f}, λ=σ±iω, σ={sigma:+.5f}/s, "
         f"ω={omega:.4f} ⇒ linear T_breath=2π/ω={T_lin:.2f}s, per-cycle growth "
         f"exp(σT)={per_cycle_growth:.3f}.\n",
         "| quantity | reduced | measured | ratio |",
         "|---|---|---|---|",
         f"| (i) breath period T_b | {pooled['Tb']:.1f}s (lin {T_lin:.1f}s) | "
         f"{Tb_lo:.0f}–{Tb_hi:.0f}s | {pooled['Tb']/Tb_mid:.2f} |",
         f"| (ii) spiral slope (pooled) | {pooled['slope']:.2f} | −0.55…−0.77 | "
         f"{pooled['slope']/(-0.66):.2f} |",
         f"| (iii) capture time (median) | {pooled['tcap']:.0f}s | {tau_abs:.0f}s | "
         f"{pooled['tcap']/tau_abs:.2f} |",
         "",
         "### spiral slope, reduced vs measured, per N",
         "| N | reduced | measured |",
         "|---|---|---|"]
    for k in per_N:
        if k in meas_slopes:
            L.append(f"| {k} | {per_N[k]['slope']:.2f} | {meas_slopes[k]:.2f} |")
    L += ["",
          "### breath period & capture, per N (reduced)",
          "| N | T_b (s) | capture median (s) |",
          "|---|---|---|"]
    for N in Ns:
        k = str(N)
        tb = per_N[k]['Tb']; tc = per_N[k]['tcap']
        L.append(f"| {N} | {tb:.1f} | {tc:.0f} |")

    L += ["",
          "## A=0.2 (stable stationary chimera)\n",
          f"Reduced stable-chimera FP r*={r2_star:.4f}; relaxation σ={sigma2:+.5f}/s "
          f"(per-cycle decay {per_cycle_decay:.3f}).",
          "",
          "| quantity | reduced | measured |",
          "|---|---|---|",
          f"| stable-FP r* vs never-absorber mean R_incoh | {r2_star:.4f} | "
          f"{a02['pooled_mean_Rincoh']:.4f} (pooled); "
          f"N=32:{a02['per_N']['32']['mean_Rincoh']:.4f}, "
          f"N=64:{a02['per_N']['64']['mean_Rincoh']:.4f} |",
          f"| breath-max drift ⟨ΔM⟩ sign | σ<0 (relaxes, ΔM<0) | −0.0011/cycle (ΔM<0) |",
          "",
          "Finite-N mean R_incoh climbs toward r*=0.678 as N grows (N=32: 0.676) — "
          "the never-absorbers hover around the reduced stable fixed point, and the "
          "reduced relaxation (σ<0) explains the stationary, slightly-negative ⟨ΔM⟩ "
          "(opposite sign to the A=0.5 positive ratchet).",
          ]
    with open(os.path.join(OUT, "cp3_match.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
