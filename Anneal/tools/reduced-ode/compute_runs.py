"""Shared step — integrate the 3-D reduced system for every A=0.5 campaign run,
mapping each logged seed to reduced ICs (rho1=Rsync0, rho2=Rincoh0, psi=dphi0,
all precomputed in transient-tests CP2). One integration per run yields the
reduced capture time, breath period, spiral-out slope, and rho1=1 transverse
stability. Writes reduced_results/reduced_runs.jsonl, consumed by CP3 and CP4.

Deterministic; the only inputs are the read-only feature rows and the config.
Usage: python3 tools/reduced-ode/compute_runs.py
"""
import json
import os
import time

import reduced_core as rc

CFG = rc.load_config()
rc.set_config(CFG)
OUT = os.path.join(rc.ROOT, CFG["output_dir"])
os.makedirs(OUT, exist_ok=True)
BETA = CFG["beta_ours"]
A = 0.5
T_MAX = CFG["cp4"]["t_max"]
THETA = CFG["boundary"]["theta"]


def main():
    feats = os.path.join(rc.ROOT, CFG["cp4"]["input_features"])
    rows = [json.loads(l) for l in open(feats) if l.strip()]
    rows = [r for r in rows if r["A"] == A]
    p = rc.Params(A=A, beta=BETA)
    t0 = time.time()
    out_path = os.path.join(OUT, "reduced_runs.jsonl")
    n = 0
    with open(out_path, "w") as f:
        for r in rows:
            state0 = [r["Rsync0"], r["Rincoh0"], r["dphi0"]]
            res = rc.reduced_run_3d(state0, p, T_MAX, CFG, theta=THETA)
            rec = {
                "N": r["N"],
                "seed": r["seed"],
                "A": A,
                "beta": BETA,
                # measured
                "t_abs_meas": r["t_abs"],
                "abs_censored_meas": r["abs_censored"],
                # ICs
                "Rincoh0": r["Rincoh0"],
                "Rsync0": r["Rsync0"],
                "dphi0": r["dphi0"],
                "absdphi0": r["absdphi0"],
                # reduced
                "t_capture": res["t_capture"],
                "captured": res["captured"],
                "breath_period": res["breath_period"],
                "spiral_slope": res["spiral_slope"],
                "n_peaks": res["n_peaks"],
                "n_sub_theta_peaks": res["n_sub_theta_peaks"],
                "rho1_dev_max": res["rho1_dev_max"],
                "transverse_rate_mean": res["transverse_rate_mean"],
                "transverse_attracting_frac": res["transverse_attracting_frac"],
            }
            f.write(json.dumps(rec) + "\n")
            n += 1
            if n % 200 == 0:
                print(f"  {n}/{len(rows)} ({time.time() - t0:.0f}s)")
    print(f"reduced_runs.jsonl: {n} runs in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
