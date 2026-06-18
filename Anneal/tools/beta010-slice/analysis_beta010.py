"""beta=0.10 robustness-slice analysis — replicates the paper pipeline exactly.

Reuses (via AST import, no reimplementation) the canonical stat kernels and
pipeline functions from:
  - tools/absorption-recampaign/analysis.py  (survival/exp-MLE, censored-Weibull,
    geometric per-pass, breath-phase Rayleigh, peak detection)
  - tools/transient-tests/analysis.py         (cluster-bootstrap ⟨ΔM⟩, ensemble M_k)

Operates on the beta=0.10 campaign + phase traces + M_k extraction. Prints every
headline number to stdout and writes absorption_results/beta010_results.json.

Run: python3 tools/beta010-slice/analysis_beta010.py
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import linregress

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from _reuse import load_funcs  # noqa: E402

ABS_PY = ROOT / "tools/absorption-recampaign/analysis.py"
TRANS_PY = ROOT / "tools/transient-tests/analysis.py"
CFG = json.loads((ROOT / "tools/absorption-recampaign/absorption.beta010.config.json").read_text())
TRANS_CFG = json.loads((ROOT / "tools/transient-tests/transient.config.json").read_text())

CAMPAIGN = ROOT / "absorption_results/absorption_campaign_beta010.jsonl"
PHASE_TRACES_PATH = ROOT / "absorption_results/phase_traces_beta010.jsonl"
OUT = ROOT / "absorption_results"
THETA = CFG["graze_criterion"]["theta"]

A_POST = 0.5
A_DEPTH = 0.444
A_STABLE = 0.2

# --------------------------------------------------------------------------- #
# Reuse the absorption pipeline functions, wiring beta=0.10 data into globals.
# --------------------------------------------------------------------------- #
A = load_funcs(ABS_PY)
# constants the reused functions close over
A["BR"] = CFG["breath"]
A["MIN_CYCLES"] = CFG["breath"]["minCyclesForPhase"]
A["THETA"] = THETA
A["W"] = CFG["graze_criterion"]["W"]
A["T_V"] = CFG["absorption_criterion"]["T_v"]
A["REC_THRESH"] = CFG["absorption_criterion"]["recoveryThreshold"]
A["REC_WIN"] = CFG["absorption_criterion"]["recoveryWindowSec"]
A["PHASE_TRACES"] = PHASE_TRACES_PATH
A["LAB"] = {A_POST: "A_post", A_DEPTH: "A_depth", A_STABLE: "A_stable"}

# data globals, in the same order the module computes them
A["camp"] = A["load_jsonl"](CAMPAIGN)
A["graze"] = A["by_point"](A["camp"], "t_graze", "graze_censored")
A["absb"] = A["by_point"](A["camp"], "t_abs", "abs_censored")
A["A_vals"] = sorted({a for a, _ in A["graze"]}, reverse=True)

surv = A["survival_table"]()
A["surv"] = surv
curves = A["plateau_curves"]()
weib = A["weibull_table"]()
geo = A["geometric_table"]()
phase_agg, phase_pooled, phase_pooled05, tb_cross = A["phase_clustering"]()

# --------------------------------------------------------------------------- #
# Reuse the transient ratchet functions for the post-homoclinic M_k sequences.
# --------------------------------------------------------------------------- #
T = load_funcs(TRANS_PY)
T["RNG"] = np.random.default_rng(TRANS_CFG["cp2"]["cv_seed"])
ALLOW_INV = TRANS_CFG["cp1"]["monotone_allowed_inversions"]


def read_jsonl(p):
    out = []
    for line in Path(p).read_text().splitlines():
        s = line.strip()
        if s:
            out.append(json.loads(s))
    return out


def ratchet_for(a_tag, Ns):
    runs_all = read_jsonl(ROOT / f"transient_results/cp1_mk_beta010_A{a_tag}.jsonl")
    rows = []
    for N in Ns:
        runs = [r for r in runs_all if r["N"] == N]
        mks = [np.asarray(r["Mk"], float) for r in runs]
        testable = [m for m in mks if len(m) >= 3]
        ratchet = [m for m in testable if T["count_inversions"](m) <= ALLOW_INV]
        strict = [m for m in testable if T["count_inversions"](m) == 0]
        incs = [np.diff(m) for m in mks if len(m) >= 2]
        grand, lo, hi = T["cluster_bootstrap_mean"](incs)
        med_cycles = float(np.median([r["cycles"] for r in runs])) if runs else float("nan")
        # spiral-out slope distribution over sub-theta approach peaks (verbatim logic)
        slopes = []
        for m in mks:
            sub_idx = np.where(m < THETA)[0]
            if len(sub_idx) < 3:
                continue
            y = np.log(THETA - m[sub_idx])
            slopes.append(float(linregress(sub_idx.astype(float), y).slope))
        rows.append(dict(
            N=N, runs=len(runs), testable=len(testable),
            ratchet_frac=(len(ratchet) / len(testable) if testable else float("nan")),
            strict_frac=(len(strict) / len(testable) if testable else float("nan")),
            mean_increment=grand, ci_lo=lo, ci_hi=hi, median_cycles=med_cycles,
            spiral_slope_median=(float(np.median(slopes)) if slopes else float("nan")),
            spiral_n_fit=len(slopes),
        ))
    return rows


ratchet_post = ratchet_for("0p5", [8, 16, 32, 64])
ratchet_depth = ratchet_for("0p444", [8, 16, 32, 64])

# --------------------------------------------------------------------------- #
# Assemble + print
# --------------------------------------------------------------------------- #
def tau_plateau(A_val):
    rows = [r for r in surv if r["A"] == A_val]
    rows.sort(key=lambda r: r["N"])
    return rows


def jget(d):
    return json.loads(json.dumps(d, default=float))


results = {
    "beta": 0.10,
    "campaign": CFG["campaign"],
    "n_runs": len(A["camp"]),
    "survival": jget(surv),
    "plateau_fits": jget(curves),
    "weibull": jget(weib),
    "geometric": jget(geo),
    "ratchet_A_post": jget(ratchet_post),
    "ratchet_A_depth": jget(ratchet_depth),
    "phase": {
        "pooled_all": {k: phase_pooled[k] for k in ["n", "mean_phase", "Rbar", "z", "p"]},
        "pooled_post": {k: phase_pooled05[k] for k in ["n", "mean_phase", "Rbar", "z", "p"]},
        "pooled_post_phis": [float(x) for x in np.asarray(phase_pooled05["phis"]).ravel()],
        "per_point": {f"A{a}_N{n}": {kk: vv for kk, vv in v.items() if kk != "phis"}
                      for (a, n), v in phase_agg.items()},
        "tb_crosscheck": tb_cross,
    },
}

print("=" * 78)
print(f"beta=0.10 ROBUSTNESS SLICE  ({CFG['campaign']}, {len(A['camp'])} runs)")
print("=" * 78)

for label, A_val in [("A_post (0.5)", A_POST), ("A_depth (0.444)", A_DEPTH), ("A_stable (0.2)", A_STABLE)]:
    rows = tau_plateau(A_val)
    if not rows:
        continue
    print(f"\n--- tau_abs(N), censored fraction @ {label} ---")
    for r in rows:
        print(f"  N={r['N']:2d}: tau_abs={r['tau_abs']:.1f}s "
              f"[{r['tau_abs_lo']:.1f},{r['tau_abs_hi']:.1f}]  km_abs={r['km_abs']:.1f}  "
              f"abs_cens_frac={r['abs_cens_frac']:.2f}  events={r['abs_events']}/{r['n']}")
    fit = curves[A_val]["fit_abs"]
    if fit:
        print(f"  plateau: tau(Nhi)/tau(Nlo)={fit['ratio']:.2f}  "
              f"exp-rate={fit['exp_c']:+.4e}/osc  AIC(exp)={fit['aic_exp']:.1f} AIC(pow)={fit['aic_pow']:.1f}")

print("\n--- censored-Weibull k_abs(N) (aging: k>1) ---")
for r in weib:
    if r["A"] in (A_POST, A_DEPTH):
        tag = "post" if r["A"] == A_POST else "depth"
        print(f"  [{tag}] N={r['N']:2d}: k_abs={r['k_abs']:.2f} [{r['k_abs_lo']:.2f},{r['k_abs_hi']:.2f}]  events={r['abs_events']}")

print("\n--- ratchet ⟨ΔM⟩ per cycle (cluster-bootstrap 95% CI) + spiral slope ---")
for tag, rr in [("A_post", ratchet_post), ("A_depth", ratchet_depth)]:
    print(f"  [{tag}]")
    for r in rr:
        print(f"    N={r['N']:2d}: ⟨ΔM⟩={r['mean_increment']:+.4f} [{r['ci_lo']:+.4f},{r['ci_hi']:+.4f}]  "
              f"ratchet_frac={r['ratchet_frac']:.2f} strict={r['strict_frac']:.2f}  "
              f"med_cycles={r['median_cycles']:.0f}  spiral_slope_med={r['spiral_slope_median']:+.3f} (n={r['spiral_n_fit']})")

print("\n--- breath-phase locking (true absorptions, Rayleigh) ---")
pp = phase_pooled05
print(f"  pooled post-homoclinic: n={pp['n']} Rbar={pp['Rbar']:.3f} mean_phase={pp['mean_phase']:.3f} p={pp['p']:.2e}")
pa = phase_pooled
print(f"  pooled all points:      n={pa['n']} Rbar={pa['Rbar']:.3f} mean_phase={pa['mean_phase']:.3f} p={pa['p']:.2e}")
for (a, n), v in sorted(phase_agg.items()):
    print(f"    A={a} N={n:2d}: n={v['n']} Rbar={v['Rbar']:.3f} p={v['p']:.2e}")
if tb_cross:
    print(f"  T_b crosscheck (JS vs Py): n={tb_cross['n']} median_rel_dev={tb_cross['median_rel_dev']*100:.2f}% corr={tb_cross['corr']:.4f}")

print("\n--- A_stable (0.2) control ---")
for r in tau_plateau(A_STABLE):
    print(f"  N={r['N']:2d}: abs_cens_frac={r['abs_cens_frac']:.2f} (never-absorber check)  "
          f"km_surv@2000={r['surv_abs_2000']:.2f}  events={r['abs_events']}/{r['n']}")

(OUT / "beta010_results.json").write_text(json.dumps(results, indent=2, default=float))
print(f"\nWrote {OUT / 'beta010_results.json'}")
