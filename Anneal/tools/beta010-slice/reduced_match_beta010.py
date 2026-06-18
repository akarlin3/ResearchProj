"""beta=0.10 reduced-vs-measured comparison at the post-homoclinic corner.

Reuses tools/reduced-ode/reduced_core.py (reduced_run_3d) on the SAME seed-mapped
collective ICs as section 6.4 (rho1=Rsync0, rho2=Rincoh0, psi=dphi0 — a pure
function of the seed, so the beta=0.05 cp2_features ICs apply unchanged at
beta=0.10). For each (N,seed) it integrates the 3-variable reduced flow at
beta=0.10 and compares:
  - breath period (reduced vs measured T_b),
  - pooled spiral slope (reduced vs measured M_k spiral),
  - median capture-time ratio measured/reduced (the ~x3 prolongation factor),
  - whether that prolongation factor is N-independent,
  - A_stable control: persistent-run mean R_incoh vs the reduced stable r*.

Writes absorption_results/beta010_reduced_match.json and prints the numbers.
Run: python3 tools/beta010-slice/reduced_match_beta010.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools/reduced-ode"))
import reduced_core as rc  # noqa: E402

CFG = rc.load_config()
rc.set_config(CFG)
BETA = 0.10
A_POST = 0.5
A_DEPTH = 0.444
A_STABLE = 0.2
NS = [8, 16, 32, 64]
T_MAX = 2000.0
OUT = ROOT / "absorption_results"


def read_jsonl(p):
    return [json.loads(s) for s in Path(p).read_text().splitlines() if s.strip()]


# ICs per (N, seed): pure function of the seed (A/beta-independent), reused from
# the section-6.4 feature file. Restrict to the beta=0.10 seed family 100000+.
feat = read_jsonl(ROOT / "transient_results/cp2_features.jsonl")
ics = {}
for r in feat:
    if r["N"] in NS and 100000 <= r["seed"] < 100100:
        ics[(r["N"], r["seed"])] = (r["Rsync0"], r["Rincoh0"], r["dphi0"])

# Measured t_abs and T_b at beta=0.10 (campaign rows).
camp = read_jsonl(OUT / "absorption_campaign_beta010.jsonl")
measured = {}  # (A,N,seed) -> row
for r in camp:
    measured[(r["A"], r["N"], r["seed"])] = r


def reduced_capture(A, N, seed):
    state0 = list(ics[(N, seed)])
    res = rc.reduced_run_3d(state0, rc.Params(A=A, beta=BETA), T_MAX, CFG)
    return res


def match_corner(A):
    per_N = {}
    pooled_ratio = []
    for N in NS:
        red_t, red_Tb, red_sl = [], [], []
        meas_t, meas_Tb = [], []
        ratios = []
        for (n, seed), ic in ics.items():
            if n != N:
                continue
            mrow = measured.get((A, N, seed))
            if mrow is None or mrow["abs_censored"]:
                continue
            res = reduced_capture(A, N, seed)
            if not res.get("captured") or res.get("t_capture") is None:
                continue
            rt = res["t_capture"]
            red_t.append(rt)
            if res.get("breath_period"):
                red_Tb.append(res["breath_period"])
            if res.get("spiral_slope") is not None:
                red_sl.append(res["spiral_slope"])
            meas_t.append(mrow["t_abs"])
            if mrow.get("T_b"):
                meas_Tb.append(mrow["T_b"])
            ratios.append(mrow["t_abs"] / rt if rt > 0 else np.nan)
        ratios = [x for x in ratios if np.isfinite(x)]
        pooled_ratio += ratios
        per_N[N] = dict(
            n=len(red_t),
            reduced_capture_med=float(np.median(red_t)) if red_t else float("nan"),
            measured_tabs_med=float(np.median(meas_t)) if meas_t else float("nan"),
            ratio_med=float(np.median(ratios)) if ratios else float("nan"),
            reduced_Tb_med=float(np.median(red_Tb)) if red_Tb else float("nan"),
            measured_Tb_med=float(np.median(meas_Tb)) if meas_Tb else float("nan"),
            reduced_spiral_med=float(np.median(red_sl)) if red_sl else float("nan"),
        )
    rN = [per_N[N]["ratio_med"] for N in NS if np.isfinite(per_N[N]["ratio_med"])]
    n_indep = dict(
        ratio_per_N={N: per_N[N]["ratio_med"] for N in NS},
        ratio_pooled_med=float(np.median(pooled_ratio)) if pooled_ratio else float("nan"),
        ratio_spread=(float(max(rN) - min(rN)) if rN else float("nan")),
        ratio_cv=(float(np.std(rN) / np.mean(rN)) if rN and np.mean(rN) else float("nan")),
    )
    return per_N, n_indep


post_perN, post_nindep = match_corner(A_POST)
depth_perN, depth_nindep = match_corner(A_DEPTH)

# A_stable control: reduced stable r* vs persistent-run mean R_incoh.
rstar = None
try:
    fp = rc.get_fixed_point(A_STABLE, BETA, branch="chimera")
    rstar = float(fp[0]) if fp else None
except Exception:
    rstar = None
stable = read_jsonl(OUT / "stable_level_beta010.jsonl")
stable_ctrl = {}
for N in [8, 64]:
    persist = [r["mean_Rincoh_steady"] for r in stable if r["N"] == N and r["abs_censored"]]
    stable_ctrl[N] = dict(
        n_persistent=len(persist),
        mean_Rincoh_persistent=float(np.mean(persist)) if persist else float("nan"),
        reduced_rstar=rstar,
    )

results = dict(beta=BETA, A_hc=0.35372,
               A_post=dict(per_N=post_perN, n_independence=post_nindep),
               A_depth=dict(per_N=depth_perN, n_independence=depth_nindep),
               A_stable=dict(reduced_rstar=rstar, control=stable_ctrl))
(OUT / "beta010_reduced_match.json").write_text(json.dumps(results, indent=2, default=float))

print("=" * 78)
print("beta=0.10 REDUCED vs MEASURED  (post-homoclinic corners)")
print("=" * 78)
for tag, perN, nind in [("A_post=0.5", post_perN, post_nindep), ("A_depth=0.444", depth_perN, depth_nindep)]:
    print(f"\n--- {tag} ---")
    print("  N | n  | reduced_cap | measured_tabs | ratio(meas/red) | red_Tb | meas_Tb | red_spiral")
    for N in NS:
        d = perN[N]
        print(f"  {N:2d}| {d['n']:3d}| {d['reduced_capture_med']:11.1f} | {d['measured_tabs_med']:13.1f} | "
              f"{d['ratio_med']:15.2f} | {d['reduced_Tb_med']:6.1f} | {d['measured_Tb_med']:7.1f} | {d['reduced_spiral_med']:+.3f}")
    print(f"  pooled prolongation factor = {nind['ratio_pooled_med']:.2f}  "
          f"(per-N spread={nind['ratio_spread']:.2f}, CV={nind['ratio_cv']:.3f}) "
          f"-> {'N-INDEPENDENT' if nind['ratio_cv'] < 0.15 else 'N-DEPENDENT'}")

print(f"\n--- A_stable=0.2 control: reduced stable r* = {rstar:.4f} ---")
for N in [8, 64]:
    d = stable_ctrl[N]
    print(f"  N={N:2d}: persistent={d['n_persistent']}  mean_Rincoh(steady)={d['mean_Rincoh_persistent']:.4f}  vs r*={rstar:.4f}")
print(f"\nWrote {OUT / 'beta010_reduced_match.json'}")
