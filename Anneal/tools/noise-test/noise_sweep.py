"""CP2 — finite-size noise test on the reduced flow (A=0.5, beta=0.05).

Euler-Maruyama sweep over N in {8,16,32,64} (sigma_noise = c/sqrt(N)) x
c in {0, 0.025, 0.05, 0.1, 0.2} (c=0 control; the others bracket the measured
finite-N R_incoh fluctuation scale c~0.05) x 200 noise realizations per cell.
ICs are the same seed-mapped collective coordinates used in section 6.4
(transient_results/cp2_features.jsonl, A=0.5).

Per cell reports: median capture time, prolongation ratio vs the deterministic
reference, Weibull-in-time shape k (rising-hazard constraint), and the Rayleigh p
of capture phase (breath-phase-locking constraint). Tests whether any (c,N) cell
reproduces BOTH (i) prolongation ~3 and (ii) N-independence of that factor.

Run: python3 tools/noise-test/noise_sweep.py
"""
import json
import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools/reduced-ode"))
sys.path.insert(0, str(ROOT / "tools/beta010-slice"))
sys.path.insert(0, str(HERE))
import reduced_core as rc  # noqa: E402
from em_core import em_run  # noqa: E402
from _reuse import load_funcs  # noqa: E402

CFG = rc.load_config()
rc.set_config(CFG)
A = 0.5
BETA = 0.05
NS = [8, 16, 32, 64]
CS = [0.0, 0.025, 0.05, 0.1, 0.2]
C_PHYS = 0.05
N_REAL = 200
DT = 0.01
T_MAX = 2000.0
OUT = ROOT / "noise_results"
OUT.mkdir(exist_ok=True)

ABS = load_funcs(ROOT / "tools/absorption-recampaign/analysis.py")
fit_weibull = ABS["fit_weibull"]
rayleigh = ABS["rayleigh"]

# Seed-mapped collective ICs (section 6.4), A=0.5.
feat = [json.loads(s) for s in (ROOT / "transient_results/cp2_features.jsonl").read_text().splitlines() if s.strip()]
ICS = {}
for N in NS:
    rows = sorted([r for r in feat if r["N"] == N], key=lambda r: r["seed"])[:N_REAL]
    ICS[N] = [(r["Rsync0"], r["Rincoh0"], r["dphi0"]) for r in rows]


def run_cell(args):
    ci, c, N = args
    ics = ICS[N]
    recs = []
    for j in range(N_REAL):
        ic = ics[j % len(ics)]
        seed = 7_000_000 + ci * 100_000 + N * 1000 + j
        r = em_run(list(ic), A, BETA, c, N, CFG, seed=seed, dt=DT, t_max=T_MAX)
        recs.append((r["t_capture"], r["n_cycles"], r["capture_phase"],
                     r["breath_period"], r["spiral_slope"], r["captured"]))
    return (c, N, recs)


def deterministic_ref():
    """Deterministic (DOP853) capture per N: per-IC array + median."""
    ref, per_ic = {}, {}
    for N in NS:
        caps = []
        for ic in ICS[N]:
            res = rc.reduced_run_3d(list(ic), rc.Params(A=A, beta=BETA), T_MAX, CFG)
            caps.append(res["t_capture"] if res["captured"] else np.nan)
        per_ic[N] = np.array(caps, float)
        ref[N] = float(np.nanmedian(per_ic[N])) if np.isfinite(per_ic[N]).any() else float("nan")
    return ref, per_ic


def main():
    print("Computing deterministic reference (DOP853) per N ...")
    det_ref, det_per_ic = deterministic_ref()
    for N in NS:
        print(f"  N={N:2d}: deterministic median capture = {det_ref[N]:.2f}s")

    cells = [(ci, c, N) for ci, c in enumerate(CS) for N in NS]
    print(f"\nRunning {len(cells)} cells x {N_REAL} realizations (EM dt={DT}) ...")
    with Pool(min(10, len(cells))) as pool:
        out = pool.map(run_cell, cells)

    recs_by_cell = {(c, N): recs for c, N, recs in out}
    matrix = {}
    for c, N, recs in out:
        caps = np.array([r[0] for r in recs if r[5] and r[0] is not None], float)
        phis = np.array([r[2] for r in recs if r[2] is not None], float)
        ncyc = np.array([r[1] for r in recs if r[5]], float)
        cap_frac = float(np.mean([r[5] for r in recs]))
        med = float(np.median(caps)) if len(caps) else float("nan")
        ratio = med / det_ref[N] if det_ref[N] and np.isfinite(det_ref[N]) else float("nan")
        wk = fit_weibull(np.maximum(caps, 0.05), np.ones(len(caps), int))["k"] if len(caps) >= 5 else float("nan")
        n_ray, mp, Rbar, z, p = rayleigh(phis) if len(phis) else (0, np.nan, np.nan, np.nan, np.nan)
        matrix[(c, N)] = dict(
            c=c, N=N, n_captured=int(len(caps)), capture_frac=cap_frac,
            median_capture=med, ratio_to_det=ratio,
            weibull_k=float(wk), median_cycles=float(np.median(ncyc)) if len(ncyc) else float("nan"),
            rayleigh_n=int(n_ray), rayleigh_Rbar=float(Rbar), rayleigh_p=float(p),
        )

    # ---- print full matrix ----
    print("\n" + "=" * 92)
    print("CP2 NOISE-TEST MATRIX  (A=0.5, beta=0.05; prolongation ratio = median capture / deterministic)")
    print("=" * 92)
    print(f"{'c':>6} {'N':>4} {'cap_frac':>8} {'med_cap':>8} {'ratio':>7} {'weib_k':>7} {'med_cyc':>7} {'ray_p':>10} {'Rbar':>6}")
    for c in CS:
        for N in NS:
            m = matrix[(c, N)]
            print(f"{c:>6.3f} {N:>4d} {m['capture_frac']:>8.2f} {m['median_capture']:>8.1f} "
                  f"{m['ratio_to_det']:>7.2f} {m['weibull_k']:>7.2f} {m['median_cycles']:>7.1f} "
                  f"{m['rayleigh_p']:>10.2e} {m['rayleigh_Rbar']:>6.2f}")
        print("-" * 92)

    # ---- c=0 must match deterministic (PER-IC, apples-to-apples) ----
    print("\n[gate] c=0 column vs deterministic reference (per-IC matched):")
    c0_ok = True
    for N in NS:
        recs = recs_by_cell[(0.0, N)]
        det = det_per_ic[N]
        rels = []
        for j, r in enumerate(recs):
            tc, det_tc = r[0], det[j % len(det)]
            if r[5] and tc is not None and np.isfinite(det_tc) and det_tc > 0:
                rels.append(abs(tc - det_tc) / det_tc)
        rels = np.array(rels)
        med_rel = float(np.median(rels)) if len(rels) else float("nan")
        frac_lt1 = float(np.mean(rels < 0.01)) if len(rels) else float("nan")
        ok = med_rel < 0.01
        c0_ok = c0_ok and ok
        print(f"  N={N:2d}: median per-IC rel={med_rel*100:.3f}%  frac<1%={frac_lt1:.2f}  {'OK' if ok else 'FAIL'}")

    # ---- targets (i) prolongation~3 and (ii) N-independence ----
    print("\n[targets] (i) prolongation factor ~3 ; (ii) N-independent across N")
    target_hits = []
    for c in CS:
        if c == 0.0:
            continue
        ratios = [matrix[(c, N)]["ratio_to_det"] for N in NS]
        ratios = [x for x in ratios if np.isfinite(x)]
        med_ratio = float(np.median(ratios)) if ratios else float("nan")
        cv = float(np.std(ratios) / np.mean(ratios)) if ratios and np.mean(ratios) else float("nan")
        prolong3 = 2.0 <= med_ratio <= 4.0
        nindep = cv < 0.15
        ks = [matrix[(c, N)]["weibull_k"] for N in NS]
        ps = [matrix[(c, N)]["rayleigh_p"] for N in NS]
        aging = all(k > 1.0 for k in ks if np.isfinite(k))
        locked = all((p < 0.05) for p in ps if np.isfinite(p))
        print(f"  c={c:.3f}: median ratio={med_ratio:.2f} (CV={cv:.3f})  "
              f"prolong~3={prolong3}  N-indep={nindep}  aging(k>1 all N)={aging}  phase-locked(all N)={locked}")
        if prolong3 and nindep:
            target_hits.append(c)

    verdict = (f"POSITIVE: c={target_hits} reproduce BOTH prolongation~3 and N-independence"
               if target_hits else
               "NEGATIVE: no (c) cell reproduces BOTH prolongation~3 AND N-independence "
               "with additive 1/sqrt(N) noise on the collective flow")
    print(f"\n[VERDICT] {verdict}")

    results = dict(
        A=A, beta=BETA, dt=DT, n_realizations=N_REAL, c_phys=C_PHYS,
        cs=CS, Ns=NS, deterministic_ref=det_ref,
        c0_gate_pass=bool(c0_ok),
        matrix={f"c{c}_N{N}": matrix[(c, N)] for c in CS for N in NS},
        target_hits=target_hits, verdict=verdict,
    )
    (OUT / "noise_results.json").write_text(json.dumps(results, indent=2, default=float))
    print(f"\nWrote {OUT / 'noise_results.json'}")


if __name__ == "__main__":
    main()
