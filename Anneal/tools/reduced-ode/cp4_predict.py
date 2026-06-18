"""CP4 — run-by-run lifetime prediction (the capstone).

Correlate the reduced 3-D capture time with the measured t_abs, per N and pooled
(Spearman), and ask: does the parameter-free reduced model beat the static-feature
regression? Benchmarks: |dphi0|-alone Spearman ~0.52 and collective-IC CV R^2 ~0.23
(transient-tests). Also: partial correlation (does reduced add beyond |dphi0|?),
the N-trend of where the N-dependence enters (seed-ensemble ICs vs dynamics), and
the rho1=1 transverse stability along trajectories.

Reads reduced_results/reduced_runs.jsonl. Writes cp4_correlation.{json,md},
cp4_scatter.{png,pdf}, cp4_rho_vs_N.{png,pdf}.
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

import reduced_core as rc

CFG = rc.load_config()
rc.set_config(CFG)
OUT = os.path.join(rc.ROOT, CFG["output_dir"])
MEAS = CFG["measured"]
CV_SEED = 12345
KFOLDS = 5


def kfold_r2(X, y, seed=CV_SEED, k=KFOLDS):
    """Out-of-fold R^2 of OLS y ~ [1, X]. X is (n,p). Deterministic folds."""
    n = len(y)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    folds = np.array_split(idx, k)
    sse, sst = 0.0, 0.0
    ybar = y.mean()
    for i in range(k):
        te = folds[i]
        trn = np.concatenate([folds[j] for j in range(k) if j != i])
        Xtr = np.hstack([np.ones((len(trn), 1)), X[trn]])
        Xte = np.hstack([np.ones((len(te), 1)), X[te]])
        beta, *_ = np.linalg.lstsq(Xtr, y[trn], rcond=None)
        pred = Xte @ beta
        sse += np.sum((y[te] - pred) ** 2)
        sst += np.sum((y[te] - ybar) ** 2)
    return 1.0 - sse / sst if sst > 0 else float("nan")


def partial_spearman(a, b, c):
    """Spearman partial correlation of a,b controlling for c (rank-residual)."""
    from scipy.stats import rankdata

    ra, rb, rc_ = rankdata(a), rankdata(b), rankdata(c)

    def resid(x, z):
        Z = np.vstack([np.ones_like(z), z]).T
        beta, *_ = np.linalg.lstsq(Z, x, rcond=None)
        return x - Z @ beta

    return float(np.corrcoef(resid(ra, rc_), resid(rb, rc_))[0, 1])


def main():
    rows = [json.loads(l) for l in open(os.path.join(OUT, "reduced_runs.jsonl"))]
    Ns = sorted(set(r["N"] for r in rows))

    def arr(rs, k):
        return np.array([r[k] for r in rs], float)

    per_N = {}
    for N in Ns:
        rs = [r for r in rows if r["N"] == N and r["t_capture"] is not None]
        tc = arr(rs, "t_capture")
        ta = arr(rs, "t_abs_meas")
        ad = arr(rs, "absdphi0")
        logta = np.log(ta)
        rho_red = spearmanr(tc, ta).statistic
        rho_dphi = spearmanr(ad, ta).statistic
        psp = partial_spearman(tc, ta, ad)
        r2_red = kfold_r2(np.log(np.clip(tc, 0.1, None)).reshape(-1, 1), logta)
        r2_dphi = kfold_r2(np.column_stack([ad, ad ** 2]), logta)
        r2_both = kfold_r2(np.column_stack([np.log(np.clip(tc, 0.1, None)), ad, ad ** 2]), logta)
        per_N[str(N)] = {
            "n": len(rs),
            "spearman_reduced_vs_tabs": float(rho_red),
            "spearman_absdphi0_vs_tabs": float(rho_dphi),
            "partial_spearman_reduced_given_absdphi0": psp,
            "cvR2_reduced": float(r2_red),
            "cvR2_absdphi0_quad": float(r2_dphi),
            "cvR2_both": float(r2_both),
            "reduced_capture_median": float(np.median(tc)),
            "ic_mean_Rincoh0": float(np.mean(arr(rs, "Rincoh0"))),
            "ic_mean_absdphi0": float(np.mean(ad)),
            "transverse_attracting_frac_mean": float(
                np.mean(arr(rs, "transverse_attracting_frac"))
            ),
            "rho1_dev_max_median": float(np.median(arr(rs, "rho1_dev_max"))),
        }

    # pooled
    rs = [r for r in rows if r["t_capture"] is not None]
    tc = arr(rs, "t_capture"); ta = arr(rs, "t_abs_meas"); ad = arr(rs, "absdphi0")
    logta = np.log(ta)
    pooled = {
        "n": len(rs),
        "spearman_reduced_vs_tabs": float(spearmanr(tc, ta).statistic),
        "spearman_absdphi0_vs_tabs": float(spearmanr(ad, ta).statistic),
        "partial_spearman_reduced_given_absdphi0": partial_spearman(tc, ta, ad),
        "cvR2_reduced": float(kfold_r2(np.log(np.clip(tc, 0.1, None)).reshape(-1, 1), logta)),
        "cvR2_absdphi0_quad": float(kfold_r2(np.column_stack([ad, ad ** 2]), logta)),
        "cvR2_both": float(kfold_r2(np.column_stack([np.log(np.clip(tc, 0.1, None)), ad, ad ** 2]), logta)),
    }

    out = {
        "checkpoint": "CP4",
        "benchmarks": {
            "absdphi0_spearman": MEAS["dphi0_spearman"],
            "collective_ic_R2": MEAS["collective_ic_R2"],
            "d0_null_rho": MEAS["d0_null_rho"],
        },
        "per_N": per_N,
        "pooled": pooled,
        "system": "3-D reduced (rho1=Rsync0, rho2=Rincoh0, psi=dphi0); not 2-D",
    }
    with open(os.path.join(OUT, "cp4_correlation.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)

    # ---- scatter predicted (reduced) vs actual (measured), colored by N ----
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    cmap = plt.get_cmap("viridis")
    for i, N in enumerate(Ns):
        rs = [r for r in rows if r["N"] == N and r["t_capture"] is not None]
        ax.scatter(arr(rs, "t_capture"), arr(rs, "t_abs_meas"),
                   s=12, alpha=0.5, color=cmap(i / max(1, len(Ns) - 1)),
                   label=f"N={N}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("reduced 3-D capture time (s)")
    ax.set_ylabel("measured t_abs (s)")
    ax.set_title(f"Reduced-predicted vs measured lifetime (A=0.5, β=0.05)\n"
                 f"pooled Spearman ρ={pooled['spearman_reduced_vs_tabs']:.3f} "
                 f"(|Δφ₀| benchmark ρ≈{MEAS['dphi0_spearman']})")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "cp4_scatter.png"), dpi=140)
    fig.savefig(os.path.join(OUT, "cp4_scatter.pdf"))
    plt.close(fig)

    # ---- rho vs N + N-trend of where N enters ----
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    rr = [per_N[str(N)]["spearman_reduced_vs_tabs"] for N in Ns]
    rd = [abs(per_N[str(N)]["spearman_absdphi0_vs_tabs"]) for N in Ns]
    rp = [per_N[str(N)]["partial_spearman_reduced_given_absdphi0"] for N in Ns]
    ax.plot(Ns, rr, "o-", label="reduced t_capture")
    ax.plot(Ns, rd, "s--", label="|Δφ₀| (static)")
    ax.plot(Ns, rp, "^:", label="partial (reduced | |Δφ₀|)")
    ax.axhline(MEAS["dphi0_spearman"], color="grey", ls=":", lw=1,
               label="benchmark 0.52")
    ax.axhline(MEAS["d0_null_rho"], color="red", ls=":", lw=1, label="D₀ null")
    ax.set_xlabel("N"); ax.set_ylabel("Spearman ρ (vs measured t_abs)")
    ax.set_title("predictive power vs N")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    ax = axes[1]
    tcm = [per_N[str(N)]["reduced_capture_median"] for N in Ns]
    icr = [per_N[str(N)]["ic_mean_Rincoh0"] for N in Ns]
    ica = [per_N[str(N)]["ic_mean_absdphi0"] for N in Ns]
    ax.plot(Ns, tcm, "o-", color="C3", label="reduced capture median (s)")
    ax.set_xlabel("N"); ax.set_ylabel("reduced capture median (s)", color="C3")
    ax2 = ax.twinx()
    ax2.plot(Ns, icr, "s--", color="C0", label="mean R_incoh0 (IC)")
    ax2.plot(Ns, ica, "^--", color="C2", label="mean |Δφ₀| (IC)")
    ax2.set_ylabel("seed-ensemble IC stats")
    ax.axhline(MEAS["tau_abs_s"], color="grey", ls=":", label="measured τ_abs 139s")
    ax.set_title("where N enters: seed-ensemble ICs vs dynamics")
    l1, la1 = ax.get_legend_handles_labels()
    l2, la2 = ax2.get_legend_handles_labels()
    ax.legend(l1 + l2, la1 + la2, fontsize=7, loc="center right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "cp4_rho_vs_N.png"), dpi=140)
    fig.savefig(os.path.join(OUT, "cp4_rho_vs_N.pdf"))
    plt.close(fig)

    # ---- markdown ----
    L = ["# CP4 — run-by-run lifetime prediction\n",
         "Reduced **3-D** system integrated per run from (ρ₁=R_sync0, ρ₂=R_incoh0, "
         "ψ=Δφ₀). Capture = first sustained min(ρ₁,ρ₂)>θ (no-recovery). "
         "Correlated with measured t_abs.\n",
         f"**Pooled (n={pooled['n']}):** Spearman(reduced, t_abs)="
         f"**{pooled['spearman_reduced_vs_tabs']:.3f}**, "
         f"|Δφ₀| static={pooled['spearman_absdphi0_vs_tabs']:.3f}, "
         f"partial(reduced | |Δφ₀|)={pooled['partial_spearman_reduced_given_absdphi0']:.3f}.",
         f"CV R²: reduced={pooled['cvR2_reduced']:.3f}, "
         f"|Δφ₀|+quad={pooled['cvR2_absdphi0_quad']:.3f}, "
         f"both={pooled['cvR2_both']:.3f}  (collective-IC benchmark R²≈"
         f"{MEAS['collective_ic_R2']}).\n",
         "| N | n | ρ(reduced) | ρ(\\|Δφ₀\\|) | partial(red\\|Δφ₀) | R²(red) | "
         "R²(\\|Δφ₀\\|q) | R²(both) | cap med (s) | rho1 dev med |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for N in Ns:
        v = per_N[str(N)]
        L.append(
            f"| {N} | {v['n']} | {v['spearman_reduced_vs_tabs']:+.3f} | "
            f"{v['spearman_absdphi0_vs_tabs']:+.3f} | "
            f"{v['partial_spearman_reduced_given_absdphi0']:+.3f} | "
            f"{v['cvR2_reduced']:.3f} | {v['cvR2_absdphi0_quad']:.3f} | "
            f"{v['cvR2_both']:.3f} | {v['reduced_capture_median']:.0f} | "
            f"{v['rho1_dev_max_median']:.3f} |"
        )
    L += ["",
          "## Transverse stability of ρ₁=1",
          "Per-run max deviation of ρ₁ from 1 stays small (median column above); "
          "the ρ₁=1 manifold's local transverse rate is attracting for a majority "
          "of the trajectory at every N (transverse_attracting_frac mean: "
          + ", ".join(f"N={N}:{per_N[str(N)]['transverse_attracting_frac_mean']:.2f}"
                      for N in Ns) + ").",
          "",
          "## Where the N-dependence enters",
          "The reduced **dynamics** (A,β) are N-independent; the only N-channel is "
          "the **seed-ensemble IC** distribution. Mean R_incoh0 and mean |Δφ₀| of "
          "the canonical seed shift with N (see cp4_rho_vs_N right panel), and the "
          "reduced capture median tracks that shift — it does NOT reproduce the flat "
          "measured τ_abs(N) plateau. So the finite-N τ-plateau and the growing "
          "spiral-rate-with-N both live in the **IC ensemble + finite-N fluctuations**, "
          "not in the reduced collective dynamics.",
          ]
    with open(os.path.join(OUT, "cp4_correlation.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nfigures: cp4_scatter.{png,pdf}, cp4_rho_vs_N.{png,pdf}")


if __name__ == "__main__":
    main()
