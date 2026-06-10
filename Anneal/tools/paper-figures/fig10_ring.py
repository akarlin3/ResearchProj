"""Figure 10 — ring chimera collapse: survival inversion, aging shape, and its
saturating thermodynamic limit.

Three panels: (a) Kaplan-Meier survival by N at beta=0.130 (median lifetime
decreases with N in the near-boundary regime); (b) censored-Weibull shape
k_hat(N) per beta with profile-likelihood intervals and the AICc-selected
bounded fits k_inf - a/N; (c) bootstrapped k_inf versus beta, all intervals
excluding the memoryless k=1 and rising toward the existence boundary.

Pure rendering of committed analysis outputs — no integration, no refitting:
reads anneal-hazard/results/{cp4_fits.json,cp_fits_N192.json,cp_fits_N256.json,
ensemble.csv,ensemble_N192.csv,ensemble_N256.csv} and
anneal-hazard/results/extrapolation/cpB_n5_fits.json.
Writes paper_figures/fig10.{pdf,png}.
Run: python3 tools/paper-figures/fig10_ring.py
"""
import csv
import json
import os

import numpy as np

import style

style.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

ROOT = style.ROOT
P = style.PALETTE
RES = os.path.join(ROOT, "anneal-hazard", "results")

BETAS = ["0.110", "0.115", "0.120", "0.125", "0.130"]
NS = [32, 64, 128, 192, 256]
# one palette colour per beta, cool -> warm toward the existence boundary
BCOL = {"0.110": P["blue"], "0.115": P["skyblue"], "0.120": P["green"],
        "0.125": P["orange"], "0.130": P["vermillion"]}


def cell_key(beta, N):
    # cp4_fits.json uses trimmed beta strings (b0.11, b0.115, ...); the N=192/256
    # files use the zero-padded form (b0.110, ...).
    short = f"b{float(beta):g}_N{N}"
    return short


def load_cells():
    cells = {}
    for fname in ("cp4_fits.json", "cp_fits_N192.json", "cp_fits_N256.json"):
        with open(os.path.join(RES, fname)) as f:
            d = json.load(f)
        for key, v in d.items():
            if key.startswith("_"):
                continue
            beta = f"{v['beta']:.3f}" if "beta" in v else None
            if beta is None:  # cp4_fits.json has no beta field; parse the key
                b, n = key[1:].split("_N")
                beta = f"{float(b):.3f}"
                N = int(n)
            else:
                N = v["N"]
            cells[(beta, N)] = v
    return cells


def km_curve(tau, event, t_max=12000.0):
    """Kaplan-Meier estimator (right-censored). Returns step (t, S) arrays."""
    order = np.argsort(tau)
    tau, event = np.asarray(tau)[order], np.asarray(event)[order]
    n = len(tau)
    t_pts, s_pts = [0.0], [1.0]
    S = 1.0
    at_risk = n
    i = 0
    while i < n:
        t = tau[i]
        d = c = 0
        while i < n and tau[i] == t:
            if event[i] == 1:
                d += 1
            else:
                c += 1
            i += 1
        if d > 0:
            S *= 1.0 - d / at_risk
            t_pts.append(t)
            s_pts.append(S)
        at_risk -= d + c
    return np.asarray(t_pts), np.asarray(s_pts)


def load_ensemble_b013():
    """(tau, event) rows at beta=0.130 keyed by N, from the three ensemble files."""
    out = {N: ([], []) for N in NS}
    for fname in ("ensemble.csv", "ensemble_N192.csv", "ensemble_N256.csv"):
        with open(os.path.join(RES, fname)) as f:
            for row in csv.DictReader(f):
                if abs(float(row["beta"]) - 0.130) > 1e-9:
                    continue
                N = int(row["N"])
                out[N][0].append(float(row["tau"]))
                out[N][1].append(int(float(row["event"])))
    return {N: v for N, v in out.items() if v[0]}


cells = load_cells()
with open(os.path.join(RES, "extrapolation", "cpB_n5_fits.json")) as f:
    extrap = json.load(f)

fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))

# (a) KM survival at beta = 0.130 by N
ax = axes[0]
ens = load_ensemble_b013()
ncol = {32: P["blue"], 64: P["skyblue"], 128: P["green"], 192: P["orange"],
        256: P["vermillion"]}
for N in NS:
    if N not in ens:
        continue
    t, S = km_curve(*ens[N])
    ax.step(t, S, where="post", color=ncol[N], label=f"N = {N}")
ax.set_xscale("log")
ax.set_xlim(50, 6000)
ax.set_ylim(0, 1.02)
ax.set_xlabel("time  (model units)")
ax.set_ylabel(r"survival $\hat S(t)$")
ax.set_title(r"(a) Ring survival, $\beta = 0.130$")
ax.legend(fontsize=8)

# (b) k_hat(N) per beta with bounded fits
ax = axes[1]
for beta in BETAS:
    ks, los, his = [], [], []
    for N in NS:
        v = cells[(beta, N)]
        ks.append(v["weibull_k"])
        los.append(v["weibull_k"] - v["weibull_k_lo"])
        his.append(v["weibull_k_hi"] - v["weibull_k"])
    ax.errorbar(NS, ks, yerr=[los, his], color=BCOL[beta], marker="o", ms=3.5,
                ls="none", capsize=2, label=rf"$\beta = {beta}$")
    b2 = extrap[beta]["bounded2"]
    Ngrid = np.geomspace(28, 320, 200)
    ax.plot(Ngrid, b2["k_inf"] - b2["a"] / Ngrid, color=BCOL[beta], lw=1.0,
            alpha=0.7)
ax.axhline(1.0, color=P["black"], ls=":", lw=0.8)
ax.set_xscale("log", base=2)
ax.set_xticks(NS)
ax.set_xticklabels(NS)
ax.set_xlabel("ring size N")
ax.set_ylabel(r"Weibull shape $\hat k$")
ax.set_title("(b) Aging shape: rise, then saturation")
ax.legend(fontsize=7.5, ncol=2, loc="lower right")

# (c) bootstrapped k_inf vs beta
ax = axes[2]
bx = [float(b) for b in BETAS]
kinf = [extrap[b]["boot_kinf_bounded2"]["median"] for b in BETAS]
lo = [extrap[b]["boot_kinf_bounded2"]["median"]
      - extrap[b]["boot_kinf_bounded2"]["lo"] for b in BETAS]
hi = [extrap[b]["boot_kinf_bounded2"]["hi"]
      - extrap[b]["boot_kinf_bounded2"]["median"] for b in BETAS]
ax.errorbar(bx, kinf, yerr=[lo, hi], color=P["purple"], marker="s", ms=4,
            capsize=3, lw=1.2)
ax.axhline(1.0, color=P["black"], ls=":", lw=0.8, label=r"$k = 1$ (memoryless)")
ax.axvspan(0.13, 0.14, color=P["grey"], alpha=0.25, lw=0,
           label=r"$\beta_c \approx 0.13$–$0.14$")
ax.set_xlim(0.1075, 0.1355)
ax.set_xlabel(r"phase-lag parameter $\beta$")
ax.set_ylabel(r"limiting shape $k_\infty$")
ax.set_title(r"(c) Bootstrap $k_\infty$: $> 1$, rising toward $\beta_c$")
ax.legend(fontsize=8, loc="lower right")

fig.tight_layout()
paths = style.savefig(fig, "fig10")
print("wrote", *paths)
