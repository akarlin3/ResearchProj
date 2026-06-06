"""Appendix figure A — beta=0.10 robustness slice summary.

Three panels: (a) absorption lifetime tau_abs(N) at the post-homoclinic corners
(A=0.5 and the depth-matched A=0.444) showing the flat plateau; (b) censored-
Weibull aging shape k_abs(N) > 1; (c) breath-phase rose of true absorptions
(pooled post-homoclinic) with the Rayleigh result.

Reads absorption_results/beta010_results.json. Writes paper_figures/figA.{pdf,png}.
Run: python3 tools/paper-figures/figA.py
"""
import json
import os

import numpy as np

import style

style.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

ROOT = style.ROOT
R = json.load(open(os.path.join(ROOT, "absorption_results/beta010_results.json")))
P = style.PALETTE

NS = [8, 16, 32, 64]


def rows_for(A):
    return sorted([r for r in R["survival"] if abs(r["A"] - A) < 1e-6], key=lambda r: r["N"])


def weib_for(A):
    return sorted([r for r in R["weibull"] if abs(r["A"] - A) < 1e-6], key=lambda r: r["N"])


fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))

# (a) tau_abs(N) plateau
ax = axes[0]
for A, col, mk, lab in [(0.5, P["blue"], "o", "A = 0.5"),
                        (0.444, P["vermillion"], "s", "A = 0.444 (depth-matched)")]:
    rr = rows_for(A)
    Ns = [r["N"] for r in rr]
    tau = [r["tau_abs"] for r in rr]
    lo = [r["tau_abs"] - r["tau_abs_lo"] for r in rr]
    hi = [r["tau_abs_hi"] - r["tau_abs"] for r in rr]
    ax.errorbar(Ns, tau, yerr=[lo, hi], color=col, marker=mk, capsize=2, label=lab)
ax.set_xscale("log", base=2)
ax.set_xticks(NS); ax.set_xticklabels(NS)
ax.set_xlabel("oscillators per population N")
ax.set_ylabel(r"$\hat\tau_{\mathrm{abs}}$  (s)")
ax.set_title("(a) Absorption lifetime plateau")
ax.set_ylim(0, max(70, ax.get_ylim()[1]))
ax.legend()

# (b) Weibull k_abs(N)
ax = axes[1]
for A, col, mk, lab in [(0.5, P["blue"], "o", "A = 0.5"),
                        (0.444, P["vermillion"], "s", "A = 0.444")]:
    ww = weib_for(A)
    Ns = [r["N"] for r in ww]
    k = [r["k_abs"] for r in ww]
    lo = [r["k_abs"] - r["k_abs_lo"] for r in ww]
    hi = [r["k_abs_hi"] - r["k_abs"] for r in ww]
    ax.errorbar(Ns, k, yerr=[lo, hi], color=col, marker=mk, capsize=2, label=lab)
ax.axhline(1.0, color=P["black"], ls=":", lw=0.8, label="k = 1 (memoryless)")
ax.set_xscale("log", base=2)
ax.set_xticks(NS); ax.set_xticklabels(NS)
ax.set_xlabel("oscillators per population N")
ax.set_ylabel(r"Weibull shape $k_{\mathrm{abs}}$")
ax.set_title("(b) Aging: increasing hazard")
ax.legend()

# (c) breath-phase rose (pooled post-homoclinic true absorptions)
ax = axes[2]
ax.remove()
ax = fig.add_subplot(1, 3, 3, projection="polar")
phis = np.asarray(R["phase"]["pooled_post_phis"], float)
pp = R["phase"]["pooled_post"]
nb = 16
counts, edges = np.histogram(phis % (2 * np.pi), bins=nb, range=(0, 2 * np.pi))
width = 2 * np.pi / nb
ax.bar(edges[:-1] + width / 2, counts, width=width, color=P["purple"], alpha=0.8,
       edgecolor="white", linewidth=0.5)
# mean resultant vector
mp, Rbar = pp["mean_phase"], pp["Rbar"]
ax.annotate("", xy=(mp, Rbar * counts.max()), xytext=(0, 0),
            arrowprops=dict(color=P["black"], width=1.2, headwidth=6))
ax.set_theta_zero_location("E")
ax.set_yticklabels([])
ax.set_title(f"(c) Breath phase of absorptions\n"
             rf"$\bar R$={Rbar:.2f}, Rayleigh $p$={pp['p']:.1e} (n={pp['n']})", pad=14)

fig.suptitle(r"$\beta = 0.10$ robustness slice (post-homoclinic corner; $A_{hc}(0.10)=0.354$)",
             y=1.04, fontsize=11)
fig.tight_layout()
paths = style.savefig(fig, "figA")
print("wrote", *paths)
