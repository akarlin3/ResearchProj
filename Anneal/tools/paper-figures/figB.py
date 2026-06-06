"""Appendix figure B — finite-size noise test on the reduced flow.

Prolongation factor (median capture / deterministic) vs noise amplitude c, one
line per N. Reference lines: deterministic (ratio = 1, no prolongation) and the
measured finite-N prolongation the candidate must reproduce (~x3). Shows that
additive 1/sqrt(N) noise reaches the target ratio only at small N (N-dependent
prolongation) -- the clean negative.

Reads noise_results/noise_results.json. Writes paper_figures/figB.{pdf,png}.
Run: python3 tools/paper-figures/figB.py
"""
import json
import os

import numpy as np

import style

style.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

ROOT = style.ROOT
R = json.load(open(os.path.join(ROOT, "noise_results/noise_results.json")))
P = style.PALETTE

CS = R["cs"]
NS = R["Ns"]
M = R["matrix"]
MEAS_PROLONG = 3.2  # measured finite-N prolongation at beta=0.05 (139/43)

colors = {8: P["blue"], 16: P["green"], 32: P["orange"], 64: P["vermillion"]}
markers = {8: "o", 16: "s", 32: "^", 64: "D"}

fig, ax = plt.subplots(figsize=(6.4, 4.4))
x = np.array(CS)
for N in NS:
    y = [M[f"c{c}_N{N}"]["ratio_to_det"] for c in CS]
    ax.plot(x, y, color=colors[N], marker=markers[N], label=f"N = {N}")

ax.axhline(1.0, color=P["black"], ls="--", lw=0.9, label="deterministic (no prolongation)")
ax.axhline(MEAS_PROLONG, color=P["purple"], ls=":", lw=1.2,
           label=rf"measured finite-$N$ target ($\approx${MEAS_PROLONG:g})")

# annotate the measured physical noise scale
cphys = R["c_phys"]
ax.axvline(cphys, color=P["grey"], ls="-", lw=0.7)
ax.annotate(rf"measured scale $c\approx{cphys:g}$", xy=(cphys, ax.get_ylim()[1]),
            xytext=(cphys + 0.005, 2.7), fontsize=7, color=P["grey"])

ax.set_xlabel(r"noise amplitude $c$  (with $\sigma_{\mathrm{noise}}=c/\sqrt{N}$)")
ax.set_ylabel(r"prolongation factor  (median capture / deterministic)")
ax.set_title("Additive $1/\\sqrt{N}$ noise: prolongation is N-dependent")
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7,
          frameon=True, borderaxespad=0.0)
ax.set_ylim(0.5, max(3.4, max(M[f"c{c}_N8"]["ratio_to_det"] for c in CS) + 0.2))

fig.tight_layout()
paths = style.savefig(fig, "figB", tight=True)
print("wrote", *paths)
