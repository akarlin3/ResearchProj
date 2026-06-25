#!/usr/bin/env python3
"""Regenerate the Procrustes manuscript figures from the seeded gate results.

Reads ../../results/phase{1,3}_*.json and writes two PDFs into this directory:

  fig_distinctness.pdf -- the conditional gap of D by true-D* regime (all / well-ID /
      strict-lo / high-D*), with bootstrap CIs. Shows the gap intensifying inside the
      companion study's identifiable region ("trust D") and vanishing in its own
      high-D* wall -- the opposite-to-identifiability gradient.
  fig_robustness.pdf   -- the well-identified gap vs SNR, with the SNR=25 honest
      boundary marked, plus the b-scheme bars.

Non-fatal in the build: if this fails, the committed figures/*.pdf are used.
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "..", "results")


def _load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


def fig_distinctness(p1, p3):
    st = p1["families"]["stretched"]
    h16 = p3["headline16"]
    # 16-seed for all/well/strict-lo; 8-seed phase1 for the high-D* wall contrast.
    rows = [
        ("all\nvoxels", h16["gap_all"], "#7f7f7f"),
        ("well-ID\n(bottom-2 D*)\n= identifiable", h16["gap_wellid"], "#3b6ea5"),
        ("strict-lo\n(bottom D*)\ndeepest", h16["gap_lo"], "#2a5783"),
        ("high-D*\n(top tercile)\n= Gauge wall", st["gap_hi"], "#c0504d"),
    ]
    labels = [r[0] for r in rows]
    pts = [r[1]["point"] for r in rows]
    los = [r[1]["point"] - r[1]["lo"] for r in rows]
    his = [r[1]["hi"] - r[1]["point"] for r in rows]
    colors = [r[2] for r in rows]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = range(len(rows))
    ax.bar(x, pts, yerr=[los, his], capsize=5, color=colors, width=0.62,
           edgecolor="black", linewidth=0.6)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("conditional gap of $D$\n(placebo $-$ worst-departure coverage)")
    ax.set_title("Misspecification breaks $D$-coverage INSIDE the identifiable region\n"
                 "(the gap intensifies where the companion says \"trust $D$\", "
                 "and vanishes in its own high-$D^{*}$ wall)", fontsize=9.5)
    for xi, p in zip(x, pts):
        ax.annotate(f"{p:+.3f}", (xi, p), textcoords="offset points",
                    xytext=(0, 6 if p >= 0 else -12), ha="center", fontsize=8)
    fig.tight_layout()
    out = os.path.join(HERE, "fig_distinctness.pdf")
    fig.savefig(out); plt.close(fig)
    print(f"[figures] wrote {out}")


def fig_robustness(p3):
    snr = p3["snr"]
    snrs = sorted(int(k) for k in snr)
    well = [snr[str(s)]["gap_wellid"]["point"] for s in snrs]
    lo = [snr[str(s)]["gap_wellid"]["point"] - snr[str(s)]["gap_wellid"]["lo"] for s in snrs]
    hi = [snr[str(s)]["gap_wellid"]["hi"] - snr[str(s)]["gap_wellid"]["point"] for s in snrs]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    a1.errorbar(snrs, well, yerr=[lo, hi], marker="o", capsize=4, color="#3b6ea5", lw=1.6)
    a1.axhline(0, color="black", lw=0.8, ls="--")
    a1.axvspan(20, 30, color="#c0504d", alpha=0.12)
    a1.annotate("SNR=25:\nseparation\nfails (noise-\ndominated)", (25, min(well)),
                fontsize=7.5, ha="center", color="#8a2a25",
                xytext=(25, max(well) * 0.45), arrowprops=dict(arrowstyle="->", color="#8a2a25"))
    a1.set_xlabel("Rician SNR")
    a1.set_ylabel("well-identified gap of $D$")
    a1.set_title("Separation intensifies with SNR;\nfails only at SNR=25 (honest boundary)",
                 fontsize=9)

    bs = p3["bscheme"]
    names = ["default-22pt", "clinical-sparse-10pt", "lowb-poor-8pt"]
    short = ["default\n22-pt", "clinical\nsparse 10-pt", "low-b-poor\n8-pt"]
    pts = [bs[n]["gap_wellid"]["point"] for n in names]
    blo = [bs[n]["gap_wellid"]["point"] - bs[n]["gap_wellid"]["lo"] for n in names]
    bhi = [bs[n]["gap_wellid"]["hi"] - bs[n]["gap_wellid"]["point"] for n in names]
    a2.bar(range(3), pts, yerr=[blo, bhi], capsize=5, color="#3b6ea5",
           edgecolor="black", linewidth=0.6, width=0.6)
    a2.axhline(0, color="black", lw=0.8)
    a2.set_xticks(range(3)); a2.set_xticklabels(short, fontsize=8)
    a2.set_ylabel("well-identified gap of $D$")
    a2.set_title("Robust across acquisition schemes\n(all survive; weakest when low-b is sparse)",
                 fontsize=9)
    fig.tight_layout()
    out = os.path.join(HERE, "fig_robustness.pdf")
    fig.savefig(out); plt.close(fig)
    print(f"[figures] wrote {out}")


def main():
    p1 = _load("phase1_gateB.json")
    p3 = _load("phase3_gateD.json")
    fig_distinctness(p1, p3)
    fig_robustness(p3)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
