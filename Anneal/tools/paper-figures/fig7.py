"""
Figure 7 — (β, A) stability diagram with operating corners (Abrams 2008 Fig. 4
style, rendered fresh from our reduced-ODE machinery).

Draws the saddle-node, Hopf, and homoclinic curves; shades and labels the four
regions (no chimera / stable chimera / breathing chimera / post-homoclinic);
marks the Takens–Bogdanov point where they meet; and marks BOTH operating
corners ((β=0.05, A=0.2) stable-chimera never-absorbers; (β=0.05, A=0.5)
post-homoclinic transient — this work). Reads the cached curves
(paper_figures/fig7_curves.json from fig7_curves.py); writes fig7.{pdf,png} and
fig7_caption.txt. Deterministic.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import style

ROOT = style.ROOT
CFG = style.load_config()
OUT = os.path.join(ROOT, CFG["output_dir"])
F7 = CFG["fig7"]


def load_curves():
    with open(os.path.join(OUT, "fig7_curves.json")) as f:
        return json.load(f)


def main():
    style.apply_style()
    c = load_curves()
    beta = np.array(c["beta_curve"])
    A_SN = np.array(c["A_SN_series"])
    A_H = np.array(c["A_H_series"])
    hc = c["homoclinic"]
    hcb = np.array([p["beta"] for p in hc])
    hcA = np.array([p["A_hc"] for p in hc])
    ax_cfg = F7["axis"]
    tb = c["tb_point"]

    fig, ax = plt.subplots(figsize=(6.6, 5.2))

    # homoclinic interpolated onto the dense β grid, restricted to its support
    hc_lo, hc_hi = hcb.min(), hcb.max()
    in_hc = (beta >= hc_lo) & (beta <= hc_hi)
    A_hc_dense = np.interp(beta, hcb, hcA)

    Amax = ax_cfg["A_max"]
    pal = style.PALETTE

    # ---- region shading (light fills) ----
    # no chimera (below saddle-node)
    ax.fill_between(beta, 0.0, A_SN, color=pal["grey"], alpha=0.16, lw=0, zorder=0)
    # stable stationary chimera (SN..Hopf)
    ax.fill_between(beta, A_SN, A_H, color=pal["green"], alpha=0.16, lw=0, zorder=0)
    # breathing chimera (Hopf..homoclinic) — only where homoclinic is defined
    ax.fill_between(beta, A_H, A_hc_dense, where=in_hc, color=pal["skyblue"],
                    alpha=0.20, lw=0, zorder=0)
    # post-homoclinic / sync only (above homoclinic)
    ax.fill_between(beta, A_hc_dense, Amax, where=in_hc, color=pal["orange"],
                    alpha=0.13, lw=0, zorder=0)

    # ---- curves ----
    # saddle-node: series as the curve, numeric as check marks
    ax.plot(beta, A_SN, color=style.ROLES["sn"], lw=1.8, zorder=4, label="saddle-node (Eq. 17)")
    snn = c["sn_numeric"]
    ax.plot([p["beta"] for p in snn], [p["A"] for p in snn], ls="none", marker="o",
            mfc="none", mec=style.ROLES["sn"], mew=1.1, ms=5, zorder=5)
    # Hopf
    ax.plot(beta, A_H, color=style.ROLES["hopf"], lw=1.8, ls="--", zorder=4,
            label="Hopf (Eq. 18)")
    hnn = c["hopf_numeric"]
    ax.plot([p["beta"] for p in hnn], [p["A"] for p in hnn], ls="none", marker="s",
            mfc="none", mec=style.ROLES["hopf"], mew=1.1, ms=5, zorder=5)
    # homoclinic: numeric bisection points + connecting curve
    ax.plot(hcb, hcA, color=style.ROLES["homoclinic"], lw=1.8, ls="-.", zorder=4,
            label="homoclinic (bisection)")
    ax.plot(hcb, hcA, ls="none", marker="D", mfc="none", mec=style.ROLES["homoclinic"],
            mew=1.1, ms=4.5, zorder=5)

    # ---- Takens–Bogdanov point ----
    ax.plot([tb["beta"]], [tb["A"]], marker="o", ms=8, mfc=style.ROLES["tb"],
            mec="white", mew=1.0, zorder=6)
    ax.annotate(f"Takens–Bogdanov\n({tb['beta']}, {tb['A']})",
                xy=(tb["beta"], tb["A"]), xytext=(tb["beta"] - 0.012, tb["A"] + 0.085),
                ha="center", fontsize=7.5,
                arrowprops=dict(arrowstyle="->", color=style.ROLES["tb"], lw=0.9))

    # ---- region labels ----
    ax.text(0.11, 0.055, "no chimera\n(below saddle-node)", ha="center", va="center",
            fontsize=8, color="#555555")
    ax.text(0.045, 0.215, "stable\nchimera", ha="center", va="center", fontsize=8,
            color=pal["green"])
    ax.text(0.105, 0.30, "breathing\nchimera", ha="center", va="center", fontsize=8,
            color=pal["blue"])
    ax.text(0.07, 0.47, "post-homoclinic  (sync only)", ha="center", va="center",
            fontsize=8.5, color=pal["vermillion"])

    # ---- operating corners ----
    corner_styles = {
        "stable": dict(marker="*", ms=15, color=style.ROLES["corner_stable"]),
        "transient": dict(marker="X", ms=11, color=style.ROLES["corner_transient"]),
    }
    for corner in c["corners"]:
        is_transient = abs(corner["A"] - 0.5) < 1e-9
        st = corner_styles["transient" if is_transient else "stable"]
        ax.plot([corner["beta"]], [corner["A"]], ls="none", mec="white", mew=1.0,
                zorder=7, **st)
        if is_transient:
            ax.annotate(corner["label"], xy=(corner["beta"], corner["A"]),
                        xytext=(corner["beta"] + 0.015, corner["A"] + 0.022),
                        fontsize=7.5, color=st["color"], ha="left",
                        arrowprops=dict(arrowstyle="->", color=st["color"], lw=0.9))
        else:
            ax.annotate(corner["label"], xy=(corner["beta"], corner["A"]),
                        xytext=(corner["beta"] + 0.015, corner["A"] - 0.045),
                        fontsize=7.5, color=st["color"], ha="left",
                        arrowprops=dict(arrowstyle="->", color=st["color"], lw=0.9))

    ax.set_xlim(ax_cfg["beta_min"], ax_cfg["beta_max"])
    ax.set_ylim(ax_cfg["A_min"], Amax)
    ax.set_xlabel(r"phase lag  $\beta$")
    ax.set_ylabel(r"coupling asymmetry  $A$")
    ax.set_title("Two-population chimera stability in the (β, A) plane", fontsize=9.5)

    # legend: curves + corner markers
    handles = [
        Line2D([], [], color=style.ROLES["sn"], lw=1.8, label="saddle-node (Eq. 17)"),
        Line2D([], [], color=style.ROLES["hopf"], lw=1.8, ls="--", label="Hopf (Eq. 18)"),
        Line2D([], [], color=style.ROLES["homoclinic"], lw=1.8, ls="-.",
               label="homoclinic (bisection)"),
        Line2D([], [], ls="none", marker="o", mfc="none", mec=pal["grey"],
               label="numeric check"),
        Line2D([], [], ls="none", **corner_styles["stable"], mec="white",
               label="A=0.2: stable chimera"),
        Line2D([], [], ls="none", **corner_styles["transient"], mec="white",
               label="A=0.5: post-homoclinic (this work)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7, framealpha=0.92,
              frameon=True, edgecolor="none")

    paths = style.savefig(fig, "fig7")
    plt.close(fig)

    sc = c["selfcheck"]
    dev = c["deviations"]
    caption = (
        f"FIG. 7. Stability of the two-population Sakaguchi–Kuramoto chimera in the "
        f"(β, A) plane, computed from the reduced order-parameter ODEs (Abrams et al., "
        f"2008): the saddle-node (Eq. 17) and Hopf (Eq. 18) series curves (numeric fold "
        f"and trace-zero locations overplotted as open markers; max deviation "
        f"{dev['A_SN_max']:.1e} and {dev['A_H_max']:.1e}) and the homoclinic curve "
        f"located by escape-to-sync bisection bound four regions — no chimera, stable "
        f"stationary chimera, breathing chimera, and post-homoclinic (sync globally "
        f"attracting) — meeting at the Takens–Bogdanov point ({c['tb_point']['beta']}, "
        f"{c['tb_point']['A']}). The two operating corners both sit at β=0.05: A=0.2 "
        f"(star) lies in the stable-chimera band (the never-absorbers), while A=0.5 (cross) "
        f"lies above the homoclinic (A_hc={sc['A_hc']:.4f}) in the post-homoclinic region, "
        f"where the chimera is only a transient and sync is globally attracting — the "
        f"regime studied in this work."
    )
    with open(os.path.join(OUT, "fig7_caption.txt"), "w") as f:
        f.write(caption + "\n")

    print("wrote:", *[os.path.relpath(p, ROOT) for p in paths])
    print("       " + os.path.relpath(os.path.join(OUT, "fig7_caption.txt"), ROOT))
    print("self-check:", sc["verdict"])


if __name__ == "__main__":
    main()
