"""Generate the railing-first manuscript figures DIRECTLY from frozen run outputs.

Every value plotted is read from a run artifact -- nothing is hand-entered:
  * Gnomon/results/reproduction.json          (K1 real railing; K2 quantile; K3 flow)
  * Gnomon/handoff/conditional_coverage.json  (R1 per-D*-tercile coverage, both conventions)
  * Sextant/results/railing_results.json       (cross-cohort generalization)

Outputs (PDF + PNG) into figures/manuscript/:
  fig1_railing_cohorts        K1 + Sextant generalization (railing rate by cohort, 95% CI)
  fig2_conditional_coverage   R1 honest-CRLB tercile coverage + floored-convention overlay
  fig3_resolution             K2 marginal quantile fix + K3 flow-vs-railed-NLLS
  graphical_abstract          compact railing-by-cohort panel, 50x60 mm

Run:  KMP_DUPLICATE_LIB_OK=TRUE conda run -n proteus python make_railing_figures.py
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                      # redraft worktree root
GNO = ROOT / "Gnomon"
SEX = ROOT / "Sextant"
OUT = HERE / "figures" / "manuscript"
OUT.mkdir(parents=True, exist_ok=True)

rep = json.loads((GNO / "results" / "reproduction.json").read_text())
cond = json.loads((GNO / "handoff" / "conditional_coverage.json").read_text())
sx = json.loads((SEX / "results" / "railing_results.json").read_text())

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                     "legend.fontsize": 7, "figure.dpi": 200, "savefig.bbox": "tight"})

C_TRUTH = "#1b4965"   # Gnomon source-of-truth
C_SX = "#5fa8d3"      # Sextant generalization
C_REF = "#bc4749"     # reference / threshold lines


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=200)
    plt.close(fig)
    print("wrote", stem)


# ---------------------------------------------------------------- Fig 1: railing by cohort
def fig1():
    real = rep["real"]
    g_pt = real["rail_tol_1e-3"]["rate"]
    g_lo, g_hi = real["rail_tol_1e-3"]["ci"]
    sxc = {c["name"]: c for c in sx["cohorts"]}

    def row(label, src, pt, lo, hi, color):
        return dict(label=label, src=src, pt=pt, lo=lo, hi=hi, color=color)

    rows = [row("OSIPI abdomen, homogeneous ROI\n(n=1932; rail_tol 1e-3)", "Gnomon (reference)",
               g_pt, g_lo, g_hi, C_TRUTH)]
    for key, lab in [("abdomen_homogeneous", "OSIPI abdomen, high-SNR subset\n(n=1618)"),
                     ("abdomen_full", "OSIPI abdomen, full ROI\n(n=19,652)"),
                     ("lihc_liver_4b_0_50_500_800", "TCGA-LIHC liver, 4-b clean\n(b 0/50/500/800)"),
                     ("lihc_liver_3b_50_400_800", "TCGA-LIHC liver, 3-b sparse\n(b 50/400/800)")]:
        ci = sxc[key]["bootstrap_ci"]
        rows.append(row(lab, "Sextant (generalization)", ci["point"], ci["lo"], ci["hi"], C_SX))

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    y = list(range(len(rows)))[::-1]
    for yi, r in zip(y, rows):
        ax.errorbar(r["pt"] * 100, yi,
                    xerr=[[(r["pt"] - r["lo"]) * 100], [(r["hi"] - r["pt"]) * 100]],
                    fmt="o", color=r["color"], ecolor=r["color"], capsize=3, ms=6)
        ax.text(r["pt"] * 100, yi + 0.18, f"{r['pt']*100:.1f}%", ha="center", va="bottom",
                fontsize=7, color=r["color"])
    ax.axvline(sx["meta"]["fashion_reported_homogeneous"] * 100, ls="--", lw=1, color=C_REF,
               label="prior report 54.7%")
    ax.axvline(30, ls=":", lw=1, color="gray", label="pre-registered replication floor (30%)")
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=6.6)
    ax.set_xlabel(r"NLLS $D^{*}$ railing rate (% of voxels, 95% CI)")
    ax.set_xlim(0, 85)
    ax.set_title("Boundary-railing of $D^{*}$ reproduces across independent in-vivo cohorts")
    # legend: source colours
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color=C_TRUTH, ls="", label="Gnomon (clean-room reference)"),
               Line2D([0], [0], marker="o", color=C_SX, ls="", label="Sextant (cross-cohort)"),
               Line2D([0], [0], ls="--", color=C_REF, label="prior report 54.7%"),
               Line2D([0], [0], ls=":", color="gray", label="replication floor 30%")]
    ax.legend(handles=handles, loc="lower left", framealpha=0.9)
    save(fig, "fig1_railing_cohorts")


# ------------------------------------------------- Fig 2: conditional coverage (R1)
def fig2():
    cc = cond["conditional_coverage"]
    terc = ["low_Dstar", "mid_Dstar", "high_Dstar"]
    tlab = ["low $D^{*}$", "mid $D^{*}$", "high $D^{*}$"]

    def series(est, conv):
        node = cc[est][conv]
        return ([node[t]["coverage"] for t in terc],
                [[node[t]["coverage"] - node[t]["ci"][0] for t in terc],
                 [node[t]["ci"][1] - node[t]["coverage"] for t in terc]])

    lap_h, lap_h_e = series("Laplace_SD", "honest")
    mc_h, mc_h_e = series("MCMC_SD", "honest")
    qu_h, qu_h_e = series("MCMC_quantile_recommended", "honest")
    lap_f, lap_f_e = series("Laplace_SD", "floored")

    import numpy as np
    x = np.arange(3)
    w = 0.22
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.bar(x - 1.5 * w, lap_h, w, yerr=lap_h_e, capsize=2, color="#a8c69f", label="Laplace SD (honest CRLB)")
    ax.bar(x - 0.5 * w, mc_h, w, yerr=mc_h_e, capsize=2, color="#4d908e", label="MCMC SD (honest CRLB)")
    ax.bar(x + 0.5 * w, qu_h, w, yerr=qu_h_e, capsize=2, color="#1b4965", label="MCMC quantile (recommended)")
    ax.bar(x + 1.5 * w, lap_f, w, yerr=lap_f_e, capsize=2, color="#d9d9d9", hatch="//",
           edgecolor="#888", label="Laplace SD (floored convention)")
    ax.axhline(0.95, ls="--", lw=1, color=C_REF, label="nominal 0.95")
    ax.set_xticks(x)
    ax.set_xticklabels(tlab)
    ax.set_ylabel("central-95% $D^{*}$ coverage")
    ax.set_ylim(0, 1.05)
    ax.set_title("The Gaussian failure is conditional: it concentrates in the high-$D^{*}$ tercile")
    ax.legend(loc="lower left", ncol=1, framealpha=0.9)
    save(fig, "fig2_conditional_coverage")


# ------------------------------------------------- Fig 3: resolution (K2 + K3)
def fig3():
    syn = rep["synthetic"]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.8, 3.2))

    # Panel A: marginal D* coverage by construction (pooled, honest) + quantile fix
    lap_pool = cond["conditional_coverage"]["Laplace_SD"]["honest"]["pooled"]
    mc_pool = cond["conditional_coverage"]["MCMC_SD"]["honest"]["pooled"]
    qd = syn["T3c_mcmc_quantile_Dstar"]
    names = ["Laplace\nSD", "MCMC\nSD", "MCMC\nquantile"]
    pts = [lap_pool["coverage"], mc_pool["coverage"], qd["coverage"]]
    errs = [[pts[0] - lap_pool["ci"][0], pts[1] - mc_pool["ci"][0], pts[2] - qd["ci"][0]],
            [lap_pool["ci"][1] - pts[0], mc_pool["ci"][1] - pts[1], qd["ci"][1] - pts[2]]]
    colors = ["#a8c69f", "#4d908e", "#1b4965"]
    axA.bar(range(3), pts, yerr=errs, capsize=3, color=colors)
    axA.axhline(0.95, ls="--", lw=1, color=C_REF)
    axA.set_xticks(range(3)); axA.set_xticklabels(names)
    axA.set_ylim(0, 1.05); axA.set_ylabel("marginal $D^{*}$ coverage")
    axA.set_title("(A) Interval shape restores\nmarginal coverage")
    for i, p in enumerate(pts):
        axA.text(i, p + 0.02, f"{p:.2f}", ha="center", fontsize=7)

    # Panel B: flow vs railed NLLS (coverage, ECE, sharpness)
    fl = syn["T4_flow_vs_railed_nlls"]["flow"]
    nl = syn["T4_flow_vs_railed_nlls"]["nlls_railed"]
    metrics = ["coverage", "ECE", "sharpness"]
    flv = [fl["coverage"], fl["ece"], fl["sharpness"]]
    nlv = [nl["coverage"], nl["ece"], nl["sharpness"]]
    import numpy as np
    x = np.arange(3); w = 0.36
    axB.bar(x - w / 2, nlv, w, color=C_REF, label="railed NLLS")
    axB.bar(x + w / 2, flv, w, color="#1b4965", label="amortized flow (NPE)")
    axB.axhline(0.95, ls="--", lw=0.8, color="gray")
    axB.set_xticks(x); axB.set_xticklabels(metrics)
    axB.set_title("(B) Amortized posterior out-calibrates\nand out-sharpens the railed fit")
    axB.set_ylabel("$D^{*}$ metric value")
    axB.legend(loc="upper right", fontsize=6.5)
    for xi, (a, b) in enumerate(zip(nlv, flv)):
        axB.text(xi - w / 2, a + 0.01, f"{a:.2f}", ha="center", fontsize=6)
        axB.text(xi + w / 2, b + 0.01, f"{b:.2f}", ha="center", fontsize=6)
    save(fig, "fig3_resolution")


# ------------------------------------------------- Graphical abstract (50x60 mm)
def graphical_abstract():
    MM = 1 / 25.4
    real = rep["real"]
    sxc = {c["name"]: c for c in sx["cohorts"]}
    rows = [("OSIPI abdomen", real["rail_tol_1e-3"]["rate"], real["rail_tol_1e-3"]["ci"], C_TRUTH),
            ("OSIPI full ROI", sxc["abdomen_full"]["bootstrap_ci"]["point"],
             [sxc["abdomen_full"]["bootstrap_ci"]["lo"], sxc["abdomen_full"]["bootstrap_ci"]["hi"]], C_SX),
            ("Liver 4-b", sxc["lihc_liver_4b_0_50_500_800"]["bootstrap_ci"]["point"],
             [sxc["lihc_liver_4b_0_50_500_800"]["bootstrap_ci"]["lo"],
              sxc["lihc_liver_4b_0_50_500_800"]["bootstrap_ci"]["hi"]], C_SX),
            ("Liver 3-b", sxc["lihc_liver_3b_50_400_800"]["bootstrap_ci"]["point"],
             [sxc["lihc_liver_3b_50_400_800"]["bootstrap_ci"]["lo"],
              sxc["lihc_liver_3b_50_400_800"]["bootstrap_ci"]["hi"]], C_SX)]
    fig, ax = plt.subplots(figsize=(60 * MM, 60 * MM))
    y = list(range(len(rows)))[::-1]
    for yi, (lab, pt, ci, col) in zip(y, rows):
        ax.barh(yi, pt * 100, color=col, height=0.6)
        ax.errorbar(pt * 100, yi, xerr=[[(pt - ci[0]) * 100], [(ci[1] - pt) * 100]],
                    fmt="none", ecolor="black", capsize=1.5, lw=0.8)
        ax.text(2, yi, lab, va="center", ha="left", fontsize=5.2, color="white")
        ax.text(pt * 100 + 1.5, yi, f"{pt*100:.0f}%", va="center", ha="left", fontsize=5.2)
    ax.set_yticks([])
    ax.set_xlim(0, 90)
    ax.set_xlabel("$D^{*}$ railing rate (%)", fontsize=6)
    ax.tick_params(axis="x", labelsize=5)
    ax.set_title("$D^{*}$ rails across cohorts", fontsize=6.5)
    fig.savefig(OUT / "graphical_abstract.pdf", bbox_inches="tight")
    fig.savefig(OUT / "graphical_abstract.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote graphical_abstract (50x60 mm target)")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    graphical_abstract()
    print("all figures ->", OUT)
