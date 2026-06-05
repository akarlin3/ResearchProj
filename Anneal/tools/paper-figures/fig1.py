"""
Figure 1 — annotated example trace of the A=0.5 chimera collapse.

Plots min(R1,R2) for the deterministically regenerated example run (R1, R2 faint
behind it), the graze/recovery thresholds, and annotates the phenomenology the
paper defines: breathing (T_b bracket), the first long graze (= the engineering
"collapse", t_graze) with its >=5 s supra-theta excursion and the sub-0.80
recovery dip where the chimera reforms, and the final crossing with no recovery
(absorption, t_abs). Small ticks mark where the shipped supervisor's 2-s detector
would have fired, showing it over-triggers on grazes that self-heal.

Reads paper_figures/fig1_trace.json (written by fig1_trace.mjs). Writes
paper_figures/fig1.{pdf,png} and fig1_caption.txt. Deterministic.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch  # noqa: F401  (kept for parity)

import style

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = style.ROOT
CFG = style.load_config()
OUT = os.path.join(ROOT, CFG["output_dir"])


def load_trace():
    with open(os.path.join(OUT, "fig1_trace.json")) as f:
        return json.load(f)


def deepest_recovery_between(R, lo, hi, recThresh):
    """Index of the deepest sub-recThresh local minimum in [lo, hi)."""
    best_i, best_v = None, recThresh
    seg = R[lo:hi]
    for k in range(1, len(seg) - 1):
        v = seg[k]
        if v < seg[k - 1] and v <= seg[k + 1] and v < best_v:
            best_v, best_i = v, lo + k
    return best_i, best_v


def main():
    style.apply_style()
    d = load_trace()
    m = d["meta"]
    sdt = m["sampleDt"]
    theta = m["theta"]
    recThr = m["recThresh"]
    R1 = np.array(d["R1"])
    R2 = np.array(d["R2"])
    Rin = np.array(d["R_incoh"])
    t = np.arange(len(Rin)) * sdt

    lab = d["labels"]
    t_graze = lab["t_graze"]
    t_abs = lab["t_abs"]
    gi = lab["grazeIndex"]
    ai = lab["absIndex"]
    Tb = d["breath"]["T_b"]
    peaks = d["breath"]["peaks"]

    # crop to the event window
    p = CFG["fig1"]["plot"]
    x0 = -p["xpad_left_s"]
    x1 = t_abs + p["xpad_right_s"]

    fig, ax = plt.subplots(figsize=(7.6, 3.7))

    # faint per-population order parameters behind the headline min
    if p["show_R1R2_faint"]:
        ax.plot(t, R1, color=style.ROLES["faint"], lw=0.7, alpha=0.55, zorder=1)
        ax.plot(t, R2, color=style.ROLES["faint"], lw=0.7, alpha=0.55, zorder=1,
                label=r"$R_1,\,R_2$")
    ax.plot(t, Rin, color=style.ROLES["primary"], lw=1.5, zorder=3,
            label=r"$\min(R_1,R_2)$")

    # thresholds
    ax.axhline(theta, color=style.ROLES["graze_thresh"], lw=1.0, ls="--", zorder=2)
    ax.axhline(recThr, color=style.ROLES["recover_thresh"], lw=1.0, ls=":", zorder=2)
    ax.text(x0 + 2, theta + 0.006, r"graze threshold $\theta=0.85$",
            color=style.ROLES["graze_thresh"], va="bottom", ha="left", fontsize=7.5)
    ax.text(x0 + 2, recThr - 0.008, r"recovery threshold $=0.80$",
            color=style.ROLES["recover_thresh"], va="top", ha="left", fontsize=7.5)

    # (i) breathing T_b bracket between the first two clean breath peaks
    if Tb and len(peaks) >= 2:
        pa, pb = peaks[0]["t"], peaks[1]["t"]
        ybr = 0.965
        ax.annotate("", xy=(pa, ybr), xytext=(pb, ybr),
                    arrowprops=dict(arrowstyle="<->", color=style.ROLES["black"], lw=1.0))
        ax.text(0.5 * (pa + pb), ybr + 0.012, r"$T_b\approx%.0f$ s" % Tb,
                ha="center", va="bottom", fontsize=8)
        for pk in peaks[:3]:
            ax.plot(pk["t"], np.interp(pk["t"], t, Rin), marker="v", ms=4,
                    color=style.ROLES["black"], zorder=4)
    ax.text(0.30 * t_graze, 0.06, "breathing\nchimera", ha="center", va="bottom",
            fontsize=8, color=style.ROLES["primary"])

    # (ii) graze event: shade the first supra-theta excursion, mark the reform dip
    ri, rv = deepest_recovery_between(list(Rin), gi, ai, recThr)
    # shade graze excursion from the crossing to the deep reform dip
    if ri is not None:
        ax.axvspan(t_graze, ri * sdt, color=style.ROLES["graze_event"], alpha=0.10, zorder=0)
        ax.annotate("graze $\\to$ chimera reforms\n(dips $<0.80$)",
                    xy=(ri * sdt, rv), xytext=(ri * sdt + 14, 0.40),
                    fontsize=8, color=style.ROLES["graze_event"],
                    arrowprops=dict(arrowstyle="->", color=style.ROLES["graze_event"], lw=1.0))

    # (iv) distinct markers for t_graze (engineering collapse) and t_abs
    ax.axvline(t_graze, color=style.ROLES["graze_event"], lw=1.2, ls="-", zorder=2)
    ax.plot([t_graze], [1.005], marker="v", ms=8, color=style.ROLES["graze_event"],
            clip_on=False, zorder=5)
    ax.text(t_graze - 3, 1.018, r"$t_{\rm graze}=%.0f$ s" % t_graze, ha="right",
            va="bottom", fontsize=7.5, color=style.ROLES["graze_event"])

    ax.axvline(t_abs, color=style.ROLES["absorption"], lw=1.4, ls="-", zorder=2)
    ax.plot([t_abs], [1.005], marker="*", ms=11, color=style.ROLES["absorption"],
            clip_on=False, zorder=5)
    ax.text(t_abs + 3, 1.018, r"$t_{\rm abs}=%.0f$ s" % t_abs, ha="left",
            va="bottom", fontsize=7.5, color=style.ROLES["absorption"])

    # (iii) absorption annotation
    ax.annotate("absorption\n(no recovery)", xy=(t_abs, theta),
                xytext=(t_abs - 6, 0.50), ha="right", fontsize=8,
                color=style.ROLES["absorption"],
                arrowprops=dict(arrowstyle="->", color=style.ROLES["absorption"], lw=1.0))
    ax.text(t_abs + 0.55 * p["xpad_right_s"], 0.52, "absorbed\n(sync)", ha="center",
            va="center", fontsize=8, color=style.ROLES["absorption"])

    # (v) supervisor 2-s detector firings (small ticks along the bottom)
    if p["show_supervisor_ticks"]:
        fy = 0.015
        ft = [f["fire_t"] for f in d["supervisor_firings"]]
        ax.plot(ft, [fy] * len(ft), marker="|", ms=9, ls="none",
                color=style.ROLES["supervisor"], mew=1.4, zorder=4,
                label="shipped 2-s detector fires")

    ax.set_xlim(x0, x1)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("time  (s)")
    ax.set_ylabel(r"order parameter")
    ax.set_title(f"A={m['A']}, β={m['beta']}, N={m['N']} (seed {m['seed']}): "
                 f"{lab['n_grazes_before_abs']} recovered grazes, then absorption",
                 fontsize=9, pad=20)
    ax.legend(loc="lower right", ncol=1, fontsize=7.5, framealpha=0.9,
              frameon=True, edgecolor="none")

    paths = style.savefig(fig, "fig1")
    plt.close(fig)

    # caption: two sentences with the run's own numbers
    caption = (
        f"FIG. 1. Engineering collapse versus absorption in a finite-N two-population "
        f"Sakaguchi–Kuramoto chimera (A={m['A']}, β={m['beta']}, N={m['N']}; seed "
        f"{m['seed']}, regenerated deterministically). A graze is a θ=0.85 crossing of "
        f"min(R₁,R₂) sustained ≥{m['W']:.0f} s that the chimera survives—it "
        f"reforms, min(R₁,R₂) dipping back below the recovery threshold 0.80—as at "
        f"the first long graze t_graze={t_graze:.0f} s (the shipped supervisor's 2-s "
        f"detector, ticks, fires here and at every subsequent graze); absorption is the "
        f"first such crossing that is NOT followed by recovery within T_v={m['T_v']:.0f} s, "
        f"reached only at t_abs={t_abs:.0f} s, after which both populations remain locked "
        f"(sync). This run grazes {lab['n_grazes_before_abs']} times—each an over-trigger of "
        f"the engineering detector—before the chimera is finally absorbed, "
        f"{t_abs - t_graze:.0f} s after the first graze."
    )
    with open(os.path.join(OUT, "fig1_caption.txt"), "w") as f:
        f.write(caption + "\n")

    print("wrote:", *[os.path.relpath(p, ROOT) for p in paths])
    print("       " + os.path.relpath(os.path.join(OUT, "fig1_caption.txt"), ROOT))
    print(f"determinism: t_graze={t_graze} t_abs={t_abs} (campaign: 83.2 / 221.9)")


if __name__ == "__main__":
    main()
