#!/usr/bin/env python
"""Manuscript consistency gate -- every load-bearing number in levy.tex traces to a seeded
results file; none is typed by hand.

build_numbers()  loads results/RESULTS_CP0.json, RESULTS_CP1.json, RESULTS_CP2.json (all
                 produced by the seeded experiments/run_cp{0,1,2}.py --full) and returns a
                 {macro: (value, why)} dict.
write_numbers()  emits numbers.tex (\\newcommand{\\num...}{...}), AUTO-GENERATED.
verify()         (a) every \\num* used in levy.tex is defined, and (b) load-bearing asserts on
                 the actual numbers (the wall lands in band, the orders are degenerate, etc.).

Exit 0 = gate green; nonzero = a macro is undefined or an assert failed (the script names it).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
TEX = HERE / "levy.tex"
NUMBERS = HERE / "numbers.tex"


def _load(name):
    return json.loads((RESULTS / name).read_text())


def _sci(x, sig=1):
    """Format x as a LaTeX scientific-notation BODY (used inside $...$): '1.7\\times10^{8}'."""
    from math import floor, log10
    if x == 0:
        return "0"
    e = int(floor(log10(abs(x))))
    m = x / (10 ** e)
    return f"{m:.{sig}f}\\times10^{{{e}}}"


def build_numbers():
    cp0 = _load("RESULTS_CP0.json")
    cp1 = _load("RESULTS_CP1.json")
    cp2 = _load("RESULTS_CP2.json")

    h0 = cp0["headline"]
    h1 = cp1["headline"]
    g1 = cp1["grid"]
    br = cp1["two_dt_break"]
    nbp = cp1["nb_persistence"]
    w4 = cp2["wall_vs_alpha_nb4"]
    band = cp0["headline"]["band"]

    nums = {
        # ---- shared regime ----
        "BandLo": (f"{int(band[0])}", "clinical-SNR band lower edge (Polders 2011)"),
        "BandHi": (f"{int(band[1])}", "clinical-SNR band upper edge (Polders 2011)"),
        "CVthr": (f"{h0['cv_threshold']:.2f}", "pre-registered relative-CRLB wall threshold cv_alpha"),
        # ---- CP0 single-order wall ----
        "CPzeroNb": (f"{int(h0['n_b'])}", "CP0 headline number of b-values (sparse clinical)"),
        "CPzeroBmax": (f"{int(h0['b_max'])}", "CP0 headline b_max (s/mm^2)"),
        "CPzeroAlpha": (f"{h0['alpha']:.2f}", "CP0 headline stretched-exp alpha"),
        "CPzeroD": (f"{h0['D']*1e3:.1f}", "CP0 tissue D (x10^-3 mm^2/s)"),
        "WallAnalytic": (f"{h0['wall_snr_analytic']:.1f}", "CP0 analytic CRLB wall SNR*"),
        "WallEmp": (f"{h0['wall_snr_emp']:.1f}", "CP0 empirical bootstrap-MLE wall SNR*"),
        "WallCIlo": (f"{h0['wall_ci'][0]:.1f}", "CP0 empirical wall SNR* 95% CI low"),
        "WallCIhi": (f"{h0['wall_ci'][1]:.1f}", "CP0 empirical wall SNR* 95% CI high"),
        "RhoAD": (f"{cp0['rho_alpha_D_nb4_bmax3000_snr30']:+.2f}",
                  "CP0 alpha-D correlation (n_b=4,b_max=3000,SNR=30)"),
        # ---- CP1 joint (alpha,beta) degeneracy ----
        "RhoABmedian": (f"{g1['rho_median']:.3f}", "median |rho_alpha_beta| over physiological grid"),
        "RhoABmin": (f"{g1['rho_min']:.3f}", "min |rho_alpha_beta| over grid"),
        "RhoABmax": (f"{g1['rho_max']:.3f}", "max |rho_alpha_beta| over grid"),
        "CondMedian": (_sci(g1["cond_median"]), "median joint-FIM condition number (sci body)"),
        "DegenRho": (f"{cp1['verdict']['degen_rho_threshold']:.2f}",
                     "pre-registered degeneracy threshold on median |rho_alpha_beta|"),
        "HeadAlpha": (f"{h1['alpha']:.2f}", "CP1 headline alpha"),
        "HeadBeta": (f"{h1['beta']:.2f}", "CP1 headline beta"),
        "HeadRhoAB": (f"{h1['rho_ab']:+.3f}", "CP1 headline analytic rho_alpha_beta"),
        "HeadCond": (_sci(h1["cond"]), "CP1 headline FIM condition number (sci body)"),
        "HeadCVa": (f"{h1['cv_alpha']:.2f}", "CP1 headline relative CRLB cv_alpha"),
        "HeadCVb": (f"{h1['cv_beta']:.2f}", "CP1 headline relative CRLB cv_beta"),
        "BootAlphaLo": (f"{h1['alpha_ci'][0]:.2f}", "CP1 bootstrap alpha_hat 95% CI low"),
        "BootAlphaHi": (f"{h1['alpha_ci'][1]:.2f}", "CP1 bootstrap alpha_hat 95% CI high"),
        "BootBetaLo": (f"{h1['beta_ci'][0]:.2f}", "CP1 bootstrap beta_hat 95% CI low"),
        "BootBetaHi": (f"{h1['beta_ci'][1]:.2f}", "CP1 bootstrap beta_hat 95% CI high"),
        "Nboot": (f"{int(h1['n_boot'])}", "CP1 bootstrap replicates"),
        "NbRhoFour": (f"{nbp['rho'][0]:.3f}", "|rho_alpha_beta| at n_b=4 (single dt)"),
        "NbRhoSixteen": (f"{nbp['rho'][-1]:.3f}", "|rho_alpha_beta| at n_b=16 (single dt)"),
        "BreakRhoSingle": (f"{br['rho_single']:.3f}", "|rho| 16 b one diffusion time"),
        "BreakRhoTwo": (f"{br['rho_two']:.3f}", "|rho| 8+8 b two diffusion times"),
        "BreakCVaSingle": (f"{br['cv_alpha_single']:.2f}", "cv_alpha single diffusion time"),
        "BreakCVaTwo": (f"{br['cv_alpha_two']:.2f}", "cv_alpha two diffusion times"),
        # ---- CP2 across-alpha robustness ----
        "AlphaPhysLo": (f"{w4['alpha'][0]:.2f}", "physiological alpha range low"),
        "AlphaPhysHi": (f"{w4['alpha'][-1]:.2f}", "physiological alpha range high"),
        "WallAlphaLo": (f"{w4['wall_min']:.1f}", "min wall SNR* across alpha (n_b=4)"),
        "WallAlphaHi": (f"{w4['wall_max']:.1f}", "max wall SNR* across alpha (n_b=4)"),
    }
    return nums, cp0, cp1, cp2


def write_numbers(nums):
    lines = [
        "% AUTO-GENERATED by paper/consistency.py from seeded results -- DO NOT EDIT BY HAND.",
        "% Every \\num* macro traces to results/RESULTS_CP{0,1,2}.json, produced by the seeded",
        "% experiments/run_cp{0,1,2}.py --full. Regenerate: python paper/consistency.py.",
        "",
    ]
    for k, (v, why) in nums.items():
        lines.append(f"\\newcommand{{\\num{k}}}{{{v}}}  % {why}")
    NUMBERS.write_text("\n".join(lines) + "\n")


def verify(nums, cp0, cp1, cp2):
    ok = True
    # (a) traceability: every \num* used in the manuscript is defined
    if TEX.exists():
        used = set(re.findall(r"\\num([A-Za-z]+)\b", TEX.read_text()))
        undefined = used - set(nums)
        print(f"  macros: used={len(used)} defined={len(nums)} undefined={sorted(undefined)}")
        if undefined:
            ok = False
    else:
        print("  (levy.tex not present yet -- skipping traceability check)")

    # (b) load-bearing asserts on the actual seeded numbers
    h0 = cp0["headline"]
    lo, hi = h0["band"]
    asserts = [
        ("CP0 analytic wall in band", lo <= h0["wall_snr_analytic"] <= hi),
        ("CP0 empirical wall in band", lo <= h0["wall_snr_emp"] <= hi),
        ("CP0 empirical wall CI ordered", h0["wall_ci"][0] <= h0["wall_ci"][1]),
        ("CP0 wall stands (not refuted)", cp0["verdict"]["wall_stands"] and not cp0["verdict"]["refuted"]),
        ("CP1 degenerate verdict matches threshold",
         cp1["verdict"]["degenerate"] == (cp1["grid"]["rho_median"] >= cp1["verdict"]["degen_rho_threshold"])),
        ("CP1 median |rho_ab| >= 0.9 (degenerate)", cp1["grid"]["rho_median"] >= 0.90),
        ("CP1 degeneracy persists with n_b (n_b=16 still >= 0.9)", cp1["nb_persistence"]["rho"][-1] >= 0.90),
        ("CP1 two-dt breaks degeneracy (rho drops)",
         cp1["two_dt_break"]["rho_two"] < cp1["two_dt_break"]["rho_single"]),
        ("CP2 wall robust across alpha (all in band)",
         cp2["verdict"]["wall_robust_across_alpha"] and not cp2["verdict"]["refuted_across_alpha"]),
        ("CP2 wall range across alpha inside band",
         lo <= cp2["wall_vs_alpha_nb4"]["wall_min"] and cp2["wall_vs_alpha_nb4"]["wall_max"] <= hi),
    ]
    for name, cond in asserts:
        print(f"  assert {name}: {'PASS' if cond else 'FAIL'}")
        ok = ok and bool(cond)
    return ok


def main():
    nums, cp0, cp1, cp2 = build_numbers()
    write_numbers(nums)
    print(f"wrote {NUMBERS.name} ({len(nums)} macros)")
    ok = verify(nums, cp0, cp1, cp2)
    print("consistency gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
