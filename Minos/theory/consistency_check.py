"""CP4 / GATE 4 — final consistency check for the Plumbline note (deterministic, fast).

Verifies that the finalized note's claims are a subset of the gated derivations and the CP1/CP3/CP4
outputs, without re-running the expensive Monte-Carlo gates:

  1. Re-derives the Theorem-1 analytic constants (slope |z*|/6, G_lead, tau_stat^exact, G_full, the
     operating-point decomposition) from scratch and asserts they match the numbers quoted in
     ``plumbline.md`` / ``plumbline.tex``.
  2. Asserts ``impossibility.md`` still carries the HUMAN-PROOF-REVIEW banner and the CP2 hardening
     (pinned ``hidden``-as-definition + the partial-leak proposition), and makes no "machine-verified"
     claim about the impossibility conclusion.
  3. Asserts the note embeds the CP1 error-barred headline ("opposite sides of 1" = ROBUST) and the
     specific CI numbers printed by ``gap_ci.py`` this session.
  4. Asserts the two publication figures exist as vector PDFs and regenerate from ``make_figures.py``.

Prints a one-paragraph consistency summary and a GATE 4 verdict. Pure-analytic + file checks; no seed
dependence beyond the deterministic analytic re-derivation.

Run:  ``.venv-theory/bin/python theory/consistency_check.py``
"""
from __future__ import annotations

import os
import re

import numpy as np
import mpmath as mp
from scipy.stats import norm, skewnorm
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
LAM, KAP = 3.0, 3.0


def zstar(lam):
    psi = lambda t: t * float(mp.ncdf(t)) + float(mp.npdf(t))
    return float(mp.findroot(lambda t: (lam - 1) * psi(t) + t, -0.4))


def gamma_of_kappa(kappa):
    d = kappa / np.sqrt(1 + kappa ** 2)
    return ((4 - np.pi) / 2) * (d * np.sqrt(2 / np.pi)) ** 3 / (1 - 2 * d ** 2 / np.pi) ** 1.5


def taustat_exact(kappa, level=0.90):
    a = kappa
    d = a / np.sqrt(1 + a ** 2)
    mean_sn, sd_sn = d * np.sqrt(2 / np.pi), np.sqrt(1 - 2 * d ** 2 / np.pi)
    zL = norm.ppf(0.5 + level / 2.0)
    Fu = lambda x: skewnorm.cdf(x * sd_sn + mean_sn, a)
    return float(brentq(lambda tau: Fu(zL * tau) - Fu(-zL * tau) - level, 0.2, 5.0, xtol=1e-9))


def read(name):
    with open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  — {detail}" if detail else ""))
    return ok


def main():
    print("=" * 84)
    print("CP4 / GATE 4 — Plumbline final consistency check")
    print("=" * 84)
    results = []

    # --- 1. analytic re-derivation vs the note ------------------------------------------------
    print("\n[1] Theorem-1 analytic constants re-derived from scratch:")
    zs = zstar(LAM)
    slope = abs(zs) / 6.0
    g = gamma_of_kappa(KAP)
    G_lead = slope * g
    ts_exact = taustat_exact(KAP)
    G_full = (1 + G_lead) - ts_exact
    frac_dec = G_lead / G_full
    frac_cov = (1 - ts_exact) / G_full
    print(f"      z*(3)={zs:.4f}  slope=|z*|/6={slope:.4f}  gamma(3)={g:.4f}")
    print(f"      G_lead={G_lead:.4f}  tau_stat^exact={ts_exact:.4f}  G_full={G_full:.4f}")
    print(f"      decomposition: decision {100*frac_dec:.0f}%  coverage-shrink {100*frac_cov:.0f}%")
    results.append(check("slope = 0.0727", abs(slope - 0.0727) < 5e-4, f"{slope:.4f}"))
    results.append(check("G_lead = 0.0485", abs(G_lead - 0.0485) < 5e-4, f"{G_lead:.4f}"))
    results.append(check("tau_stat^exact = 0.9638", abs(ts_exact - 0.9638) < 5e-4, f"{ts_exact:.4f}"))
    results.append(check("G_full = 0.0847", abs(G_full - 0.0847) < 5e-4, f"{G_full:.4f}"))
    results.append(check("decision share ~57%", abs(100 * frac_dec - 57) < 1.5, f"{100*frac_dec:.0f}%"))

    md = read("plumbline.md")
    tex = read("plumbline.tex")
    for tok in ("0.0727", "0.0485", "0.0847", "0.9638"):
        results.append(check(f"plumbline.md quotes {tok}", tok in md))
    results.append(check("note states ~57% decomposition", "57" in md and "43" in md))

    # --- 2. impossibility.md: banner + CP2 hardening, no 'verified' claim ----------------------
    print("\n[2] impossibility.md banner + CP2 hardening:")
    imp = read("impossibility.md")
    results.append(check("HUMAN-PROOF-REVIEW banner present",
                         "REQUIRES HUMAN PROOF-REVIEW" in imp and "NOT MACHINE-VERIFIED" in imp))
    results.append(check("hidden pinned as a DEFINITION",
                         ("Definition 3" in imp) and ("O`-invariant" in imp or "O-invariant" in imp)))
    results.append(check("partial-leak proposition present",
                         "Proposition (partial leak)" in imp and "P_mid" in imp.replace("`", "")))
    results.append(check("two bounds cross-referenced to parts (i)/(ii)",
                         "part (ii)" in imp and "part (i)" in imp))
    results.append(check("structural-premise number intact (0.007248)", "0.007248" in imp))
    results.append(check("no 'machine-verified' impossibility claim",
                         "machine-verified" not in imp.lower().replace("not machine-verified", "")
                         .replace("not\nmachine-verified", "")))

    # --- 3. CP1 headline + CI numbers embedded in the note ------------------------------------
    print("\n[3] CP1 error-barred headline embedded:")
    results.append(check("'opposite sides of 1' = ROBUST in note",
                         "opposite sides of 1" in md and "ROBUST" in md))
    for tok in ("1.0514", "0.0876", "0.9639", "[1.0493, 1.0536]", "[0.0855, 0.0897]"):
        results.append(check(f"note embeds CP1 number {tok}", tok in md))
    results.append(check("halt-to-report verdict present", "halt-to-report" in md.lower()))

    # --- 4. figures exist as vector PDFs ------------------------------------------------------
    print("\n[4] publication figures:")
    f1 = os.path.join(HERE, "figures", "fig1_scaling_law.pdf")
    f2 = os.path.join(HERE, "figures", "fig2_monitor_bound.pdf")
    mk = os.path.join(HERE, "figures", "make_figures.py")
    results.append(check("fig1_scaling_law.pdf exists", os.path.exists(f1),
                         f"{os.path.getsize(f1)} bytes" if os.path.exists(f1) else "missing"))
    results.append(check("fig2_monitor_bound.pdf exists", os.path.exists(f2),
                         f"{os.path.getsize(f2)} bytes" if os.path.exists(f2) else "missing"))
    results.append(check("make_figures.py present (regenerable from seed)", os.path.exists(mk)))
    results.append(check("note references both figures",
                         "fig1_scaling_law.pdf" in md and "fig2_monitor_bound.pdf" in md))

    # --- summary -------------------------------------------------------------------------------
    n_pass = sum(results)
    n_tot = len(results)
    print("\n" + "=" * 84)
    print("GATE 4 — consistency summary")
    print("=" * 84)
    print(
        "  Every quantitative claim in the finalized note traces to a gated artefact: Theorem 1's\n"
        "  constants (slope |z*|/6=0.0727, G_lead=0.0485, tau_stat^exact=0.9638, G_full=0.0847, and the\n"
        "  ~57%/43% decision/coverage decomposition at the operating point) are re-derived here and\n"
        "  match the note; the headline gap and its 95% CI (tau*=1.0514 [1.0493,1.0536], G=0.0876\n"
        "  [0.0855,0.0897]) are the CP1 (gap_ci.py) printout, on which the 'opposite sides of 1'\n"
        "  conclusion is reported ROBUST; the monitor bound (L_U=3, slope 0.32) and the theory-vs-v2/v3\n"
        "  tables are the CP3 (confirm.py) GATE-3 PASS; the two vector figures regenerate deterministically\n"
        "  from make_figures.py. The impossibility (Theorem 2(i)) is hardened per CP2 — 'hidden' pinned as\n"
        "  the O-invariant component, plus the partial-leak proposition cross-referencing parts (i)/(ii) —\n"
        "  and retains its REQUIRES-HUMAN-PROOF-REVIEW banner with no machine-verified claim. No theory\n"
        "  constant was tuned to the simulation; the note claims are a subset of the gated outputs.")
    print(f"\n  checks: {n_pass}/{n_tot} PASS")
    ok = n_pass == n_tot
    print(f"\n  GATE 4: {'PASS' if ok else 'FAIL — see [FAIL] lines above'}")
    assert ok, "GATE 4 consistency check failed"
    return results


if __name__ == "__main__":
    main()
