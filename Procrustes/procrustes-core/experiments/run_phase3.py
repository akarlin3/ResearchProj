"""PHASE 3 / GATE D -- robustness envelope (honest boundaries).

Stress the heavy-tail separation across the realistic acquisition/analysis
envelope and report, with bootstrap CIs, where it holds and where it breaks:

  * SNR sweep        : Rician SNR in {25, 35, 50, 75, 100},
  * sample-size sweep: n voxels in {300, 500, 1000, 2000},
  * b-scheme sweep   : Lattice default (22-pt) vs a sparse clinical 10-pt scheme
                       vs a low-b-poor 8-pt scheme,
  * more seeds       : a 16-seed confirmation of the default-config headline.

For each condition we report the overall conditional gap and the well-identified-
D* gap (Gauge's identifiable region) with two-level cluster-bootstrap 95% CIs.
A condition is marked SURVIVES iff the well-ID gap's CI lower bound > 0; any
condition that fails is reported as an honest boundary, not buried.

Usage:  python experiments/run_phase3.py [n_seeds]   (default 8; default-config
        headline is always confirmed at 16 seeds)
"""
from __future__ import annotations

import dataclasses
import json
import os
import sys

import numpy as np

from procrustes import (ProcrustesConfig, STRETCHED, separation_detail,
                        cluster_bootstrap_gap, GLOBAL_SEED)

_RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "results")

# Alternative acquisition schemes (s/mm^2).
B_DEFAULT = None  # Lattice 22-point default
B_CLINICAL_SPARSE = (0, 10, 20, 40, 60, 80, 150, 300, 500, 800)   # 10-pt clinical
B_LOWB_POOR = (0, 20, 50, 100, 200, 400, 600, 800)                # 8-pt, sparse low-b


def envelope(cfg, seeds):
    """Bootstrap gap (ALL) and well-ID gap for STRETCHED under one config."""
    det = [separation_detail(STRETCHED, cfg, seed=s) for s in seeds]
    return {"gap_all": cluster_bootstrap_gap(det, "all"),
            "gap_wellid": cluster_bootstrap_gap(det, "wellid"),
            "gap_lo": cluster_bootstrap_gap(det, "lo")}


def _survives(g):
    return bool(g["gap_wellid"]["lo"] > 0.0)


def _row(label, g):
    ga, gw = g["gap_all"], g["gap_wellid"]
    flag = "SURVIVES" if _survives(g) else "** FAILS **"
    return (f"  {label:<22}: gap(all) {ga['point']:+.3f} [{ga['lo']:+.3f},{ga['hi']:+.3f}]"
            f"   well-ID {gw['point']:+.3f} [{gw['lo']:+.3f},{gw['hi']:+.3f}]   {flag}")


def main(n_seeds: int = 8):
    base = ProcrustesConfig()
    seeds = [GLOBAL_SEED + 1009 * i for i in range(n_seeds)]
    seeds16 = [GLOBAL_SEED + 1009 * i for i in range(16)]
    print(f"# PHASE 3 / GATE D -- robustness envelope (stretched-exp, heavy-tail)")
    print(f"# base n={base.n}, SNR={base.snr}; {n_seeds} seeds/condition; "
          f"well-ID = Gauge identifiable region (bottom-2 D* terciles)\n")
    payload = {"config": {"n_seeds": n_seeds, "base_n": base.n, "base_snr": base.snr},
               "snr": {}, "n": {}, "bscheme": {}, "headline16": {}}

    # ---- default-config headline at 16 seeds ---------------------------------
    print("=" * 96)
    print("DEFAULT CONFIG headline (16 seeds): the load-bearing gaps with bootstrap CIs")
    print("=" * 96)
    g16 = envelope(base, seeds16)
    print(_row("SNR50 n1000 (default)", g16))
    payload["headline16"] = g16
    print()

    # ---- SNR sweep -----------------------------------------------------------
    print("=" * 96)
    print("SNR SWEEP (n=1000)")
    print("=" * 96)
    for snr in (25, 35, 50, 75, 100):
        g = envelope(dataclasses.replace(base, snr=float(snr)), seeds)
        print(_row(f"SNR={snr}", g))
        payload["snr"][str(snr)] = g
    print()

    # ---- sample-size sweep ---------------------------------------------------
    print("=" * 96)
    print("SAMPLE-SIZE SWEEP (SNR=50)")
    print("=" * 96)
    for n in (300, 500, 1000, 2000):
        g = envelope(dataclasses.replace(base, n=int(n)), seeds)
        print(_row(f"n={n}", g))
        payload["n"][str(n)] = g
    print()

    # ---- b-scheme sweep ------------------------------------------------------
    print("=" * 96)
    print("B-SCHEME SWEEP (SNR=50, n=1000)")
    print("=" * 96)
    for name, bv in (("default-22pt", B_DEFAULT),
                     ("clinical-sparse-10pt", B_CLINICAL_SPARSE),
                     ("lowb-poor-8pt", B_LOWB_POOR)):
        g = envelope(dataclasses.replace(base, bvalues=bv), seeds)
        print(_row(name, g))
        payload["bscheme"][name] = g
    print()

    # ---- envelope verdict ----------------------------------------------------
    def collect(d): return [(k, v) for k, v in d.items()]
    all_conditions = (collect(payload["snr"]) + collect(payload["n"])
                      + collect(payload["bscheme"]))
    failed = [k for k, v in all_conditions if not _survives(v)]
    survived = [k for k, v in all_conditions if _survives(v)]
    print("#" * 96)
    print("GATE D VERDICT -- robustness envelope:")
    print(f"  SURVIVES (well-ID gap CI-lo > 0) in {len(survived)}/{len(all_conditions)} conditions.")
    if failed:
        print(f"  HONEST BOUNDARY -- separation does NOT survive in: {failed}")
    else:
        print(f"  No boundary found within the swept envelope; separation robust throughout.")
    print(f"  Default headline (16 seeds): conditional gap {g16['gap_all']['point']:+.3f} "
          f"[{g16['gap_all']['lo']:+.3f},{g16['gap_all']['hi']:+.3f}]; "
          f"well-ID gap {g16['gap_wellid']['point']:+.3f} "
          f"[{g16['gap_wellid']['lo']:+.3f},{g16['gap_wellid']['hi']:+.3f}].")
    print("#" * 96)
    payload["gate_d"] = {"survived": survived, "failed": failed,
                         "n_conditions": len(all_conditions)}

    os.makedirs(_RESULTS, exist_ok=True)
    with open(os.path.join(_RESULTS, "phase3_gateD.json"), "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"\n[wrote results/phase3_gateD.json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 8))
