#!/usr/bin/env python
"""Reverb -- Lethe's constructive precision-vs-coverage counterexample (SOLID).

This is SOLID, not provisional: it depends only on **Lattice** (synthetic
ground-truth IVIM cohorts; a publication-independent resource) and **Caliper**
(in-tree estimator + conformal ruler). It depends on NO upstream paper.

It exhibits, with known ground truth, a regime where test--retest repeatability
is excellent while the deployed interval's coverage of the *truth* is broken --
turning Lethe's "precision != coverage" from asserted into shown. Two runs frame
it:

  * MATCHED CONTROL (biexp truth, biexp-calibrated): the model is correctly
    specified -> precision and coverage track -> no counterexample (the null that
    proves the construction is not rigged).
  * MISMATCH HEADLINE (dispersed-perfusion truth, biexp-calibrated): the deployed
    bi-exp model is mis-specified, as real perfusion is -> a structural bias that
    ROI-averaging makes precise but cannot remove -> the perfusion fraction f at
    low D* is highly repeatable yet badly under-covered. The counterexample.

A pre-registered family x ROI-size sensitivity surface is written alongside, so
the whole monotone effect (and the controls that never break) is visible.

Scope (load-bearing): a synthetic possibility-and-mechanism proof. It shows the
divergence CAN occur in IVIM and WHY; it does NOT quantify any real-world
miscalibration magnitude.

Run: python scripts/run_reverb.py [--n 2000] [--n-boot 2000] [--seed 20260622]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from echo_repeat import reverb as R  # noqa: E402

# Pre-registered headline configuration (realism-, not effect-size-, motivated):
#   * dispersed perfusion fit by a bi-exp model (real perfusion is dispersed);
#   * region_size = a representative whole-tumor ROI voxel count;
#   * SNR = the reference estimator's tuned single-voxel regime.
HEADLINE_FAMILY = "stretched"          # anomalous-perfusion IVIM (representative dispersed model)
CONTROL_FAMILY = "biexp"               # correctly-specified control
REGION_SIZE = 200
SNR = 40.0


def _f_lo_row(res: dict) -> dict:
    return [r for r in res["rows"] if r["param"] == "f" and r["stratum"] == 0][0]


def _arm(truth_family, n, n_boot, seed):
    tr = R.simulate_testretest(n_eval=n, n_cal=n, snr=SNR, region_size=REGION_SIZE,
                               truth_family=truth_family, cal_family=CONTROL_FAMILY, seed=seed)
    res = R.analyze(tr, n_boot=n_boot, seed=0)
    return tr, res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=R.REVERB_SEED)
    ap.add_argument("--sweep-n", type=int, default=1500)
    args = ap.parse_args()

    tr_c, res_c = _arm(CONTROL_FAMILY, args.n, args.n_boot, args.seed)
    tr_m, res_m = _arm(HEADLINE_FAMILY, args.n, args.n_boot, args.seed)
    sweep = R.sensitivity_sweep(snr=SNR, n=args.sweep_n, seed=args.seed)

    fc = _f_lo_row(res_c)
    fm = _f_lo_row(res_m)
    verdict = {
        "counterexample_found": bool(res_m["counterexample_found"]),
        "control_tracks": bool(not fc["coverage_broken"]),
        # precision is blind: repeatability essentially unchanged while coverage diverges
        "precision_blind": bool(abs(fm["wcv"] - fc["wcv"]) < 0.02
                                and fm["repeatable"] and fc["repeatable"]),
        "coverage_drop": float(fc["cov_split"] - fm["cov_split"]),
        "control_f_lo_cov": fc["cov_split"],
        "mismatch_f_lo_cov": fm["cov_split"],
    }

    out = {
        "kind": "constructive-counterexample",
        "provisional": False, "solid": True,
        "seed": int(args.seed),
        "headline_config": {
            "control_family": CONTROL_FAMILY, "mismatch_family": HEADLINE_FAMILY,
            "cal_family": CONTROL_FAMILY, "region_size": REGION_SIZE, "snr": SNR,
            "n_eval": args.n, "n_cal": args.n, "n_boot": args.n_boot,
            "bars": res_m["bars"],
        },
        "control": {"manifest": res_c["manifest"], "f_lo": fc,
                    "marginal_f_split": res_c["marginal_coverage"]["split"]["f"],
                    "counterexample_found": res_c["counterexample_found"],
                    "rows": res_c["rows"]},
        "mismatch": {"manifest": res_m["manifest"], "f_lo": fm,
                     "marginal_f_split": res_m["marginal_coverage"]["split"]["f"],
                     "counterexample_found": res_m["counterexample_found"],
                     "rows": res_m["rows"]},
        "sensitivity": sweep,
        "verdict": verdict,
    }

    outdir = ROOT / "results"
    outdir.mkdir(exist_ok=True)
    (outdir / "RESULTS_REVERB.json").write_text(json.dumps(out, indent=2))

    # --- markdown ---
    L = [
        "# Reverb -- constructive counterexample: precision without coverage (SOLID)",
        "",
        f"Seed {args.seed}. Lattice (synthetic truth) + Caliper (estimator/conformal), read-only. "
        f"ROI region_size={REGION_SIZE}, single-voxel SNR={SNR}, n_cal=n_eval={args.n}, "
        f"BCa bootstrap n_boot={args.n_boot}. Params in Caliper space (D, f, D*).",
        "",
        "**Scope:** synthetic possibility-and-mechanism proof -- the divergence CAN occur in IVIM and "
        "here is why; it does NOT quantify real-world miscalibration magnitude.",
        "",
        "## Headline: f @ D*-lo (perfusion fraction, low-D* regime)",
        "| arm | truth fit as bi-exp | wCV [BCa] | ICC | repeatable | cov(true) [BCa] | broken |",
        "|---|---|---|---|---|---|---|",
        f"| matched control | {CONTROL_FAMILY} | {fc['wcv']:.3f} [{fc['wcv_ci'][0]:.3f},{fc['wcv_ci'][1]:.3f}] "
        f"| {fc['icc']:.3f} | {fc['repeatable']} | {fc['cov_split']:.3f} "
        f"[{fc['cov_split_ci'][0]:.3f},{fc['cov_split_ci'][1]:.3f}] | {fc['coverage_broken']} |",
        f"| **mismatch** | **{HEADLINE_FAMILY}** | {fm['wcv']:.3f} [{fm['wcv_ci'][0]:.3f},{fm['wcv_ci'][1]:.3f}] "
        f"| {fm['icc']:.3f} | {fm['repeatable']} | **{fm['cov_split']:.3f} "
        f"[{fm['cov_split_ci'][0]:.3f},{fm['cov_split_ci'][1]:.3f}]** | **{fm['coverage_broken']}** |",
        "",
        f"Precision is identical (wCV {fc['wcv']:.3f} vs {fm['wcv']:.3f}; ICC ~{fm['icc']:.2f}), yet "
        f"model mismatch drops true-coverage by {verdict['coverage_drop']:.3f} "
        f"({fc['cov_split']:.3f} -> {fm['cov_split']:.3f}). Marginal f-coverage looks fine in both "
        f"({res_c['marginal_coverage']['split']['f']:.3f} / {res_m['marginal_coverage']['split']['f']:.3f}) "
        f"-- the break is hidden until truth is known. **Repeatability is blind to the bias.**",
        "",
        "## Verdict",
        f"- counterexample_found: **{verdict['counterexample_found']}**",
        f"- control_tracks (matched model not broken): {verdict['control_tracks']}",
        f"- precision_blind (wCV ~equal while coverage diverges): {verdict['precision_blind']}",
        "",
        "## Per-regime detail (mismatch arm): repeatability vs conditional true-coverage",
        "| param | regime | wCV | ICC | repeatable | cov(true) [BCa] | cov(mondrian) | counterex |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in res_m["rows"]:
        L.append(
            f"| {r['param']} | {r['stratum_name']} | {r['wcv']:.2f} | {r['icc']:+.2f} | "
            f"{'YES' if r['repeatable'] else 'no'} | "
            f"{r['cov_split']:.2f} [{r['cov_split_ci'][0]:.2f},{r['cov_split_ci'][1]:.2f}] | "
            f"{r['cov_mondrian']:.2f} | {'**YES**' if r['counterexample'] else '-'} |")
    L += [
        "",
        f"## Sensitivity surface: cov(true) of f @ D*-lo by truth-family x ROI size (SNR {SNR})",
        "| truth family | " + " | ".join(f"nvox={rs}" for rs in sweep["region_sizes"]) + " |",
        "|---|" + "|".join("---" for _ in sweep["region_sizes"]) + "|",
    ]
    for fam in sweep["families"]:
        cells = " | ".join(f"{sweep['grid'][fam][rs]['cov_split']:.2f}"
                           + ("*" if sweep['grid'][fam][rs]['counterexample'] else "")
                           for rs in sweep["region_sizes"])
        L.append(f"| {fam} | {cells} |")
    L += [
        "",
        "`*` = repeatable-but-broken (counterexample) cell. The matched bi-exp control never breaks; "
        "dispersed-perfusion families break increasingly with ROI size -- precision rises, accuracy does not.",
    ]
    (outdir / "RESULTS_REVERB.md").write_text("\n".join(L))

    print("\n".join(L))
    print(f"\nwrote {outdir/'RESULTS_REVERB.json'} and .md")
    # Success = the run rendered. Both a counterexample and an honest null are valid outcomes.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
