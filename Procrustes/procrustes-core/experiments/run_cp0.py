"""CP0 separation driver -- the headline numbers, transcribed into RESULTS.md.

Runs the three-part separation + observable diagnostic for all three departure
families across several seeds and prints across-seed means with 95% CIs.

    python experiments/run_cp0.py            # default 8 seeds
    python experiments/run_cp0.py 4          # fewer seeds (faster)
"""
from __future__ import annotations

import sys

import numpy as np

from procrustes import ProcrustesConfig, FAMILIES, run_separation
from procrustes.config import STRETCHED, LOGNORMAL, TRIEXP
from procrustes.conformal import split_indices
from procrustes.generators import cohort_at, fit_cohort
from procrustes.diagnostic import diagnose
from procrustes.seeding import make_rng, GLOBAL_SEED

ORDER = [STRETCHED, LOGNORMAL, TRIEXP]


def ci(vals):
    a = np.asarray(vals, float)
    m = a.mean()
    se = a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0
    return m, m - 1.96 * se, m + 1.96 * se


def diagnostic_for(fam, cfg, seed):
    """Recompute fits once to score the observable diagnostic for one seed."""
    cohorts = {v: cohort_at(fam.lattice_family, fam.knob, v, cfg) for v in (fam.limit, fam.values[-1] if fam.values[-1] != fam.limit else fam.values[0])}
    worst = fam.values[-1] if fam.values[-1] != fam.limit else fam.values[0]
    truth = cohorts[fam.limit].params[:, cfg.target_index].copy()
    fits = {v: fit_cohort(cohorts[v], cfg.snr) for v in cohorts}
    _, test = split_indices(cfg.n, make_rng(seed))
    return diagnose(fits, truth, fam.limit, worst, test)


def main(n_seeds: int = 8):
    cfg = ProcrustesConfig()
    seeds = [GLOBAL_SEED + 1009 * i for i in range(n_seeds)]
    print(f"# Procrustes CP0 separation  (n={cfg.n}, SNR={cfg.snr}, alpha={cfg.alpha}, "
          f"{n_seeds} seeds, target={cfg.target_param})")
    print(f"# nominal coverage = {cfg.nominal:.2f}; well-ID = bottom-2 D* terciles\n")

    for fam in ORDER:
        res = [run_separation(fam, cfg, seed=s) for s in seeds]
        diag = [diagnostic_for(fam, cfg, s) for s in seeds]
        print("=" * 78)
        print(f"FAMILY {fam.label}   [pre-registered expectation: {fam.expect.upper()}]")
        print("=" * 78)
        for name, vals in [
            ("marginal coverage (departure-blind)", [r.marginal for r in res]),
            ("conditional GAP (placebo - worst)", [r.gap for r in res]),
            ("well-identified-D* GAP", [r.wellid_gap for r in res]),
            ("bias growth ratio |worst|/|limit|", [r.bias_ratio for r in res]),
            ("diagnostic AUC (best residual stat)", [d["auc_best"] for d in diag]),
        ]:
            m, lo, hi = ci(vals)
            print(f"  {name:<38}: {m:7.3f}  [{lo:.3f}, {hi:.3f}]")
        print("  conditional coverage by departure (mean):")
        for v in fam.values:
            cm, clo, chi = ci([r.cond[v] for r in res])
            wm, wlo, whi = ci([r.cond_wellid[v] for r in res])
            print(f"     {fam.knob}={v:>4}: cond={cm:.3f}[{clo:.3f},{chi:.3f}]  "
                  f"well-ID-D*={wm:.3f}[{wlo:.3f},{whi:.3f}]")
        print()

    print("Ref: Gauge naive-transfer misspecification monitor AUC ~ 0.501 (hidden channel).")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
