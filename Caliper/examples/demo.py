"""Caliper end-to-end demo (one command, fixed seeds).

    synthetic IVIM  ->  MAF posterior  ->  split-conformal  ->  calibration ruler

Run:  python examples/demo.py

Requires the estimator extra (torch):  pip install -e ".[estimator]"

Every number printed here is produced by this run -- the README quotes this
script's output verbatim. The scenario is a realistic *deployment shift*: the
MAF is trained/calibrated for fitting at high SNR but evaluated at a lower SNR,
which is where model-based IVIM UQ is known to be over-confident.
"""
from __future__ import annotations

import numpy as np

from caliper import metrics as M
from caliper.conformal import SplitConformalQuantile
from caliper.estimator_maf import MAFPosterior
from caliper.forward import PARAM_NAMES, synthetic_cohort

# --- configuration (fixed seeds for full reproducibility) ------------------- #
TRAIN_SNR, DEPLOY_SNR = 60.0, 25.0
LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
ALPHA = 0.10  # 90% central intervals


def main() -> None:
    print("Caliper demo -- synthetic IVIM -> MAF -> conformal -> ruler\n")

    # 1. synthetic, PHI-free cohorts (in-repo generator, seeded)
    train = synthetic_cohort(n=8000, snr=TRAIN_SNR, seed=0)  # fit the flow
    cal = synthetic_cohort(n=2000, snr=DEPLOY_SNR, seed=1)   # conformal split
    test = synthetic_cohort(n=4000, snr=DEPLOY_SNR, seed=2)  # evaluation
    print(f"cohorts: train(SNR{TRAIN_SNR:.0f}, n={len(train)}), "
          f"cal(SNR{DEPLOY_SNR:.0f}, n={len(cal)}), "
          f"test(SNR{DEPLOY_SNR:.0f}, n={len(test)}); "
          f"b-values={train.bvalues.size}")

    # 2. train the MAF posterior
    est = MAFPosterior(n_bvalues=train.bvalues.size, epochs=60, seed=0)
    est.fit(train.signals, train.params)
    print(f"MAF trained: final NLL/dim = {est.history_[-1]:.4f}\n")

    # 3. predict posterior quantiles
    q_cal = est.predict_quantiles(cal.signals, LEVELS)
    q_test = est.predict_quantiles(test.signals, LEVELS)

    # 4a. score the RAW MAF (conditioning conditional coverage on true params)
    raw = M.score_quantiles(test.params, q_test, LEVELS, alpha=ALPHA,
                            param_names=PARAM_NAMES, conditioning=test.params)
    print(M.format_scorecard(raw, title="RAW MAF (held-out synthetic)"))
    print()

    # 4b. conformalize on the calibration split, then re-score
    cq = SplitConformalQuantile(LEVELS).calibrate(q_cal, cal.params)
    q_corr = cq.apply(q_test)
    cor = M.score_quantiles(test.params, q_corr, LEVELS, alpha=ALPHA,
                            param_names=PARAM_NAMES, conditioning=test.params)
    print(M.format_scorecard(cor, title="MAF + SPLIT-CONFORMAL"))

    # 5. headline summary
    print("\n=== headline: marginal coverage (nominal 0.900) ===")
    print(f"{'param':>6} {'raw cover':>10} {'conf cover':>11} "
          f"{'raw|gap|':>9} {'conf|gap|':>10}")
    for r, c in zip(raw, cor):
        print(f"{r.name:>6} {r.coverage:>10.3f} {c.coverage:>11.3f} "
              f"{abs(r.coverage_gap):>9.3f} {abs(c.coverage_gap):>10.3f}")

    print("\nHonest caveat -- conditional coverage on the high-D* tercile:")
    print(f"  raw   D* g2 (high) = {raw[2].conditional[2]:.3f}")
    print(f"  conf  D* g2 (high) = {cor[2].conditional[2]:.3f}  "
          "(still < nominal: identifiability limit, not fixed by conformal)")


if __name__ == "__main__":
    main()
