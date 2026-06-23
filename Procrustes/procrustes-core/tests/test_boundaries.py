"""Pre-registered BOUNDARY gates -- the wedge is scoped, not universal.

These pin the honest limits discovered while hardening the probe, so a future
over-claim (e.g. "all misspecification breaks D") trips a gate.
"""
import numpy as np
import pytest

from procrustes import ProcrustesConfig, TRIEXP, STRETCHED, run_separation
from procrustes.generators import make_cohort

CFG = ProcrustesConfig(n=500)


def _cohort(family, extra):
    return make_cohort(family=family, n=CFG.n, snr=CFG.snr, seed=CFG.seed,
                       prior=CFG.prior, noise=CFG.noise, extra=extra)


def test_triexp_is_null():
    # BOUNDARY: a *faster* tri-exp third pool decays away from the high-b tissue
    # slope, so it does NOT bias D -- conditional coverage stays flat. The wedge
    # is mechanism-specific (high-b aliasing), not generic misspecification.
    res = run_separation(TRIEXP, CFG, seed=CFG.seed)
    assert res.gap < CFG.null_gap_max, f"tri-exp gap {res.gap:.3f} should be ~null"
    assert res.bias_ratio < 2.0


def test_bi_exp_limit_is_exact_continuity():
    # The placebo is real: at the bi-exp limit (beta=1) the stretched generator
    # reduces EXACTLY to the bi-exponential cohort (same base params, same seed).
    stretched_lim = _cohort(STRETCHED.lattice_family, {STRETCHED.knob: STRETCHED.limit})
    biexp = _cohort("biexp", None)
    max_resid = float(np.max(np.abs(stretched_lim.signals_clean - biexp.signals_clean)))
    assert max_resid < 1e-9, f"continuity residual {max_resid:.2e} (should be ~0)"
