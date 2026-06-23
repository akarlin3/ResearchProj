"""GATE -- the observable misspecification diagnostic (the risk-bearing half).

Pins two pre-registered facts: the residual-structure diagnostic DETECTS the
heavy-tail channel that breaks D-coverage (beating Gauge's AUC~0.5 hidden-channel
monitor), while pure dispersion stays LESS detectable -- the honest scoping.
"""
import pytest

from procrustes import ProcrustesConfig, STRETCHED, LOGNORMAL
from procrustes.conformal import split_indices
from procrustes.generators import cohort_at, fit_cohort
from procrustes.diagnostic import diagnose
from procrustes.seeding import make_rng

CFG = ProcrustesConfig(n=500)


def _diag(fam):
    worst = fam.values[-1] if fam.values[-1] != fam.limit else fam.values[0]
    cohorts = {v: cohort_at(fam.lattice_family, fam.knob, v, CFG) for v in (fam.limit, worst)}
    truth = cohorts[fam.limit].params[:, CFG.target_index].copy()
    fits = {v: fit_cohort(cohorts[v], CFG.snr) for v in cohorts}
    _, test = split_indices(CFG.n, make_rng(CFG.seed))
    return diagnose(fits, truth, fam.limit, worst, test)


def test_stretched_is_detectable():
    # Beats chance / Gauge's hidden-channel monitor (AUC ~ 0.501).
    d = _diag(STRETCHED)
    assert d["auc_best"] > CFG.diagnostic_auc_min - 0.02, f"AUC {d['auc_best']:.3f}"


def test_dispersion_is_a_harder_channel():
    # BOUNDARY (honest scoping): pure dispersion is less detectable than the
    # heavy-tail channel -- consistent with Gauge's hidden-channel finding.
    assert _diag(LOGNORMAL)["auc_best"] < _diag(STRETCHED)["auc_best"]
