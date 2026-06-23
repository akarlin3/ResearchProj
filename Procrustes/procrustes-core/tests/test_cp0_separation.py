"""GATE CP0 -- the three-part separation on the heavy-tail (stretched) family.

Each test pins one pre-registered claim; a regression that kills the wedge trips
the matching gate.  Fast config (small N, single seed); the publication numbers
live in experiments/run_cp0.py + RESULTS.md.
"""
import pytest

from procrustes import ProcrustesConfig, STRETCHED, run_separation

CFG = ProcrustesConfig(n=500)


@pytest.fixture(scope="module")
def res():
    return run_separation(STRETCHED, CFG, seed=CFG.seed)


def test_marginal_holds(res):
    # (a) departure-blind conformal keeps ~nominal MARGINAL coverage despite the
    # misspecified bi-exp fit (Lei 2018) -- this is the trap, not the contribution.
    assert abs(res.marginal - CFG.nominal) <= 0.05


def test_conditional_fails(res):
    # (b) + REFUTE R1: conditional coverage of D must degrade along the latent
    # departure axis. gap ~ 0 would mean misspecification does not break it -> dead.
    assert res.gap > CFG.break_gap_min, f"conditional gap {res.gap:.3f} too small (R1)"


def test_distinct_from_gauge(res):
    # (c) + REFUTE R2: the failure must SURVIVE inside the well-identified D*
    # subset, where Gauge's high-D* wall does not operate ("trust D"). If it
    # vanishes there, it was Gauge's identifiability wall all along -> dead.
    assert res.wellid_gap > CFG.wellid_gap_min, f"well-ID gap {res.wellid_gap:.3f} (R2)"


def test_placebo_not_broken(res):
    # REFUTE R3: the bi-exp-limit (beta=1) placebo row must NOT itself under-cover
    # -- otherwise the binned failure is an artefact, not misspecification.
    assert res.cond[res.limit] >= CFG.nominal - 0.04


def test_bias_is_the_mechanism(res):
    # Signed D-bias must grow with departure: heavy-tail perfusion aliases into
    # the high-b tissue slope. (Variance-only inflation would not break it here.)
    assert res.bias_ratio > 2.0


def test_within_departure_recalibration_restores(res):
    # The failure is CONDITIONAL on the (unobservable) departure: an oracle that
    # recalibrates within the departure level recovers ~nominal coverage.
    assert abs(res.within[res.worst] - CFG.nominal) <= 0.06
