"""Unit tests for echo_repeat.reverb -- the constructive precision-vs-coverage
counterexample (CP1 harness reproducibility + Lattice read-only reuse, and the
CP2 counterexample-or-null gate).

These need the in-tree siblings Lattice (synthetic DRO) and Caliper (estimator +
conformal). They are skipped automatically if either is absent (e.g. when Lethe
is extracted as a standalone repo)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from echo_repeat import _paths, reverb as R  # noqa: E402

_HAVE_SIBLINGS = _paths.lattice_available() and _paths.caliper_available()
pytestmark = pytest.mark.skipif(not _HAVE_SIBLINGS,
                                reason="needs in-tree Lattice + Caliper siblings")


# --------------------------------------------------------------------------- #
# CP1 harness: reproducibility, test/retest structure, Lattice read-only reuse
# --------------------------------------------------------------------------- #
def test_reproducible_from_seed():
    a = R.simulate_testretest(n_eval=200, n_cal=200, snr=40.0, seed=R.REVERB_SEED)
    b = R.simulate_testretest(n_eval=200, n_cal=200, snr=40.0, seed=R.REVERB_SEED)
    assert np.array_equal(a.point_test, b.point_test)
    assert np.array_equal(a.q_test, b.q_test)
    assert np.array_equal(a.y_eval, b.y_eval)


def test_test_and_retest_share_truth_but_differ():
    tr = R.simulate_testretest(n_eval=200, n_cal=200, snr=40.0, seed=R.REVERB_SEED)
    # two acquisitions of one truth: estimates differ (independent noise)...
    assert not np.array_equal(tr.point_test, tr.point_retest)
    # ...but it is the SAME truth (one cohort, shared y_eval)
    assert tr.point_test.shape == tr.point_retest.shape == tr.y_eval.shape


def test_roi_acquire_reuses_lattice_generator_not_a_copy():
    # region_size=1 must be EXACTLY Lattice's Rician generator (no reimplementation)
    _paths.add_lattice()
    from lattice import generators as G
    clean = np.linspace(0.2, 1.0, 11)[None, :].repeat(7, axis=0)
    seedrng = lambda: np.random.default_rng(12345)
    via_reverb = R.roi_acquire(clean, 30.0, 1, seedrng())
    via_lattice = G.add_rician_noise(clean, 30.0, seedrng())
    assert np.array_equal(via_reverb, via_lattice)


def test_roi_acquire_reduces_spread_by_sqrt_n():
    # ROI-averaging must cut the random spread ~ sqrt(region_size) (the precision boost)
    _paths.add_lattice()
    clean = np.full((400, 11), 0.6)
    s1 = R.roi_acquire(clean, 25.0, 1, np.random.default_rng(1))
    s100 = R.roi_acquire(clean, 25.0, 100, np.random.default_rng(1))
    # std across regions shrinks markedly with averaging
    assert s100.std() < 0.25 * s1.std()


def test_truth_in_caliper_space_units_and_columns():
    _paths.add_lattice()
    import lattice
    c = lattice.make_cohort("biexp", n=500, snr=40, seed=R.REVERB_SEED, bvalues=np.array([0.0, 50, 800]))
    y = R.truth_in_caliper_space(c)
    D, f, Dstar = y[:, 0], y[:, 1], y[:, 2]
    # Caliper space: D,D* in 1e-3 mm^2/s (so D ~ [0.5,3], D* ~ [10,100]); f in [0.05,0.40]
    assert 0.4 < D.min() and D.max() < 3.2
    assert 9.0 < Dstar.min() and Dstar.max() < 101.0
    assert 0.04 < f.min() and f.max() < 0.41


# --------------------------------------------------------------------------- #
# Repeatability metric sanity
# --------------------------------------------------------------------------- #
def test_wcv_zero_and_icc_one_on_identical():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert R.within_subject_cv(x, x) == pytest.approx(0.0, abs=1e-12)
    assert R.icc_oneway(x, x) == pytest.approx(1.0, abs=1e-9)


def test_wcv_positive_under_noise():
    rng = np.random.default_rng(0)
    truth = rng.uniform(1.0, 2.0, 2000)
    a = truth + rng.normal(0, 0.1, 2000)
    b = truth + rng.normal(0, 0.1, 2000)
    assert R.within_subject_cv(a, b) > 0.02
    assert R.repeatability_pass(R.within_subject_cv(a, b), R.icc_oneway(a, b))


# --------------------------------------------------------------------------- #
# CP2 gate: the counterexample (mismatch) vs the matched control (null)
# --------------------------------------------------------------------------- #
def test_matched_control_tracks_no_counterexample():
    # correctly-specified model: precision and coverage track -> f@D*-lo not broken
    tr = R.simulate_testretest(n_eval=1500, n_cal=1500, snr=40.0, region_size=200,
                               truth_family="biexp", cal_family="biexp", seed=R.REVERB_SEED)
    cell = R.cell_point(tr, "f", 0)
    assert cell["repeatable"]
    assert not cell["coverage_broken"]
    assert not cell["counterexample"]


def test_model_mismatch_yields_precise_but_broken():
    # dispersed perfusion fit as bi-exp: f@D*-lo stays precise yet coverage breaks
    tr = R.simulate_testretest(n_eval=1500, n_cal=1500, snr=40.0, region_size=200,
                               truth_family="stretched", cal_family="biexp", seed=R.REVERB_SEED)
    cell = R.cell_point(tr, "f", 0)
    assert cell["repeatable"], cell          # precision intact
    assert cell["coverage_broken"], cell     # accuracy broken
    assert cell["counterexample"], cell


def test_precision_blind_to_bias():
    # the crux: repeatability is essentially identical control vs mismatch, while
    # coverage diverges -- repeatability cannot see the bias
    common = dict(n_eval=1500, n_cal=1500, snr=40.0, region_size=200, cal_family="biexp",
                  seed=R.REVERB_SEED)
    ctrl = R.cell_point(R.simulate_testretest(truth_family="biexp", **common), "f", 0)
    mis = R.cell_point(R.simulate_testretest(truth_family="stretched", **common), "f", 0)
    assert abs(ctrl["wcv"] - mis["wcv"]) < 0.02            # precision ~ unchanged
    assert ctrl["cov_split"] - mis["cov_split"] > 0.10     # coverage clearly drops


def test_analyze_verdict_and_marginal_looks_fine():
    tr = R.simulate_testretest(n_eval=1500, n_cal=1500, snr=40.0, region_size=200,
                               truth_family="stretched", cal_family="biexp", seed=R.REVERB_SEED)
    res = R.analyze(tr, n_boot=300, seed=0)
    assert res["counterexample_found"]
    # split-conformal restores ~nominal MARGINAL coverage -> "looks fine globally"
    assert 0.80 < res["marginal_coverage"]["split"]["f"] < 0.97
