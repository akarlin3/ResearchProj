"""Reverb -- Lethe's *constructive* counterexample: precision without coverage.

Lethe's central claim is asserted on real ACRIN-6698 data, where it can only be
*argued*: repeatability is precision, not coverage, so a ground-truth-free
validator can certify that an interval is *sized* to scan--rescan noise without
certifying that it *covers the truth*. With no per-voxel ground truth in vivo,
the gap between the two cannot be exhibited -- only reasoned about.

Reverb closes that gap with a controlled synthetic demonstration. It draws a
**Lattice** cohort with known ground truth ``(D, D*, f)``, simulates a same-truth
test--retest pair (two independent noise draws of one clean signal), fits both
with **Caliper**'s reference IVIM estimator, deploys a conformal interval, and
then -- because the truth is known -- measures *both*:

  * **repeatability** (precision): test--retest agreement of the point estimate
    (within-subject CV and ICC), and
  * **coverage** (accuracy): whether the deployed interval contains the *true*
    parameter, marginally and *conditionally* per true-D* regime.

The construction exhibits a regime -- expected at high D*, the identifiability
wall Gauge and Fashion both localise -- where repeatability looks acceptable
while conditional truth-coverage is broken. That is precision != coverage turned
from *asserted* into *shown*.

Scope (load-bearing). This is a synthetic *possibility-and-mechanism* proof: it
shows the divergence **can** occur in IVIM and **why**. It does **not** quantify
the magnitude of any real-world miscalibration. Read every number here as "this
divergence can occur, and here is the mechanism," never as "this is how
miscalibrated real IVIM is."

Reuse posture (guardrail). Reverb is a thin test--retest + analysis layer. The
synthetic cohort/generator is **Lattice** (read-only); the estimator, conformal
correction, and coverage ruler are **Caliper** (read-only); the bootstrap is
``echo_repeat.statistic`` (in-package). Reverb duplicates none of them.

Parameter-space note. Lattice's ground truth is ordered ``(D, D*, f)`` with
``D, D*`` in mm^2/s; Caliper's estimator and conformal layer use ``(D, f, D*)``
with ``D, D*`` in 1e-3 mm^2/s. Everything in this module works in **Caliper
space** -- :func:`truth_in_caliper_space` converts Lattice ground truth into it
once, so estimates, truth, intervals, and strata are all column- and
unit-consistent.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import _paths
from . import statistic as st

# --- reproducibility constants (no wall-clock; everything seeds from these) ----
REVERB_SEED = 20260622          # date-stamped default seed
_OFFSET_CAL = 0                 # calibration cohort truth stream
_OFFSET_EVAL = 555              # eval (test/retest) cohort truth stream
_OFFSET_CAL_NOISE = 1300        # calibration single-draw noise stream
_OFFSET_TEST_NOISE = 2100       # eval "test" acquisition noise stream
_OFFSET_RETEST_NOISE = 2200     # eval "retest" acquisition noise stream (independent of test)

# Caliper conventions (see module docstring).
PARAM_NAMES = ("D", "f", "Dstar")
DSTAR_COL = 2                   # stratify on true D* (Caliper column order)
LEVELS = (0.05, 0.25, 0.50, 0.75, 0.95)   # symmetric; central 0.90 at alpha=0.10
ALPHA = 0.10
STRATUM_NAMES = {0: "D*-lo", 1: "D*-mid", 2: "D*-hi"}

# Pre-registered bars (stated before any run; not tuned to the result).
WCV_BAR = 0.20      # within-subject CV <= 0.20 => "acceptable repeatability" (QIBA-style)
ICC_BAR = 0.75      # ICC >= 0.75 => "good" agreement (Koo & Li 2016)
COVERAGE_BROKEN = 0.65   # conditional truth-coverage point < 0.65 => "broken" (Lethe's own bar)
NOMINAL = 1.0 - ALPHA    # 0.90


# ===========================================================================
# Lattice <-> Caliper bridges (read-only use of Lattice)
# ===========================================================================
def truth_in_caliper_space(cohort) -> np.ndarray:
    """Convert a Lattice cohort's ground truth ``(D, D*, f)`` mm^2/s to Caliper
    ``(D, f, D*)`` in 1e-3 mm^2/s. Returns ``(n, 3)``."""
    D = cohort.params[:, 0] * 1e3
    Dstar = cohort.params[:, 1] * 1e3
    f = cohort.params[:, 2]
    return np.column_stack([D, f, Dstar])


def roi_acquire(clean, snr, region_size, rng):
    """One region-level acquisition: ROI-mean of ``region_size`` noisy voxels.

    Each region is homogeneous (its ``region_size`` voxels share one clean
    signal ``clean[i]``); the acquisition draws ``region_size`` *independent*
    Rician voxels and averages them, exactly Lethe's whole-tumor ROI-mean
    reduction. ROI-averaging cuts the random measurement spread by ~``sqrt(
    region_size)`` (the effective-SNR boost Lethe's real-data model uses) while
    leaving any *systematic* bias (Rician floor, model mismatch) untouched -- the
    decoupling of precision from accuracy this counterexample turns on.

    ``region_size=1`` reduces to a single per-voxel Rician draw. The Rician
    generator is Lattice's (``add_rician_noise``) -- reused, not reimplemented.
    """
    _paths.add_lattice()
    from lattice import generators as G
    clean = np.asarray(clean, dtype=float)
    if region_size <= 1:
        return G.add_rician_noise(clean, snr, rng)
    n, n_b = clean.shape
    tiled = np.broadcast_to(clean[:, None, :], (n, region_size, n_b))
    voxels = G.add_rician_noise(tiled, snr, rng)   # (n, region_size, n_b)
    return voxels.mean(axis=1)                      # ROI-mean -> (n, n_b)


# ===========================================================================
# CP1 -- the test--retest harness on Lattice
# ===========================================================================
@dataclass
class TestRetest:
    """One synthetic test--retest experiment on a Lattice cohort (Caliper space).

    Calibration split (single acquisition) + an evaluation split acquired twice
    (test + retest) from one shared ground truth. All point estimates / quantiles
    are produced by Caliper's reference estimator; ``y_*`` are Lattice ground
    truth converted to Caliper space.
    """
    bvalues: np.ndarray
    levels: np.ndarray
    snr: float
    seed: int
    region_size: int
    cal_family: str           # the believed/modeled family the interval is calibrated on
    truth_family: str         # the deployment "reality" the test/retest are drawn from
    # calibration split
    q_cal: np.ndarray         # (n_cal, P, L) raw quantiles
    y_cal: np.ndarray         # (n_cal, P) truth
    # evaluation split (same truth, two acquisitions)
    y_eval: np.ndarray        # (n_eval, P) truth
    point_test: np.ndarray    # (n_eval, P)
    point_retest: np.ndarray  # (n_eval, P)
    q_test: np.ndarray        # (n_eval, P, L) raw quantiles
    q_retest: np.ndarray      # (n_eval, P, L) raw quantiles
    dstar_true: np.ndarray    # (n_eval,) true D* for stratification
    param_names: tuple = PARAM_NAMES

    def __len__(self) -> int:
        return self.y_eval.shape[0]

    @property
    def model_mismatch(self) -> bool:
        """True when deployment reality differs from the calibrated model."""
        return self.cal_family != self.truth_family

    def manifest(self) -> dict:
        return {
            "n_cal": int(self.y_cal.shape[0]),
            "n_eval": int(self.y_eval.shape[0]),
            "region_size": int(self.region_size),
            "cal_family": self.cal_family,
            "truth_family": self.truth_family,
            "model_mismatch": self.model_mismatch,
            "bvalues": self.bvalues.tolist(),
            "levels": self.levels.tolist(),
            "snr": float(self.snr),
            "seed": int(self.seed),
            "param_names": list(self.param_names),
            "param_units": "Caliper space: (D, f, D*); D,D* in 1e-3 mm^2/s",
        }


def simulate_testretest(
    n_eval: int = 2000,
    n_cal: int = 2000,
    snr: float = 40.0,
    region_size: int = 1,
    truth_family: str = "biexp",
    cal_family: str = "biexp",
    seed: int = REVERB_SEED,
    bvalues: np.ndarray | None = None,
    levels=LEVELS,
    estimator_kwargs: dict | None = None,
) -> TestRetest:
    """Build a reproducible region-level test--retest experiment on Lattice (CP1).

    Steps (each reuses Lattice/Caliper read-only):
      1. Draw a calibration cohort from the *believed* model ``cal_family``
         (default ``biexp``) and an evaluation cohort from deployment *reality*
         ``truth_family``; both carry known ground truth on the estimator's
         b-schedule.
      2. Calibration split: one region-level acquisition -> Caliper quantiles +
         truth. The conformal interval will be calibrated here, so it is blind to
         any ``cal_family != truth_family`` mismatch -- exactly Lethe's posture of
         calibrating on a synthetic model and deploying on reality.
      3. Evaluation split: TWO independent region-level acquisitions of the *same*
         clean truth (test, retest) -> two point fits + two quantile blocks.

    With ``truth_family == cal_family`` the model is correctly specified (the
    control, which tracks); with a dispersed ``truth_family`` fit by the bi-exp
    estimator, a structural mismatch bias appears that ROI-averaging makes precise
    but cannot remove -- the counterexample. Fully determined by ``seed`` (no
    wall-clock); identical inputs reproduce bit-for-bit.
    """
    _paths.add_lattice()
    _paths.add_caliper()
    import lattice
    from caliper.estimator_reference import ReferenceIVIMEstimator
    from caliper.forward import DEFAULT_BVALUES as CALIPER_BVALUES

    levels = np.asarray(levels, dtype=float)
    b = CALIPER_BVALUES.copy() if bvalues is None else np.asarray(bvalues, dtype=float)
    est = ReferenceIVIMEstimator(bvalues=b, **(estimator_kwargs or {}))

    # --- calibration cohort (believed model): one region-level acquisition ----
    cal = lattice.make_cohort(cal_family, n=n_cal, snr=snr, seed=seed + _OFFSET_CAL,
                              prior="realistic", noise="none", bvalues=b)
    rng_cal = np.random.default_rng(seed + _OFFSET_CAL_NOISE)
    sig_cal = roi_acquire(cal.signals_clean, snr, region_size, rng_cal)
    q_cal = est.predict_quantiles(sig_cal, levels)
    y_cal = truth_in_caliper_space(cal)

    # --- evaluation cohort (deployment reality): same truth acquired twice ----
    ev = lattice.make_cohort(truth_family, n=n_eval, snr=snr, seed=seed + _OFFSET_EVAL,
                             prior="realistic", noise="none", bvalues=b)
    rng_test = np.random.default_rng(seed + _OFFSET_TEST_NOISE)
    rng_retest = np.random.default_rng(seed + _OFFSET_RETEST_NOISE)
    sig_test = roi_acquire(ev.signals_clean, snr, region_size, rng_test)
    sig_retest = roi_acquire(ev.signals_clean, snr, region_size, rng_retest)
    y_eval = truth_in_caliper_space(ev)

    return TestRetest(
        bvalues=b, levels=levels, snr=float(snr), seed=int(seed),
        region_size=int(region_size), cal_family=cal_family, truth_family=truth_family,
        q_cal=q_cal, y_cal=y_cal,
        y_eval=y_eval,
        point_test=est.predict_point(sig_test),
        point_retest=est.predict_point(sig_retest),
        q_test=est.predict_quantiles(sig_test, levels),
        q_retest=est.predict_quantiles(sig_retest, levels),
        dstar_true=y_eval[:, DSTAR_COL],
    )


# ===========================================================================
# Repeatability (precision) metrics -- bias-blind by construction
# ===========================================================================
def within_subject_cv(a, b) -> float:
    """QIBA-style within-subject coefficient of variation from paired repeats.

    ``wCV = sqrt( mean_i[ (d_i / m_i)^2 / 2 ] )`` with ``d_i = a_i - b_i`` and
    ``m_i = (a_i + b_i)/2``. Dimensionless; small = repeatable.
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = 0.5 * (a + b)
    ok = np.abs(m) > 0
    d = (a[ok] - b[ok]) / m[ok]
    return float(np.sqrt(np.mean(d ** 2) / 2.0)) if d.size else float("nan")


def icc_oneway(a, b) -> float:
    """One-way random-effects ICC(1,1) agreement for paired repeats.

    ``ICC = (MSB - MSW) / (MSB + MSW)`` with between/within mean squares for
    ``k=2`` repeats. 1 = perfect agreement, 0 = none; scale-relative, so it is
    the natural companion to wCV.
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    n = a.size
    if n < 2:
        return float("nan")
    means = 0.5 * (a + b)
    msw = float(np.mean(((a - b) ** 2) / 2.0))       # within-subject MS
    msb = 2.0 * float(np.var(means, ddof=1))          # between-subject MS (k=2)
    denom = msb + msw
    return float((msb - msw) / denom) if denom > 0 else float("nan")


def repeatability_pass(wcv: float, icc: float) -> bool:
    """Pre-registered repeatability bar: ``wCV <= WCV_BAR`` AND ``ICC >= ICC_BAR``."""
    return bool(np.isfinite(wcv) and np.isfinite(icc) and wcv <= WCV_BAR and icc >= ICC_BAR)


# ===========================================================================
# CP2 -- the counterexample: repeatability vs coverage, per regime
# ===========================================================================
def _estimator_bounds(bvalues):
    """The clip bounds Caliper's reference estimator applies to (D, f, D*).

    Railing onto a bound is the cleanest *precise-but-biased* mechanism (a bound
    is the same every draw -> repeatable -> but wrong if truth is interior); it
    is exactly Fashion's boundary-railing diagnostic. Bounds mirror
    ``estimator_reference.predict_point``.
    """
    return {"D": (0.1, 4.0), "f": (0.001, 0.95), "Dstar": (1.0, 300.0)}


def _rail_fraction(values, name) -> float:
    lo, hi = _estimator_bounds(None)[name]
    v = np.asarray(values, float)
    at = np.isclose(v, lo, rtol=1e-3, atol=1e-9) | np.isclose(v, hi, rtol=1e-3, atol=1e-9)
    return float(np.mean(at)) if v.size else float("nan")


def conformal_intervals(tr: TestRetest):
    """Deploy three intervals on the test acquisition (Caliper, read-only).

    Returns ``(intervals, groups_eval)`` where ``intervals`` maps method ->
    ``(lo, hi)`` each ``(n_eval, P)`` central ``(1-ALPHA)`` interval:

      * ``raw``      -- the estimator's own over-confident quantiles (no fix);
      * ``split``    -- split-conformal/CQR: restores *marginal* coverage (the
        "looks fine globally" certificate a ground-truth-free pipeline issues);
      * ``mondrian`` -- group-conditional CQR on true-D* terciles: the oracle
        fix that restores conditional coverage *only by spending width*.

    ``groups_eval`` are the true-D* terciles used as the conditional lens.
    """
    _paths.add_caliper()
    from caliper.conformal import SplitConformalQuantile, MondrianConformalQuantile
    from caliper.metrics import central_interval, tercile_groups

    levels = tr.levels
    g_cal = tercile_groups(tr.y_cal[:, DSTAR_COL])
    g_eval = tercile_groups(tr.dstar_true)

    lo_raw, hi_raw = central_interval(tr.q_test, levels, ALPHA)

    q_split = SplitConformalQuantile(levels).calibrate(tr.q_cal, tr.y_cal).apply(tr.q_test)
    lo_s, hi_s = central_interval(q_split, levels, ALPHA)

    q_mond = (MondrianConformalQuantile(levels)
              .calibrate(tr.q_cal, tr.y_cal, g_cal)
              .apply(tr.q_test, g_eval))
    lo_m, hi_m = central_interval(q_mond, levels, ALPHA)

    intervals = {
        "raw": (lo_raw, hi_raw),
        "split": (lo_s, hi_s),
        "mondrian": (lo_m, hi_m),
    }
    return intervals, g_eval


def _bca(stat_fn, n, n_boot, seed):
    point, lo, hi = st.bca_ci(stat_fn, n, n_boot=n_boot, alpha=0.05, seed=seed)
    return point, (lo, hi)


def analyze(tr: TestRetest, n_boot: int = 2000, seed: int = 0) -> dict:
    """Compute the per-regime repeatability-vs-coverage table and the verdict.

    For each parameter and each true-D* tercile, report repeatability (wCV, ICC,
    with BCa CIs and the pre-registered pass), the marginally-conformal interval's
    *conditional* truth-coverage (with a BCa CI), the Mondrian fix's coverage and
    width, and the mechanism localisers (mean true error, boundary-railing
    fraction). The counterexample is any (parameter, regime) where repeatability
    PASSES while split-conformal conditional coverage is BROKEN.
    """
    _paths.add_caliper()
    from caliper.metrics import empirical_coverage

    intervals, g_eval = conformal_intervals(tr)
    lo_s, hi_s = intervals["split"]
    lo_m, hi_m = intervals["mondrian"]
    lo_r, hi_r = intervals["raw"]
    strata = sorted(np.unique(g_eval).tolist())

    rows = []
    counterexamples = []
    for p, pname in enumerate(tr.param_names):
        for k, g in enumerate(strata):
            idx = np.where(g_eval == g)[0]
            n_g = idx.size
            a = tr.point_test[idx, p]; bvec = tr.point_retest[idx, p]
            yt = tr.y_eval[idx, p]
            # per-parameter column slices of each interval, restricted to this stratum
            lo_sp = lo_s[idx, p]; hi_sp = hi_s[idx, p]
            lo_rp = lo_r[idx, p]; hi_rp = hi_r[idx, p]
            lo_mp = lo_m[idx, p]; hi_mp = hi_m[idx, p]

            # repeatability (precision) + BCa CIs
            wcv, wcv_ci = _bca(lambda ii: within_subject_cv(a[ii], bvec[ii]),
                               n_g, n_boot, seed + 10 * p + g)
            icc, icc_ci = _bca(lambda ii: icc_oneway(a[ii], bvec[ii]),
                               n_g, n_boot, seed + 1000 + 10 * p + g)
            repeatable = repeatability_pass(wcv, icc)

            # coverage (accuracy vs known truth), split-conformal, + BCa CI
            cov_split, cov_split_ci = _bca(
                lambda ii: empirical_coverage(yt[ii], lo_sp[ii], hi_sp[ii]),
                n_g, n_boot, seed + 2000 + 10 * p + g)
            cov_raw = float(empirical_coverage(yt, lo_rp, hi_rp))
            cov_mond = float(empirical_coverage(yt, lo_mp, hi_mp))
            broken = bool(cov_split < COVERAGE_BROKEN and cov_split_ci[1] < NOMINAL)

            is_ce = bool(repeatable and broken)
            row = {
                "param": pname, "stratum": int(g), "stratum_name": STRATUM_NAMES.get(int(g), f"g{g}"),
                "n": int(n_g),
                "wcv": wcv, "wcv_ci": list(wcv_ci),
                "icc": icc, "icc_ci": list(icc_ci),
                "repeatable": repeatable,
                "cov_raw": cov_raw,
                "cov_split": cov_split, "cov_split_ci": list(cov_split_ci),
                "cov_mondrian": cov_mond,
                "width_split": float(np.mean(hi_s[idx, p] - lo_s[idx, p])),
                "width_mondrian": float(np.mean(hi_m[idx, p] - lo_m[idx, p])),
                "coverage_broken": broken,
                "true_abs_err": float(np.mean(np.abs(a - yt))),
                "rail_frac": _rail_fraction(a, pname),
                "counterexample": is_ce,
            }
            rows.append(row)
            if is_ce:
                counterexamples.append((pname, int(g)))

    # marginal coverage per method (the "looks fine globally" check)
    marginal = {}
    for meth, (lo, hi) in intervals.items():
        marginal[meth] = {
            pname: float(empirical_coverage(tr.y_eval[:, p], lo[:, p], hi[:, p]))
            for p, pname in enumerate(tr.param_names)
        }

    # "track together" null check: do coverage-broken regimes also fail repeatability?
    broken_rows = [r for r in rows if r["coverage_broken"]]
    track_together = bool(broken_rows) and all(not r["repeatable"] for r in broken_rows)

    return {
        "manifest": tr.manifest(),
        "bars": {"wcv_bar": WCV_BAR, "icc_bar": ICC_BAR,
                 "coverage_broken": COVERAGE_BROKEN, "nominal": NOMINAL},
        "n_boot": int(n_boot),
        "rows": rows,
        "marginal_coverage": marginal,
        "counterexamples": counterexamples,
        "counterexample_found": bool(counterexamples),
        "track_together_null": track_together,
    }


def cell_point(tr: TestRetest, param: str = "f", stratum: int = 0) -> dict:
    """No-bootstrap repeatability + split-conformal coverage for one (param, stratum).

    The lightweight cell used by the sensitivity sweep; mirrors one row of
    :func:`analyze` without the BCa bootstrap.
    """
    _paths.add_caliper()
    from caliper.metrics import empirical_coverage
    p = list(tr.param_names).index(param)
    intervals, g_eval = conformal_intervals(tr)
    lo_s, hi_s = intervals["split"]
    m = g_eval == stratum
    a = tr.point_test[m, p]; bvec = tr.point_retest[m, p]; t = tr.y_eval[m, p]
    wcv = within_subject_cv(a, bvec); icc = icc_oneway(a, bvec)
    cov = float(empirical_coverage(t, lo_s[m, p], hi_s[m, p]))
    return {"param": param, "stratum": stratum, "n": int(m.sum()),
            "wcv": wcv, "icc": icc, "repeatable": repeatability_pass(wcv, icc),
            "cov_split": cov, "coverage_broken": bool(cov < COVERAGE_BROKEN),
            "counterexample": bool(repeatability_pass(wcv, icc) and cov < COVERAGE_BROKEN)}


def sensitivity_sweep(
    families=("biexp", "dispersion_gamma", "dispersion_lognormal", "stretched", "triexp"),
    region_sizes=(50, 100, 200, 400),
    snr: float = 40.0,
    n: int = 1500,
    seed: int = REVERB_SEED,
    param: str = "f",
    stratum: int = 0,
) -> dict:
    """Pre-registered honesty surface for the counterexample cell (default f@D*-lo).

    Sweeps deployment-truth family (matched ``biexp`` control + dispersed
    realities) against ROI size at fixed SNR, calibrating each on ``biexp``.
    Reports the cell's wCV and split-conformal coverage so the reader sees the
    whole monotone surface -- the matched control that never breaks, and the
    dispersed families that break increasingly with ROI size. No cell is hidden.
    """
    grid = {}
    for fam in families:
        grid[fam] = {}
        for rs in region_sizes:
            tr = simulate_testretest(n_eval=n, n_cal=n, snr=snr, region_size=rs,
                                     truth_family=fam, cal_family="biexp", seed=seed)
            c = cell_point(tr, param=param, stratum=stratum)
            grid[fam][int(rs)] = {"wcv": c["wcv"], "cov_split": c["cov_split"],
                                  "repeatable": c["repeatable"], "broken": c["coverage_broken"],
                                  "counterexample": c["counterexample"]}
    return {"param": param, "stratum": stratum, "stratum_name": STRATUM_NAMES.get(stratum),
            "snr": snr, "n": n, "families": list(families),
            "region_sizes": list(region_sizes), "grid": grid}


def format_report(result: dict) -> str:
    """Render the analysis as a human-readable markdown report."""
    m = result["manifest"]; bars = result["bars"]
    L = [
        "# Reverb -- constructive counterexample: precision without coverage",
        "",
        f"Lattice biexp cohort (read-only), Caliper estimator + conformal (read-only). "
        f"n_cal={m['n_cal']}, n_eval={m['n_eval']}, SNR={m['snr']}, seed={m['seed']}, "
        f"b-values={m['bvalues']}.",
        f"Params in Caliper space (D, f, D*); D,D* in 1e-3 mm^2/s. "
        f"Bars: wCV<={bars['wcv_bar']}, ICC>={bars['icc_bar']} => repeatable; "
        f"split-conformal conditional coverage <{bars['coverage_broken']} (CI hi < {bars['nominal']}) => broken.",
        "",
        "## Marginal coverage (the 'looks fine globally' certificate)",
        "| method | D | f | D* |",
        "|---|---|---|---|",
    ]
    for meth in ("raw", "split", "mondrian"):
        mc = result["marginal_coverage"][meth]
        L.append(f"| {meth} | {mc['D']:.3f} | {mc['f']:.3f} | {mc['Dstar']:.3f} |")
    L += [
        "",
        "## Per-regime: repeatability (precision) vs conditional truth-coverage (accuracy)",
        "| param | regime | n | wCV [CI] | ICC [CI] | repeatable | cov(split) [CI] | cov(mondrian) | width split->mond | true|err| | rail | COUNTEREX |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in result["rows"]:
        L.append(
            f"| {r['param']} | {r['stratum_name']} | {r['n']} | "
            f"{r['wcv']:.2f} [{r['wcv_ci'][0]:.2f},{r['wcv_ci'][1]:.2f}] | "
            f"{r['icc']:+.2f} [{r['icc_ci'][0]:+.2f},{r['icc_ci'][1]:+.2f}] | "
            f"{'YES' if r['repeatable'] else 'no'} | "
            f"{r['cov_split']:.2f} [{r['cov_split_ci'][0]:.2f},{r['cov_split_ci'][1]:.2f}] | "
            f"{r['cov_mondrian']:.2f} | "
            f"{r['width_split']:.3g}->{r['width_mondrian']:.3g} | "
            f"{r['true_abs_err']:.3g} | {r['rail_frac']:.0%} | "
            f"{'**YES**' if r['counterexample'] else '-'} |"
        )
    L.append("")
    if result["counterexample_found"]:
        ces = ", ".join(f"{p} @ {STRATUM_NAMES.get(g, g)}" for p, g in result["counterexamples"])
        L.append(f"**Counterexample(s) found:** {ces} -- repeatability clears its bar while "
                 f"split-conformal conditional truth-coverage is broken. Precision != coverage, shown.")
    elif result["track_together_null"]:
        L.append("**Honest null:** every broken-coverage regime *also* fails repeatability "
                 "(precision and coverage track together). This would weaken Lethe's claim and is "
                 "reported as the finding, not tuned away.")
    else:
        L.append("**Honest null:** no regime pairs acceptable repeatability with broken coverage. "
                 "Reported as the finding.")
    return "\n".join(L)
