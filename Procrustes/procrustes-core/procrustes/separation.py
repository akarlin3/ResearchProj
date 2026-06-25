"""The CP0 separation: marginal holds, conditional fails, distinct from Gauge.

For one departure family we observe the SAME voxel population under increasing
departure from bi-exponential and ask three questions about the conformal
intervals for the *well-identified* tissue-diffusion map D:

  (a) MARGINAL holds      -- a single departure-blind conformal radius gives
                             ~nominal pooled coverage (confirms Lei 2018; NOT the
                             contribution).
  (b) CONDITIONAL fails   -- coverage degrades monotonically along the latent
                             departure axis; an oracle that knows the departure
                             recalibrates it away (within-departure radius).
  (c) DISTINCT from Gauge -- the failure persists, undiminished, inside the
                             *well-identified* D* subset (bottom-2 terciles),
                             where Gauge's high-D* identifiability wall does not
                             operate and Gauge's triage rule says "trust D".

The signed bias of D-hat, growing with departure, is the mechanism: heavy-tailed
perfusion misspecification aliases into the high-b tissue slope.

Pre-registered REFUTE (any => wedge dead for this family):
  R1  conditional gap ~ 0 (no departure-driven failure),
  R2  gap vanishes inside the well-identified D* subset (it was Gauge all along),
  R3  the bi-exp-limit placebo row is itself broken (artefact, not misspec).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .conformal import conformal_radius, coverage, split_indices
from .generators import cohort_at, fit_cohort
from .seeding import make_rng


@dataclass
class SeparationResult:
    family: str
    knob: str
    values: tuple
    limit: float
    marginal: float                 # (a) pooled departure-blind coverage
    cond: dict                      # value -> conditional coverage (global radius)
    cond_wellid: dict               # value -> coverage in well-identified D* subset
    within: dict                    # value -> within-departure recalibrated coverage
    bias: dict                      # value -> signed mean bias of D-hat
    radius: float
    cond_lo: dict = None            # value -> coverage in strict bottom-D* tercile
    cond_hi: dict = None            # value -> coverage in top-D* tercile (Gauge's wall)
    n_wellid: int = 0               # voxels in the well-ID (bottom-2-tercile) subset
    n_lo: int = 0                   # voxels in the strict bottom tercile

    @property
    def worst(self) -> float:
        return self.values[-1] if self.values[-1] != self.limit else self.values[0]

    @property
    def gap(self) -> float:
        return self.cond[self.limit] - self.cond[self.worst]

    @property
    def wellid_gap(self) -> float:
        return self.cond_wellid[self.limit] - self.cond_wellid[self.worst]

    @property
    def lo_gap(self) -> float:
        """Conditional gap in the STRICT bottom-D* tercile (deepest in Gauge's
        identifiable region; the strongest form of refute R2)."""
        return self.cond_lo[self.limit] - self.cond_lo[self.worst]

    @property
    def hi_gap(self) -> float:
        """Conditional gap in the TOP-D* tercile -- Gauge's own ill-posed wall."""
        return self.cond_hi[self.limit] - self.cond_hi[self.worst]

    @property
    def bias_ratio(self) -> float:
        return abs(self.bias[self.worst]) / (abs(self.bias[self.limit]) + 1e-9)


def run_separation(family_cfg, cfg, seed: int | None = None) -> SeparationResult:
    """Run the three-part CP0 separation for one departure family."""
    seed = cfg.seed if seed is None else seed
    fam, knob, values, limit = (family_cfg.lattice_family, family_cfg.knob,
                                family_cfg.values, family_cfg.limit)
    ti = cfg.target_index

    cohorts = {v: cohort_at(fam, knob, v, cfg) for v in values}
    ref = cohorts[limit]
    truth = ref.params[:, ti].copy()
    Dstar = ref.Dstar.copy()
    # same base parameters across departure is the property the placebo relies on
    for v in values:
        assert np.allclose(cohorts[v].params, ref.params), "base-param invariance violated"

    fits = {v: fit_cohort(cohorts[v], cfg.snr) for v in values}
    Dhat = {v: fits[v]["Dhat"] for v in values}

    calib, test = split_indices(cfg.n, make_rng(seed))

    # (a) departure-blind global conformal radius -> marginal coverage
    radius = conformal_radius(
        np.concatenate([np.abs(Dhat[v][calib] - truth[calib]) for v in values]), cfg.alpha)
    marginal = coverage(
        np.concatenate([np.abs(Dhat[v][test] - truth[test]) for v in values]), radius)

    # (b) conditional coverage by departure + signed bias
    cond, bias, within = {}, {}, {}
    for v in values:
        cond[v] = coverage(np.abs(Dhat[v][test] - truth[test]), radius)
        bias[v] = float(np.mean(Dhat[v][test] - truth[test]))
        rv = conformal_radius(np.abs(Dhat[v][calib] - truth[calib]), cfg.alpha)
        within[v] = coverage(np.abs(Dhat[v][test] - truth[test]), rv)

    # (c) well-identified D* subset: does the gap SURVIVE inside Gauge's
    # identifiable region?  Use Gauge's EXACT partition of the latent axis --
    # true-D* terciles via np.quantile([1/3, 2/3]) + np.digitize (see
    # Gauge/gauge/conditional.py:157-159 and conditional_attack._regime_from_true):
    #   regime 0 = lo D*, 1 = mid D*, 2 = hi D* (Gauge's ill-posed wall).
    # Gauge's identifiable region (where it says "trust D") is the bottom-2
    # terciles {0, 1}; the strict bottom tercile {0} is the deepest interior of
    # it (the strongest refute-R2 test).
    edges = np.quantile(Dstar, [1.0 / 3.0, 2.0 / 3.0])
    regime = np.digitize(Dstar, edges)              # 0=lo, 1=mid, 2=hi (Gauge)
    tmask = np.zeros(cfg.n, bool); tmask[test] = True
    sub = tmask & (regime <= 1)                      # well-ID = bottom-2 terciles
    sub_lo = tmask & (regime == 0)                   # strict bottom tercile
    sub_hi = tmask & (regime == 2)                   # Gauge's high-D* wall
    cond_wellid = {v: coverage(np.abs(Dhat[v][sub] - truth[sub]), radius) for v in values}
    cond_lo = {v: coverage(np.abs(Dhat[v][sub_lo] - truth[sub_lo]), radius) for v in values}
    cond_hi = {v: coverage(np.abs(Dhat[v][sub_hi] - truth[sub_hi]), radius) for v in values}

    return SeparationResult(fam, knob, values, limit, marginal,
                            cond, cond_wellid, within, bias, radius,
                            cond_lo=cond_lo, cond_hi=cond_hi,
                            n_wellid=int(sub.sum()), n_lo=int(sub_lo.sum()))


def separation_detail(family_cfg, cfg, seed: int | None = None) -> dict:
    """Per-voxel material for bootstrap CIs of the separation gaps.

    Same computation as :func:`run_separation` but returns the raw TEST-voxel
    absolute errors at the placebo (limit) and worst-departure knobs, the fixed
    departure-blind conformal radius, and the Gauge-exact D*-regime masks
    (well-ID = bottom-2 terciles; lo = strict bottom tercile; hi = Gauge's
    high-D* wall) -- all aligned to the test-voxel order so a paired voxel
    bootstrap can resample them jointly.
    """
    seed = cfg.seed if seed is None else seed
    fam, knob, values, limit = (family_cfg.lattice_family, family_cfg.knob,
                                family_cfg.values, family_cfg.limit)
    ti = cfg.target_index
    worst = values[-1] if values[-1] != limit else values[0]

    cohorts = {v: cohort_at(fam, knob, v, cfg) for v in (limit, worst)}
    ref = cohorts[limit]
    truth = ref.params[:, ti].copy()
    Dstar = ref.Dstar.copy()
    fits = {v: fit_cohort(cohorts[v], cfg.snr) for v in (limit, worst)}
    Dhat = {v: fits[v]["Dhat"] for v in (limit, worst)}

    calib, test = split_indices(cfg.n, make_rng(seed))
    # departure-blind radius from BOTH knobs' calibration residuals (marginal recipe)
    radius = conformal_radius(
        np.concatenate([np.abs(Dhat[v][calib] - truth[calib]) for v in (limit, worst)]),
        cfg.alpha)

    signed_limit = Dhat[limit][test] - truth[test]
    signed_worst = Dhat[worst][test] - truth[test]
    err_limit = np.abs(signed_limit)
    err_worst = np.abs(signed_worst)
    # signed-bias mechanism: heavy-tail perfusion aliases into the high-b tissue
    # slope and biases D-hat; the bias grows from placebo to worst departure.
    bias_limit = float(np.mean(signed_limit))
    bias_worst = float(np.mean(signed_worst))
    bias_ratio = abs(bias_worst) / (abs(bias_limit) + 1e-9)
    # marginal coverage under the SAME departure-blind radius (placebo+worst
    # pooled) -- so the marginal claim and the gaps share one radius recipe.
    marginal = coverage(np.concatenate([err_limit, err_worst]), radius)

    edges = np.quantile(Dstar, [1.0 / 3.0, 2.0 / 3.0])
    regime = np.digitize(Dstar, edges)
    reg_test = regime[test]
    return {
        "family": fam, "knob": knob, "limit": limit, "worst": worst,
        "radius": float(radius), "marginal": float(marginal),
        "bias_limit": bias_limit, "bias_worst": bias_worst, "bias_ratio": bias_ratio,
        "err_limit": err_limit, "err_worst": err_worst,
        "wellid": reg_test <= 1, "lo": reg_test == 0, "hi": reg_test == 2,
        "all": np.ones(len(test), bool),
        "n_test": len(test),
        "n_wellid": int((reg_test <= 1).sum()), "n_lo": int((reg_test == 0).sum()),
    }
