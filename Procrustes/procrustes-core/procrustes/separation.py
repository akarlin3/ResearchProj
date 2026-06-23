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

    # (c) well-identified D* subset (bottom-2 terciles): does the gap survive?
    cut = np.quantile(Dstar, cfg.wellid_quantile)
    tmask = np.zeros(cfg.n, bool); tmask[test] = True
    sub = tmask & (Dstar <= cut)
    cond_wellid = {v: coverage(np.abs(Dhat[v][sub] - truth[sub]), radius) for v in values}

    return SeparationResult(fam, knob, values, limit, marginal,
                            cond, cond_wellid, within, bias, radius)
