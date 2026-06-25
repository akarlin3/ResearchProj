"""Frozen configuration for the Procrustes CP0 separation experiment.

A ``DepartureFamily`` names one non-bi-exponential Lattice generator, the knob
that tunes its departure from bi-exponential, and the value of that knob at the
exact bi-exponential limit (the built-in placebo / continuity point).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .seeding import GLOBAL_SEED


@dataclass(frozen=True)
class DepartureFamily:
    """A non-bi-exp Lattice family + its departure axis."""
    lattice_family: str          # key into lattice.FAMILIES
    knob: str                    # the `extra` parameter that tunes departure
    limit: float                 # knob value at the exact bi-exp limit (placebo)
    values: tuple                # ascending-departure sweep, values[?]==limit
    label: str                   # human label
    expect: str                  # pre-registered: "break" | "weak" | "null"


# The three departure families probed at the gate. The *mechanism* (validated in
# POSITIONING.md) is high-b aliasing: only departures that push perfusion signal
# into the high-b regime -- where the tissue D is read -- bias D and break its
# conditional coverage. Stretched-exp heavy tail => break; pure dispersion =>
# weak; a *faster* tri-exp third pool decays away from high-b => null.
STRETCHED = DepartureFamily("stretched", "beta", 1.0,
                            (1.0, 0.9, 0.8, 0.7, 0.6, 0.5), "stretched-exp (beta)", "break")
LOGNORMAL = DepartureFamily("dispersion_lognormal", "cv", 0.0,
                            (0.0, 0.2, 0.4, 0.6, 0.8), "log-normal dispersion (cv)", "weak")
TRIEXP = DepartureFamily("triexp", "g", 0.0,
                         (0.0, 0.1, 0.2, 0.3, 0.4), "tri-exponential third pool (g)", "null")

FAMILIES = {f.lattice_family: f for f in (STRETCHED, LOGNORMAL, TRIEXP)}


@dataclass(frozen=True)
class ProcrustesConfig:
    n: int = 1000               # voxels per cohort
    snr: float = 50.0           # Rician SNR
    alpha: float = 0.10         # target miscoverage => 90% intervals
    prior: str = "realistic"    # Lattice parameter prior
    noise: str = "rician"
    seed: int = GLOBAL_SEED
    target_param: str = "D"     # the *well-identified* tissue-diffusion map
    target_index: int = 0       # column of (D, Dstar, f)
    wellid_quantile: float = 2 / 3   # bottom-2-tercile D* = well-identified regime
    bvalues: tuple = None       # acquisition b-scheme; None => Lattice default (22-pt)

    # pre-registered gate thresholds (see POSITIONING.md, refute conditions)
    marginal_tol: float = 0.03       # |marginal - (1-alpha)| must stay within this
    break_gap_min: float = 0.06      # CP0-b: a "break" family's conditional gap
    wellid_gap_min: float = 0.03     # CP0-c: gap must survive in well-ID-D* subset
    null_gap_max: float = 0.03       # boundary: a "null" family's gap stays below this
    diagnostic_auc_min: float = 0.60 # diagnostic must clear chance for a "break" family

    @property
    def nominal(self) -> float:
        return 1.0 - self.alpha
