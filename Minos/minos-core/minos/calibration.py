"""The decision–calibration gap (Headline B).

Two scales for the reported error bar, computed by **independent** code paths:

* ``tau_stat`` — *statistical* calibration: the scale giving nominal central-interval
  coverage of the truth. Root-find on coverage; reads only :mod:`diagnostics` (the error
  law and the reported intervals). **It never touches the utility.**
* ``tau_star`` — *decision* calibration: ``argmax_tau E[U_posterior(tau)]``. Dense grid +
  local refine; reads only :mod:`voi` (the decision/utility core). **It never touches
  coverage.**

The gap ``G = tau_star - tau_stat`` (and ratio ``tau_star/tau_stat``) is the quantity the
neighbours (DCA/net-benefit, ISPOR VoI, ARCliDS) do not compute. Math: DESIGN_B Section 3.

This module also promotes the v1 trust-gate footnote to a first-class number: the
break-even shift ``delta_be`` where ``VoTG(delta) = 0`` (DESIGN_B Section 5 / GATE 3).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .config import MinosConfig
from .generative import BaseDraws

# Default search windows / reference level (documented modelling choices, DESIGN_B Sec 3).
TAU_STAT_BRACKET = (0.2, 5.0)
TAU_STAR_BOUNDS = (0.5, 2.5)
TAU_STAR_GRID_STEP = 0.05
TAU_STAR_FIT_HALFWIDTH = 0.35   # parabola-fit window half-width (in tau) around the peak
FLAT_EU_TOL = 1e-4              # EU(tau) span below this => decision is tau-insensitive
L_REF = 0.90


# --------------------------------------------------------------------------------------
# tau_stat — statistical calibration (coverage root-find; utility never read)
# --------------------------------------------------------------------------------------
def coverage_at(base: BaseDraws, cfg: MinosConfig, tau: float, *, level: float = L_REF,
                delta: float = 0.0, shift=False) -> float:
    """Central-``level`` interval coverage of ``N(mu, (tau*s)^2)`` against the truth.

    Thin wrapper over :func:`minos.diagnostics.central_interval_coverage` — this function
    (and ``tau_stat``) imports *only* the diagnostics path, never the decision core.
    """
    from .diagnostics import central_interval_coverage

    return central_interval_coverage(base, cfg, level, tau=tau, delta=delta, shift=shift)


def tau_stat(base: BaseDraws, cfg: MinosConfig, *, level: float = L_REF,
             delta: float = 0.0, shift=False, bracket=TAU_STAT_BRACKET) -> float:
    """Scale at which the reported interval has nominal coverage: root of ``C(tau)-level``.

    Coverage is monotone increasing in ``tau`` (wider interval covers more), so the root
    is unique; we bracket it on ``[0.2, 5]`` (coverage -> 0 / 1 at the ends).
    """
    lo, hi = bracket
    f = lambda t: coverage_at(base, cfg, t, level=level, delta=delta, shift=shift) - level
    flo, fhi = f(lo), f(hi)
    assert flo < 0 < fhi, f"coverage root not bracketed on {bracket}: f(lo)={flo}, f(hi)={fhi}"
    return float(brentq(f, lo, hi, xtol=1e-4, rtol=1e-6))


# --------------------------------------------------------------------------------------
# tau_star — decision calibration (grid + local refine; coverage never read)
# --------------------------------------------------------------------------------------
def posterior_eu(base: BaseDraws, cfg: MinosConfig, tau: float, *,
                 delta: float = 0.0, shift=False) -> float:
    """``E[U(a*(N(mu,(tau s)^2)), theta)]`` — the decision-side objective for ``tau_star``.

    Thin wrapper over :func:`minos.voi.expected_utility` — this function (and
    ``tau_star``) imports *only* the decision/utility path, never coverage.
    """
    from .voi import expected_utility

    return expected_utility("posterior", base, cfg, tau=tau, delta=delta, shift=shift)


def tau_star(base: BaseDraws, cfg: MinosConfig, *, delta: float = 0.0, shift=False,
             bounds=TAU_STAR_BOUNDS, step: float = TAU_STAR_GRID_STEP,
             halfwidth: float = TAU_STAR_FIT_HALFWIDTH):
    """``argmax_tau E[U_posterior(tau)]`` by a dense grid + local **parabola-fit** refine.

    Returns ``(tau, eu)``. CRN makes ``EU(tau)`` a smooth deterministic function of ``tau`` for a
    given sample, but near a *flat* optimum (the well-specified ``kappa=0`` case) a raw argmax
    chases finite-sample bumps. Fitting a quadratic to ``EU(tau)`` over a window around the grid
    peak and taking its vertex is robust to those bumps: it estimates the optimum from the curve's
    shape, not one noisy point. The vertex is clamped to the fit window; a flat curve (the
    ``FLAT_EU_TOL`` guard, e.g. symmetric utility lambda=1, where the escalate boundary is
    tau-independent) or a non-concave fit returns ``tau*=1`` (no warranted deviation).

    Note: ``tau*`` is only well-identified where the decision actually depends on ``tau`` — report
    centres near a threshold. Far from any threshold the EU curve is degenerately flat and ``tau*``
    is unidentified; the gap sweep keeps ``rho`` in that resolved regime (DESIGN_B Section 5).
    """
    lo, hi = bounds
    grid = np.round(np.arange(lo, hi + step / 2, step), 4)
    eus = np.array([posterior_eu(base, cfg, t, delta=delta, shift=shift) for t in grid])
    one_eu = float(posterior_eu(base, cfg, 1.0, delta=delta, shift=shift))
    if float(np.ptp(eus)) < FLAT_EU_TOL:
        return 1.0, one_eu
    i = int(np.argmax(eus))
    win = np.abs(grid - grid[i]) <= halfwidth + 1e-9
    x, y = grid[win], eus[win]
    if x.size >= 3:
        a, b, c = np.polyfit(x, y, 2)
        if a < 0:  # concave -> interior maximum at the vertex
            vtx = -b / (2.0 * a)
            if x.min() <= vtx <= x.max():
                return float(vtx), float(np.polyval([a, b, c], vtx))
    return 1.0, one_eu


# --------------------------------------------------------------------------------------
# the gap
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class GapResult:
    tau_stat: float
    tau_star: float
    gap: float          # tau_star - tau_stat
    ratio: float        # tau_star / tau_stat
    level: float
    eu_star: float

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return (f"tau_stat={self.tau_stat:.4f}  tau_star={self.tau_star:.4f}  "
                f"G={self.gap:+.4f}  ratio={self.ratio:.4f}  (L={self.level})")


def gap(base: BaseDraws, cfg: MinosConfig, *, level: float = L_REF,
        delta: float = 0.0, shift=False) -> GapResult:
    """Compute ``(tau_stat, tau_star, G, ratio)`` via the two independent paths."""
    ts = tau_stat(base, cfg, level=level, delta=delta, shift=shift)
    tstar, eu = tau_star(base, cfg, delta=delta, shift=shift)
    return GapResult(tau_stat=ts, tau_star=tstar, gap=tstar - ts,
                     ratio=tstar / ts, level=level, eu_star=eu)


# --------------------------------------------------------------------------------------
# break-even shift for the trust-gate (CP3 — promotes the v1 footnote to a number)
# --------------------------------------------------------------------------------------
def break_even_shift(base: BaseDraws, cfg: MinosConfig, *, tau: float = 1.0,
                     bracket=(0.0, 1.5)) -> float:
    """Shift ``delta_be`` at which ``VoTG(delta) = 0`` (gate stops being a net cost).

    ``VoTG(0) < 0`` (the fixed false-positive cost with no corruption to repair) and
    ``VoTG(delta) > 0`` once the shift is severe enough; the crossing is ``delta_be``.
    """
    from .gate import votg

    lo, hi = bracket
    f = lambda d: votg(base, cfg, delta=d, tau=tau)
    flo, fhi = f(lo), f(hi)
    assert flo < 0 < fhi, f"VoTG root not bracketed on {bracket}: f(lo)={flo}, f(hi)={fhi}"
    return float(brentq(f, lo, hi, xtol=1e-4, rtol=1e-6))
