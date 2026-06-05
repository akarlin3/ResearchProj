"""
Reduced-ODE core for the two-population Sakaguchi-Kuramoto chimera.

All dynamical equations are transcribed VERBATIM from
  Abrams, Mirollo, Strogatz & Wiley, PRL 101, 084103 (2008),
  arXiv:0806.0594v2 (26 Aug 2008).
Equation numbers are cited at each function. Conventions match the shipped
finite-N integrator (tools/chimera-campaign/integrator.mjs), verified in the
campaign CP0:  mu=(1+A)/2, nu=(1-A)/2, alpha=pi/2-beta, omega=0.

The 2-D system is the rho1==1 invariant-manifold restriction (Eq. 12); the
3-D system is the full rotational-symmetry reduction (Eq. 10). We DERIVE the
3-D RHS and confirm it collapses to the 2-D RHS at rho1=1 (see test_core.py),
so a single sign error cannot hide in either.

Nothing here is fitted to the finite-N data; this module is pure dynamics.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
# Repo root = two levels up from tools/reduced-ode (input/output paths in the
# config are relative to it), independent of the caller's cwd.
ROOT = os.path.dirname(os.path.dirname(HERE))


def load_config(path=None):
    path = path or os.path.join(HERE, "reduced.config.json")
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    """Coupling parameters in the paper's conventions."""

    A: float
    beta: float

    @property
    def mu(self) -> float:
        return (1.0 + self.A) / 2.0

    @property
    def nu(self) -> float:
        return (1.0 - self.A) / 2.0

    @property
    def alpha(self) -> float:
        return math.pi / 2.0 - self.beta


# ---------------------------------------------------------------------------
# Right-hand sides
# ---------------------------------------------------------------------------
def rhs_2d(state, p: Params):
    """2-D invariant-manifold system, Eq. (12).  state=(r, psi), r=rho2.

      rdot   = ((1 - r^2)/2) [ mu r cos a + nu cos(psi - a) ]
      psidot = ((1 + r^2)/(2 r)) [ mu r sin a - nu sin(psi - a) ]
                 - mu sin a - nu r sin(psi + a)
    """
    r, psi = state
    mu, nu, a = p.mu, p.nu, p.alpha
    rdot = 0.5 * (1.0 - r * r) * (mu * r * math.cos(a) + nu * math.cos(psi - a))
    psidot = (
        (1.0 + r * r) / (2.0 * r) * (mu * r * math.sin(a) - nu * math.sin(psi - a))
        - mu * math.sin(a)
        - nu * r * math.sin(psi + a)
    )
    return np.array([rdot, psidot])


def rhs_3d(state, p: Params):
    """Full 3-D reduction (rho1, rho2, psi) of Eq. (10), psi = phi1 - phi2.

    Derived from Eq. (10) per population (a_sigma = rho_sigma e^{-i phi_sigma}),
    using cos(phi2-phi1-a)=cos(psi+a), cos(phi1-phi2-a)=cos(psi-a),
    sin(phi1-phi2+a)=sin(psi+a), sin(phi2-phi1+a)=sin(a-psi):

      rho1dot = ((1 - rho1^2)/2)[ mu rho1 cos a + nu rho2 cos(psi + a) ]
      rho2dot = ((1 - rho2^2)/2)[ mu rho2 cos a + nu rho1 cos(psi - a) ]
      phi1dot = -((1 + rho1^2)/(2 rho1))[ mu rho1 sin a + nu rho2 sin(psi + a) ]
      phi2dot = -((1 + rho2^2)/(2 rho2))[ mu rho2 sin a + nu rho1 sin(a - psi) ]
      psidot  = phi1dot - phi2dot

    At rho1==1 this reduces exactly to rhs_2d (verified in test_core.py).
    """
    r1, r2, psi = state
    mu, nu, a = p.mu, p.nu, p.alpha
    ca, sa = math.cos(a), math.sin(a)
    r1dot = 0.5 * (1.0 - r1 * r1) * (mu * r1 * ca + nu * r2 * math.cos(psi + a))
    r2dot = 0.5 * (1.0 - r2 * r2) * (mu * r2 * ca + nu * r1 * math.cos(psi - a))
    phi1dot = -(1.0 + r1 * r1) / (2.0 * r1) * (mu * r1 * sa + nu * r2 * math.sin(psi + a))
    phi2dot = -(1.0 + r2 * r2) / (2.0 * r2) * (mu * r2 * sa + nu * r1 * math.sin(a - psi))
    return np.array([r1dot, r2dot, phi1dot - phi2dot])


def transverse_rate(state, p: Params):
    """d(rho1dot)/d(rho1) evaluated at rho1=1 along a trajectory: the local
    transverse (Lyapunov-type) growth rate of the rho1=1 manifold.

    rho1dot = ((1-rho1^2)/2) G,  G = mu rho1 cos a + nu rho2 cos(psi+a).
    d/drho1 = -rho1 G + ((1-rho1^2)/2)(mu cos a).  At rho1=1: -G|_{rho1=1}.
    Negative => manifold attracting transversally.
    """
    _, r2, psi = state
    mu, nu, a = p.mu, p.nu, p.alpha
    G = mu * math.cos(a) + nu * r2 * math.cos(psi + a)
    return -G


# ---------------------------------------------------------------------------
# Stationary chimeras, Eqs. (13)-(14)
# ---------------------------------------------------------------------------
def fixed_point_family(beta, psi_grid):
    """Stationary-chimera family parametrised by psi, Eqs. (13)-(14).

      r = sqrt[ sin(2b+psi) / ( sin(2b-psi) + 2 sin psi ) ]          (14)
      A = [ sin(b+psi) + r sin b ] / [ sin(b+psi) - r sin b ]        (13)

    Returns arrays (psi, r, A) keeping only entries with a real r in (0,1)
    and finite A.
    """
    b = beta
    out_psi, out_r, out_A = [], [], []
    for psi in psi_grid:
        num = math.sin(2 * b + psi)
        den = math.sin(2 * b - psi) + 2 * math.sin(psi)
        if den == 0:
            continue
        rad = num / den
        if rad <= 0:
            continue
        r = math.sqrt(rad)
        if not (0.0 < r < 1.0):
            continue
        dA = math.sin(b + psi) - r * math.sin(b)
        if dA == 0:
            continue
        A = (math.sin(b + psi) + r * math.sin(b)) / dA
        out_psi.append(psi)
        out_r.append(r)
        out_A.append(A)
    return np.array(out_psi), np.array(out_r), np.array(out_A)


# ---------------------------------------------------------------------------
# Jacobian of the 2-D system (analytic-free, central differences)
# ---------------------------------------------------------------------------
def jacobian_2d(state, p: Params, h=1e-6):
    J = np.zeros((2, 2))
    for j in range(2):
        sp = np.array(state, dtype=float)
        sm = np.array(state, dtype=float)
        sp[j] += h
        sm[j] -= h
        J[:, j] = (rhs_2d(sp, p) - rhs_2d(sm, p)) / (2 * h)
    return J


def fp_stability(r, A, beta):
    """det, trace, eigenvalues of the 2-D Jacobian at fixed point (r,psi*),
    where psi* is recovered from r,A via Eq. (13). Returns dict."""
    p = Params(A=A, beta=beta)
    # psi from Eq (13): A (sin(b+psi) - r sinb) = sin(b+psi) + r sinb
    # => sin(b+psi) (A-1) = r sinb (A+1) => sin(b+psi) = r sinb (A+1)/(A-1)
    b = beta
    s = r * math.sin(b) * (A + 1.0) / (A - 1.0)
    s = max(-1.0, min(1.0, s))
    # two psi branches; pick the one that actually zeroes the RHS best
    best = None
    for psi in (math.asin(s) - b, math.pi - math.asin(s) - b):
        res = np.linalg.norm(rhs_2d([r, psi], p))
        if best is None or res < best[0]:
            best = (res, psi)
    psi = best[1]
    J = jacobian_2d([r, psi], p)
    ev = np.linalg.eigvals(J)
    return {
        "psi": psi,
        "rhs_resid": best[0],
        "det": float(np.linalg.det(J)),
        "trace": float(np.trace(J)),
        "eig": ev,
        "stable": bool(np.all(ev.real < 0)),
    }


# ---------------------------------------------------------------------------
# Bifurcation series, Eqs. (17)-(18)
# ---------------------------------------------------------------------------
def A_SN_series(beta, cfg):
    c = cfg["series"]["A_SN_coeffs"]
    b = beta
    return c["b1"] * b + c["b2"] * b**2 + c["b3"] * b**3 + c["b4"] * b**4 + c["b5"] * b**5


def A_H_series(beta, cfg):
    s = cfg["series"]
    b = beta
    return s["A_H_const"] + s["A_H_b2"] * b**2 + s["A_H_b4"] * b**4


# ---------------------------------------------------------------------------
# Numerical bifurcation location along the Eqs (13)-(14) family
# ---------------------------------------------------------------------------
def _det_trace_at_A(A, beta, branch="chimera"):
    """det(J), trace(J) of the 2-D Jacobian at the stationary chimera with this
    A on the chosen r-branch. The family is multivalued in A (a fold born at the
    saddle-node); the STABLE chimera sits on the lower-r branch (a spiral), the
    SADDLE on the upper-r branch. branch in {'chimera'|'lower', 'saddle'|'upper'}.
    """
    branch = {"chimera": "lower", "saddle": "upper"}.get(branch, branch)
    # Solve for r at fixed A by sweeping the family and interpolating: easier to
    # invert via the parametric family.
    psi_grid = np.linspace(-1.4, 1.4, 20001)
    psi, r, Aarr = fixed_point_family(beta, psi_grid)
    # candidate indices where A crosses the target
    rs = []
    for i in range(len(Aarr) - 1):
        if (Aarr[i] - A) * (Aarr[i + 1] - A) <= 0 and Aarr[i] != Aarr[i + 1]:
            frac = (A - Aarr[i]) / (Aarr[i + 1] - Aarr[i])
            rs.append(r[i] + frac * (r[i + 1] - r[i]))
    if not rs:
        return None
    rs = sorted(set(round(x, 9) for x in rs))
    rpick = rs[-1] if branch == "upper" else rs[0]
    st = fp_stability(rpick, A, beta)
    return st


def locate_saddle_node(beta, a_lo=0.001, a_hi=0.6):
    """SN = fold of the stationary-chimera family in A: the minimum A reached by
    the family (dA/dpsi = 0). Found by scanning the parametric family."""
    psi_grid = np.linspace(-1.4, 1.4, 200001)
    psi, r, A = fixed_point_family(beta, psi_grid)
    m = (A > a_lo) & (A < a_hi) & (r > 0) & (r < 1)
    if not np.any(m):
        return None
    Am = A[m]
    i = int(np.argmin(Am))
    return float(Am[i]), float(r[m][i]), float(psi[m][i])


def get_fixed_point(A, beta, branch="chimera"):
    """Return (r, psi, stability-dict) for the stationary chimera at this (A,beta)
    on the chosen r-branch, or None if the family does not reach this A.
    branch 'chimera' (lower r, stable spiral) or 'saddle' (upper r)."""
    branch = {"chimera": "lower", "saddle": "upper"}.get(branch, branch)
    st = _det_trace_at_A(A, beta, branch=branch)
    if st is None:
        return None
    # recover r from the family for this branch
    psi_grid = np.linspace(-1.4, 1.4, 200001)
    psi, r, Aarr = fixed_point_family(beta, psi_grid)
    rs = []
    for i in range(len(Aarr) - 1):
        if (Aarr[i] - A) * (Aarr[i + 1] - A) <= 0 and Aarr[i] != Aarr[i + 1]:
            frac = (A - Aarr[i]) / (Aarr[i + 1] - Aarr[i])
            rs.append(r[i] + frac * (r[i + 1] - r[i]))
    if not rs:
        return None
    rs = sorted(set(round(x, 9) for x in rs))
    rpick = rs[-1] if branch == "upper" else rs[0]
    return rpick, st["psi"], st


def locate_hopf(beta, a_lo=0.05, a_hi=0.6):
    """Hopf along the lower-r (stable-chimera) branch: trace(J)=0 with det(J)>0."""
    psi_grid = np.linspace(-1.4, 1.4, 60001)
    psi, r, A = fixed_point_family(beta, psi_grid)
    order = np.argsort(A)
    A_s, r_s = A[order], r[order]

    def tr(Aval):
        st = _det_trace_at_A(Aval, beta, branch="chimera")
        return st["trace"] if st else float("nan")

    # scan for a sign change of trace on the upper branch
    grid = np.linspace(max(a_lo, A_s.min() + 1e-4), min(a_hi, A_s.max() - 1e-4), 400)
    prev_A, prev_t = None, None
    for Aval in grid:
        t = tr(Aval)
        if prev_t is not None and np.isfinite(t) and np.isfinite(prev_t) and prev_t * t < 0:
            try:
                Ah = brentq(tr, prev_A, Aval, xtol=1e-8)
                st = _det_trace_at_A(Ah, beta, branch="chimera")
                if st and st["det"] > 0:
                    return float(Ah), float(st["psi"])
            except ValueError:
                pass
        prev_A, prev_t = Aval, t
    return None


# ---------------------------------------------------------------------------
# Homoclinic location (period divergence / collision with the saddle)
# ---------------------------------------------------------------------------
def _escaped_to_sync(A, beta, t_max, sync_r, npts=8000):
    """Integrate the 2-D system from just outside the chimera fixed point.
    Returns (escaped, info). escaped=True if r reaches sync (r>sync_r) — the
    breathing cycle has been destroyed (post-homoclinic, sync attracting).
    escaped=False if r stays bounded below sync (a persistent breathing cycle)."""
    fp = get_fixed_point(A, beta, branch="chimera")
    if fp is None:
        return None, {"reason": "no chimera FP"}
    r0, psi0, _ = fp
    p = Params(A=A, beta=beta)
    # perturb slightly outward (toward larger r, the saddle side)
    state0 = [min(0.995, r0 + 1e-3), psi0]

    def sync_event(t, y):
        return y[0] - sync_r

    sync_event.terminal = True
    sync_event.direction = 1
    t_eval = np.linspace(0, t_max, npts)
    sol = integrate_2d(state0, p, t_max, _CFG, t_eval=t_eval, events=sync_event)
    escaped = sol.t_events is not None and len(sol.t_events[0]) > 0
    rmax = float(np.max(sol.y[0]))
    info = {
        "escaped": bool(escaped),
        "t_escape": float(sol.t_events[0][0]) if escaped else None,
        "r_max": rmax,
        "r_final": float(sol.y[0][-1]),
    }
    return bool(escaped), info


def locate_homoclinic(beta, cfg):
    """Bisection on A for the homoclinic: the breathing cycle's period diverges
    as it collides with the saddle. Below A_hc a trajectory from the chimera FP
    settles on a bounded breathing cycle; above A_hc it escapes to sync (r->1).
    Returns dict with the bracketing interval and per-A escape diagnostics."""
    h = cfg["homoclinic"]
    a_lo, a_hi = h["A_lo"], h["A_hi"]
    t_max, sync_r, tol = h["t_max"], h["sync_r"], h["bisect_tol"]
    trace = []
    lo_esc, lo_info = _escaped_to_sync(a_lo, beta, t_max, sync_r)
    hi_esc, hi_info = _escaped_to_sync(a_hi, beta, t_max, sync_r)
    trace.append({"A": a_lo, **(lo_info or {})})
    trace.append({"A": a_hi, **(hi_info or {})})
    if lo_esc is None or hi_esc is None:
        return {"error": "FP missing at bracket", "trace": trace}
    if lo_esc or not hi_esc:
        # bracket invalid (lo should be bounded, hi should escape)
        return {
            "error": "bracket does not straddle homoclinic",
            "a_lo_escaped": lo_esc,
            "a_hi_escaped": hi_esc,
            "trace": trace,
        }
    while a_hi - a_lo > tol:
        mid = 0.5 * (a_lo + a_hi)
        esc, info = _escaped_to_sync(mid, beta, t_max, sync_r)
        trace.append({"A": mid, **(info or {})})
        if esc:
            a_hi = mid
        else:
            a_lo = mid
    return {
        "A_hc_bracket": [a_lo, a_hi],
        "A_hc": 0.5 * (a_lo + a_hi),
        "width": a_hi - a_lo,
        "trace": sorted(trace, key=lambda d: d["A"]),
    }


# ---------------------------------------------------------------------------
# Integration helpers
# ---------------------------------------------------------------------------
def integrate_2d(state0, p: Params, t_max, cfg, t_eval=None, events=None):
    ig = cfg["integrator"]
    sol = solve_ivp(
        lambda t, y: rhs_2d(y, p),
        (0.0, t_max),
        np.asarray(state0, float),
        method=ig["method"],
        rtol=ig["rtol"],
        atol=ig["atol"],
        max_step=ig["max_step_frac"] * t_max if ig.get("max_step_frac") else np.inf,
        t_eval=t_eval,
        events=events,
        dense_output=False,
    )
    return sol


def solve_ivp_branch(state0, p: Params, t_signed, npts):
    """Integrate the 2-D flow over signed time (negative = backward, for stable
    manifolds). Returns an (npts,2) array of (r, psi), clipped to r in (0,1.001)."""
    ig = _CFG["integrator"]
    t_eval = np.linspace(0, t_signed, npts)
    sol = solve_ivp(
        lambda t, y: rhs_2d(y, p),
        (0.0, t_signed),
        np.asarray(state0, float),
        method=ig["method"],
        rtol=ig["rtol"],
        atol=ig["atol"],
        t_eval=t_eval,
    )
    arr = sol.y.T.copy()
    # keep only physical r
    good = (arr[:, 0] > 0) & (arr[:, 0] < 1.001)
    return arr[good]


def integrate_3d(state0, p: Params, t_max, cfg, t_eval=None, events=None):
    ig = cfg["integrator"]
    sol = solve_ivp(
        lambda t, y: rhs_3d(y, p),
        (0.0, t_max),
        np.asarray(state0, float),
        method=ig["method"],
        rtol=ig["rtol"],
        atol=ig["atol"],
        max_step=ig["max_step_frac"] * t_max if ig.get("max_step_frac") else np.inf,
        t_eval=t_eval,
        events=events,
    )
    return sol


# ---------------------------------------------------------------------------
# Capture time in the 3-D system (maps the finite-N absorption label)
# ---------------------------------------------------------------------------
def capture_time_3d(state0, p: Params, t_max, cfg, npts=None):
    """First-passage 'capture' of the 3-D reduced flow, mapping the finite-N
    absorption-grade label: capture = first time min(rho1,rho2) crosses
    theta=0.85 AND never dips back below recoveryThreshold for recoveryWindowSec
    thereafter (no-recovery). In the deterministic reduced flow at A=0.5 the
    crossing is effectively permanent, so this reduces to the first sustained
    crossing. Also returns the transverse-stability summary of the rho1=1 manifold
    along the trajectory.

    Returns dict: captured, t_capture (or None=censored), m_max, transverse_*.
    """
    b = cfg["boundary"]
    theta, recThr, recWin = b["theta"], b["recoveryThreshold"], b["recoveryWindowSec"]
    sample_dt = cfg["timescale"]["sampleStride"]
    if npts is None:
        npts = int(round(t_max / sample_dt)) + 1
    t_eval = np.linspace(0.0, t_max, npts)
    sol = integrate_3d(state0, p, t_max, cfg, t_eval=t_eval)
    r1, r2 = sol.y[0], sol.y[1]
    m = np.minimum(r1, r2)
    t = sol.t
    dt = t[1] - t[0] if len(t) > 1 else sample_dt
    win = max(1, int(round(recWin / dt)))
    t_capture = None
    for i in range(len(m)):
        if m[i] > theta:
            # confirm no recovery (dip below recThr) within recWin after i
            j1 = min(len(m), i + win + 1)
            if np.all(m[i:j1] >= recThr):
                t_capture = float(t[i])
                break
    # transverse stability of rho1=1 along the trajectory
    tr_rates = np.array(
        [transverse_rate([r1[k], r2[k], sol.y[2][k]], p) for k in range(len(t))]
    )
    return {
        "captured": t_capture is not None,
        "t_capture": t_capture,
        "m_max": float(np.max(m)),
        "rho1_min": float(np.min(r1)),
        "rho1_dev_max": float(np.max(np.abs(1.0 - r1))),
        "transverse_rate_mean": float(np.mean(tr_rates)),
        "transverse_rate_max": float(np.max(tr_rates)),
        "transverse_attracting_frac": float(np.mean(tr_rates < 0)),
    }


def reduced_run_3d(state0, p: Params, t_max, cfg, theta=None):
    """Single 3-D integration yielding BOTH the capture time and the breath/spiral
    metrics, with an early-stop event (min(rho1,rho2) >= 0.97, safely past theta
    and the no-recovery band) so capturing runs integrate only to ~capture, not
    t_max. Used by CP3 and CP4.

    Capture maps the finite-N absorption-grade label: first time min(rho1,rho2)
    crosses theta AND stays >= recoveryThreshold for recoveryWindowSec after
    (no-recovery). Breath period = mean spacing of min-maxima; spiral_slope = OLS
    of log(theta - peak) vs cycle index over SUB-theta peaks (the finite-N
    estimator). Also returns the transverse stability of rho1=1 along the path.
    """
    b = cfg["boundary"]
    theta = b["theta"] if theta is None else theta
    recThr, recWin = b["recoveryThreshold"], b["recoveryWindowSec"]
    sample_dt = cfg["timescale"]["sampleStride"]
    npts = int(round(t_max / sample_dt)) + 1
    t_eval = np.linspace(0.0, t_max, npts)

    stop_level = 0.97

    def reach(t, y):
        return min(y[0], y[1]) - stop_level

    reach.terminal = True
    reach.direction = 1
    sol = integrate_3d(state0, p, t_max, cfg, t_eval=t_eval, events=reach)
    r1, r2, psi = sol.y[0], sol.y[1], sol.y[2]
    m = np.minimum(r1, r2)
    t = sol.t
    dt = t[1] - t[0] if len(t) > 1 else sample_dt
    win = max(1, int(round(recWin / dt)))

    # capture = first sustained theta crossing (no recovery within recWin)
    t_capture = None
    cap_idx = None
    for i in range(len(m)):
        if m[i] > theta:
            j1 = min(len(m), i + win + 1)
            if np.all(m[i:j1] >= recThr):
                t_capture = float(t[i])
                cap_idx = i
                break

    # breath/spiral over the pre-capture transient
    hi = cap_idx if cap_idx is not None else len(m)
    mpre = m[: hi + 1] if hi < len(m) else m
    tpre = t[: hi + 1] if hi < len(m) else t
    pk = [
        i for i in range(1, len(mpre) - 1)
        if mpre[i] > mpre[i - 1] and mpre[i] >= mpre[i + 1]
    ]
    ptimes = tpre[pk]
    pvals = mpre[pk]
    period = float(np.mean(np.diff(ptimes))) if len(ptimes) >= 2 else None
    sub = pvals < theta
    n_sub = int(np.sum(sub))
    slope = None
    if n_sub >= 3:
        k = np.arange(len(pvals))[sub]
        y = np.log(np.clip(theta - pvals[sub], 1e-9, None))
        A_ = np.vstack([k, np.ones_like(k)]).T
        slope = float(np.linalg.lstsq(A_, y, rcond=None)[0][0])

    tr_rates = np.array(
        [transverse_rate([r1[k], r2[k], psi[k]], p) for k in range(len(t))]
    )
    return {
        "captured": t_capture is not None,
        "t_capture": t_capture,
        "m_max": float(np.max(m)),
        "rho1_dev_max": float(np.max(np.abs(1.0 - r1))),
        "breath_period": period,
        "n_peaks": len(pk),
        "n_sub_theta_peaks": n_sub,
        "spiral_slope": slope,
        "transverse_rate_mean": float(np.mean(tr_rates)),
        "transverse_attracting_frac": float(np.mean(tr_rates < 0)),
        "peak_vals": pvals.tolist(),
        "peak_times": ptimes.tolist(),
    }


def breath_metrics_3d(state0, p: Params, t_max, cfg, theta, npts=None):
    """Per-run breath/spiral metrics from the 3-D transient: successive maxima of
    min(rho1,rho2), their period (mean spacing), and the spiral-out slope
    = linear fit of log(theta - m_max_k) vs cycle index k over SUB-theta maxima
    (matching the finite-N estimator). Returns dict (period, slope, n_peaks...)."""
    sample_dt = cfg["timescale"]["sampleStride"]
    if npts is None:
        npts = int(round(t_max / sample_dt)) + 1
    t_eval = np.linspace(0.0, t_max, npts)
    sol = integrate_3d(state0, p, t_max, cfg, t_eval=t_eval)
    m = np.minimum(sol.y[0], sol.y[1])
    t = sol.t
    # local maxima
    pk = [i for i in range(1, len(m) - 1) if m[i] > m[i - 1] and m[i] >= m[i + 1]]
    ptimes = t[pk]
    pvals = m[pk]
    period = float(np.mean(np.diff(ptimes))) if len(ptimes) >= 2 else None
    # spiral slope over sub-theta peaks
    sub = pvals < theta
    slope = None
    n_sub = int(np.sum(sub))
    if n_sub >= 3:
        k = np.arange(len(pvals))[sub]
        y = np.log(theta - pvals[sub])
        A_ = np.vstack([k, np.ones_like(k)]).T
        slope = float(np.linalg.lstsq(A_, y, rcond=None)[0][0])
    return {
        "period": period,
        "n_peaks": len(pk),
        "n_sub_theta_peaks": n_sub,
        "spiral_slope": slope,
        "peak_vals": pvals.tolist(),
        "peak_times": ptimes.tolist(),
    }


# ---------------------------------------------------------------------------
# Limit-cycle classification & period (2-D)
# ---------------------------------------------------------------------------
def classify_2d(p: Params, r0, psi0, t_settle, t_window, const_tol=1e-4, npts=20000):
    """Integrate, discard transient, then classify the steady state of r(t):
    'fixed' (r ~ constant), 'cycle' (bounded oscillation, returns period), or
    'sync' (r -> 1). Period from mean spacing of r-maxima."""
    t_eval = np.linspace(0.0, t_settle + t_window, npts)
    sol = integrate_2d([r0, psi0], p, t_settle + t_window, _CFG, t_eval=t_eval)
    t, r = sol.t, sol.y[0]
    mask = t >= t_settle
    tw, rw = t[mask], r[mask]
    if rw.size < 10:
        return {"kind": "error"}
    if rw.min() > 0.999:
        return {"kind": "sync", "r_final": float(rw[-1])}
    amp = rw.max() - rw.min()
    if amp < const_tol:
        return {"kind": "fixed", "r_star": float(np.mean(rw)), "amp": float(amp)}
    # local maxima
    peaks = [
        i
        for i in range(1, len(rw) - 1)
        if rw[i] > rw[i - 1] and rw[i] >= rw[i + 1]
    ]
    period = None
    if len(peaks) >= 2:
        dts = np.diff(tw[peaks])
        # use the later half (settled)
        dts = dts[len(dts) // 2 :] if len(dts) > 2 else dts
        period = float(np.mean(dts))
    return {
        "kind": "cycle",
        "amp": float(amp),
        "r_min": float(rw.min()),
        "r_max": float(rw.max()),
        "period": period,
        "n_peaks": len(peaks),
    }


# Module-level config handle (set by load at import for integrate_* convenience).
_CFG = load_config()


def set_config(cfg):
    global _CFG
    _CFG = cfg
