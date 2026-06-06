"""Euler-Maruyama integration of the 3-variable reduced flow with finite-size
additive noise sigma_noise = c/sqrt(N) on the collective variables (rho1, rho2,
psi).

Reuses tools/reduced-ode/reduced_core.py for the deterministic drift (rhs_3d) and
Params. The capture / breath-period / spiral-slope detection replicates
reduced_run_3d EXACTLY (same theta / recovery rule, same sub-theta spiral fit, on
min(rho1,rho2) sampled at sampleStride), so c=0 reproduces the deterministic
capture time.

Noise model: dX = f(X) dt + sigma_noise dW, sigma_noise = c/sqrt(N) (the SDE
diffusion coefficient), applied to all three collective variables. rho1, rho2 are
kept in [0,1] by REFLECTION at the boundaries (documented; clipping gives nearly
identical medians but reflection conserves probability). psi evolves unbounded
(it enters the drift only through trig functions).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reduced-ode"))
import reduced_core as rc  # noqa: E402


def _reflect(x):
    # reflect into [0,1]
    if x > 1.0:
        return 2.0 - x
    if x < 0.0:
        return -x
    return x


def _detect(m, dt_sample, theta, recThr, recWin):
    """Capture time + breath period + spiral slope from a sampled min(rho1,rho2)
    series — verbatim port of reduced_run_3d's detection."""
    win = max(1, int(round(recWin / dt_sample)))
    t = np.arange(len(m)) * dt_sample
    t_capture, cap_idx = None, None
    for i in range(len(m)):
        if m[i] > theta:
            j1 = min(len(m), i + win + 1)
            if np.all(m[i:j1] >= recThr):
                t_capture = float(t[i])
                cap_idx = i
                break
    hi = cap_idx if cap_idx is not None else len(m)
    mpre = m[: hi + 1] if hi < len(m) else m
    tpre = t[: hi + 1] if hi < len(m) else t
    pk = [i for i in range(1, len(mpre) - 1)
          if mpre[i] > mpre[i - 1] and mpre[i] >= mpre[i + 1]]
    ptimes = tpre[pk]
    pvals = mpre[pk]
    period = float(np.mean(np.diff(ptimes))) if len(ptimes) >= 2 else None
    sub = pvals < theta
    n_sub = int(np.sum(sub))
    slope = None
    if n_sub >= 3:
        k = np.arange(len(pvals))[sub]
        y = np.log(np.clip(theta - pvals[sub], 1e-9, None))
        Amat = np.vstack([k, np.ones_like(k)]).T
        slope = float(np.linalg.lstsq(Amat, y, rcond=None)[0][0])
    # capture phase relative to the preceding breath peak (for phase-locking test)
    capture_phase = None
    if t_capture is not None and period and len(ptimes) >= 1:
        frac = ((t_capture - ptimes[-1]) / period) % 1.0
        capture_phase = float(2 * np.pi * frac)
    return dict(captured=t_capture is not None, t_capture=t_capture,
                breath_period=period, n_peaks=len(pk), n_sub_theta_peaks=n_sub,
                spiral_slope=slope, n_cycles=max(0, len(pk) - 1),
                capture_phase=capture_phase, peak_vals=pvals.tolist())


def em_run(state0, A, beta, c, N, cfg, seed, dt=0.05, t_max=2000.0):
    """One Euler-Maruyama realization. Returns the same detection dict as
    reduced_run_3d (capture time, breath period, spiral slope, cycle count)."""
    p = rc.Params(A=A, beta=beta)
    b = cfg["boundary"]
    theta, recThr, recWin = b["theta"], b["recoveryThreshold"], b["recoveryWindowSec"]
    dt_sample = cfg["timescale"]["sampleStride"]
    sample_every = max(1, int(round(dt_sample / dt)))
    sigma = c / np.sqrt(N)
    sqdt = np.sqrt(dt)
    rng = np.random.default_rng(seed)
    n_steps = int(round(t_max / dt))

    x = np.array(state0, float)
    m_samples = [min(x[0], x[1])]
    stop_level = 0.97
    for step in range(1, n_steps + 1):
        f = rc.rhs_3d(x, p)
        if sigma > 0:
            xi = rng.standard_normal(3)
            x = x + dt * np.asarray(f) + sigma * sqdt * xi
        else:
            x = x + dt * np.asarray(f)
        x[0] = _reflect(x[0])
        x[1] = _reflect(x[1])
        if step % sample_every == 0:
            mm = min(x[0], x[1])
            m_samples.append(mm)
            # early stop once safely captured (mirrors reduced_run_3d's event)
            if mm >= stop_level and len(m_samples) > int(round(recWin / dt_sample)):
                break
    return _detect(np.asarray(m_samples), dt_sample, theta, recThr, recWin)
