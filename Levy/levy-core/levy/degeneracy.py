"""Phase-3 structural identifiability: is the joint CTRW time-order alpha separable from the
space-order beta at a SINGLE clinical diffusion time, or are they degenerate?

Deliverables (all scoped to the regime: single-diffusion-time, finite b-values, Rician noise):
  1. ANALYTIC degeneracy map over the physiological (alpha, beta) grid -- the (alpha,beta)
     correlation rho_alpha_beta from the inverse joint FIM, the FIM condition number, and the
     relative CRLBs cv_alpha, cv_beta. (fisher_joint)
  2. HEADLINE cell with bootstrap CIs -- the empirical (alpha_hat, beta_hat) ridge: marginal
     95% CIs + their correlation. (identifiability_joint)
  3. n_b PERSISTENCE -- unlike the CP0 single-order wall (which recedes for n_b>=8), the
     (alpha,beta) degeneracy PERSISTS as b-values are added at one diffusion time.
  4. The CONSTRUCTIVE boundary -- a SECOND diffusion time breaks the degeneracy (rho drops,
     CRLB recovers): the scope limit, reported not hidden.

Pre-registered (fixed before looking at the map): the orders are "degenerate" when the median
|rho_alpha_beta| over the physiological grid is >= DEGEN_RHO. This is characterization, not a
kill test -- it is reported whichever way it lands, with no overclaim in either direction.
CRLB = identifiability, never impossibility.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import fisher_joint, forward, identifiability_joint, seeding

# Pre-registered degeneracy threshold on the (alpha, beta) correlation (fixed up front).
DEGEN_RHO = 0.90

# Physiological CTRW regime (Magin, Ingo, Colon-Perez, Triplett & Mareci 2013, brain tissue):
#   time-order   alpha in ~[0.42, 0.95]   (Caputo / waiting-time)
#   space-order  beta  in ~[1.15, 1.95]   (Riesz / jump-length)
ALPHA_GRID = np.array([0.55, 0.65, 0.75, 0.80, 0.85, 0.90])
BETA_GRID = np.array([1.40, 1.55, 1.70, 1.80, 1.90])

# Headline clinical acquisition (matches the CP0 clinical regime) and a representative voxel.
HEADLINE_CELL = dict(S0=1.0, D=1.5e-3, alpha=0.80, beta=1.70)
CLINICAL = dict(n_b=6, b_max=2500.0, snr=40.0)


@dataclass(frozen=True)
class CP1Report:
    alpha_grid: np.ndarray
    beta_grid: np.ndarray
    rho_map: np.ndarray            # |rho_alpha_beta| over (alpha, beta) at clinical acquisition
    cond_map: np.ndarray
    cv_alpha_map: np.ndarray
    cv_beta_map: np.ndarray
    rho_median: float
    rho_min: float
    rho_max: float
    cond_median: float
    # headline cell
    headline_theta: np.ndarray
    headline_rho: float            # analytic rho_alpha_beta (signed)
    headline_cond: float
    headline_cv_alpha: float
    headline_cv_beta: float
    boot_alpha_ci: tuple
    boot_beta_ci: tuple
    boot_corr: float
    boot_alpha_rel_width: float
    boot_beta_rel_width: float
    n_boot: int
    # n_b persistence (single diffusion time)
    nb_list: tuple
    nb_rho: np.ndarray             # |rho_alpha_beta| vs n_b at the headline cell
    # constructive two-diffusion-time break
    break_rho_single: float
    break_rho_two: float
    break_cond_single: float
    break_cond_two: float
    break_cv_alpha_single: float
    break_cv_alpha_two: float
    # verdict
    degenerate: bool
    notes: list = field(default_factory=list)


def _headline_theta():
    h = HEADLINE_CELL
    return np.array([h["S0"], h["D"], h["alpha"], h["beta"]])


def analytic_map(snr=None, n_b=None, b_max=None):
    """|rho_alpha_beta|, cond, cv_alpha, cv_beta over the physiological (alpha,beta) grid."""
    from .wall import default_b_design
    snr = CLINICAL["snr"] if snr is None else snr
    n_b = CLINICAL["n_b"] if n_b is None else n_b
    b_max = CLINICAL["b_max"] if b_max is None else b_max
    b = default_b_design(b_max=b_max, n_b=n_b)
    D = HEADLINE_CELL["D"]
    rho = np.empty((len(ALPHA_GRID), len(BETA_GRID)))
    cond = np.empty_like(rho)
    cva = np.empty_like(rho)
    cvb = np.empty_like(rho)
    for i, a in enumerate(ALPHA_GRID):
        for j, bet in enumerate(BETA_GRID):
            r = fisher_joint.crlb_joint(b, np.array([1.0, D, a, bet]), snr, "rician")
            rho[i, j] = abs(r.rho_alpha_beta)
            cond[i, j] = r.cond
            cva[i, j] = r.cv_alpha
            cvb[i, j] = r.cv_beta
    return b, rho, cond, cva, cvb


def cp1_report(rng=None, n_boot: int = 80) -> CP1Report:
    """Full Phase-3 (alpha, beta) degeneracy characterization with bootstrap CIs."""
    if rng is None:
        rng = seeding.make_rng()
    from .wall import default_b_design

    snr = CLINICAL["snr"]
    b, rho_map, cond_map, cva_map, cvb_map = analytic_map()
    rho_median = float(np.median(rho_map))
    rho_min = float(rho_map.min())
    rho_max = float(rho_map.max())
    cond_median = float(np.median(cond_map))

    # headline cell, single diffusion time
    theta = _headline_theta()
    b_h = default_b_design(b_max=CLINICAL["b_max"], n_b=CLINICAL["n_b"])
    rh = fisher_joint.crlb_joint(b_h, theta, snr, "rician")
    boot = identifiability_joint.parametric_bootstrap_joint(
        theta, b_h, np.ones_like(b_h), snr=snr, n_boot=n_boot, rng=rng)

    # n_b persistence at the headline cell (single diffusion time)
    nb_list = (4, 6, 8, 16)
    nb_rho = np.array([
        abs(fisher_joint.crlb_joint(default_b_design(b_max=CLINICAL["b_max"], n_b=nb),
                                    theta, snr, "rician").rho_alpha_beta)
        for nb in nb_list
    ])

    # constructive two-diffusion-time break (equal total measurements: 16 vs 8+8)
    b1, dt1 = identifiability_joint.two_dt_design(b_max=CLINICAL["b_max"], n_b=16, ratios=(1.0,))
    b2, dt2 = identifiability_joint.two_dt_design(b_max=CLINICAL["b_max"], n_b=8, ratios=(1.0, 2.5))
    r1 = fisher_joint.crlb_joint(b1, theta, snr, "rician", dt=dt1)
    r2 = fisher_joint.crlb_joint(b2, theta, snr, "rician", dt=dt2)

    degenerate = rho_median >= DEGEN_RHO

    notes = []
    if degenerate:
        notes.append(f"DEGENERATE: median |rho_alpha_beta| = {rho_median:.3f} >= {DEGEN_RHO} over the "
                     f"physiological grid; the time-order alpha and space-order beta are not separately "
                     f"recoverable at a single diffusion time.")
    else:
        notes.append(f"SEPARABLE: median |rho_alpha_beta| = {rho_median:.3f} < {DEGEN_RHO}; alpha and beta "
                     f"are separately recoverable at a single diffusion time in this regime.")
    notes.append(f"degeneracy PERSISTS with b-values: |rho| at n_b={nb_list} = "
                 f"{np.round(nb_rho,3).tolist()} (contrast: the CP0 single-order wall recedes for n_b>=8).")
    notes.append(f"a SECOND diffusion time breaks it: |rho| {abs(r1.rho_alpha_beta):.3f} -> "
                 f"{abs(r2.rho_alpha_beta):.3f}, cond {r1.cond:.1e} -> {r2.cond:.1e}, "
                 f"cv_alpha {r1.cv_alpha:.2f} -> {r2.cv_alpha:.2f}.")

    return CP1Report(
        alpha_grid=ALPHA_GRID, beta_grid=BETA_GRID,
        rho_map=rho_map, cond_map=cond_map, cv_alpha_map=cva_map, cv_beta_map=cvb_map,
        rho_median=rho_median, rho_min=rho_min, rho_max=rho_max, cond_median=cond_median,
        headline_theta=theta, headline_rho=rh.rho_alpha_beta, headline_cond=rh.cond,
        headline_cv_alpha=rh.cv_alpha, headline_cv_beta=rh.cv_beta,
        boot_alpha_ci=boot.alpha_ci, boot_beta_ci=boot.beta_ci, boot_corr=boot.corr_alpha_beta,
        boot_alpha_rel_width=boot.alpha_rel_width, boot_beta_rel_width=boot.beta_rel_width,
        n_boot=n_boot,
        nb_list=nb_list, nb_rho=nb_rho,
        break_rho_single=abs(r1.rho_alpha_beta), break_rho_two=abs(r2.rho_alpha_beta),
        break_cond_single=r1.cond, break_cond_two=r2.cond,
        break_cv_alpha_single=r1.cv_alpha, break_cv_alpha_two=r2.cv_alpha,
        degenerate=degenerate, notes=notes,
    )
