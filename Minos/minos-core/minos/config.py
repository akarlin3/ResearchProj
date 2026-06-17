"""Frozen configuration for the Minos-Core synthetic decision model.

All numbers that define the toy live here. The dataclass is frozen so a run cannot
silently mutate its own parameters mid-experiment.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

from .seeding import GLOBAL_SEED


@dataclass(frozen=True)
class MinosConfig:
    # --- utility (Section 1) ---
    t1: float = 0.0          # spare | treat threshold
    t2: float = 2.0          # treat | escalate threshold
    k_under: float = 2.0     # under-treatment slope (penalised more)
    k_over: float = 1.0      # over-treatment slope

    # --- latent prior (Section 2) ---
    # latent_mode="mixture" (v1 default): 3-component Gaussian mixture, symmetric about
    #   the decision midpoint (t1+t2)/2 = 1.0 so the calibrated reported posterior
    #   (tau=1) is the decision-optimal error bar (DESIGN 6.2). This is BASELINE_V1.
    # latent_mode="gaussian" (v2 gap sweep): single Gaussian N(theta_mean, theta_std^2),
    #   used so the threshold-proximity axis rho = dist(E[theta], nearest threshold)/s
    #   is well-defined (DESIGN_B Section 4).
    latent_mode: str = "mixture"
    mix_weights: Tuple[float, ...] = (0.35, 0.30, 0.35)
    mix_means: Tuple[float, ...] = (-1.0, 1.0, 3.0)
    mix_stds: Tuple[float, ...] = (1.0, 1.0, 1.0)
    theta_mean: float = 1.0  # gaussian-latent mean (latent_mode="gaussian")
    theta_std: float = 1.0   # gaussian-latent std  (latent_mode="gaussian")

    # --- measurement (Section 2) ---
    s: float = 0.5           # intrinsic measurement spread
    tau: float = 1.0         # calibration knob (1 = calibrated)

    # --- misspecification (DESIGN_B Section 2) ---
    # The reported posterior is the Gaussian moment-match N(mu, (tau*s)^2) of a SKEWED
    # true error d = mu - theta. kappa is the skew-normal shape; kappa=0 -> d ~ N(0,s^2)
    # exactly -> v1 (well-specified). kappa>0 -> heavy under-treatment tail.
    kappa: float = 0.0       # report-error skew shape (0 = Gaussian = v1)
    family: str = "skewnorm" # misspecification family; "gaussian" forces kappa-free N(0,s^2)
    # posterior_centric=False (v1): theta is the prior draw, mu = theta + error (forward
    #   measurement; theta|mu carries prior shrinkage). True (v2 gap): the latent draw is
    #   the report centre mu and the truth is theta = mu + s*u, so E[theta|mu]=mu exactly
    #   and tau*=1 at kappa=0 for ANY latent/utility -- the gap is then pure skew effect.
    posterior_centric: bool = False

    # --- shift (Section 3) ---
    delta: float = 0.0       # distribution-shift magnitude
    alpha: float = 0.5       # noise-inflation gain (overconfidence under shift)
    beta: float = 5.0        # downward-bias gain (systematic under-treatment)

    # --- trust-gate (Section 5) ---
    w_train_mean: float = 0.0
    w_train_std: float = 1.0
    q_gate: float = 0.995    # in-distribution quantile -> threshold (0.5% nominal FPR)

    # --- estimation ---
    n_voxels: int = 200_000
    seed: int = GLOBAL_SEED

    def __post_init__(self) -> None:
        assert self.t1 < self.t2, "need t1 < t2"
        # v2 relaxes the v1 strict ``>`` to ``>=`` so the symmetric-utility slice lambda=1
        # (k_under == k_over) is a legal config; the v1 default (2, 1) still satisfies it.
        assert self.k_under >= self.k_over > 0, "under-treatment must cost at least as much"
        assert self.s > 0 and self.tau > 0
        assert len(self.mix_weights) == len(self.mix_means) == len(self.mix_stds)
        assert abs(sum(self.mix_weights) - 1.0) < 1e-9, "mixture weights must sum to 1"
        assert 0.0 < self.q_gate < 1.0
        assert self.n_voxels > 0
        assert self.latent_mode in ("mixture", "gaussian")
        assert self.family in ("skewnorm", "gaussian")
        assert self.kappa >= 0.0 and self.theta_std > 0

    def replace(self, **kw) -> "MinosConfig":
        """Return a copy with overrides (frozen dataclasses need this helper)."""
        from dataclasses import replace as _replace

        return _replace(self, **kw)


DEFAULT = MinosConfig()

#: v1 symmetric, well-specified baseline (mixture latent, Gaussian error). Retained
#: verbatim so the v1 result is reproducible inside v2 as the degenerate ``G ~ 0`` case.
BASELINE_V1 = MinosConfig()


#: Spare/treat threshold offset for the gap configs: pushed far below ``t2`` so the gap
#: experiment is a clean *treat vs escalate* decision around the single active threshold
#: ``t2`` (spare never fires), leaving room for the threshold-proximity (``rho``) sweep.
GAP_T1_OFFSET = 20.0


def gaussian_latent_config(rho: float, kappa: float = 0.0, lam: float = 2.0,
                           theta_std: float = 0.5, base: "MinosConfig" = DEFAULT,
                           **kw) -> "MinosConfig":
    """A single-Gaussian-latent config for the gap sweep at a target ``(rho, kappa, lam)``.

    The report centres are ``mu ~ N(theta_mean, theta_std^2)`` with ``theta_mean = t2 -
    rho*s``, so ``rho`` is the distance from ``E[mu]`` to the (single active) upper
    threshold ``t2`` in units of ``s``. ``lam`` is the under:over cost ratio, realised as
    ``k_under = lam`` with ``k_over = 1``. The spare threshold is placed ``GAP_T1_OFFSET``
    below ``t2`` so only the treat/escalate boundary is active (DESIGN_B Section 4).
    """
    theta_mean = base.t2 - rho * base.s
    return base.replace(latent_mode="gaussian", posterior_centric=True,
                        t1=base.t2 - GAP_T1_OFFSET, theta_mean=theta_mean,
                        theta_std=theta_std, kappa=kappa, k_under=lam, k_over=1.0, **kw)


#: Default *misspecified* gap config (strong skew, asymmetric utility, near the threshold).
DEFAULT_MISSPEC = gaussian_latent_config(rho=0.5, kappa=3.0, lam=3.0)
