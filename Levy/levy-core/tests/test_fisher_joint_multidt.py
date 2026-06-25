"""A second diffusion time breaks the (alpha,beta) degeneracy ANALYTICALLY: at equal total
measurement count, the two-Delta FIM is far better conditioned and |rho_alpha_beta| drops.
"""
import numpy as np

from levy import fisher_joint, identifiability_joint as ij

THETA = np.array([1.0, 1.5e-3, 0.80, 1.7])


def test_dt_argument_defaults_to_single_dt():
    # crlb_joint(b, ...) and crlb_joint(b, ..., dt=ones) must agree.
    from levy import wall
    b = wall.default_b_design(2500.0, 8)
    r1 = fisher_joint.crlb_joint(b, THETA, snr=40.0)
    r2 = fisher_joint.crlb_joint(b, THETA, snr=40.0, dt=np.ones_like(b))
    assert np.isclose(r1.rho_alpha_beta, r2.rho_alpha_beta, rtol=1e-9)


def test_second_diffusion_time_breaks_degeneracy():
    # Same total number of measurements (16): one design uses a single Delta with 16 b-values,
    # the other splits them across two diffusion times. The two-Delta design is better
    # conditioned and less (alpha,beta)-degenerate.
    b1, dt1 = ij.two_dt_design(b_max=2500.0, n_b=16, ratios=(1.0,))          # single Delta, 16 b
    b2, dt2 = ij.two_dt_design(b_max=2500.0, n_b=8, ratios=(1.0, 2.5))        # two Delta, 8+8
    r1 = fisher_joint.crlb_joint(b1, THETA, snr=40.0, dt=dt1)
    r2 = fisher_joint.crlb_joint(b2, THETA, snr=40.0, dt=dt2)
    assert abs(r2.rho_alpha_beta) < abs(r1.rho_alpha_beta), \
        f"two-dt rho={r2.rho_alpha_beta} vs single-dt rho={r1.rho_alpha_beta}"
    assert r2.cond < r1.cond, f"two-dt cond={r2.cond:.2e} vs single={r1.cond:.2e}"
    # and the relative CRLB on alpha improves with the second diffusion time
    assert r2.cv_alpha < r1.cv_alpha
