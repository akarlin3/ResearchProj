"""CP1 degeneracy report: structure, scoping, and the constructive two-diffusion-time break."""
import numpy as np

from levy import degeneracy, seeding


def test_analytic_map_shapes_and_bounds():
    b, rho, cond, cva, cvb = degeneracy.analytic_map()
    assert rho.shape == (len(degeneracy.ALPHA_GRID), len(degeneracy.BETA_GRID))
    assert np.all(rho >= 0.0) and np.all(rho <= 1.0 + 1e-9)
    assert np.all(cond >= 1.0)
    assert np.all(cva > 0) and np.all(cvb > 0)


def test_report_verdict_consistent_with_threshold():
    rep = degeneracy.cp1_report(rng=seeding.make_rng(11), n_boot=12)
    assert rep.degenerate == (rep.rho_median >= degeneracy.DEGEN_RHO)
    assert rep.rho_min <= rep.rho_median <= rep.rho_max


def test_second_diffusion_time_reduces_degeneracy_in_report():
    rep = degeneracy.cp1_report(rng=seeding.make_rng(12), n_boot=12)
    # the constructive boundary: two diffusion times lower |rho|, condition number, and cv_alpha
    assert rep.break_rho_two < rep.break_rho_single
    assert rep.break_cond_two < rep.break_cond_single
    assert rep.break_cv_alpha_two < rep.break_cv_alpha_single


def test_bootstrap_cis_present():
    rep = degeneracy.cp1_report(rng=seeding.make_rng(13), n_boot=16)
    assert rep.boot_alpha_ci[0] <= rep.boot_alpha_ci[1]
    assert rep.boot_beta_ci[0] <= rep.boot_beta_ci[1]
    assert -1.0 - 1e-9 <= rep.boot_corr <= 1.0 + 1e-9
