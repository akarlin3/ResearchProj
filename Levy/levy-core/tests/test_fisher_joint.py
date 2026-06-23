"""Joint 4-parameter Fisher/CRLB layer for theta = (S0, D, alpha, beta).

The structural-identifiability deliverable: the (alpha, beta) correlation rho_alpha_beta from
the inverse Fisher matrix (-> +-1 signals the time/space orders are degenerate at a single
diffusion time), alongside the alpha-D and beta-D trade-offs and the FIM condition number.
"""
import numpy as np

from levy import fisher, forward, fisher_joint


def test_fim_is_4x4_symmetric_psd():
    b = forward_b()
    fim = fisher_joint.fisher_matrix_joint(b, THETA, snr=40.0, model="rician")
    assert fim.shape == (4, 4)
    assert np.allclose(fim, fim.T, rtol=1e-10)
    w = np.linalg.eigvalsh(fim)
    assert np.all(w >= -1e-9), f"FIM not PSD: eigvals {w}"


def test_rician_approaches_gaussian_at_high_snr():
    b = forward_b()
    fg = fisher_joint.fisher_matrix_joint(b, THETA, snr=2000.0, model="gaussian")
    fr = fisher_joint.fisher_matrix_joint(b, THETA, snr=2000.0, model="rician")
    rel = np.linalg.norm(fg - fr) / np.linalg.norm(fg)
    assert rel < 1e-2, f"Rician FIM should approach Gaussian at high SNR; rel={rel}"


def test_s0_d_block_matches_cp0_at_alpha1_beta2_gaussian():
    # At alpha=1, beta=2 the joint (S0,D) Jacobian columns equal CP0's, so the (S0,D) 2x2
    # sub-block of the joint Gaussian FIM equals the CP0 Gaussian FIM (S0,D) sub-block.
    b = forward_b()
    theta_j = np.array([1.0, 1.5e-3, 1.0, 2.0])
    fj = fisher_joint.fisher_matrix_joint(b, theta_j, snr=50.0, model="gaussian")
    fc = fisher.fisher_matrix(b, np.array([1.0, 1.5e-3, 1.0]), snr=50.0, model="gaussian")
    iS0, iD = forward.IDX_JOINT["S0"], forward.IDX_JOINT["D"]
    jS0, jD = forward.IDX["S0"], forward.IDX["D"]
    assert np.isclose(fj[iS0, iS0], fc[jS0, jS0], rtol=1e-4)
    assert np.isclose(fj[iD, iD], fc[jD, jD], rtol=1e-3)
    assert np.isclose(fj[iS0, iD], fc[jS0, jD], rtol=1e-3)


def test_crlb_joint_diagnostics_present_and_sane():
    b = forward_b()
    r = fisher_joint.crlb_joint(b, THETA, snr=40.0, model="rician")
    # CRLBs positive
    assert r.crlb_alpha > 0 and r.crlb_beta > 0 and r.crlb_D > 0
    # relative CRLBs (cv) positive
    assert r.cv_alpha > 0 and r.cv_beta > 0
    # correlations in [-1, 1]
    for rho in (r.rho_alpha_beta, r.rho_alpha_D, r.rho_beta_D):
        assert -1.0 - 1e-9 <= rho <= 1.0 + 1e-9
    # condition number finite and >= 1
    assert np.isfinite(r.cond) and r.cond >= 1.0


def test_more_b_values_do_not_increase_crlb():
    # Adding b-values cannot reduce information: CRLB_alpha at n_b=16 <= n_b=4 (same b_max).
    from levy import wall
    b4 = wall.default_b_design(b_max=2500.0, n_b=4)
    b16 = wall.default_b_design(b_max=2500.0, n_b=16)
    a4 = fisher_joint.crlb_joint(b4, THETA, snr=40.0).crlb_alpha
    a16 = fisher_joint.crlb_joint(b16, THETA, snr=40.0).crlb_alpha
    assert a16 <= a4 * (1 + 1e-9)


# --- helpers / fixtures -------------------------------------------------------
THETA = np.array([1.0, 1.5e-3, 0.80, 1.7])   # representative physiological CTRW voxel


def forward_b():
    from levy import wall
    return wall.default_b_design(b_max=2500.0, n_b=8)
