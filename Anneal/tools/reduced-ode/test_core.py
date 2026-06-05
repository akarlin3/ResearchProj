"""Self-tests for reduced_core: the 3-D RHS must collapse to the 2-D RHS on the
rho1=1 manifold, and the Eqs (13)-(14) fixed points must zero the 2-D RHS.
Run: python3 tools/reduced-ode/test_core.py  (exit 0 = pass)."""
import math
import sys

import numpy as np

import reduced_core as rc


def test_3d_reduces_to_2d():
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(2000):
        A = rng.uniform(0.05, 0.6)
        beta = rng.uniform(0.0, 0.25)
        r = rng.uniform(0.05, 0.99)
        psi = rng.uniform(-math.pi, math.pi)
        p = rc.Params(A=A, beta=beta)
        d2 = rc.rhs_2d([r, psi], p)
        d3 = rc.rhs_3d([1.0, r, psi], p)  # rho1=1 manifold
        # 3-D returns (rho1dot, rho2dot, psidot); compare rho2dot<->rdot, psidot<->psidot
        err = max(abs(d3[1] - d2[0]), abs(d3[2] - d2[1]), abs(d3[0]))
        worst = max(worst, err)
    print(f"[1] 3D->2D max |residual| over 2000 random pts: {worst:.3e}")
    assert worst < 1e-9, "3-D RHS does not reduce to 2-D on rho1=1"


def test_fixed_points_zero_rhs():
    """Eqs (13)-(14) family must give rdot=psidot=0 in the 2-D system."""
    worst = 0.0
    for beta in (0.05, 0.1, 0.15):
        psi_grid = np.linspace(-1.4, 1.4, 4001)
        psi, r, A = rc.fixed_point_family(beta, psi_grid)
        for pj, rj, Aj in zip(psi, r, A):
            if not (0.0 < Aj < 0.8):
                continue
            p = rc.Params(A=Aj, beta=beta)
            d = rc.rhs_2d([rj, pj], p)
            worst = max(worst, float(np.linalg.norm(d)))
    print(f"[2] Eqs(13-14) fixed-point max |RHS|: {worst:.3e}")
    assert worst < 1e-9, "stationary-chimera relations do not zero the 2-D RHS"


def test_sync_is_fixed_point():
    """Full sync rho1=rho2=1, psi=0 must be a fixed point of the 3-D system."""
    worst = 0.0
    for A in (0.2, 0.5):
        for beta in (0.05, 0.1):
            d = rc.rhs_3d([1.0, 1.0, 0.0], rc.Params(A=A, beta=beta))
            worst = max(worst, float(np.linalg.norm(d)))
    print(f"[3] sync state (1,1,0) max |RHS|: {worst:.3e}")
    assert worst < 1e-12


if __name__ == "__main__":
    test_3d_reduces_to_2d()
    test_fixed_points_zero_rhs()
    test_sync_is_fixed_point()
    print("ALL CORE TESTS PASS")
    sys.exit(0)
