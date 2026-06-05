"""Phase-portrait utility for the 2-D invariant-manifold system (Eq. 12),
paper Fig. 3 style: vector field, fixed points (stable chimera spiral + saddle),
the saddle's stable/unstable manifolds, and representative trajectories from our
canonical-seed neighbourhood. Used by CP2 and the report. PNG+PDF, deterministic.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import reduced_core as rc


def _manifolds(saddle_state, p, t=400.0, eps=1e-4, npts=4000):
    """Stable (backward) and unstable (forward) manifolds of the saddle."""
    J = rc.jacobian_2d(saddle_state, p)
    ev, evec = np.linalg.eig(J)
    branches = {"unstable": [], "stable": []}
    for k in range(2):
        v = np.real(evec[:, k])
        v = v / np.linalg.norm(v)
        lam = np.real(ev[k])
        sign_t = +1.0 if lam > 0 else -1.0  # forward for unstable, backward for stable
        key = "unstable" if lam > 0 else "stable"
        for s in (+1.0, -1.0):
            st0 = np.array(saddle_state) + s * eps * v
            te = np.linspace(0, sign_t * t, npts)
            sol = rc.solve_ivp_branch(st0, p, sign_t * t, npts)
            branches[key].append(sol)
    return branches


def make_portrait(A, beta, seed_ics, out_base, cfg, title=None, traj_t=400.0):
    p = rc.Params(A=A, beta=beta)
    fig, ax = plt.subplots(figsize=(7.2, 5.4))

    # vector field (direction only)
    psis = np.linspace(-1.2, 1.2, 34)
    rs = np.linspace(0.02, 0.999, 30)
    PS, RS = np.meshgrid(psis, rs)
    U = np.zeros_like(PS)
    V = np.zeros_like(RS)
    for i in range(PS.shape[0]):
        for j in range(PS.shape[1]):
            d = rc.rhs_2d([RS[i, j], PS[i, j]], p)  # (rdot, psidot)
            V[i, j], U[i, j] = d[0], d[1]
    mag = np.hypot(U, V) + 1e-12
    ax.quiver(PS, RS, U / mag, V / mag, mag, cmap="Greys", alpha=0.45,
              scale=38, width=0.0022, pivot="mid")

    # sync line r=1
    ax.axhline(1.0, color="#888", lw=1.0, ls=":", label="sync (r=1)")

    # fixed points
    chim = rc.get_fixed_point(A, beta, branch="chimera")
    sad = rc.get_fixed_point(A, beta, branch="saddle")
    if chim is not None:
        r_c, psi_c, st_c = chim
        col = "#1a9850" if st_c["stable"] else "#d73027"
        lab = "chimera FP (stable)" if st_c["stable"] else "chimera FP (unstable spiral)"
        ax.plot(psi_c, r_c, "o", ms=11, mfc=col, mec="k", mew=1.2, label=lab, zorder=6)
    if sad is not None:
        r_s, psi_s, st_s = sad
        ax.plot(psi_s, r_s, "X", ms=12, mfc="#4575b4", mec="k", mew=1.2,
                label="saddle", zorder=6)
        # manifolds
        try:
            br = _manifolds([r_s, psi_s], p, t=traj_t)
            for k, (key, color, ls) in enumerate(
                [("unstable", "#d73027", "-"), ("stable", "#4575b4", "--")]
            ):
                for m, seg in enumerate(br[key]):
                    ax.plot(seg[:, 1], seg[:, 0], color=color, ls=ls, lw=1.4,
                            alpha=0.9, zorder=4,
                            label=(f"{key} manifold" if m == 0 else None))
        except Exception as e:  # pragma: no cover
            print("manifold warn:", e)

    # representative trajectories from canonical-seed neighbourhood
    for k, (r0, psi0) in enumerate(seed_ics):
        sol = rc.integrate_2d([min(0.999, r0), psi0], p, traj_t, cfg,
                              t_eval=np.linspace(0, traj_t, 3000))
        ax.plot(sol.y[1], sol.y[0], color="#542788", lw=1.0, alpha=0.65,
                zorder=5, label=("canonical-seed trajectory" if k == 0 else None))
        ax.plot(psi0, r0, ".", color="#542788", ms=7, zorder=5)

    ax.set_xlabel(r"$\psi = \phi_1-\phi_2$")
    ax.set_ylabel(r"$r=\rho_2$ (incoherent-population coherence)")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(0.0, 1.02)
    ax.set_title(title or f"2-D phase portrait  A={A}, β={beta}")
    ax.legend(loc="lower left", fontsize=7, framealpha=0.9, ncol=2)
    fig.tight_layout()
    fig.savefig(out_base + ".png", dpi=140)
    fig.savefig(out_base + ".pdf")
    plt.close(fig)
    return {
        "A": A,
        "beta": beta,
        "chimera_fp": None if chim is None else {"r": chim[0], "psi": chim[1],
                                                 "stable": chim[2]["stable"]},
        "saddle_fp": None if sad is None else {"r": sad[0], "psi": sad[1]},
    }
