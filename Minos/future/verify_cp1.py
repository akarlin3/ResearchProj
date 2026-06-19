"""CP1 gate -- the theory half reproduces when driven from ``future/`` by read-only import.

Three independent checks, all green => CP1 gate passes:

  (1) WIRING       every dependency path (minos-core, theory, Fashion, Gauge) resolves and
                   imports, through ``future/_paths.py`` only -- nothing copied or edited.

  (2) LIVE IMPORT  wire minos-core, import the validated v1-v3 package, and *recompute* the
                   canonical decision-calibration gap at the default cell
                   (kappa=3, lambda=3, rho=0.5). Asserts the published v2/v3 sign structure:
                   tau_stat < 1 < tau_star, so G = tau_star - tau_stat > 0. (tau_star is a
                   shallow optimum and noisy at finite N; the tight numeric agreement with
                   the leading-order theory value 1.0485 is GATE 3's job, re-run in check 3.)

  (3) GATE REPLAY  run the three *unmodified* theory gate scripts as subprocesses and assert
                   each prints its 'GATE n ... PASS' line. This is the read-only reproduction
                   of the whole theory half from future/.

Run:  <proteus python> Minos/future/verify_cp1.py        # add --full for full-N GATE 3
Exit 0 = CP1 gate green; nonzero = a check failed (the script names which).
"""
from __future__ import annotations

import os
import subprocess
import sys
import time

import _paths  # noqa: E402  (same directory)

FAST = "--full" not in sys.argv
DEFAULT_CELL = dict(kappa=3.0, lam=3.0, rho=0.5)
N_LIVE = 300_000
LIVE_SEEDS = (11, 12, 13)
THEORY_VALUE_TAUSTAR = 1.0485  # leading-order Theorem 1, default cell (documentation anchor)


def _hr(title: str) -> None:
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78)


def check_wiring() -> None:
    _hr("CP1 check 1/3 -- dependency wiring (read-only, via future/_paths.py)")
    paths = _paths.add_all()
    for k, v in paths.items():
        print(f"  {k:<12} -> {v}")
    # Each package must import from where _paths put it.
    import minos  # noqa: F401
    from minos import calibration, decision, monitor, generative, correction, voi  # noqa: F401
    print(f"  import minos      OK  ({minos.__file__})")
    import uq  # noqa: F401  (Fashion -- PROVISIONAL)
    from uq import calib as _fc  # noqa: F401
    print(f"  import uq (Fashion) OK  -- PROVISIONAL  ({uq.__file__})")
    import gauge  # noqa: F401  (Gauge -- PROVISIONAL)
    from gauge import conformal as _gc, monitor as _gm  # noqa: F401
    assert hasattr(_gc, "split_conformal") and hasattr(_gm, "DeploymentMonitor")
    print(f"  import gauge        OK  -- PROVISIONAL  ({gauge.__file__})")
    print("  WIRING: PASS")


def check_live_import() -> None:
    _hr("CP1 check 2/3 -- live recompute of the canonical gap (imported minos package)")
    _paths.add_minos_core()
    from minos.config import MinosConfig, gaussian_latent_config
    from minos.calibration import gap
    from minos.generative import make_population
    from minos.seeding import make_rng

    cfg = gaussian_latent_config(
        rho=DEFAULT_CELL["rho"], kappa=DEFAULT_CELL["kappa"],
        lam=DEFAULT_CELL["lam"], base=MinosConfig(n_voxels=N_LIVE),
    )
    ts, tstar, g = [], [], []
    t0 = time.time()
    for sd in LIVE_SEEDS:
        base = make_population(cfg.replace(seed=sd), make_rng(sd))
        gr = gap(base, cfg)
        ts.append(gr.tau_stat); tstar.append(gr.tau_star); g.append(gr.gap)
    import statistics as st
    ts_m, tstar_m, g_m = st.mean(ts), st.mean(tstar), st.mean(g)
    print(f"  default cell kappa={DEFAULT_CELL['kappa']} lambda={DEFAULT_CELL['lam']} "
          f"rho={DEFAULT_CELL['rho']};  N={N_LIVE} x {len(LIVE_SEEDS)} seeds  "
          f"({time.time()-t0:.1f}s)")
    print(f"  tau_stat = {ts_m:.4f}   tau_star = {tstar_m:.4f}   G = {g_m:+.4f}")
    print(f"  (leading-order theory tau_star = {THEORY_VALUE_TAUSTAR} -- exact agreement is GATE 3)")
    # The structural prediction of Theorem 1: skew + asymmetric cost => positive gap.
    assert tstar_m > ts_m, f"expected tau_star > tau_stat, got {tstar_m} <= {ts_m}"
    assert g_m > 0.0, f"expected positive gap, got G={g_m}"
    assert ts_m < 1.0, f"expected tau_stat < 1 (coverage shrink), got {ts_m}"
    assert tstar_m > 1.0, f"expected tau_star > 1 (decision widening), got {tstar_m}"
    print("  LIVE IMPORT: PASS  (imported package reproduces the v2/v3 gap structure)")


def check_gate_replay() -> None:
    _hr(f"CP1 check 3/3 -- theory gate replay (subprocess, {'FAST' if FAST else 'FULL'})")
    theory_dir = str(_paths.THEORY)
    scripts = ["gap_scaling.py", "detectability.py", "confirm.py"]
    env = dict(os.environ)
    if FAST:
        env["MINOS_FAST"] = "1"
    for s in scripts:
        path = os.path.join(theory_dir, s)
        t0 = time.time()
        r = subprocess.run([sys.executable, path], capture_output=True, text=True, env=env)
        out = r.stdout + r.stderr
        passed = (r.returncode == 0) and ("PASS" in out)
        tag = "PASS" if passed else "FAIL"
        print(f"  {s:<22} -> {tag}  (rc={r.returncode}, {time.time()-t0:.1f}s)")
        if not passed:
            print("  ---- captured tail ----")
            print("\n".join(out.splitlines()[-15:]))
            raise SystemExit(f"GATE REPLAY FAILED on {s}")
    print("  GATE REPLAY: PASS  (all three theory gates reproduce from future/)")


def main() -> int:
    print("CP1 verification -- Minos future/ wired to the theory half (read-only)")
    check_wiring()
    check_live_import()
    check_gate_replay()
    _hr("CP1 GATE: PASS")
    print("  theory half reproduces via read-only import from future/;")
    print("  Fashion/Gauge wiring resolves (PROVISIONAL, exercised at CP2/CP3).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
