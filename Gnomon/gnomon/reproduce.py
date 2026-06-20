"""CP3 reproduction driver: run the rebuild, compare to the frozen manifest. [CP3]

The one-command verdict. It runs the clean rebuild end-to-end --

    1. synthetic headline cohorts (Lattice)  -> NLLS, Laplace, MCMC, MAF
    2. real OSIPI abdomen ROI                 -> NLLS railing rate
    3. score everything with gnomon.metrics; bootstrap-CI every load-bearing number

-- then compares each result to the pinned target in :mod:`gnomon.manifest` using
the **frozen** tolerances, and emits:

    * a per-target table: claimed | rebuilt | CI | PASS/FAIL
    * the overall verdict:
        - REPRODUCES  -> Fashion's numbers are real; rejection was presentational;
                         Gnomon becomes the clean reference + complete methods.
        - DOES NOT REPRODUCE -> a precise divergence report (which numbers, by how
                         much, likely cause). This is a HARD HALT (CP3).

Results are written to ``results/reproduction.json`` (committed) so the verdict is
auditable. No number in this file is hard-coded; all come from the live run.
"""
from __future__ import annotations


def run(*args, **kwargs):  # [CP3]
    raise NotImplementedError("Gnomon CP3: reproduce-or-refute driver + verdict.")


if __name__ == "__main__":  # pragma: no cover
    run()
