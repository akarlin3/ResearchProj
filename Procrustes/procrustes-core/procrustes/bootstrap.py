"""Bootstrap confidence intervals for the load-bearing separation gaps.

The headline numbers (conditional gap, well-identified-D* gap, diagnostic AUC)
are reported with bootstrap CIs so the claims are falsifiable, not point
assertions. Two estimators are provided:

  * ``across_seed_ci`` -- normal-approximation CI over the per-seed point
    estimates (the house convention; captures cohort + split + noise variability
    between seeds). Primary for >= 8 seeds.
  * ``cluster_bootstrap_gap`` -- a two-level (seed, voxel) cluster bootstrap of a
    coverage gap on a chosen D*-regime mask: resample seeds with replacement,
    then resample test voxels within each sampled seed (paired across the
    placebo and worst-departure knobs), recompute the across-seed mean gap. This
    captures BOTH the between-seed and the within-seed sampling uncertainty.

Both are seeded (no wall-clock) and therefore reproducible.
"""
from __future__ import annotations

import numpy as np

from .seeding import GLOBAL_SEED


def across_seed_ci(vals) -> tuple:
    """(mean, lo, hi) normal-approx 95% CI over per-seed point estimates."""
    a = np.asarray(vals, float)
    m = float(a.mean())
    se = float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0
    return m, m - 1.96 * se, m + 1.96 * se


def _gap_on_mask(err_limit, err_worst, radius, mask):
    """Coverage gap (placebo - worst) on a boolean voxel mask, fixed radius."""
    m = mask
    if not m.any():
        return np.nan
    cov_lim = float((err_limit[m] <= radius).mean())
    cov_wrst = float((err_worst[m] <= radius).mean())
    return cov_lim - cov_wrst


def cluster_bootstrap_gap(details, mask_key: str, n_boot: int = 2000,
                          seed: int = GLOBAL_SEED) -> dict:
    """Two-level cluster bootstrap of a coverage gap across seeds.

    ``details`` is a list of dicts from ``separation.separation_detail`` (one per
    seed). ``mask_key`` in {"all", "wellid", "lo", "hi"} selects the D*-regime
    subset. Returns the point estimate (across-seed mean) plus the percentile
    95% CI from ``n_boot`` cluster-bootstrap replicates.
    """
    rng = np.random.default_rng(seed)
    n_seeds = len(details)

    point = float(np.nanmean([
        _gap_on_mask(d["err_limit"], d["err_worst"], d["radius"], d[mask_key])
        for d in details]))

    stats = np.empty(n_boot)
    for b in range(n_boot):
        sd_idx = rng.integers(0, n_seeds, size=n_seeds)   # resample seeds
        seed_gaps = []
        for si in sd_idx:
            d = details[si]
            n = d["n_test"]
            vi = rng.integers(0, n, size=n)               # resample voxels (paired)
            mask = d[mask_key][vi]
            if not mask.any():
                continue
            el = d["err_limit"][vi][mask]
            ew = d["err_worst"][vi][mask]
            r = d["radius"]
            seed_gaps.append(float((el <= r).mean()) - float((ew <= r).mean()))
        if seed_gaps:
            stats[b] = np.mean(seed_gaps)
        else:
            stats[b] = np.nan
    lo, hi = np.nanpercentile(stats, [2.5, 97.5])
    return {"point": point, "lo": float(lo), "hi": float(hi),
            "mask": mask_key, "n_boot": n_boot, "n_seeds": n_seeds}
