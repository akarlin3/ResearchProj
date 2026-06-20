"""Calibration ruler: coverage / ECE / sharpness, re-derived independently. [CP2]

These are standard, published quantities (Gneiting & Raftery 2007; quantile-
regression calibration), re-derived here from their definitions -- NOT imported from
``caliper.metrics`` (the forbidden module). For a central interval [lo, hi] at
nominal level ``1 - alpha`` over ground truth ``y``:

* **coverage** = mean( lo <= y <= hi )                    (target: vs nominal)
* **ECE**      = mean over a level grid of |empirical_coverage(level) - level|
* **sharpness**= mean interval width (must accompany coverage -- wide intervals
                 cover cheaply)
* plus pinball loss and the interval score, as proper-scoring cross-checks.

A ``conditional_coverage`` by D*-tercile probe is included to expose the high-D*
under-coverage. Everything is numpy-only and documented to its reference. This file
is the independent re-derivation that lets CP3 score the rebuild without touching
Fashion's or Caliper's ruler.
"""
from __future__ import annotations


def score_quantiles(*args, **kwargs):  # [CP2]
    raise NotImplementedError(
        "Gnomon CP2: independent coverage/ECE/sharpness ruler (numpy-only)."
    )
