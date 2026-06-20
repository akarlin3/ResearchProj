"""MAF amortized posterior (neural posterior estimation).  [CP2 — spec below]

Clean-room masked-autoregressive-flow posterior q(theta | signal), trained by
neural posterior estimation on simulated (theta, signal) pairs from
:mod:`gnomon.forward` over the manifest priors. This is the "amortized posterior"
arm: one network gives a per-voxel posterior at inference time, from which Gnomon
reads quantile intervals to score coverage/ECE/sharpness.

Target T4 (directional): on a held-out synthetic set the flow has lower D* ECE and
is sharper than the boundary-railed NLLS baseline, at coverage >= the baseline --
because a flow captures the skewed D* posterior shape that a Gaussian cannot.

CP2 documents the FULL training spec (architecture/transforms, simulation budget,
prior, optimizer, epochs/early-stopping, seed) -- the training completeness Fashion
was flagged for. Requires the optional ``[flow]`` extra (torch); the rest of Gnomon
runs without it.
"""
from __future__ import annotations


class MAFPosterior:  # [CP2]
    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "Gnomon CP2: clean-room MAF/NPE posterior (needs the [flow] extra)."
        )
