"""Procrustes -- misspecification-aliasing of a calibrated error bar.

Thesis (scoped & pre-registered, see POSITIONING.md): fitting a bi-exponential
IVIM model on NON-bi-exponential truth leaves distribution-free *marginal*
coverage intact, yet breaks the *conditional* coverage of the well-identified
tissue-diffusion map D along the latent misspecification axis -- a failure that
survives inside the well-identified D* regime and is therefore mechanistically
distinct from Gauge's within-model identifiability wall.  Heavy-tailed perfusion
departures (stretched-exp) alias into the high-b tissue slope and do the damage;
a faster tri-exp third pool does not (null); pure dispersion is a hidden channel.

Throughline: misspecification is invisible from inside the model -- another
observable/hidden instance.  Ground truth is the Lattice DRO; fully synthetic.
"""
from __future__ import annotations

from .config import (
    ProcrustesConfig,
    DepartureFamily,
    FAMILIES,
    STRETCHED,
    LOGNORMAL,
    TRIEXP,
)
from .seeding import GLOBAL_SEED, make_rng
from .separation import SeparationResult, run_separation, separation_detail
from .diagnostic import diagnose, auc, spearman, STATS
from .bootstrap import across_seed_ci, cluster_bootstrap_gap

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "ProcrustesConfig",
    "DepartureFamily",
    "FAMILIES",
    "STRETCHED",
    "LOGNORMAL",
    "TRIEXP",
    "GLOBAL_SEED",
    "make_rng",
    "SeparationResult",
    "run_separation",
    "separation_detail",
    "diagnose",
    "auc",
    "spearman",
    "STATS",
    "across_seed_ci",
    "cluster_bootstrap_gap",
]
