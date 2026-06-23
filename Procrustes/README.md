# Procrustes

*Misspecification-aliasing of a calibrated error bar* — fitting a bi-exponential
IVIM model on non-bi-exponential truth keeps **marginal** coverage but breaks the
**conditional** coverage of the *well-identified* tissue-diffusion map D, distinct
from Gauge's within-model identifiability wall.

- [`procrustes-core/`](procrustes-core/) — the clean-room core: CP0 separation,
  boundary gates, and the observable misspecification diagnostic. Ground truth is
  the Lattice DRO (seed-generated, no data files).
- [`procrustes-core/POSITIONING.md`](procrustes-core/POSITIONING.md) — the novelty
  gate record (vs Gauge, Lei 2018, Barber 2021, Wang–Tamir–Bush 2026, IVIM model
  selection).

Status: **scaffold (CP0)**. Venue `[confirm venue]`.
