# Roadmap — deferred work

Caliper ships only the **un-gated** calibration tooling: the model-agnostic
ruler (`caliper.metrics`), the MAF posterior estimator (`caliper.estimator_maf`),
the split-conformal / CQR wrapper (`caliper.conformal`), the synthetic data
generator (`caliper.forward`), tests, CI, and a one-command demo.

## Deliberately NOT implemented here

The following are **deferred** and gated on a separate, unpublished paper (the
Minos work). **None of it is implemented in this repository** — this note exists
only to record the boundary:

- **Value-of-information (VoI)** scoring — quantifying the expected utility gain
  from additional acquisition or a tighter posterior.
- **Decision-gap** analysis — the gap between calibrated-uncertainty-aware
  decisions and an oracle.
- **Deployment validity monitor** — online detection of calibration drift /
  distribution shift in deployment.
- A **citable, JOSS-style packaged release** with the above modules.

These are intentionally out of scope. Do not add them to this repository until
the gating paper is published and the work is explicitly authorized.

## Possible un-gated extensions (safe to consider)

These stay within the current scope and could be added without the gated paper:

- Additional estimators behind the same `predict_quantiles` contract
  (e.g. Bayesian least-squares, dropout ensembles) for ruler comparisons.
- Mondrian / group-conditional conformal to *partially* mitigate the D\*
  conditional-coverage gap (with honest reporting that the identifiability
  limit remains).
- More noise models (e.g. non-central chi for multi-coil) in `caliper.forward`.
