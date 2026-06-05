# paper-figures — Chaos (AIP) manuscript Figures 1 & 7

Produces two publication figures by **composing existing merged work** — the
traced runs from the absorption re-campaign (`absorption_results/`) and the
reduced-ODE bifurcation machinery (`tools/reduced-ode/`). New code only here;
all prior `*_results/` are **read-only**; everything is deterministic and
regenerable from the committed `figures.config.json`.

A shared style module (`style.py`) holds the rcParams and one colorblind-safe
palette (Okabe–Ito / Wong) so the whole figure set has a consistent look and can
be regenerated through it.

## Figures

| Fig | What | Outputs |
| --- | --- | --- |
| **1** | Annotated example trace of the A=0.5 collapse: `min(R₁,R₂)` with the graze/recovery thresholds; T_b breathing bracket; the first long graze (`t_graze`, engineering "collapse") with its ≥5 s supra-θ excursion and sub-0.80 reform dip; the final no-recovery crossing (`t_abs`, absorption); and ticks where the shipped 2-s supervisor would fire. | `paper_figures/fig1.{pdf,png}`, `fig1_caption.txt` |
| **7** | (β, A) stability diagram (Abrams 2008 Fig. 4 style, from our own machinery): saddle-node (Eq. 17) and Hopf (Eq. 18) series curves + numeric check marks; homoclinic curve by escape-to-sync bisection; the four regions; the Takens–Bogdanov point; and both operating corners. | `paper_figures/fig7.{pdf,png}`, `fig7_caption.txt` |

## Run

```
pip install numpy scipy matplotlib   # if not present
python3 tools/paper-figures/run_all.py            # full pipeline (~5 min; homoclinic dominates)
python3 tools/paper-figures/run_all.py --skip-curves   # reuse cached fig7_curves.json
```

Step by step: `cp0_select.py` (audit + pin-checked Fig.1 run selection) →
`node fig1_trace.mjs > paper_figures/fig1_trace.json` (deterministic trace) →
`fig1.py` → `fig7_curves.py` (SN/Hopf/homoclinic, ~4–5 min) → `fig7.py` →
`make_report.py`.

## Files

| File | Role |
| --- | --- |
| `figures.config.json` | Single committed config: palette, Fig.1 run (pinned seed) + criterion, Fig.7 β grids + thresholds + corners. |
| `style.py` | Shared rcParams + palette + `savefig` (PDF+PNG @ 300 dpi). |
| `cp0_select.py` | CP0 light audit; deterministic Fig.1 candidate scoring; asserts the pinned seed. |
| `fig1_trace.mjs` | Regenerates the chosen run (R₁, R₂, min, sync) from the seed via the shipped-identical RK4; drives the absorption labeler; replays the supervisor detector. |
| `fig1.py` | Plots Figure 1 + caption. |
| `fig7_curves.py` | Computes & caches the bifurcation curves; β=0.05 homoclinic self-check (must reproduce 0.4096). |
| `fig7.py` | Plots Figure 7 + caption. |
| `make_report.py` | `paper_figures/FIGURES_REPORT.md`. |
| `run_all.py` | Driver. |

## Provenance / determinism

- `fig1_trace.mjs` re-integrates from the logged seed with the **shipped-identical**
  RK4 (`tools/chimera-campaign/integrator.mjs`) and the single-source-of-truth
  labeler (`tools/absorption-recampaign/labeling.mjs`), so `t_graze`/`t_abs`
  reproduce the campaign **bit-for-bit** (verified in the report).
- `fig7_curves.py` calls the reduced-ODE core (`tools/reduced-ode/reduced_core.py`)
  for Eqs 17/18 series, the numeric SN/Hopf locators, and the homoclinic
  bisection; the β=0.05 point reproduces **A_hc = 0.4096** within bracket
  tolerance (self-check, PASS).
- All paths are config-relative; no randomness beyond the logged seeds.
