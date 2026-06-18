# Two-system reconstruction of the NLD manuscript — changelog

`death_of_a_chimera_nld.tex` rebuilt from the mean-field-only paper into the full
two-system manuscript **"Chimera Collapse Ages: Topology-Dependent Finite-Size
Scaling in Mean-Field and Ring Oscillator Systems"** (Springer `sn-jnl`), folding
in the C1/C2 critical-review results as the validated form of the ring and
mean-field claims.

## Recovery audit (what was recovered vs. reconstructed)

**Nothing was recovered from git history — a prior two-system LaTeX source does
not exist.** Pickaxe searches (`git log --all -S`) for "Chimera Collapse on a
Nonlocal Ring", "Topology-Dependent Finite-Size Scaling", and "Collapse Ages"
across all files in all local and remote refs returned zero hits; the complete
historical inventory of `.tex`/`.pdf` paths contains only the three
`death_of_a_chimera*` sources, the supplement, and analysis figures. No compiled
"Chimera Collapse Ages" PDF exists in the repository, in any worktree, or as an
untracked file. **All ring-paper text is therefore reconstruction, not
recovery**, drawn exclusively from committed, validated sources:

| Content | Source |
|---|---|
| Ring model, event criterion, campaign design (Sec. 2.2) | `anneal-hazard/PREREGISTRATION.md`, `anneal-hazard/SUMMARY.md`, `anneal-hazard/src/ring_fast.py` (ρ_k definition) |
| Sec. 7 structured-hazard results (15/15 cells) | `anneal-hazard/results/cp4_fits.{csv,json}`, `SUMMARY.md` |
| Sec. 7 object validation (C2/CP2) | `critical_review/cp2_validation/cp2_summary.json`, `paper/critical_review_supplement.tex` §2 |
| Sec. 7 criterion robustness (C2/CP3) | `critical_review/cp3_criterion/cp3_criterion.json`, supplement §2 |
| Sec. 7 dwell/descent decoupling | `anneal-hazard/SUMMARY.md` (mechanism block) |
| Appendix C extrapolation | `anneal-hazard/results/extrapolation/{VERDICT.md, CPB_n5_REPORT.txt, cpB_n5_fits.json, k_vs_N_table.csv}`, `analysis/cp_b_n5.py` (bootstrap method) |
| Sec. 4 C1 stratification insert | `critical_review/cp1_stratified/cp1_stratified.json`, supplement §1 and §3 |
| Two-system abstract / intro / discussion framing | newly drafted (no prior text exists); ring sentences follow the C2-validated framing of supplement §3 |

Where the supplement provided prose (C1, C2), the manuscript text follows it
closely — it is the committed, post-review-validated wording.

## Every inserted number and its source

- **C1 (Sec. 4.3):** 84/84 size×stratum cells with k>1 and CI excluding 1; k∈[2.00,5.00] under the reduced-prediction (t_capture) binning; LRT p ≤ 1e-7; CV 0.22–0.56 tracking the Weibull CV–k relation at correlation 0.996; 0% censoring (all 1,400 runs) — `critical_review/cp1_stratified/cp1_stratified.json` (verified to match `critical_review/CHANGELOG.md` exactly).
- **Ring campaign (Sec. 2.2, 7.1–7.2):** 5 β × {32,64,128} × 300 = 4,500 runs + 5 β × {192,256} × 300 = 3,000 runs (7,500 total); k̂ range 1.13–1.68; LRT p 2e-2…5e-27; hazard rise factors 1.46–2.27; runs-test p 1e-27…1e-58; ε_std robustness ≲0.07; censoring 0.8% overall, 9.3% worst cell (β=0.110, N=32), 0% at N≥64 — `anneal-hazard/results/cp4_fits.{csv,json}`, `cp_fits_N192.json`, `cp_fits_N256.json`, `SUMMARY.md`. All 25 (β,N) cells re-verified directly from the JSONs: every CI excludes 1, every LRT p<0.05.
- **Lifetime inversion:** median τ at β=0.130: 872→422→404→349→345 for N=32→256 — `k_vs_N_table.csv` + `cp_fits_N192.json`.
- **C2 object validation (Sec. 7.3):** 0/600 τ mismatches; canonical 83.7%/89.3%; multi-head 14.3%/8.0%; degenerate 2.0%/2.7%; chimera-death 98.0%/97.3%; k̂ 1.22[1.12,1.32]→1.24[1.14,1.35] (β=0.110) and 1.47[1.35,1.58]→1.48[1.36,1.60] (β=0.130); contrast-gate sweep [0.06,0.12] — `critical_review/cp2_validation/cp2_summary.json`, supplement Table 2.
- **C2 criterion robustness (Sec. 7.4, Fig. C1):** τ̃(256)/τ̃(32) = 0.34–0.40 under all four criteria; ρ_std criteria 0% censoring at every N; mean-coh artifact: 6.7% twisted-state (ρ_k≈0.50) mis-censoring at N=256, corrected k̂(N)=1.26,1.44,1.62,1.54,1.57; KM log-cumulative-hazard slope 1.5–2.1 — `critical_review/cp3_criterion/cp3_criterion.json`, supplement §2.
- **Dwell (Sec. 7.5):** mean 52.8, median 50.5, CV 0.72; Spearman |ρ| ≤ 0.14 — `anneal-hazard/SUMMARY.md` / `cp4_fits.csv`.
- **Extrapolation (Sec. 7.6, Appendix C, Table C1):** bounded2 k∞ = 1.35[1.27,1.43], 1.49[1.40,1.56], 1.56[1.48,1.65], 1.62[1.53,1.71], 1.65[1.56,1.74]; AICc gaps +2.98, +5.42, +2.29, +2.51, +4.65; n=4→n=5 firming 2.32→2.98, 2.10→2.29, 3.87→4.65; satgen3 γ = 1.63, 1.51, 1.14, 1.48, 1.51 (agreement within 0.02–0.04); divergent-fit slopes log 0.05–0.18 / power 0.04–0.13; χ² = 8–11 (non-monotonic cells) vs 0.70 (β=0.120); bootstrap B=5,000 Gaussian resamples in the symmetrized profile σ — `cpB_n5_fits.json`, `CPB_n5_REPORT.txt`, `VERDICT.md`.

## Reconciliation

The brief's claimed k∞ table (1.35, 1.49, 1.56, 1.62, 1.65) matches the committed
`cpB_n5_fits.json` bounded2 bootstrap medians (1.3531, 1.4841, 1.5643, 1.6208,
1.6495) exactly at 2-decimal rounding. **Zero drift in any cell; no
recovered-vs-validated conflicts arose** (there is no recovered text). The
critical-review CHANGELOG numbers were independently re-verified against the raw
JSON outputs: all exact or rounding-level matches, no discrepancies.

## Figures

- **Fig. 10 (new):** ring survival/k̂(N)/k∞ summary. Generated by the new
  `tools/paper-figures/fig10_ring.py` → `paper_figures/fig10.{pdf,png}`. Pure
  rendering of committed analysis outputs (`cp4_fits.json`, `cp_fits_N192.json`,
  `cp_fits_N256.json`, the three ensemble CSVs, `cpB_n5_fits.json`) — **no
  campaign was re-run and no fit recomputed**, so no Zenodo config-set bump is
  required; the prior PNG-only versions (`cpA/cpB/cpB_n5_k_vs_N.png`) remain
  untouched in `anneal-hazard/results/extrapolation/`.
- **Fig. 11 (new):** pre-death spatial fields — included directly from the
  committed `paper/critical_review_figs/cp2_object.pdf` (the strongest new C2
  figure: it shows what the criterion kills).
- **Fig. C1 (new):** criterion comparison — committed
  `paper/critical_review_figs/cp3_criterion.pdf`.
- Existing figures fig1–fig9, figA, figB unchanged. Numbering is monotonic in
  order of first citation: main 1–11, appendix A1, B1, C1. No figure could not
  be produced; nothing required regeneration from configs.

## Carryover fixes

- Section 7's system reference points at the ring equation `eq:ring` (Eq. 2),
  not the reduced-flow equation (`eq:threeD`).
- Language tempering was already applied in the NLD conversion (no instance of
  "immortal"/"never lets go"/"folklore" remained); one new occurrence inherited
  from the supplement ("spurious immortal tail") was rephrased to "spurious
  never-dying tail".
- Springer-compliant AI disclosure confirmed present in Acknowledgments
  (drafting+code assistance, author sole and accountable, AI not an author).
- **Build fix:** the shipped `sn-jnl.cls` is **v0.1 (2019)** — not v2.x as the
  NLD conversion changelog stated — and references `manyfoot` and `xcolor`
  commands without loading them, so the manuscript did not compile on current
  TeX distributions. Two `\RequirePackage` lines (manyfoot, xcolor) were added
  before `\documentclass` with an explanatory comment. No other class change.

## Build verification (tectonic, sn-jnl/sn-mathphys Numbered)

- Compiles clean: 25 pages, all references and citations resolve (no undefined),
  no missing figures.
- Two-system title and abstract; Sec. 2.2 (ring definition, own labeled
  equation); Sec. 7 (ring collapse, C2-validated, no deferrals); Appendix C
  (extrapolation + Table C1 + Fig. C1); cross-topology Discussion ("aging
  transfers, scaling does not") and Conclusion; keywords (+ "nonlocal
  coupling") and Declarations present.
- The single overfull-hbox warning (16.9 pt, Table 3 / `tab:match`) is
  byte-identical in the pre-edit manuscript — pre-existing, not introduced.

## Not done (by design)

- No campaigns re-run; primary data untouched.
- Nothing committed or pushed — left for review per the brief.
