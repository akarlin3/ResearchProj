# Echo manuscript (built PASS-only, CP4)

The full manuscript (`echo.tex`, `ebgaramond` + `microtype`) and its traceability gate
(`consistency.py` → auto-generated `numbers.tex`) are assembled at **CP4, only if CP3
returns PASS**. On a Lethe-fold or Reverb-fold the deliverable is instead the honest-
limitation note / synthetic-harness spec described in `../PROMOTION.md`.

Build (once `echo.tex` exists):
```bash
bash build.sh        # tectonic preferred; falls back to pdflatex x2
```

Discipline (mirrors `Minos/future/paper/`):
- every `\num*` macro in `echo.tex` is defined in `numbers.tex`, which `consistency.py`
  regenerates from the seeded `results/RESULTS_*.json`;
- every Fashion/Minos/Gauge-dependent number carries the PROVISIONAL marker;
- `consistency.py` is the CP4 gate (zero undefined macros + internal-consistency asserts).
