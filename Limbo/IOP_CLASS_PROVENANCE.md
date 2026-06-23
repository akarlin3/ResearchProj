# Vendored IOP class files — provenance

`limbo.tex` targets **Physics in Medicine & Biology** (Topical Review) and is typeset with
IOP Publishing's official journal class. The class is **not on CTAN** and **cannot be fetched
by tectonic**, so the required files are vendored here so the manuscript builds self-contained.

| file | role | license |
|---|---|---|
| `iopjournal.cls` | IOP Publishing journal article class (2025) | LaTeX Project Public License (LPPL) 1.3c — redistribution/modification permitted |
| `orcid.pdf` | ORCID iD logo used by the class's `\orcid{}` command | IOP author-template asset, shipped with the class |

- **Source:** IOP Publishing official LaTeX template,
  `https://publishingsupport.iopscience.iop.org/wp-content/uploads/2025/07/ioplatextemplate.zip`
  (linked from `https://publishingsupport.iopscience.iop.org/questions/latex-template/`).
- **Retrieved:** 2026-06-22.
- **Integrity (sha256):**
  - `iopjournal.cls`: `0fa5404542b8ba606b9c63a55cd818733cf259cb93ea34198a4bdd8d665f4b94`
  - `orcid.pdf`: `7276a964b7439d1a69b55f63568a0f4f2126874712d3850d242dde3e2394047f`
- **Why vendored, not fetched:** `iopart`/`iopjournal` are distributed by IOP for author use and
  are absent from the CTAN bundle tectonic resolves against (verified: `ctan.org/pkg/iopart` →
  "Not found"; a minimal `\documentclass{iopjournal}` doc fails under tectonic until the `.cls`
  is local). The class header carries the LPPL, which permits redistribution.

The header text written by `\articletype{Topical Review}` shows the class's placeholder
journal branding ("Journal Name"); IOP sets the real running head at production. This is expected
for an author-prepared submission and does not affect content.
