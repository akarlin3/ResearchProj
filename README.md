# ResearchProj

Consolidated monorepo bringing together my `proj*` research repositories. Each
project is preserved in its own subdirectory **with full commit history**.

| Subfolder | Origin repo | History |
|-----------|-------------|---------|
| `Forge/` | projForge | full history |
| `Gauge/` | projGauge | full history |
| `Minos/` | projMinos | full history |
| `Ouroboros/` | projOuroboros | full history |
| `Proteus/` | projProteus | full history |
| `Fashion/` | projFashion | fork — **only my own 21 commits**; upstream (`OSIPI/TF2.4_IVIM-MRI_CodeCollection`) history re-rooted to a single fork-point snapshot |

Each subdirectory's history was rewritten with `git-filter-repo` and combined with
`git merge --allow-unrelated-histories`, so `git log -- <Subfolder>/` shows that
project's original commits, authors, and dates.
