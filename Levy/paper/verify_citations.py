#!/usr/bin/env python
"""Citation gate for the Levy manuscript: every \\cite key resolves to a \\bibitem, every
\\bibitem is cited (no orphans), and the forward-cited Ouroboros tooling is present.

Run:  <proteus python> paper/verify_citations.py
Exit 0 = all citations resolve and Ouroboros is forward-cited; nonzero names the problem.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE / "levy.tex"


def main() -> int:
    src = TEX.read_text()
    cited = set()
    for m in re.findall(r"\\cite[tp]?\{([^}]*)\}", src):
        cited.update(k.strip() for k in m.split(",") if k.strip())
    defined = set(re.findall(r"\\bibitem\{([^}]+)\}", src))

    print(f"  cited keys   ({len(cited)}): {sorted(cited)}")
    print(f"  bibitem keys ({len(defined)}): {sorted(defined)}")

    ok = True
    unresolved = cited - defined
    if unresolved:
        print(f"  FAIL: cited but undefined: {sorted(unresolved)}")
        ok = False
    else:
        print("  all \\cite keys resolve to a \\bibitem: PASS")

    orphans = defined - cited
    if orphans:
        print(f"  FAIL: defined but never cited (orphans): {sorted(orphans)}")
        ok = False
    else:
        print("  no orphan bibitems: PASS")

    # Ouroboros must be forward-cited (tooling provenance), with the finalization marker present.
    if "ouroboros" not in defined:
        print("  FAIL: Ouroboros not in the bibliography (must be forward-cited)")
        ok = False
    elif "ouroboros" not in cited:
        print("  FAIL: Ouroboros bibitem present but never \\cite'd in the text")
        ok = False
    elif "FORWARD-CITE-ouroboros" not in src:
        print("  FAIL: Ouroboros finalization marker (% FORWARD-CITE-ouroboros) missing")
        ok = False
    else:
        print("  Ouroboros forward-cited with finalization marker: PASS")

    print("citation gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
