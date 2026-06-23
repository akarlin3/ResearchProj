#!/usr/bin/env python3
"""Citation verification gate for the Procrustes manuscript.

House discipline: no fabricated DOIs; in-program companions are forward-cited as
unpublished (no invented DOI) with a finalization checklist; every other DOI must
be well-formed and (when a network is available) resolve.

Checks
------
1. Every @unpublished entry (Gauge, Lattice) carries NO doi field and a note that
   states its status -- the honest forward-cite. Emits a finalization checklist
   line for each: "on publication, set the DOI in release_config.json and swap the
   bib entry".
2. Every entry WITH a doi has a well-formed DOI (10.NNNN/...). If the network is
   reachable, each DOI is resolved via https://doi.org/<doi> (HEAD, short timeout);
   unreachable network -> reported as UNVERIFIED (non-fatal), bad/4xx -> FAIL.
3. No DOI on an entry whose note marks it unverified/TODO (guards against a
   fabricated DOI slipping in next to a TODO-AUTHORS marker).
4. Every \\cite/\\citep key used in procrustes.tex is defined in refs.bib.

Exit 0 = all structural checks pass (DOI network resolution is advisory); 1 = a
structural problem (malformed/blocked DOI, undefined key, mis-formed forward-cite).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIB = HERE / "refs.bib"
TEX = HERE / "procrustes.tex"

# in-program companions that MUST be forward-cited (unpublished, no DOI)
FORWARD_CITE = {"gauge", "lattice"}
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")


def parse_bib(text):
    """Minimal BibTeX parser -> {key: {"type":..., "doi":..., "note":..., "raw":...}}."""
    entries = {}
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,]+),", text):
        etype, key = m.group(1).lower(), m.group(2).strip()
        start = m.end()
        depth, i = 1, m.start(0) + text[m.start(0):].index("{") + 1
        # find matching closing brace from the entry's opening brace
        j = i
        while j < len(text) and depth:
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
            j += 1
        body = text[start:j - 1]
        doi = None
        dm = re.search(r"doi\s*=\s*[{\"]([^}\"]+)[}\"]", body, re.I)
        if dm:
            doi = dm.group(1).strip()
        nm = re.search(r"note\s*=\s*[{\"](.+?)[}\"]\s*,?\s*\n", body, re.I | re.S)
        note = (nm.group(1) if nm else "").strip()
        entries[key] = {"type": etype, "doi": doi, "note": note}
    return entries


def try_resolve(doi):
    """Return 'OK' / 'FAIL(code)' / 'UNVERIFIED(no-net)'. Advisory; short timeout."""
    import urllib.request
    import urllib.error
    url = "https://doi.org/" + doi
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=6) as r:
            return "OK" if r.status < 400 else f"FAIL({r.status})"
    except urllib.error.HTTPError as e:
        # doi.org often 302->publisher; many publishers 403 a HEAD bot -> not fatal
        return "OK" if e.code in (301, 302, 303, 403, 405) else f"FAIL({e.code})"
    except Exception:
        return "UNVERIFIED(no-net)"


def main():
    print("Procrustes citation verification gate")
    print("=" * 64)
    if not BIB.exists():
        print(f"FAIL: missing {BIB.name}", file=sys.stderr)
        return 1
    entries = parse_bib(BIB.read_text())
    fails, checklist = [], []

    # (1) forward-cite companions: unpublished, no DOI, status note
    for key in sorted(FORWARD_CITE):
        if key not in entries:
            fails.append(f"forward-cite companion '{key}' missing from refs.bib")
            continue
        e = entries[key]
        if e["type"] != "unpublished":
            fails.append(f"'{key}' must be @unpublished (is @{e['type']})")
        if e["doi"]:
            fails.append(f"'{key}' is unpublished but carries a DOI (no invented DOI allowed)")
        if not e["note"]:
            fails.append(f"'{key}' must carry a status note (honest forward-cite)")
        checklist.append(f"  [ ] {key}: on publication, set its DOI in release_config.json "
                         f"and swap refs.bib @unpublished -> published entry")

    # (2)/(3) DOI well-formedness + (advisory) resolution
    print("  DOI checks (resolution is advisory):")
    for key, e in sorted(entries.items()):
        if e["doi"] is None:
            continue
        if not DOI_RE.match(e["doi"]):
            fails.append(f"'{key}' has a malformed DOI: {e['doi']!r}")
            continue
        if re.search(r"todo|verify|tbd|xxxx", e["note"], re.I):
            fails.append(f"'{key}' has a DOI next to an unverified/TODO note "
                         f"(possible fabrication): {e['doi']}")
        status = try_resolve(e["doi"])
        if status.startswith("FAIL"):
            fails.append(f"'{key}' DOI does not resolve: {e['doi']} -> {status}")
        print(f"    {key:16s} {e['doi']:34s} {status}")

    # (4) every cite key in the manuscript is defined
    if TEX.exists():
        used = set(re.findall(r"\\cite[a-z]*\{([^}]+)\}", TEX.read_text()))
        used = {k.strip() for grp in used for k in grp.split(",")}
        undefined = sorted(used - set(entries))
        if undefined:
            fails.append(f"cite keys used but not in refs.bib: {undefined}")

    print("\n  Finalization checklist (forward-cited, unpublished):")
    for line in checklist:
        print(line)

    print("=" * 64)
    if fails:
        for f in fails:
            print(f"  CITE FAIL: {f}", file=sys.stderr)
        print("CITATION GATE: FAIL", file=sys.stderr)
        return 1
    print("CITATION GATE: PASS (structural; DOI resolution advisory where network allows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
