"""S4 — Catalytic geometry gate: Ser-His-Asp triad + oxyanion hole.

Given folded models (or control crystal structures), test for a serine-hydrolase
catalytic machine: a Ser-His-Asp(/Glu) triad in hydrogen-bonding geometry plus a
backbone-amide oxyanion hole. This is an UNSUPERVISED detector — it does not know
where the active site is; it enumerates every geometrically valid triad. Geometry
windows come from config (s4_geometry); nothing is hardcoded.

A triad PASSES only if both triad hydrogen bonds are within tolerance AND a backbone
oxyanion hole is present. The control-recovery harness (tests/test_s4_geometry.py and
the --recover CLI flag) then checks, on known controls, that the blind detector finds
the documented catalytic serine — and correctly finds nothing at the mutated S165A
site of 6THS.

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.s4_geometry \
        --pdb structures/6EQE.pdb --out data/processed/s4_triads.json
    PYTHONPATH=src python -m proteus.s4_geometry --recover   # run control checks
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

from proteus.utils import (
    DEFAULT_CONFIG,
    REPO,
    atom_coord,
    atom_coords,
    backbone_amide_donors,
    euclidean,
    load_structure,
    protein_atoms,
    residue_iter,
)

# Carboxylate oxygen names per acidic residue (His ND1 hydrogen-bonds to one of these).
_ACID_O = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}


def _unit_sphere(n: int = 200) -> np.ndarray:
    """Deterministic ~uniform unit vectors (Fibonacci sphere) — fixed sampling, no RNG,
    so oxyanion-hole detection is reproducible without consuming the global seed."""
    i = np.arange(n) + 0.5
    y = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - y * y, 0.0, None))
    theta = np.pi * (3.0 - np.sqrt(5.0)) * np.arange(n)
    return np.column_stack((np.cos(theta) * r, y, np.sin(theta) * r))


_SPHERE = _unit_sphere(200)

# Documented catalytic serines for the positive controls (supervised recovery check).
# 6THS is the S165A inactivated trap: its catalytic serine is mutated, so recovery
# of a triad at 165 is EXPECTED to fail — that is a correct null, not a bug.
KNOWN_SITES = {
    "6EQE": {"role": "positive", "name": "IsPETase", "ser": 160, "his": 237, "acid": 206},
    "4EB0": {"role": "positive", "name": "LCC_WT", "ser": 165, "his": 242, "acid": 210},
    "6THS": {"role": "trap", "name": "LCC_ICCG_S165A", "ser": 165, "expect_triad": False},
}


def _residue_tables(protein):
    """Collect the atoms S4 needs from each residue: Ser OG; His NE2/ND1; acid O's."""
    sers, hiss, acids = [], [], []
    for chain, rid, rname, sub in residue_iter(protein):
        if rname == "SER":
            og = atom_coord(sub, "OG")
            if og is not None:
                sers.append({"chain": chain, "res_id": rid, "OG": og})
        elif rname == "HIS":
            ne2, nd1 = atom_coord(sub, "NE2"), atom_coord(sub, "ND1")
            if ne2 is not None and nd1 is not None:
                hiss.append({"chain": chain, "res_id": rid, "NE2": ne2, "ND1": nd1})
        elif rname in _ACID_O:
            os_ = atom_coords(sub, _ACID_O[rname])
            if os_:
                acids.append({"chain": chain, "res_id": rid, "res_name": rname, "O": os_})
    return sers, hiss, acids


def _oxyanion_hole(og, ne2, donors, probe_min, probe_max, max_dist, min_donors, ser_res_id):
    """Find the oxyanion hole near Ser OG.

    The oxyanion (substrate carbonyl / tetrahedral-intermediate O-) sits in the
    hemisphere AWAY from the His, ~probe_min..probe_max A from OG. We scan that shell
    for the point with the most backbone-amide N donors within `max_dist`, requiring
    >= min_donors that share that point — so the donors genuinely converge on a single
    oxyanion site rather than merely being scattered near OG. The Ser's own amide is
    excluded; the hole is formed by other backbone NH groups (e.g. IsPETase Met161 +
    the Trp185 loop; LCC Ser165's i+1/i+2 amides)."""
    away = og - ne2
    norm = float(np.linalg.norm(away))
    if norm == 0.0 or not donors:
        return {"hole_ok": False, "donor_count": 0, "donors": [], "probe": None}
    away /= norm
    dcoords = np.array([d[3] for d in donors])
    dmask = np.array([d[1] != ser_res_id for d in donors])  # exclude Ser's own N
    radii = np.linspace(probe_min, probe_max, 5)
    best = {"hole_ok": False, "donor_count": 0, "donors": [], "probe": None}
    for u in _SPHERE:
        if float(u @ away) < 0.0:           # keep only the hemisphere opposite the His
            continue
        for r in radii:
            probe = og + u * r
            dist = np.linalg.norm(dcoords - probe, axis=1)
            sel = (dist <= max_dist) & dmask
            count = int(sel.sum())
            if count > best["donor_count"]:
                found = [{"chain": donors[j][0], "res_id": donors[j][1],
                          "res_name": donors[j][2], "dist": round(float(dist[j]), 3)}
                         for j in np.nonzero(sel)[0]]
                found.sort(key=lambda f: f["dist"])
                best = {"hole_ok": count >= min_donors, "donor_count": count,
                        "donors": found,
                        "probe": [round(float(x), 3) for x in probe]}
                if count >= min_donors and count >= 3:
                    return best  # strong hole — no need to keep searching
    return best


def detect_triads(arr, cfg: dict) -> list:
    """Enumerate every Ser-His-(Asp|Glu) triad with valid H-bond geometry + oxyanion
    hole. Returns a list of triad dicts (one per Ser-His pair, best acid kept)."""
    s4 = cfg["s4_geometry"]
    d_sh = float(s4["ser_og_his_ne2_max"])
    d_ha = float(s4["his_nd1_acid_max"])
    d_ox = float(s4["oxyanion_hole_max"])
    probe_min = float(s4["oxyanion_probe_min"])
    probe_max = float(s4["oxyanion_probe_max"])
    min_donors = int(s4["oxyanion_min_donors"])

    protein = protein_atoms(arr)
    sers, hiss, acids = _residue_tables(protein)
    donors = backbone_amide_donors(protein)

    triads = []
    for ser in sers:
        for his in hiss:
            d1 = euclidean(ser["OG"], his["NE2"])
            if d1 > d_sh:
                continue
            # best (closest) acid carboxylate O to this His ND1
            best_acid, best_d2 = None, None
            for acid in acids:
                d2 = min(euclidean(his["ND1"], o) for o in acid["O"])
                if d2 <= d_ha and (best_d2 is None or d2 < best_d2):
                    best_acid, best_d2 = acid, d2
            if best_acid is None:
                continue
            ox = _oxyanion_hole(ser["OG"], his["NE2"], donors, probe_min, probe_max,
                                d_ox, min_donors, ser["res_id"])
            triads.append({
                "ser": {"chain": ser["chain"], "res_id": ser["res_id"]},
                "his": {"chain": his["chain"], "res_id": his["res_id"]},
                "acid": {"chain": best_acid["chain"], "res_id": best_acid["res_id"],
                         "res_name": best_acid["res_name"]},
                "ser_og_his_ne2": round(d1, 3),
                "his_nd1_acid": round(best_d2, 3),
                "oxyanion": ox,
                "triad_geometry_ok": True,
                "passes": bool(ox["hole_ok"]),
            })
    # rank: passing first, then tightest combined H-bond geometry
    triads.sort(key=lambda t: (not t["passes"],
                               t["ser_og_his_ne2"] + t["his_nd1_acid"]))
    return triads


def analyze_model(pdb_path: str, cfg: dict) -> dict:
    """Run triad detection on model 1 of a PDB. Returns a per-model result with the
    full triad list, a triad_found flag, and the best (catalytic) triad if any."""
    arr = load_structure(pdb_path)
    triads = detect_triads(arr, cfg)
    passing = [t for t in triads if t["passes"]]
    best = passing[0] if passing else (triads[0] if triads else None)
    return {
        "model": os.path.basename(pdb_path),
        "n_triads": len(triads),
        "n_passing": len(passing),
        "triad_found": len(passing) > 0,
        "best": best,
        "triads": triads,
    }


def _catalytic_ser(result: dict, res_id: int):
    for t in result["triads"]:
        if t["ser"]["res_id"] == res_id:
            return t
    return None


def recovery_report(cfg: dict, struct_dir: str) -> list:
    """For each known control, report whether blind detection recovered the documented
    catalytic serine (and that 6THS correctly yields no triad there)."""
    lines = []
    for pdb_id, site in KNOWN_SITES.items():
        path = os.path.join(struct_dir, f"{pdb_id}.pdb")
        if not os.path.exists(path):
            lines.append({"pdb": pdb_id, **site, "status": "missing"})
            continue
        res = analyze_model(path, cfg)
        t = _catalytic_ser(res, site["ser"])
        rec = {"pdb": pdb_id, "name": site["name"], "role": site["role"],
               "known_ser": site["ser"], "triad_found_model": res["triad_found"]}
        if site.get("expect_triad") is False:
            # the trap: success == NO triad at the mutated serine
            rec["expected"] = "no triad at mutated S165A"
            rec["recovered_site_triad"] = t is not None
            rec["status"] = "EXPECTED-NULL-OK" if t is None else "UNEXPECTED-TRIAD"
        else:
            rec["recovered_site_triad"] = t is not None
            if t is not None:
                rec.update(ser_og_his_ne2=t["ser_og_his_ne2"],
                           his_nd1_acid=t["his_nd1_acid"],
                           his=t["his"]["res_id"], acid=t["acid"]["res_id"],
                           oxyanion_donors=t["oxyanion"]["donor_count"],
                           passes=t["passes"])
            rec["status"] = "RECOVERED" if (t is not None and t["passes"]) else "MISSED"
        lines.append(rec)
    return lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", action="append", default=None,
                    help="PDB file(s) to analyze (repeatable)")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "processed", "s4_triads.json"),
                    help="path to write per-model triad results")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--recover", action="store_true",
                    help="run the control-recovery check on KNOWN_SITES and print it")
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"))
    args = ap.parse_args(argv)

    from proteus.utils import load_config  # noqa: PLC0415
    cfg = load_config(args.config)

    if args.recover:
        rep = recovery_report(cfg, args.struct_dir)
        print("S4 RECOVERY (blind detection vs documented catalytic serine):")
        for r in rep:
            print(f"  {r['pdb']} {r.get('name','')} [{r.get('role','')}]: {r['status']}", end="")
            if r["status"] == "RECOVERED":
                print(f"  Ser{r['known_ser']} OG..NE2={r['ser_og_his_ne2']}A "
                      f"ND1..acid={r['his_nd1_acid']}A his={r['his']} acid={r['acid']} "
                      f"oxN={r['oxyanion_donors']}")
            elif r["status"].startswith("EXPECTED-NULL"):
                print(f"  (S165A trap: no triad at 165 — model triad_found="
                      f"{r['triad_found_model']})")
            else:
                print()
        return 0

    pdbs = args.pdb or [os.path.join(args.struct_dir, f"{p}.pdb")
                        for p in ("6EQE", "4EB0", "6THS", "1TCA")]
    pdbs = [p for p in pdbs if os.path.exists(p)]
    if not pdbs:
        print("no input PDBs found — fetch controls or pass --pdb", file=sys.stderr)
        return 2

    results = {os.path.basename(p): analyze_model(p, cfg) for p in pdbs}
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
        fh.write("\n")
    for name, r in results.items():
        b = r["best"]
        tag = (f"Ser{b['ser']['res_id']} (OG..NE2={b['ser_og_his_ne2']}A)"
               if b else "no triad")
        print(f"[s4] {name}: triad_found={r['triad_found']} "
              f"n_passing={r['n_passing']} best={tag}")
    print(f"[s4] results -> {os.path.relpath(args.out, os.getcwd())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
