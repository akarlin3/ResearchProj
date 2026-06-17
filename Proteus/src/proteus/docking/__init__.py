"""Docking (P4) — AutoDock Vina wrappers (CPU, LOCAL).

Dock a PET-mimic ligand (default BHET) into the catalytic cleft of the screen's
PETase-like hits and rank by binding affinity. The search box is centred on the
active site found by S4 (the catalytic Ser OG) — so docking probes the SAME site
the geometry/cleft gates scored, not an arbitrary pocket. GPU docking
(GNINA / DiffDock) and Chai-1 cofolding are GCE burst targets, not here.

The Vina call is INJECTABLE (a `scorer` callable), exactly like the S3 fold
backend: the orchestration — active-site box placement, candidate iteration,
affinity ranking, top-N selection, output — is exercised by tests with a
deterministic fake scorer on a CPU host, while the real AutoDock Vina backend
(`vina_scorer`, lazy-imported) runs on the M4. The Vina search is seeded from the
global `random_seed` (never hardcoded).

A candidate is dockable only if S4 finds a catalytic Ser to centre the box on;
ranking is by affinity (more negative kcal/mol = tighter).

Local usage, from the repo root (after the S4/S5 screen):
    PYTHONPATH=src python -m proteus.docking \
        --candidates data/processed/s4s5_candidates.json \
        --models structures/folded \
        --out data/processed/docking
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from proteus.s4_geometry import analyze_model
from proteus.utils import (
    DEFAULT_CONFIG,
    REPO,
    atom_coord,
    load_config,
    load_structure,
    protein_atoms,
    residue_iter,
)


def _ser_og_coord(pdb_path: str, ser_res_id: int):
    """Coordinate of the catalytic Ser OG (the box centre), or None if absent."""
    arr = protein_atoms(load_structure(pdb_path))
    for _chain, rid, rname, sub in residue_iter(arr):
        if rid == ser_res_id and rname == "SER":
            return atom_coord(sub, "OG")
    return None


def box_center_from_model(pdb_path: str, cfg: dict):
    """Find the active-site box centre. Runs S4; centres on the catalytic Ser OG of
    the best passing triad. Returns (center_xyz | None, ser_res_id | None)."""
    s4 = analyze_model(pdb_path, cfg)
    if not s4["triad_found"] or s4["best"] is None:
        return None, None
    ser = s4["best"]["ser"]["res_id"]
    og = _ser_og_coord(pdb_path, ser)
    return (None if og is None else np.asarray(og, dtype=float)), ser


def dock_model(pdb_path: str, cfg: dict, scorer, cand_id: str | None = None,
               mean_plddt=None) -> dict:
    """Dock one model: place the box on its catalytic Ser OG and score with `scorer`.

    `scorer(receptor, center, box_size, seed, exhaustiveness, n_poses, ligand)`
    returns {"affinity": float(kcal/mol), "n_poses": int}. Returns a candidate record.
    """
    dk = cfg["docking"]
    rec = {"id": cand_id or os.path.basename(pdb_path), "mean_plddt": mean_plddt,
           "catalytic_ser": None, "box_center": None, "affinity": None,
           "n_poses": 0, "docked": False}

    center, ser = box_center_from_model(pdb_path, cfg)
    if center is None:
        rec["stage_failed"] = "no_catalytic_site"  # S4 found no triad to dock into
        return rec
    rec["catalytic_ser"] = ser
    rec["box_center"] = [round(float(x), 3) for x in center]

    box_size = [float(x) for x in dk.get("box_size", [20.0, 20.0, 20.0])]
    seed = int(cfg["random_seed"])
    try:
        out = scorer(receptor=pdb_path, center=center, box_size=box_size, seed=seed,
                     exhaustiveness=int(dk.get("exhaustiveness", 8)),
                     n_poses=int(dk.get("n_poses", 9)),
                     ligand=dk.get("ligand_pdbqt") or dk.get("ligand"))
    except Exception as exc:  # noqa: BLE001 — a single docking failure must not kill the batch
        rec["stage_failed"] = f"dock_error: {exc}"
        return rec
    rec["affinity"] = round(float(out["affinity"]), 3)
    rec["n_poses"] = int(out.get("n_poses", 1))
    rec["docked"] = True
    return rec


def _inputs_from_candidates(candidates_json: str, models_dir: str,
                            hits_only: bool = True) -> list[dict]:
    """Read screen output (s4s5_candidates.json) and resolve PDB paths under
    models_dir. By default docks only the PETase-like hits."""
    with open(candidates_json) as fh:
        summary = json.load(fh)
    out = []
    for c in summary.get("candidates", []):
        if hits_only and not c.get("petase_like_hit"):
            continue
        pdb = os.path.join(models_dir, f"{c['id']}.pdb")
        if os.path.exists(pdb):
            out.append({"id": c["id"], "pdb": pdb, "mean_plddt": c.get("mean_plddt")})
    return out


def _inputs_from_dir(models_dir: str) -> list[dict]:
    return [{"id": os.path.splitext(fn)[0], "pdb": os.path.join(models_dir, fn),
             "mean_plddt": None}
            for fn in sorted(os.listdir(models_dir)) if fn.endswith(".pdb")]


def dock_models(inputs: list[dict], cfg: dict, scorer) -> dict:
    """Dock every input model, rank by affinity (ascending = tightest first), keep
    the configured top_n. Returns a ranked summary."""
    candidates = [dock_model(i["pdb"], cfg, scorer, cand_id=i["id"],
                             mean_plddt=i.get("mean_plddt")) for i in inputs]
    docked = [c for c in candidates if c["docked"]]
    docked.sort(key=lambda c: c["affinity"])  # more negative kcal/mol = tighter binder
    for rank, c in enumerate(docked, 1):
        c["rank"] = rank
    top_n = cfg.get("docking", {}).get("top_n")
    kept = docked[: int(top_n)] if top_n else docked
    return {
        "ligand": cfg.get("docking", {}).get("ligand"),
        "box_size": cfg.get("docking", {}).get("box_size"),
        "random_seed": cfg.get("random_seed"),
        "n_input": len(inputs),
        "n_docked": len(docked),
        "n_failed": len(candidates) - len(docked),
        "top_n": top_n,
        "ranking": [c["id"] for c in docked],
        "kept_ids": [c["id"] for c in kept],
        "best_affinity": docked[0]["affinity"] if docked else None,
        "candidates": candidates,
    }


def write_outputs(summary: dict, out_prefix: str) -> tuple[str, str]:
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    js = out_prefix + ".json"
    with open(js, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    cols = ["rank", "id", "mean_plddt", "catalytic_ser", "affinity", "n_poses",
            "box_center", "docked", "kept"]
    kept = set(summary["kept_ids"])
    csv_path = out_prefix + ".csv"
    by_id = {c["id"]: c for c in summary["candidates"]}
    order = summary["ranking"] + [c["id"] for c in summary["candidates"]
                                  if c["id"] not in summary["ranking"]]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for cid in order:
            c = by_id[cid]
            w.writerow({
                "rank": c.get("rank", ""), "id": c["id"],
                "mean_plddt": c.get("mean_plddt"), "catalytic_ser": c["catalytic_ser"],
                "affinity": c["affinity"], "n_poses": c["n_poses"],
                "box_center": c["box_center"], "docked": c["docked"],
                "kept": cid in kept,
            })
    return js, csv_path


# --------------------------------------------------------------------------- #
# Real AutoDock Vina backend (lazy — only built when actually docking on the Mac)
# --------------------------------------------------------------------------- #
def prepare_receptor_pdbqt(pdb_path: str, out_pdbqt: str, obabel_bin: str = "obabel",
                           ph: float = 7.4) -> str:
    """Convert a receptor PDB to a rigid AutoDock PDBQT (Open Babel: add polar H at
    the given pH, assign Gasteiger charges + AutoDock atom types). Meeko/ADFR are
    equivalent alternatives; Open Babel is the dependency-light default. Returns
    `out_pdbqt`."""
    if shutil.which(obabel_bin) is None:
        raise FileNotFoundError(
            f"'{obabel_bin}' not on PATH — install Open Babel (or prep the receptor "
            "PDBQT with Meeko/ADFR) to dock a folded .pdb model")
    os.makedirs(os.path.dirname(os.path.abspath(out_pdbqt)), exist_ok=True)
    proc = subprocess.run(
        [obabel_bin, pdb_path, "-O", out_pdbqt, "-xr", "-p", str(ph)],
        capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.exists(out_pdbqt):
        raise RuntimeError(f"obabel receptor prep failed (rc={proc.returncode}): "
                           f"{proc.stderr[-400:]}")
    return out_pdbqt


def vina_scorer(sf_name: str = "vina", obabel_bin: str = "obabel"):
    """Build the real AutoDock Vina scorer callable. Imports `vina` lazily (it is a
    LOCAL-only dependency). A receptor passed as a `.pdb` (e.g. a folded model) is
    auto-prepped to PDBQT via Open Babel; a `.pdbqt` is used as-is. The ligand must
    be a prepared PDBQT (the committed PET-mimic, controls/ligands/bhet.pdbqt, or a
    Meeko prep). The box is centred on the catalytic Ser OG."""
    from vina import Vina  # noqa: PLC0415

    def _score(receptor, center, box_size, seed, exhaustiveness, n_poses, ligand):
        if not (isinstance(ligand, str) and ligand.endswith(".pdbqt") and os.path.exists(ligand)):
            raise RuntimeError(
                "vina_scorer needs a prepared ligand PDBQT (docking.ligand_pdbqt); "
                "use controls/ligands/bhet.pdbqt or prep the PET-mimic with Meeko")
        tmp = None
        if receptor.endswith((".pdb", ".ent")):
            tmp = tempfile.mkdtemp(prefix="dock_rec_")
            receptor_pdbqt = prepare_receptor_pdbqt(
                receptor, os.path.join(tmp, "receptor.pdbqt"), obabel_bin)
        else:
            receptor_pdbqt = receptor
        try:
            v = Vina(sf_name=sf_name, seed=int(seed), verbosity=0)
            v.set_receptor(receptor_pdbqt)
            v.set_ligand_from_file(ligand)
            v.compute_vina_maps(center=[float(x) for x in center],
                                box_size=[float(x) for x in box_size])
            v.dock(exhaustiveness=int(exhaustiveness), n_poses=int(n_poses))
            energies = v.energies(n_poses=int(n_poses))
            best = float(min(e[0] for e in energies)) if len(energies) else float(v.score()[0])
            return {"affinity": best, "n_poses": int(len(energies)) or 1}
        finally:
            if tmp:
                shutil.rmtree(tmp, ignore_errors=True)

    return _score


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", default=None,
                    help="screen output JSON (s4s5_candidates.json); docks the hits")
    ap.add_argument("--models", default=os.path.join(REPO, "structures", "folded"),
                    help="dir holding the candidate PDBs")
    ap.add_argument("--all", action="store_true",
                    help="with --candidates, dock every screened model, not only hits")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "processed", "docking"),
                    help="output path prefix (.json + .csv)")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if args.candidates:
        if not os.path.exists(args.candidates):
            print(f"candidates JSON not found: {args.candidates}", file=sys.stderr)
            return 2
        inputs = _inputs_from_candidates(args.candidates, args.models, hits_only=not args.all)
    else:
        if not os.path.isdir(args.models):
            print(f"models dir not found: {args.models}", file=sys.stderr)
            return 2
        inputs = _inputs_from_dir(args.models)
    if not inputs:
        print("no models to dock (no hits / no PDBs found)", file=sys.stderr)
        return 1

    try:
        scorer = vina_scorer()
    except Exception as exc:  # noqa: BLE001
        print(f"AutoDock Vina not available ({exc}). Install `vina` (+ a prepared "
              "ligand PDBQT) to dock locally.", file=sys.stderr)
        return 3

    summary = dock_models(inputs, cfg, scorer)
    js, csv_path = write_outputs(summary, args.out)
    print(f"[dock] {summary['n_docked']}/{summary['n_input']} docked "
          f"(best affinity {summary['best_affinity']} kcal/mol); "
          f"kept top {len(summary['kept_ids'])}.")
    print(f"[dock] -> {os.path.relpath(csv_path, os.getcwd())} ; "
          f"{os.path.relpath(js, os.getcwd())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
