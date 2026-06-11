"""Utils — shared helpers (config loading, seeding, PDB IO + geometry).

Every stage reads config/proteus.yaml the same way (load_config) and seeds RNGs
from the single global random_seed (get_seed). Structure helpers wrap biotite so
S4 (triad geometry) and S5 (cleft metrics) share one PDB-parsing/geometry path.
"""
from __future__ import annotations

import os
from typing import Iterable, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
DEFAULT_CONFIG = os.path.join(REPO, "config", "proteus.yaml")


# --------------------------------------------------------------------------- #
# Config + seeding
# --------------------------------------------------------------------------- #
def load_config(path: str = DEFAULT_CONFIG) -> dict:
    """Load config/proteus.yaml into a dict (PyYAML). Raises if the file or
    PyYAML is missing — the calibration stages cannot run without real config."""
    import yaml  # noqa: PLC0415 - import lazily so a bare `import proteus.utils` is cheap
    with open(path) as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg


def get_seed(cfg: dict) -> int:
    """Return the single global random_seed every stochastic stage must read."""
    if "random_seed" not in cfg:
        raise KeyError("random_seed missing from config")
    return int(cfg["random_seed"])


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def euclidean(a, b) -> float:
    """Euclidean distance between two 3-vectors (or any equal-length vectors)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.linalg.norm(a - b))


# --------------------------------------------------------------------------- #
# PDB IO + atom selection (biotite)
# --------------------------------------------------------------------------- #
def load_structure(path: str, model: int = 1):
    """Load one model of a PDB file as a biotite AtomArray."""
    from biotite.structure.io.pdb import PDBFile  # noqa: PLC0415
    return PDBFile.read(path).get_structure(model=model)


def protein_atoms(arr):
    """Subset of `arr` that are standard amino-acid atoms (drops waters/ligands/ions)."""
    import biotite.structure as struc  # noqa: PLC0415
    return arr[struc.filter_amino_acids(arr)]


def residue_iter(arr) -> Iterable[tuple]:
    """Yield (chain_id, res_id, res_name, residue_atom_array) for each residue."""
    import biotite.structure as struc  # noqa: PLC0415
    starts = struc.get_residue_starts(arr)
    bounds = list(starts) + [arr.array_length()]
    for i in range(len(starts)):
        sub = arr[bounds[i]:bounds[i + 1]]
        yield str(sub.chain_id[0]), int(sub.res_id[0]), str(sub.res_name[0]), sub


def atom_coord(residue_arr, atom_name: str) -> Optional[np.ndarray]:
    """Coord of a named atom within a single-residue AtomArray, or None if absent."""
    mask = residue_arr.atom_name == atom_name
    if not mask.any():
        return None
    return np.asarray(residue_arr.coord[mask][0], dtype=float)


def atom_coords(residue_arr, atom_names: Iterable[str]) -> list:
    """Coords of any of `atom_names` present in a single-residue AtomArray."""
    names = list(set(atom_names))
    mask = np.isin(residue_arr.atom_name, names)
    return [np.asarray(c, dtype=float) for c in residue_arr.coord[mask]]


def backbone_amide_donors(arr) -> list:
    """List of (chain_id, res_id, res_name, N_coord) for every backbone amide N.

    Backbone N-H groups are the oxyanion-hole hydrogen-bond donors S4 looks for.
    """
    mask = arr.atom_name == "N"
    out = []
    for chain, rid, rname, coord in zip(
        arr.chain_id[mask], arr.res_id[mask], arr.res_name[mask], arr.coord[mask]
    ):
        out.append((str(chain), int(rid), str(rname), np.asarray(coord, dtype=float)))
    return out


def per_atom_sasa(arr, probe_radius: float = 1.4) -> np.ndarray:
    """Per-atom solvent-accessible surface area (Shrake-Rupley) for `arr`.

    Returns an array aligned to `arr` (NaN for atoms biotite excludes). Compute on
    a protein-only array for clean values.
    """
    import biotite.structure as struc  # noqa: PLC0415
    return struc.sasa(arr, probe_radius=probe_radius)
