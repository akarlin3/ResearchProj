"""Atlas retrieval front-end tests (proteus.atlas + proteus.atlas_screen).

Pure-logic units run with no network and no fpocket (HTTP is monkeypatched at
`proteus.atlas._request`). The screen-wiring integration test reuses the validated
S4/S5 path and skips without fpocket + the control structures, exactly like the
calibration/screen tests.
"""
from __future__ import annotations

import json
import os
import shutil

import pytest

import proteus.atlas as atlas
from proteus.atlas import (
    _parse_search_results,
    foldclass_union,
    mean_plddt,
    query_set_hash,
    run_retrieval,
)
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(REPO, "structures")
CONFIG = os.path.join(REPO, "config", "proteus.yaml")

CONTROL_IDS = ["6EQE", "4EB0", "6THS", "1TCA", "1EA5", "1CRL", "1EVQ", "8B4U",
               "4WFI", "4CG1"]


def _controls_present() -> bool:
    return all(os.path.exists(os.path.join(STRUCT, f"{p}.pdb")) for p in CONTROL_IDS)


def _cfg():
    cfg = load_config(CONFIG)
    assert "atlas" in cfg, "config must carry the atlas block"
    return cfg


# --------------------------------------------------------------------------- #
# Result parsing (tolerant to field-name variants)
# --------------------------------------------------------------------------- #
def test_parse_search_results_extracts_accession_bits_eval_tm():
    payload = {"results": [{"db": "highquality_clust30", "alignments": [
        {"target": "MGYP000911143812", "score": 240.0, "eval": 1.2e-9, "prob": 0.98},
        {"target": "MGYP000059025561 extra", "bits": 51.0, "evalue": 0.004,
         "tmScore": 0.6},
        {"target": "not-an-mgyp", "score": 999},          # dropped (no MGYP token)
    ]}]}
    hits = _parse_search_results(payload)
    accs = {h["accession"]: h for h in hits}
    assert set(accs) == {"MGYP000911143812", "MGYP000059025561"}
    assert accs["MGYP000911143812"]["bits"] == 240.0
    assert accs["MGYP000911143812"]["tm"] == 0.98
    assert accs["MGYP000059025561"]["bits"] == 51.0
    assert accs["MGYP000059025561"]["evalue"] == 0.004


def test_parse_search_results_handles_nested_alignment_lists():
    payload = {"results": [{"alignments": [[
        {"target": "MGYP000000000042", "score": 70}]]}]}
    hits = _parse_search_results(payload)
    assert len(hits) == 1 and hits[0]["accession"] == "MGYP000000000042"


def test_parse_search_results_empty_is_empty():
    assert _parse_search_results({}) == []
    assert _parse_search_results({"results": []}) == []


# --------------------------------------------------------------------------- #
# Fold-class union query set + reproducibility hash
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not _controls_present(), reason="control structures not fetched")
def test_foldclass_union_is_broad_and_excludes_trap():
    cfg = _cfg()
    members = foldclass_union(cfg, STRUCT)
    ids = {m["id"] for m in members}
    # Broad fold class: PETases AND non-PETase hydrolases (unseeded).
    assert {"IsPETase", "LCC_WT"} <= ids          # PETases present
    assert {"CalB", "Est2", "AChE_Tc"} & ids      # non-PETase hydrolases present
    # The configured trap exclusion (LCC_ICCG / 6THS) is dropped.
    assert "LCC_ICCG" not in ids


@pytest.mark.skipif(not _controls_present(), reason="control structures not fetched")
def test_query_set_hash_is_stable_and_order_independent():
    cfg = _cfg()
    members = foldclass_union(cfg, STRUCT)
    h1 = query_set_hash(members)
    h2 = query_set_hash(list(reversed(members)))
    assert h1 == h2 and len(h1) == 16


# --------------------------------------------------------------------------- #
# mean pLDDT from a PDB B-factor column
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not _controls_present(), reason="control structures not fetched")
def test_mean_plddt_reads_ca_bfactors():
    val, n = mean_plddt(os.path.join(STRUCT, "4EB0.pdb"))
    assert n > 0 and isinstance(val, float)


def test_mean_plddt_no_atoms(tmp_path):
    p = tmp_path / "empty.pdb"
    p.write_text("HEADER something\nEND\n")
    assert mean_plddt(str(p)) == (None, 0)


# --------------------------------------------------------------------------- #
# Retrieval driver: fallback path + min_bits/hit_cap, fully offline (monkeypatched)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not _controls_present(), reason="control structures not fetched")
def test_run_retrieval_falls_back_when_hosted_search_down(tmp_path, monkeypatch):
    """Hosted search 503 -> representative fallback -> fetch (all monkeypatched).
    Verifies the driver logs the fallback, applies the cap, fetches, and tags."""
    cfg = _cfg()
    cfg["atlas"]["backend"] = "api"
    cfg["atlas"]["allow_representative_fallback"] = True
    cfg["atlas"]["pilot_sample_size"] = 3
    cfg["atlas"]["hit_cap"] = 2

    fake_pdb = b"HEADER  test\nATOM      1  CA  ALA A   1       0.0     0.0     0.0  1.00  0.80           C\nEND\n"
    lookup = b"\n".join(f"{i}\tMGYP00000000{i:04d}\t0".encode() for i in range(20))

    def fake_request(url, *, data=None, headers=None, timeout=60, max_bytes=None):
        if "searchStructure/ticket" in url:            # submit -> 503 (down)
            return 503, b"Service Temporarily Unavailable"
        if "highquality_clust30.lookup" in url:        # representative slice
            return 206, lookup
        if "fetchPredictedStructure" in url:           # structure fetch
            return 200, fake_pdb
        if "metagenomics" in url:                      # biome (no mapping)
            return 404, b""
        return 404, b""

    monkeypatch.setattr(atlas, "_request", fake_request)

    out_json = tmp_path / "atlas_hits.json"
    hits_dir = tmp_path / "atlas_hits"
    summary = run_retrieval(cfg, STRUCT, str(out_json), str(hits_dir), log=lambda *a: None)

    assert summary["discovery_mode"] == "representative_fallback"
    assert summary["fallback_note"]                       # the 503 was recorded
    assert summary["n_hits_kept"] == 2                    # hit_cap applied
    assert summary["n_fetched"] == 2
    assert all(h["biome"] == "unknown" for h in summary["hits"])  # non-blocking biome
    assert out_json.exists()
    assert len(list(hits_dir.glob("*.pdb"))) == 2


# --------------------------------------------------------------------------- #
# Screen wiring at the widened operating point (integration; needs fpocket)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(shutil.which("fpocket") is None or not _controls_present(),
                    reason="fpocket + control structures required")
def test_resolve_widened_operating_point_matches_validation_line():
    from proteus.atlas_screen import resolve_operating_point
    cfg = _cfg()
    cfg["atlas"]["operating_point"] = "widened"
    op = resolve_operating_point(cfg, STRUCT)
    assert op["line"] == "widened"
    # widened line sits below the production line and near the -1.156 validation
    # value (fpocket jitter ~±0.01); it must keep all three divergent positives.
    assert op["threshold"] < op["production_threshold"]
    assert -1.25 < op["threshold"] < -1.05
    assert set(op["widened"]["includes_recovery"]) == {"PET46", "Cut190", "TfCut2"}
