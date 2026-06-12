"""GCE burst launcher tests (proteus.launch).

The launcher only PLANS (and validates) the burst by default — `gcloud`/`gsutil`
are never invoked unless --execute is passed. So the plan builder and input
validation are pure and fully testable here; no GCP, no network.
"""
from __future__ import annotations

import copy
import json

import pytest

from proteus import launch
from proteus.utils import load_config


def _cfg(project="my-proj", bucket="gs://my-bucket", image="us-docker.pkg.dev/p/r/fold:cu124"):
    cfg = copy.deepcopy(load_config())
    cfg.setdefault("compute", {}).setdefault("gce_burst", {})
    cfg["compute"]["gce_burst"].update(project=project, bucket=bucket, image=image)
    return cfg


def _manifest(tmp_path, n=3, run_location="gce"):
    p = tmp_path / "s3_job_manifest.json"
    p.write_text(json.dumps({"run_location": run_location, "n_sequences": n,
                             "sequences": [{"id": f"c{i}"} for i in range(n)]}))
    return str(p)


def _shortlist(tmp_path, n=3):
    p = tmp_path / "s2_shortlist.fasta"
    p.write_text("".join(f">c{i}\nMKLV\n" for i in range(n)))
    return str(p)


def test_validate_accepts_coherent_inputs(tmp_path):
    v = launch.validate_inputs(_manifest(tmp_path, 3), _shortlist(tmp_path, 3))
    assert v["ok"] is True and v["errors"] == []
    assert v["manifest"]["n_sequences"] == 3


def test_validate_flags_count_mismatch(tmp_path):
    v = launch.validate_inputs(_manifest(tmp_path, 5), _shortlist(tmp_path, 3))
    assert v["ok"] is False
    assert any("shortlist has 3" in e and "manifest says 5" in e for e in v["errors"])


def test_validate_refuses_non_gce_manifest(tmp_path):
    v = launch.validate_inputs(_manifest(tmp_path, 3, run_location="vast"),
                               _shortlist(tmp_path, 3))
    assert v["ok"] is False
    assert any("run_location" in e and "gce" in e for e in v["errors"])


def test_validate_missing_files(tmp_path):
    v = launch.validate_inputs(str(tmp_path / "nope.json"), str(tmp_path / "nope.fasta"))
    assert v["ok"] is False
    assert any("manifest not found" in e for e in v["errors"])


def test_build_plan_orders_the_burst_cycle(tmp_path):
    man, fa = _manifest(tmp_path), _shortlist(tmp_path)
    plan = launch.build_plan(_cfg(), man, fa)
    names = [s["name"] for s in plan]
    assert names == ["stage_up", "create_instance", "fold", "stage_down", "delete_instance"]
    by = {s["name"]: s["cmd"] for s in plan}

    # stage_up ships ONLY the shortlist + manifest to the bucket's in/ prefix
    assert by["stage_up"][:3] == ["gsutil", "-m", "cp"]
    assert man in by["stage_up"] and fa in by["stage_up"]
    assert by["stage_up"][-1] == "gs://my-bucket/in/"
    # create: gcloud compute instances create, with project/zone/accelerator + SPOT
    assert by["create_instance"][:4] == ["gcloud", "compute", "instances", "create"]
    assert "my-proj" in by["create_instance"]
    assert "type=nvidia-l4,count=1" in by["create_instance"]
    assert "SPOT" in by["create_instance"]
    # fold: gcloud ssh runs the fold container against the staged manifest/fasta
    assert by["fold"][:3] == ["gcloud", "compute", "ssh"]
    remote = by["fold"][-1]
    assert "docker run --gpus all" in remote
    assert "us-docker.pkg.dev/p/r/fold:cu124" in remote
    assert "/data/proteus/in/s3_job_manifest.json" in remote
    assert "gs://my-bucket/in/*" in remote and "gs://my-bucket/out/" in remote
    # stage_down pulls results from the bucket to the local folded dir; delete cleans up
    assert by["stage_down"][-2] == "gs://my-bucket/out/*"
    assert by["stage_down"][-1].endswith("structures/folded/")
    assert by["delete_instance"][:4] == ["gcloud", "compute", "instances", "delete"]
    assert "--quiet" in by["delete_instance"]


def test_build_plan_uses_placeholders_when_unset(tmp_path):
    # default config has project/bucket/image blank -> clearly-flagged placeholders
    plan = launch.build_plan(load_config(), _manifest(tmp_path), _shortlist(tmp_path))
    by = {s["name"]: s["cmd"] for s in plan}
    assert any("PROJECT" in tok for tok in by["create_instance"])
    assert any("BUCKET" in tok for tok in by["stage_up"])
    assert any("IMAGE" in tok for tok in by["fold"])


def test_build_plan_on_demand_drops_spot(tmp_path):
    cfg = _cfg()
    cfg["compute"]["gce_burst"]["spot"] = False
    plan = launch.build_plan(cfg, _manifest(tmp_path), _shortlist(tmp_path))
    create = next(s["cmd"] for s in plan if s["name"] == "create_instance")
    assert "SPOT" not in create


def test_main_dry_run_does_not_execute(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(launch.subprocess, "run",
                        lambda *a, **k: pytest.fail("dry-run must not execute commands"))
    rc = launch.main(["--manifest", _manifest(tmp_path), "--shortlist", _shortlist(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "DRY-RUN" in out and "gcloud compute instances create" in out


def test_main_returns_error_on_bad_inputs(tmp_path):
    rc = launch.main(["--manifest", _manifest(tmp_path, 5),
                      "--shortlist", _shortlist(tmp_path, 3)])
    assert rc == 2
