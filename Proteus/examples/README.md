# Examples — a ready-to-fold S3 hand-off

A real worked output of the local narrowing pipeline, so you can run the GCE fold
(S3) without first standing up the full S0–S2 toolchain locally.

- **`demo_cutinases_shortlist.fasta`** — 15 reviewed cutinases (Swiss-Prot), the
  surviving shortlist from a live `fetch_corpus → corpus → S0 → S1 → S2` run on a
  mixed corpus (15 cutinases + 10 myoglobins). S2's α/β-hydrolase fold-class triage
  kept all 15 cutinases and dropped all the all-α myoglobins.
- **`demo_cutinases_s3_manifest.json`** — the matching S3 job manifest
  (`run_location: gce`, per-sequence sha256, fold params from `config/proteus.yaml`).

Fold it on GCE (after `compute.gce_burst.{bucket,image}` are set — see `gce/sync.md`):

```bash
PYTHONPATH=src python -m proteus.launch \
  --manifest examples/demo_cutinases_s3_manifest.json \
  --shortlist examples/demo_cutinases_shortlist.fasta            # dry-run plan
PYTHONPATH=src python -m proteus.launch \
  --manifest examples/demo_cutinases_s3_manifest.json \
  --shortlist examples/demo_cutinases_shortlist.fasta --execute  # up -> fold -> down -> delete
```

These are checked-in outputs (not test fixtures); regenerate your own with
`python -m proteus.pipeline` once your corpus + S2 reference DBs are configured.
