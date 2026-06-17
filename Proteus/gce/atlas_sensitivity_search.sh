#!/bin/bash
# Atlas sensitivity re-search (GCE startup-script).
#
# Re-search the PETase/cutinase query structures against the prebuilt ESMAtlas
# highquality_clust30 Foldseek DB at MAXIMUM sensitivity with RELAXED bits/e-value,
# to reach the sub-25%-identity structural tail the original sweep (-e 0.01,
# --max-seqs 300000, top-300 fetch budget) never retrieved. The original sweep was
# ALREADY -s 9.5; the new levers here are -e (relaxed 10^6x) and --max-seqs.
#
# Runs unattended on VM boot, stages result_sensitivity.m8.gz + provenance to GCS,
# and writes a STATUS marker (running|done|failed) the launcher polls. The VM has a
# 3h max-run-duration + DELETE termination as a billing backstop; the launcher also
# deletes it on STATUS=done.
#
# Downstream (partition / identity-binning / fetch / screen at -1.1587 / seqid) runs
# LOCALLY off this m8 — the live fetchPredictedStructure API + fpocket are on the Mac.
set -Eeuo pipefail

GCS="gs://projproteus-fold/sensitivity/2026-06-14"
LOG=/var/log/sensitivity.log
exec > >(tee -a "$LOG") 2>&1

# Relaxed-search parameters (the experiment's knobs).
SENS=9.5            # -s : max Foldseek sensitivity (same as the original sweep)
MAXSEQS=4000000     # --max-seqs : lift the prefilter cap far above the original 300000
EVALUE=10000        # -e : relax 10^6x vs the original 0.01 (min_bits effectively 0)
ALNTYPE=2           # 3Di+AA local alignment (no _ca C-alpha set)

echo "===== Atlas sensitivity re-search START $(date -u) ====="
echo "params: -s $SENS --max-seqs $MAXSEQS -e $EVALUE --alignment-type $ALNTYPE"

# Guarantee the Cloud SDK (gsutil) before any staging — Ubuntu GCE images usually
# ship it as a snap, but don't assume.
if ! command -v gsutil >/dev/null 2>&1; then
  echo "gsutil absent — installing Cloud SDK"
  snap install google-cloud-cli --classic 2>/dev/null || {
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -y
    apt-get install -y apt-transport-https ca-certificates gnupg curl
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
      | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
      > /etc/apt/sources.list.d/google-cloud-sdk.list
    apt-get update -y && apt-get install -y google-cloud-cli
  }
fi
command -v gsutil

on_err() {
  rc=$?
  echo "===== FAILED rc=$rc at line ${BASH_LINENO[0]} $(date -u) ====="
  gsutil cp "$LOG" "$GCS/sensitivity-search.log" || true
  printf 'failed %s rc=%s line=%s\n' "$(date -u)" "$rc" "${BASH_LINENO[0]}" \
    | gsutil cp - "$GCS/STATUS" || true
}
trap on_err ERR

printf 'running %s\n' "$(date -u)" | gsutil cp - "$GCS/STATUS"

# ---- tooling -------------------------------------------------------------- #
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y wget aria2 pigz curl coreutils
cd /opt
wget -q https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xzf foldseek-linux-avx2.tar.gz
export PATH=/opt/foldseek/bin:$PATH
FOLDSEEK_VER=$(foldseek version 2>/dev/null || echo unknown)
echo "foldseek $FOLDSEEK_VER ; nproc=$(nproc) ; mem=$(free -g | awk '/Mem:/{print $2}')G"
df -h / | tail -1

# ---- prebuilt 3Di+AA DB (skip the 88 GB _ca) ------------------------------ #
mkdir -p /data/db && cd /data/db
MAN=https://raw.githubusercontent.com/facebookresearch/esm/main/scripts/atlas/v0/highquality_clust30/foldseekdb.txt
curl -s "$MAN" | grep -v '_ca' > manifest.txt
echo "DB parts to fetch:"; cat manifest.txt
aria2c -x16 -s16 -j3 --retry-wait=5 --max-tries=8 -i manifest.txt
echo "--- DB files ---"; ls -la; du -sh /data/db
DB=/data/db/highquality_clust30
DB_ENTRIES=$(wc -l < "${DB}.index")
echo "DB entries (.index lines): $DB_ENTRIES"

# ---- query DB from the 5 PETase/cutinase anchors -------------------------- #
mkdir -p /data/q && cd /data
gsutil -m cp "$GCS/queries/"*.pdb /data/q/
ls -la /data/q
foldseek createdb /data/q/6EQE.pdb /data/q/4EB0.pdb /data/q/8B4U.pdb \
                  /data/q/4WFI.pdb /data/q/4CG1.pdb /data/querydb
echo "--- query DB headers ---"; foldseek prefixid /data/querydb_h /dev/stdout 2>/dev/null | head || true

# ---- relaxed search + convertalis (8-col, matches the original result.m8) -- #
mkdir -p /data/tmp
echo "===== foldseek search $(date -u) ====="
/usr/bin/time -v foldseek search /data/querydb "$DB" /data/result /data/tmp \
  -s "$SENS" --max-seqs "$MAXSEQS" -e "$EVALUE" --alignment-type "$ALNTYPE" \
  --threads "$(nproc)" 2>&1 | tail -40
echo "===== convertalis $(date -u) ====="
foldseek convertalis /data/querydb "$DB" /data/result /data/result_sensitivity.m8 \
  --format-output "query,target,evalue,bits,fident,alnlen,qlen,tlen" \
  --threads "$(nproc)"

cd /data
ROWS=$(wc -l < result_sensitivity.m8)
UNIQ=$(cut -f2 result_sensitivity.m8 | sort -u | wc -l)
echo "rows=$ROWS unique_targets=$UNIQ"
echo "--- per-query rows ---"; cut -f1 result_sensitivity.m8 | sort | uniq -c
echo "--- fident histogram (per alignment row, all queries) ---"
awk -F'\t' '{b=int($5*100/5)*5; h[b]++} END{for(k=0;k<=100;k+=5) printf "  [%3d-%3d%%) %d\n",k,k+5,h[k]+0}' \
  result_sensitivity.m8
echo "--- rows below the original gates (NEW retrieval) ---"
awk -F'\t' '$3>0.01{e++} $4<50{b++} END{printf "  e>0.01: %d ; bits<50: %d\n", e+0, b+0}' \
  result_sensitivity.m8

# ---- provenance + push ---------------------------------------------------- #
cat > provenance.json <<EOF
{
  "run_date_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "host": "$(hostname) $(curl -s -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/machine-type 2>/dev/null | awk -F/ '{print $NF}')",
  "foldseek_version": "$FOLDSEEK_VER",
  "db": "ESMAtlas highquality_clust30 prebuilt 3Di+AA Foldseek DB",
  "db_entries": $DB_ENTRIES,
  "queries": ["6EQE","4EB0","8B4U","4WFI","4CG1"],
  "search_params": {"s": $SENS, "max_seqs": $MAXSEQS, "e": $EVALUE, "alignment_type": $ALNTYPE},
  "original_sweep_params": {"s": 9.5, "max_seqs": 300000, "e": 0.01, "min_bits": 50},
  "result_rows": $ROWS,
  "unique_targets": $UNIQ
}
EOF
cat provenance.json
pigz -k -f result_sensitivity.m8
gsutil cp result_sensitivity.m8.gz "$GCS/"
gsutil cp provenance.json "$GCS/"
gsutil cp "$LOG" "$GCS/sensitivity-search.log"
printf 'done %s rows=%s uniq=%s\n' "$(date -u)" "$ROWS" "$UNIQ" | gsutil cp - "$GCS/STATUS"
echo "===== DONE $(date -u) — result_sensitivity.m8.gz staged to $GCS ====="
