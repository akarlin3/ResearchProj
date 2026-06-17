#!/bin/bash
# CP1 preprocessing for the sensitivity re-search — memory-safe (streaming sort+awk),
# so a 20M-row / 9.8M-target m8 never lands in a Python dict. Emits, for the relaxed
# search and the original e<=0.01 baseline:
#   - per-target aggregate (acc, best_bits, max_fident%, best_query, fident_query)
#   - identity-band counts (max-PETASE fident, the bands Figure 1 found empty)
#   - NEW-vs-original target counts per band (set difference)
#   - selected.tsv : stratified top-by-bits screen list for CP2
#
# Usage: sensitivity_cp1.sh <relaxed.m8.gz> <orig.m8.gz> <outdir> <tmpdir>
#   strata edges + per-stratum via env: STRATA="0,12,18,25" PERSTRATUM=700
set -Eeuo pipefail
REL="$1"; ORIG="$2"; OUT="$3"; TMP="$4"
mkdir -p "$OUT" "$TMP"
STRATA="${STRATA:-0,12,18,25}"
PERSTRATUM="${PERSTRATUM:-700}"
PETASE_RE='^(6EQE|4EB0|8B4U|4WFI|4CG1)'

# Stream rows -> sort by target -> awk groups consecutive same-target rows, emitting
# per-target best_bits + max_fident (over the PETASE queries) in O(1) memory.
# m8 cols: query target(.pdb.gz) evalue bits fident alnlen qlen tlen
aggregate() {  # <m8.gz> <out.tsv>
  gunzip -c "$1" \
    | awk -F'\t' -v re="$PETASE_RE" '$1 ~ re' \
    | sort -t$'\t' -k2,2 -S 512M -T "$TMP" \
    | awk -F'\t' '
        function flush(){ if(t!=""){ printf "%s\t%d\t%.4f\t%s\t%s\n", t, bb, mf*100, bq, fq } }
        { acc=$2; sub(/\..*/,"",acc); q=substr($1,1,4);
          if(acc!=t){ flush(); t=acc; bb=-1; mf=-1; bq=""; fq="" }
          b=$4+0; f=$5+0;
          if(b>bb){bb=b; bq=q}
          if(f>mf){mf=f; fq=q} }
        END{ flush() }' > "$2"
}

echo "[cp1] aggregating relaxed -> per_target_relaxed.tsv"
aggregate "$REL" "$OUT/per_target_relaxed.tsv"
echo "[cp1] aggregating original -> per_target_orig.tsv"
aggregate "$ORIG" "$OUT/per_target_orig.tsv"

# Identity-band counts (max-PETASE fident %). Band edges match Figure 1.
bandcounts() {  # <per_target.tsv>
  awk -F'\t' '{f=$3;
      if(f<20)b="<20%"; else if(f<25)b="20-25%"; else if(f<30)b="25-30%";
      else if(f<40)b="30-40%"; else if(f<60)b="40-60%"; else b=">60%";
      c[b]++ }
    END{ split("<20% 20-25% 25-30% 30-40% 40-60% >60%",o," ");
      for(i=1;i<=6;i++) printf "%s\t%d\n", o[i], c[o[i]]+0 }' "$1"
}
echo "[cp1] relaxed band counts:"; bandcounts "$OUT/per_target_relaxed.tsv" | tee "$OUT/bins_relaxed.tsv"
echo "[cp1] original band counts:"; bandcounts "$OUT/per_target_orig.tsv" | tee "$OUT/bins_orig.tsv"

# NEW targets (relaxed minus original) overall + per band.
cut -f1 "$OUT/per_target_relaxed.tsv" | sort -S 256M -T "$TMP" > "$TMP/rel_accs.txt"
cut -f1 "$OUT/per_target_orig.tsv"    | sort -S 256M -T "$TMP" > "$TMP/orig_accs.txt"
comm -23 "$TMP/rel_accs.txt" "$TMP/orig_accs.txt" > "$TMP/new_accs.txt"
NEW_TOTAL=$(wc -l < "$TMP/new_accs.txt")
echo "[cp1] NEW targets vs original: $NEW_TOTAL"
# per-band counts of the NEW targets
LC_ALL=C join -t$'\t' -1 1 -2 1 \
    <(sort -S 256M -T "$TMP" "$TMP/new_accs.txt") \
    <(sort -t$'\t' -k1,1 -S 256M -T "$TMP" "$OUT/per_target_relaxed.tsv") \
  | bandcounts /dev/stdin > "$OUT/bins_new.tsv"
echo "[cp1] NEW-by-band:"; cat "$OUT/bins_new.tsv"

# Stratified selection: top --PERSTRATUM by bits within each fident band -> selected.tsv
IFS=',' read -r -a E <<< "$STRATA"
: > "$OUT/selected.tsv"
for ((i=0; i<${#E[@]}-1; i++)); do
  lo=${E[i]}; hi=${E[i+1]}
  navail=$(awk -F'\t' -v lo="$lo" -v hi="$hi" '$3>=lo && $3<hi' "$OUT/per_target_relaxed.tsv" | wc -l)
  awk -F'\t' -v lo="$lo" -v hi="$hi" '$3>=lo && $3<hi' "$OUT/per_target_relaxed.tsv" \
    | sort -t$'\t' -k2,2 -nr -S 256M -T "$TMP" | head -n "$PERSTRATUM" >> "$OUT/selected.tsv"
  echo "[cp1] stratum ${lo}-${hi}%: $navail available, took top $PERSTRATUM by bits"
done
# de-dup selected (a target can only be in one band, but guard anyway)
sort -u -t$'\t' -k1,1 "$OUT/selected.tsv" -o "$OUT/selected.tsv"
echo "[cp1] selected $(wc -l < "$OUT/selected.tsv") targets -> $OUT/selected.tsv"

# CP1 numbers as JSON for the report.
SUBCUT=$(awk -F'\t' '$3<25' "$OUT/per_target_relaxed.tsv" | wc -l)
REL_TOT=$(wc -l < "$OUT/per_target_relaxed.tsv")
python3 - "$OUT" "$NEW_TOTAL" "$SUBCUT" "$REL_TOT" <<'PY'
import json, sys, os
out, new_total, subcut, rel_tot = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
def rd(p):
    return [{"bin": l.split("\t")[0], "n": int(l.split("\t")[1])}
            for l in open(p).read().splitlines() if l.strip()]
json.dump({"n_unique_targets": rel_tot, "n_sub_cut_total": subcut,
           "n_new_targets_vs_original": new_total,
           "bins_fident": rd(f"{out}/bins_relaxed.tsv"),
           "orig_bins_fident": rd(f"{out}/bins_orig.tsv"),
           "new_bins_fident": rd(f"{out}/bins_new.tsv")},
          open(f"{out}/cp1.json","w"), indent=2)
print("[cp1] wrote", f"{out}/cp1.json")
PY
echo "[cp1] DONE"
