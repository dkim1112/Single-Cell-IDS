#!/bin/bash
set -u
D="$HOME/mnt/Single-Cell-IDS/data/GSE194315_PBMC-01-07_processed_data_files"
OUT="$HOME/qcpbmc/pbmc_percell_qc.tsv"
for s in "$@"; do
  [ -s "$HOME/qcpbmc/done_$s" ] && { echo "skip $s"; continue; }
  zcat "$D/$s.features.tsv.gz" | awk -F'\t' '$2 ~ /^MT-/ {print NR}' > "$HOME/qcpbmc/mt_$s.txt"
  zcat "$D/$s.barcodes.tsv.gz" > "$HOME/qcpbmc/bc_$s.txt"
  zcat "$D/$s.matrix.mtx.gz" | awk -v S="$s" -v MTF="$HOME/qcpbmc/mt_$s.txt" -v BCF="$HOME/qcpbmc/bc_$s.txt" '
    BEGIN { while ((getline r < MTF) > 0) ismt[r+0]=1
            nb=0; while ((getline b < BCF) > 0) bc[++nb]=b }
    /^%/ { next }
    !hdr { hdr=1; next }
    { c=$2+0; tot[c]+=$3; det[c]++; if (ismt[$1+0]) mt[c]+=$3 }
    END { for (c=1; c<=nb; c++) { t=tot[c]+0; p=(t>0)?100.0*(mt[c]+0)/t:0
            printf "%s\t%s\t%d\t%d\t%.4f\n", S, bc[c], det[c]+0, t, p } }' >> "$OUT"
  rm -f "$HOME/qcpbmc/mt_$s.txt" "$HOME/qcpbmc/bc_$s.txt"
  touch "$HOME/qcpbmc/done_$s"
  echo "done $s"
done
