#!/bin/bash -l
# Phase 3 data fetch: Jurkat (hg19). Login node; downloads + light format conversion only. Idempotent.
#   bash scripts/00c_fetch_jurkat.sh 2>&1 | tee logs/00c_fetch_jurkat.log
#
# Libraries   GSE117832 (Chu et al. 2018, ChRO-seq; Danko-lab proseq pipeline, hg19 bigWigs):
#               GSM3309955 Jurkat nuclei PRO-seq; GSM3309956 ChRO-seq batch 1; GSM3309957 ChRO-seq
#               batch 2 (no RNase); GSM3309958 leChRO-seq (length-extension, RNase)
#             FANTOM5 Jurkat hCAGE CNhs11253 (hg19 CTSS) -> pl/mn bigWigs (minus negative)
# Labels      UW DNase Jurkat reps 1+2 (ENCODE hg19, UCSC wgEncodeUwDnase)
#             H3K27ac reps 1+2 (GSM4551183/4 = SRX8336404/5; ChIP-Atlas MACS2 q<1e-5 peaks, hg19)
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$PROJ/data"
KENT=/programs/kentUtils/bin
mkdir -p "$DATA"/{libs/raw,ref/dhs,ref/aux,ref/cage} "$PROJ/logs"
TMP="${TMPDIR:-$DATA/libs/raw}"

get() {
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then echo "skip  $dest"; return; fi
    echo "fetch $url"
    curl -fsSL --retry 5 --retry-delay 10 -o "$dest.part" "$url"
    mv "$dest.part" "$dest"
}
# re-encode a bigWig without explicit zeros (same treatment as phase 2)
clean_bw() {  # clean_bw SRC DEST
    if [[ -s "$2" ]]; then echo "skip  $2"; return; fi
    local bg; bg=$(mktemp "$TMP/clean.XXXXXX.bg")
    $KENT/bigWigToBedGraph "$1" /dev/stdout | awk '$4 != 0' | LC_ALL=C sort -k1,1 -k2,2n > "$bg"
    $KENT/bedGraphToBigWig "$bg" "$DATA/ref/hg19.chrom.sizes" "$2.part"
    mv "$2.part" "$2"; rm -f "$bg"
}

GEO=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3309nnn
declare -A LIB=(
    [Jurkat_PROseq]="GSM3309955/GSM3309955_5587_5598_24204_HGC2FBGXX_J_NUC_TTAGGC_R1"
    [Jurkat_ChROseq1]="GSM3309956/GSM3309956_5587_5598_24205_HGC2FBGXX_J_CHR_TGACCA_R1"
    [Jurkat_ChROseq2]="GSM3309957/GSM3309957_Jurkat_ChRO_NoRNase"
    [Jurkat_leChROseq]="GSM3309958/GSM3309958_Jurkat_ChRO_RNase"
)
for lib in "${!LIB[@]}"; do
    for s in plus minus; do
        short=$([[ $s == plus ]] && echo pl || echo mn)
        get "$GEO/${LIB[$lib]%%/*}/suppl/${LIB[$lib]#*/}_$s.bw" "$DATA/libs/raw/$lib.$short.geo.bw"
        clean_bw "$DATA/libs/raw/$lib.$short.geo.bw" "$DATA/libs/$lib.$short.bw"
    done
done

# FANTOM5 Jurkat CAGE (the directory listing double-encodes the file name)
F5=https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE
get "$F5/acute%2520lymphoblastic%2520leukemia%2520%2528T-ALL%2529%2520cell%2520line%253aJurkat.CNhs11253.10464-106H5.hg19.ctss.bed.gz" \
    "$DATA/ref/cage/Jurkat_CAGE_CNhs11253.hg19.ctss.bed.gz"
for s in pl mn; do
    dest="$DATA/libs/Jurkat_CAGE.$s.bw"
    [[ -s "$dest" ]] && { echo "skip  $dest"; continue; }
    bg=$(mktemp "$TMP/ctss.XXXXXX.bg")
    zcat "$DATA/ref/cage/Jurkat_CAGE_CNhs11253.hg19.ctss.bed.gz" \
        | awk -v s=$s 'BEGIN{OFS="\t"} (s=="pl" && $6=="+") {print $1,$2,$3,$5} (s=="mn" && $6=="-") {print $1,$2,$3,-$5}' \
        | awk 'NR==FNR {ok[$1]=1; next} ($1 in ok)' "$DATA/ref/hg19.chrom.sizes" - \
        | LC_ALL=C sort -k1,1 -k2,2n > "$bg"
    $KENT/bedGraphToBigWig "$bg" "$DATA/ref/hg19.chrom.sizes" "$dest.part"; mv "$dest.part" "$dest"; rm -f "$bg"
done

# labels
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC
for r in 1 2; do
    get $UCSC/wgEncodeUwDnase/wgEncodeUwDnaseJurkatPkRep$r.narrowPeak.gz "$DATA/ref/dhs/Jurkat_UW_rep$r.hg19.narrowPeak.gz"
done
CA=https://chip-atlas.dbcls.jp/data/hg19/eachData/bed05
get $CA/SRX8336404.05.bed "$DATA/ref/aux/Jurkat_H3K27ac_rep1.hg19.bed"
get $CA/SRX8336405.05.bed "$DATA/ref/aux/Jurkat_H3K27ac_rep2.hg19.bed"

for f in "$DATA"/libs/Jurkat_*.bw; do
    echo "$(basename "$f"): $($KENT/bigWigInfo "$f" | grep -E '^(min|max|mean|basesCovered)' | tr -s ' ' | tr '\n' ' ')"
done
