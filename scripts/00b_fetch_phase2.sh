#!/bin/bash -l
# Phase 2 data fetch (hg38: HCT116, A549, K562 PRO-cap). Login node; downloads + light merging only.
# Idempotent. Replicate bigWigs are kept in data/libs/raw/ and summed into data/libs/<lib>.{pl,mn}.bw
# (single-replicate libraries are copied).
#   bash scripts/00b_fetch_phase2.sh 2>&1 | tee logs/00b_fetch_phase2.log
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$PROJ/data"
source "$PROJ/scripts/site.sh"
KENT="$KENT_BIN"
ENC=https://www.encodeproject.org/files
mkdir -p "$DATA"/{libs/raw,ref/dhs,ref/aux} "$PROJ/logs"
TMP="${TMPDIR:-$DATA/libs/raw}"

get() {  # get URL DEST
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then echo "skip  $dest"; return; fi
    echo "fetch $url"
    curl -fsSL --retry 5 --retry-delay 10 -o "$dest.part" "$url"
    mv "$dest.part" "$dest"
}
enc() { get "$ENC/$1/@@download/$1.${2}" "$3"; }

# Sum replicate bigWigs and re-encode without explicit zeros (some ENCODE minus-strand files, and
# bigWigMerge output, store 0 over the whole genome). bigWigMerge drops values <= threshold, so the
# threshold is set far below 0 to keep the negative minus-strand values.
merge_bw() {  # merge_bw DEST SRC...
    local dest="$1"; shift
    if [[ -s "$dest" ]]; then echo "skip  $dest"; return; fi
    local bg; bg=$(mktemp "$TMP/merge.XXXXXX.bg")
    if [[ $# -eq 1 ]]; then
        $KENT/bigWigToBedGraph "$1" /dev/stdout | awk '$4 != 0' > "$bg"
    else
        $KENT/bigWigMerge -threshold=-1e30 "$@" /dev/stdout | awk '$4 != 0' > "$bg"
    fi
    LC_ALL=C sort -k1,1 -k2,2n "$bg" -o "$bg"
    $KENT/bedGraphToBigWig "$bg" "$DATA/ref/hg38.chrom.sizes" "$dest.part"
    mv "$dest.part" "$dest"; rm -f "$bg"
}

# ---- libraries (config/libraries.tsv, phase 2; pl/mn columns are comma-separated ENCFF accessions) ----
tail -n +2 "$PROJ/config/libraries.tsv" | while IFS=$'\t' read -r lib phase cell assay signal genome role exp pl mn notes; do
    [[ "$phase" == "2" ]] || continue
    for strand in pl mn; do
        accs=$([[ $strand == pl ]] && echo "$pl" || echo "$mn")
        srcs=()
        for a in ${accs//,/ }; do
            enc "$a" bigWig "$DATA/libs/raw/$lib.$strand.$a.bw"
            srcs+=("$DATA/libs/raw/$lib.$strand.$a.bw")
        done
        merge_bw "$DATA/libs/$lib.$strand.bw" "${srcs[@]}"
    done
done

# ---- hg38 label inputs ----
# ENCODE DNase-seq peaks (one experiment per cell; the dREG paper's Duke & UW pair has no hg38 analogue)
enc ENCFF240LRP bed.gz "$DATA/ref/dhs/HCT116_DNase.hg38.bed.gz"
enc ENCFF128ZVL bed.gz "$DATA/ref/dhs/A549_DNase.hg38.bed.gz"
enc ENCFF185XRG bed.gz "$DATA/ref/dhs/K562_DNase.hg38.bed.gz"
# ENCODE H3K27ac peaks: HCT116 ENCSR661KMA (Bernstein, replicated peaks); A549 ENCSR778NQS
# (Reddy, untreated, pseudoreplicated = ENCODE default)
enc ENCFF349LKU bed.gz "$DATA/ref/aux/HCT116_H3K27ac.hg38.narrowPeak.gz"
enc ENCFF747IZX bed.gz "$DATA/ref/aux/A549_H3K27ac.hg38.narrowPeak.gz"

# ---- convention check: phase-1 inputs are integer counts, minus strand negative ----
for f in "$DATA"/libs/{HCT116,A549,K562}_PRO*.bw; do
    echo "$(basename "$f"): $($KENT/bigWigInfo "$f" | grep -E '^(min|max|mean|basesCovered)' | tr -s ' ' | tr '\n' ' ')"
done

# ---- GM12878 PINTS-paper-style labels ----
# HiDRA (Wang 2018, GSE104001; GM12878 ATAC fragments in a STARR-seq vector), hg19 fragment-group counts
mkdir -p "$DATA/ref/hidra" "$DATA/ref/cage" "$DATA/ref/rna"
get https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104001/suppl/GSE104001_HiDRA_counts_per_fragmentgroup.txt.gz \
    "$DATA/ref/hidra/GSE104001_HiDRA_counts_per_fragmentgroup.txt.gz"
# GM12878 nAnT-iCAGE (same FANTOM5/PINTS series as K562: GSE118075), hg19 CTSS
GEO=https://ftp.ncbi.nlm.nih.gov/geo/samples
get $GEO/GSM3318nnn/GSM3318207/suppl/GSM3318207_CNhi10917_biologicalRep1-GM12878-CAGE_ACC.ctss.bed.gz "$DATA/ref/cage/GM12878_CAGE_rep1.ctss.bed.gz"
get $GEO/GSM3318nnn/GSM3318208/suppl/GSM3318208_CNhi10917_biologicalRep2-GM12878-CAGE_ATG.ctss.bed.gz "$DATA/ref/cage/GM12878_CAGE_rep2.ctss.bed.gz"
# GM12878 RSEM transcript quantifications, GRCh38 / GENCODE V24: total ENCSR000AEE, polyA ENCSR000COQ
for acc in ENCFF658DGX ENCFF462UGG; do enc $acc tsv "$DATA/ref/rna/GM12878_total_$acc.tsv"; done
for acc in ENCFF879KFK ENCFF087ZDU; do enc $acc tsv "$DATA/ref/rna/GM12878_polyA_$acc.tsv"; done
