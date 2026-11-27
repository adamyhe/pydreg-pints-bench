#!/bin/bash -l
# Phase 1 data fetch. Runs on a login node (needs internet), downloads only; ~1.5 GB total.
# Idempotent: existing files are skipped. Review before running:
#   bash scripts/00_fetch_data.sh 2>&1 | tee logs/00_fetch_data.log
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$PROJ/data"
mkdir -p "$DATA"/{libs,ref/dhs,ref/cage,ref/rna,labels/dreg_ftp,labels/pints_supp} "$PROJ/logs"

get() {  # get URL DEST [extra curl args]; ENCODE rejects browser user agents, Springer requires one
    local url="$1" dest="$2"; shift 2
    if [[ -s "$dest" ]]; then echo "skip  $dest"; return; fi
    echo "fetch $url"
    curl -fsSL --retry 5 --retry-delay 10 "$@" -o "$dest.part" "$url"
    mv "$dest.part" "$dest"
}

# ---- run-on / TSS libraries (config/libraries.tsv, phase 1) ----
tail -n +2 "$PROJ/config/libraries.tsv" | while IFS=$'\t' read -r lib phase cell assay signal genome role exp pl mn notes; do
    [[ "$phase" == "1" ]] || continue
    get "$pl" "$DATA/libs/$lib.pl.bw"
    get "$mn" "$DATA/libs/$lib.mn.bw"
done

# ---- reference genome files ----
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath
get $UCSC/hg19/bigZips/hg19.chrom.sizes                     "$DATA/ref/hg19.chrom.sizes"
get $UCSC/hg38/bigZips/hg38.chrom.sizes                     "$DATA/ref/hg38.chrom.sizes"
get $UCSC/hg19/liftOver/hg19ToHg38.over.chain.gz            "$DATA/ref/hg19ToHg38.over.chain.gz"
get $UCSC/hg38/liftOver/hg38ToHg19.over.chain.gz            "$DATA/ref/hg38ToHg19.over.chain.gz"
get https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz "$DATA/ref/hg19-blacklist.v2.bed.gz"
get https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz "$DATA/ref/hg38-blacklist.v2.bed.gz"
# GENCODE v24 = the annotation PINTS used for promoters (GRCh38)
get https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz "$DATA/ref/gencode.v24.annotation.gtf.gz"

# ---- PINTS-paper label inputs ----
# Supplementary Tables 1-5 (Table 2 = known enhancers, Table 3 = non-enhancers)
get "https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01211-7/MediaObjects/41587_2022_1211_MOESM3_ESM.xlsx" \
    "$DATA/labels/pints_supp/41587_2022_1211_MOESM3_ESM.xlsx" -A "Mozilla/5.0"
# K562 nAnT-iCAGE (SRR7641288/9 in PINTS Supp Table 1 = GSE118075 CAGE reps), hg19 CTSS
GEO=https://ftp.ncbi.nlm.nih.gov/geo/samples
get $GEO/GSM3318nnn/GSM3318239/suppl/GSM3318239_CNhi10920_biologicalRep1-K562-CAGE_ACC.ctss.bed.gz "$DATA/ref/cage/K562_CAGE_rep1.ctss.bed.gz"
get $GEO/GSM3318nnn/GSM3318240/suppl/GSM3318240_CNhi10920_biologicalRep2-K562-CAGE_ATG.ctss.bed.gz "$DATA/ref/cage/K562_CAGE_rep2.ctss.bed.gz"
# K562 RSEM transcript quantifications, GRCh38 / GENCODE V24.
#   total RNA-seq ENCSR000AEN (= PINTS libraries ENCLB057ZZZ/058ZZZ); polyA ENCSR000CPH (CSHL, whole cell)
ENC=https://www.encodeproject.org/files
for acc in ENCFF345SBQ ENCFF010QAI; do get "$ENC/$acc/@@download/$acc.tsv" "$DATA/ref/rna/K562_total_$acc.tsv"; done
for acc in ENCFF764GTR ENCFF261DSR; do get "$ENC/$acc/@@download/$acc.tsv" "$DATA/ref/rna/K562_polyA_$acc.tsv"; done

# ---- dREG-paper label inputs ----
# Test-set labels as distributed on the Danko lab hub (used by TROGDOR); hg19
for f in k562.positive.bed.gz k562.negative.bed.gz K562.positive.bed.gz GM12878.positive.bed.gz; do
    get "ftp://cbsuftp.tc.cornell.edu/danko/hub/TROGDOR/$f" "$DATA/labels/dreg_ftp/$f"
done
# DHS used for dREG negatives (Duke + UW, hg19). K562 from UCSC (the files dREG's scripts name);
# GM12878 UW = the merged file dREG's scripts name (lab storage, read-only copy, approved);
# GM12878 Duke from UCSC.
ENCDCC=$UCSC/hg19/encodeDCC
get $ENCDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz   "$DATA/ref/dhs/K562_Duke.hg19.narrowPeak.gz"
get $ENCDCC/wgEncodeUwDgf/wgEncodeUwDgfK562Pk.narrowPeak.gz                      "$DATA/ref/dhs/K562_UW.hg19.narrowPeak.gz"
get $ENCDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseGm12878Pk.narrowPeak.gz "$DATA/ref/dhs/GM12878_Duke.hg19.narrowPeak.gz"
if [[ ! -s "$DATA/ref/dhs/GM12878_UW.hg19.narrowPeak.gz" ]]; then
    gzip -c /fs/cbsudanko/storage/data/hg19/gm12878/dnase/uw.merge.narrowPeak.bed > "$DATA/ref/dhs/GM12878_UW.hg19.narrowPeak.gz"
fi
# GM12878 GRO-cap (stands in for the GRO-cap HMM calls when rebuilding GM12878 negatives)
get $GEO/GSM1480nnn/GSM1480323/suppl/GSM1480323_GM12878_GROcap_wTAP_plus.bigWig  "$DATA/ref/GM12878_GROcap.pl.bw"
get $GEO/GSM1480nnn/GSM1480323/suppl/GSM1480323_GM12878_GROcap_wTAP_minus.bigWig "$DATA/ref/GM12878_GROcap.mn.bw"
# K562 GRO-cap HMM calls used by dREG (lab storage, read-only copy, approved) -- for the
# negative-rebuild validation on K562
if [[ ! -s "$DATA/ref/K562_grocap_hmm.hg19.bed" ]]; then
    cp /fs/cbsudanko/storage/archive/DANKO-2016/train_holdout_svm/hg19.k562.new_hmm2b.post2.bed "$DATA/ref/K562_grocap_hmm.hg19.bed"
fi

# ---- dREG auxiliary labels (Wang 2019 Fig 1B; Danko 2015 Fig 1b): ENCODE hg19 Broad H3K27ac + ChromHMM ----
mkdir -p "$DATA/ref/aux"
for cell in K562 Gm12878; do
    C=${cell^^}; [[ $C == GM12878 ]] || C=K562
    get $ENCDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${cell}H3k27acStdPk.broadPeak.gz "$DATA/ref/aux/${C}_H3K27ac.hg19.broadPeak.gz"
    get $ENCDCC/wgEncodeBroadHmm/wgEncodeBroadHmm${cell}HMM.bed.gz                        "$DATA/ref/aux/${C}_ChromHMM.hg19.bed.gz"
done

# ---- pydreg model weights (so GPU nodes need no internet) ----
cd "$PROJ"
uv run python - <<'EOF'
from huggingface_hub import hf_hub_download
from pydreg import models
for fn in (models.HF_MODEL_FILENAME, models.HF_RF_MODEL_FILENAME):
    print("cached", hf_hub_download(repo_id=models.HF_REPO_ID, filename=fn))
EOF

echo "done. sizes:"; du -sh "$DATA"/libs "$DATA"/ref "$DATA"/labels
