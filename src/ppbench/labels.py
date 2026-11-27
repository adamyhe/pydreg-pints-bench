"""Build the static (tool-independent) label sets for phase 1.

PINTS-paper labels (Yao et al. 2022, Methods "ROC" + PINTS_analysis/tools_sensitivity_specificity.py):
  enhancers   Supp. Table 2, K562 rows supported by CRISPR/CRISPRi (803 elements, hg38)
  non-enh     Supp. Table 3 (7,097 elements, hg38) -> >=250 bp, merged, minus promoters
              (TSS +/-500 bp of every GENCODE v24 transcript)
  promoters   GENCODE v24 protein-coding TSS +/-500 bp with divergent K562 CAGE within +/-1 kb,
              FPKM >= 1 in both total and polyA RNA-seq, not in an overlapping-gene pair.
              Sampled per library to the number of enhancers (done at evaluation time, after the
              per-assay read filter, as PINTS did).
  Deviations (documented in BENCHMARK.md): CAGE window uses the stated intent (both strands within
  TSS +/-1 kb) rather than the code's TSS-500-anchored offsets; OverGeneDB is approximated by
  GENCODE gene-body overlaps; RNA-seq values are ENCODE RSEM FPKM (mean of 2 reps).

dREG-paper labels (Wang et al. 2019, Methods; Danko-Lab/dREG GR_submit_2018/fig-S2-pr-roc.R):
  K562 positive = Duke DHS & UW DHS & GRO-cap HMM (k562.positive.bed.gz, hub copy)
  K562 negative = not Duke DHS, not UW DHS, not GRO-cap HMM (k562.negative.bed.gz, hub copy)
  GM12878 positive = hub copy; GM12878 negative rebuilt with the same recipe, with nonzero GRO-cap
  signal standing in for the (unavailable) GM12878 GRO-cap HMM calls. The recipe is validated by
  rebuilding the K562 negatives and comparing them to the originals.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import bioframe as bf
import numpy as np
import pandas as pd

from .genome import BED3, DATA, REF, clean, genome_view, kent, liftover, merge, overlaps_any, read_bed, write_bed

LABELS = DATA / "labels"
BUILT = LABELS / "built"
SUPP = LABELS / "pints_supp" / "41587_2022_1211_MOESM3_ESM.xlsx"
GTF = REF / "gencode.v24.annotation.gtf.gz"
BW2BG = kent("bigWigToBedGraph")


# ------------------------------------------------------------------ helpers
def _supp_table(sheet: str) -> pd.DataFrame:
    raw = pd.read_excel(SUPP, sheet_name=sheet, header=None)
    hdr = raw.index[raw.iloc[:, 0].astype(str).str.strip() == "Chromosome"][0]
    df = raw.iloc[hdr + 1:].copy()
    df.columns = [str(c).strip() for c in raw.iloc[hdr]]
    df = df.dropna(subset=["Chromosome"])
    out = pd.DataFrame({"chrom": df["Chromosome"].astype(str),
                        "start": df["Start"].astype(np.int64), "end": df["End"].astype(np.int64)})
    for c in df.columns[3:]:
        if c != "nan":
            out[c] = df[c].values
    return out


def _gtf_transcripts() -> pd.DataFrame:
    gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None,
                      names=["chrom", "src", "feature", "start", "end", "score", "strand", "frame", "attr"],
                      usecols=["chrom", "feature", "start", "end", "strand", "attr"])
    gtf = gtf[gtf.feature.isin(["transcript", "gene"])].copy()
    for k in ("gene_id", "transcript_id", "gene_type", "gene_name"):
        gtf[k] = gtf.attr.str.extract(rf'{k} "([^"]+)"', expand=False)
    gtf["start"] -= 1
    return gtf.drop(columns="attr")


def _tss_windows(tx: pd.DataFrame, half: int = 500) -> pd.DataFrame:
    tss = np.where(tx.strand == "+", tx.start, tx.end - 1)
    out = tx.copy()
    out["tss"] = tss
    out["start"] = (tss - half).clip(min=0)
    out["end"] = tss + half
    return out


def _ctss_hg38(paths) -> pd.DataFrame:
    parts = [read_bed(p, names=BED3 + ["name", "count", "strand"]) for p in paths]
    df = pd.concat(parts, ignore_index=True)
    df["count"] = df["count"].astype(float)
    return liftover(df, "hg19", "hg38")


def _rna_fpkm(paths) -> pd.Series:
    tabs = [pd.read_csv(p, sep="\t", usecols=["transcript_id", "FPKM"]).set_index("transcript_id")["FPKM"]
            for p in paths]
    return pd.concat(tabs, axis=1).mean(axis=1)


def _nonzero_intervals(bw: Path, pad: int) -> pd.DataFrame:
    with tempfile.NamedTemporaryFile(suffix=".bg") as tf:
        subprocess.run([BW2BG, str(bw), tf.name], check=True)
        bg = pd.read_csv(tf.name, sep="\t", header=None, names=BED3 + ["v"])
    bg = bg[bg.v != 0]
    bg["start"] = (bg.start - pad).clip(lower=0)
    bg["end"] = bg.end + pad
    return merge(bg)


def _save(df, family, cell, name, genome, manifest):
    p = BUILT / family / cell / f"{name}.{genome}.bed.gz"
    write_bed(df, p)
    manifest[f"{family}/{cell}/{name}.{genome}"] = {"n": int(len(df)), "bp": int((df.end - df.start).sum())}
    return p


# ------------------------------------------------------------------ PINTS-paper labels
def _promoter_pool(cell: str, tx: pd.DataFrame, tx_all: pd.DataFrame) -> pd.DataFrame:
    """PINTS promoter pool (hg38): protein-coding TSS +/-500 with CAGE on both strands within TSS +/-1 kb,
    FPKM >= 1 in both total and polyA RNA-seq (mean of reps), not in an overlapping-gene pair."""
    pc = tx_all[tx_all.gene_type == "protein_coding"]
    prom = _tss_windows(pc)
    cage = _ctss_hg38(sorted((REF / "cage").glob(f"{cell}_CAGE_rep*.ctss.bed.gz")))
    win = prom[BED3].copy()
    win["start"] = (prom.tss - 1000).clip(lower=0)
    win["end"] = prom.tss + 1000
    pl = bf.count_overlaps(win, cage[cage.strand == "+"][BED3])["count"].to_numpy()
    mn = bf.count_overlaps(win, cage[cage.strand == "-"][BED3])["count"].to_numpy()
    prom = prom[(pl > 0) & (mn > 0)]
    total = _rna_fpkm(sorted((REF / "rna").glob(f"{cell}_total_*.tsv")))
    polya = _rna_fpkm(sorted((REF / "rna").glob(f"{cell}_polyA_*.tsv")))
    prom = prom[(prom.transcript_id.map(total) >= 1) & (prom.transcript_id.map(polya) >= 1)]
    genes = tx[(tx.feature == "gene") & tx.gene_type.isin(["protein_coding", "lincRNA", "antisense"])]
    g = genes[BED3 + ["gene_name"]].reset_index(drop=True)
    ov = bf.overlap(g, g, suffixes=("", "_2"))
    ov = ov[ov.gene_name != ov.gene_name_2]
    overlapping = set(ov.gene_name)
    prom = prom[~prom.gene_name.isin(overlapping)]
    prom = prom[BED3 + ["transcript_id", "gene_name", "strand"]]

    return prom


def build_pints_k562(manifest: dict) -> None:
    enh_all = _supp_table("SupplementaryTable2")
    k562 = enh_all[enh_all["Cell line"] == "K562"]
    crispr = k562[k562["Supported by"].str.contains("CRISPR")]
    nonenh = _supp_table("SupplementaryTable3")

    tx = _gtf_transcripts()
    tx_all = tx[tx.feature == "transcript"]
    prom_any = _tss_windows(tx_all)  # PINTS: promoters of every transcript, used to clean negatives

    # sanity check on the genome build of the released tables: PINTS removed promoter-overlapping
    # non-enhancers in hg38, so few released non-enhancers should touch an hg38 GENCODE promoter
    frac_prom = float(overlaps_any(nonenh, prom_any).mean())
    manifest["check/nonenh_frac_overlapping_hg38_promoters"] = round(frac_prom, 4)

    ne = nonenh[(nonenh.end - nonenh.start) >= 250]
    ne = merge(bf.sort_bedframe(ne))
    ne = ne[~overlaps_any(ne, prom_any)]

    prom = _promoter_pool("K562", tx, tx_all)

    for genome in ("hg38", "hg19"):
        def L(df):
            return clean(liftover(df, "hg38", genome), genome)
        _save(L(crispr[BED3]), "pints", "K562", "enhancers_crispr", genome, manifest)
        _save(L(k562[BED3]), "pints", "K562", "enhancers_all", genome, manifest)
        _save(L(ne), "pints", "K562", "nonenhancers", genome, manifest)
        _save(L(prom), "pints", "K562", "promoter_pool", genome, manifest)
        _save(L(prom_any[BED3]), "pints", "K562", "promoters_any_tx", genome, manifest)


# ------------------------------------------------------------------ dREG-paper labels
def _dreg_negative(dhs_files, transcribed: pd.DataFrame, positives: pd.DataFrame, genome="hg19"):
    excl = [read_bed(f)[BED3] for f in dhs_files] + [transcribed[BED3], positives[BED3]]
    excl = merge(bf.sort_bedframe(pd.concat(excl, ignore_index=True)))
    return bf.subtract(genome_view(genome), excl)[BED3]


def _jaccard(a: pd.DataFrame, b: pd.DataFrame) -> float:
    a, b = merge(bf.sort_bedframe(a)), merge(bf.sort_bedframe(b))
    inter = bf.overlap(a, b, how="inner", return_overlap=True)
    i = int((inter["overlap_end"] - inter["overlap_start"]).sum())
    u = int((a.end - a.start).sum() + (b.end - b.start).sum()) - i
    return i / u if u else float("nan")


def build_dreg(manifest: dict) -> None:
    hub = LABELS / "dreg_ftp"
    dhs = REF / "dhs"
    k_pos = merge(bf.sort_bedframe(read_bed(hub / "k562.positive.bed.gz")[BED3]))
    k_neg = read_bed(hub / "k562.negative.bed.gz")[BED3]
    k_pos_hub2 = read_bed(hub / "K562.positive.bed.gz")[BED3]
    manifest["check/k562_positive_vs_K562.positive_jaccard"] = round(_jaccard(k_pos, k_pos_hub2), 4)

    # validate the negative-rebuild recipe on K562 (original GRO-cap HMM, then GRO-cap-signal proxy)
    k_dhs = [dhs / "K562_Duke.hg19.narrowPeak.gz", dhs / "K562_UW.hg19.narrowPeak.gz"]
    hmm_f = REF / "K562_grocap_hmm.hg19.bed"  # not public; optional (config/site.env DREG_K562_HMM_BED)
    proxy = _nonzero_intervals(DATA / "libs" / "K562_GROcap.pl.bw", 250)
    proxy = merge(bf.sort_bedframe(pd.concat([proxy, _nonzero_intervals(DATA / "libs" / "K562_GROcap.mn.bw", 250)])))
    rebuilt_proxy = _dreg_negative(k_dhs, proxy, k_pos)
    if hmm_f.exists():
        rebuilt_hmm = _dreg_negative(k_dhs, read_bed(hmm_f)[BED3], k_pos)
        manifest["check/k562_negative_rebuilt_hmm_jaccard"] = round(_jaccard(k_neg, rebuilt_hmm), 4)
    manifest["check/k562_negative_rebuilt_grocap_signal_jaccard"] = round(_jaccard(k_neg, rebuilt_proxy), 4)

    g_pos = merge(bf.sort_bedframe(read_bed(hub / "GM12878.positive.bed.gz")[BED3]))
    g_dhs = [dhs / "GM12878_Duke.hg19.narrowPeak.gz", dhs / "GM12878_UW.hg19.narrowPeak.gz"]
    g_proxy = merge(bf.sort_bedframe(pd.concat([_nonzero_intervals(REF / "GM12878_GROcap.pl.bw", 250),
                                                _nonzero_intervals(REF / "GM12878_GROcap.mn.bw", 250)])))
    g_neg = _dreg_negative(g_dhs, g_proxy, g_pos)

    # Phase 1 libraries are all hg19 (the labels' native build). Negative *regions* are multi-kb
    # complements that liftOver cannot map at minMatch 0.95; phase 2 lifts them after tiling.
    g = "hg19"
    _save(clean(k_pos, g), "dreg", "K562", "positive", g, manifest)
    _save(clean(k_neg, g, "subtract"), "dreg", "K562", "negative", g, manifest)
    _save(clean(rebuilt_proxy, g, "subtract"), "dreg", "K562", "negative_rebuilt", g, manifest)
    _save(clean(g_pos, g), "dreg", "GM12878", "positive", g, manifest)
    _save(clean(g_neg, g, "subtract"), "dreg", "GM12878", "negative_rebuilt", g, manifest)


# ------------------------------------------------------------------ dREG auxiliary labels
# Wang 2019 Fig 1B: recall of "DHSs marked by H3K27ac" (Duke & UW DHS & ENCODE H3K27ac). Danko 2015
# Fig 1b: recall of ChromHMM promoters/enhancers; its tuning negatives excluded DHS and ChromHMM
# promoter/enhancer/insulator states. Neither paper paired these with negatives, so we add two:
#   quiet          genome - (any DHS, H3K27ac, ChromHMM states 1-8), tiled like the dREG negatives;
#   dhs_unmarked   Duke & UW DHSs with no H3K27ac and outside ChromHMM states 1-7 (accessible but
#                  unmarked, mostly CTCF/insulator sites): a hard specificity test, noisier as a label.
# None of these use nascent or TSS data, so they are not circular for any library.
CHROMHMM_ACTIVE = {"1_Active_Promoter", "4_Strong_Enhancer", "5_Strong_Enhancer"}
CHROMHMM_REG = {"1_Active_Promoter", "2_Weak_Promoter", "3_Poised_Promoter", "4_Strong_Enhancer",
                "5_Strong_Enhancer", "6_Weak_Enhancer", "7_Weak_Enhancer", "8_Insulator"}
DISTAL = 1000  # bp from any GENCODE TSS


def _dhs_intersect(files) -> pd.DataFrame:
    a, b = (bf.sort_bedframe(read_bed(f)[BED3]) for f in files)
    ov = bf.overlap(a, b, how="inner", return_overlap=True, suffixes=("", "_b"))
    ov = pd.DataFrame({"chrom": ov.chrom, "start": ov.overlap_start, "end": ov.overlap_end})
    return merge(bf.sort_bedframe(ov[ov.end > ov.start]))


def build_dreg_aux(manifest: dict) -> None:
    g = "hg19"
    aux = REF / "aux"
    tx = _gtf_transcripts()
    tss = _tss_windows(tx[tx.feature == "transcript"], DISTAL)[BED3].drop_duplicates()
    tss19 = merge(bf.sort_bedframe(liftover(tss, "hg38", g)))
    for cell in ("K562", "GM12878"):
        dfiles = [REF / "dhs" / f"{cell}_Duke.hg19.narrowPeak.gz", REF / "dhs" / f"{cell}_UW.hg19.narrowPeak.gz"]
        dhs_any = merge(bf.sort_bedframe(pd.concat([read_bed(f)[BED3] for f in dfiles])))
        dhs_both = _dhs_intersect(dfiles)
        k27 = merge(bf.sort_bedframe(read_bed(aux / f"{cell}_H3K27ac.hg19.broadPeak.gz")[BED3]))
        hmm = read_bed(aux / f"{cell}_ChromHMM.hg19.bed.gz", names=BED3 + ["state"])
        active = merge(bf.sort_bedframe(hmm[hmm.state.isin(CHROMHMM_ACTIVE)]))
        reg = merge(bf.sort_bedframe(hmm[hmm.state.isin(CHROMHMM_REG)]))
        reg17 = merge(bf.sort_bedframe(hmm[hmm.state.isin(CHROMHMM_REG - {"8_Insulator"})]))

        marked = dhs_both[overlaps_any(dhs_both, k27)]
        distal = marked[~overlaps_any(marked, tss19)]
        unmarked = dhs_both[~overlaps_any(dhs_both, k27) & ~overlaps_any(dhs_both, reg17)]
        quiet = bf.subtract(genome_view(g), merge(bf.sort_bedframe(pd.concat([dhs_any, k27, reg]))))[BED3]
        _save(clean(marked, g), "dreg_aux", cell, "pos_dhs_h3k27ac", g, manifest)
        _save(clean(distal, g), "dreg_aux", cell, "pos_dhs_h3k27ac_distal", g, manifest)
        _save(clean(active, g), "dreg_aux", cell, "pos_chromhmm_active", g, manifest)
        _save(clean(quiet, g, "subtract"), "dreg_aux", cell, "neg_quiet", g, manifest)
        _save(clean(unmarked, g), "dreg_aux", cell, "neg_dhs_unmarked", g, manifest)
        manifest[f"check/{cell}_dhs_both_frac_h3k27ac"] = round(len(marked) / max(len(dhs_both), 1), 4)


# ------------------------------------------------------------------ phase 2 (hg38)
# PINTS-on-TSS labels for 3' run-on libraries: positives are PINTS's default bidirectional calls on
# the matched-cell PRO-cap library; negatives follow the dREG recipe with PRO-cap in place of the
# GRO-cap HMM: genome - (DNase U PRO-cap signal +/- 250 bp U positives), tiled at evaluation.
# These labels come from a PINTS run, so they favour PINTS's element definition; they are not used
# on the PRO-cap library they were built from.
TSS_LIB = {"HCT116": "HCT116_PROcap", "A549": "A549_PROcap", "K562": "K562_PROcap", "Jurkat": "Jurkat_CAGE"}
DNASE38 = {c: REF / "dhs" / f"{c}_DNase.hg38.bed.gz" for c in ("HCT116", "A549", "K562")}
H3K27AC38 = {c: REF / "aux" / f"{c}_H3K27ac.hg38.narrowPeak.gz" for c in ("HCT116", "A549")}


def _tss19_38(half=DISTAL):
    tx = _gtf_transcripts()
    return merge(bf.sort_bedframe(_tss_windows(tx[tx.feature == "transcript"], half)[BED3].drop_duplicates()))


def _tss_pints_cell(cell, lib, g, dnase, tss, manifest) -> dict | None:
    from . import calls as C
    try:
        calls = C.load_pints_default(lib)
    except FileNotFoundError:
        print(f"skip tss_pints/{cell}: no PINTS calls for {lib}")
        return None
    pos = merge(bf.sort_bedframe(calls[BED3]))
    sig = merge(bf.sort_bedframe(pd.concat([_nonzero_intervals(DATA / "libs" / f"{lib}.{s}.bw", 250)
                                            for s in ("pl", "mn")])))
    neg = bf.subtract(genome_view(g), merge(bf.sort_bedframe(pd.concat([dnase, sig, pos]))))[BED3]
    variants = {"pints_tss": pos, "pints_tss_dhs": pos[overlaps_any(pos, dnase)],
                "pints_tss_distal": pos[~overlaps_any(pos, tss)]}
    for name, df in variants.items():
        _save(clean(df, g), "tss_pints", cell, f"pos_{name}", g, manifest)
    _save(clean(neg, g, "subtract"), "tss_pints", cell, "neg_tss_silent", g, manifest)
    manifest[f"check/{cell}_pints_tss_frac_dnase"] = round(float(overlaps_any(pos, dnase).mean()), 4)
    return variants


def build_tss_pints(manifest: dict) -> None:
    g = "hg38"
    tss38 = _tss19_38()
    for cell, lib in TSS_LIB.items():
        if cell not in DNASE38:
            continue
        dnase = merge(bf.sort_bedframe(read_bed(DNASE38[cell])[BED3]))
        variants = _tss_pints_cell(cell, lib, g, dnase, tss38, manifest)
        if variants is not None and cell == "K562":
            # K562 run-on libraries are hg19: lift the positives; negatives = dREG hub negatives (no
            # Duke/UW DHS, no GRO-cap HMM) minus the lifted positives +/- 250 bp.
            h = "hg19"
            hub_neg = read_bed(LABELS / "dreg_ftp" / "k562.negative.bed.gz")[BED3]
            lifted = {n: clean(liftover(df, g, h), h) for n, df in variants.items()}
            pad = lifted["pints_tss"].assign(start=lambda d: (d.start - 250).clip(lower=0), end=lambda d: d.end + 250)
            for name, df in lifted.items():
                _save(df, "tss_pints", cell, f"pos_{name}", h, manifest)
            _save(clean(bf.subtract(bf.sort_bedframe(hub_neg), merge(bf.sort_bedframe(pad)))[BED3], h, "subtract"),
                  "tss_pints", cell, "neg_tss_silent", h, manifest)


def build_dreg_aux_hg38(manifest: dict) -> None:
    """dreg_aux for hg38 cells: one ENCODE DNase set (no Duke & UW pair) and H3K27ac; no ChromHMM, so
    'quiet' excludes only DNase and H3K27ac, and 'dhs_unmarked' is DNase without H3K27ac."""
    g = "hg38"
    tss38 = _tss19_38()
    for cell, kf in H3K27AC38.items():
        dnase = merge(bf.sort_bedframe(read_bed(DNASE38[cell])[BED3]))
        k27 = merge(bf.sort_bedframe(read_bed(kf)[BED3]))
        marked = dnase[overlaps_any(dnase, k27)]
        _save(clean(marked, g), "dreg_aux", cell, "pos_dhs_h3k27ac", g, manifest)
        _save(clean(marked[~overlaps_any(marked, tss38)], g), "dreg_aux", cell, "pos_dhs_h3k27ac_distal", g, manifest)
        _save(clean(bf.subtract(genome_view(g), merge(bf.sort_bedframe(pd.concat([dnase, k27]))))[BED3], g, "subtract"),
              "dreg_aux", cell, "neg_quiet", g, manifest)
        _save(clean(dnase[~overlaps_any(dnase, k27)], g), "dreg_aux", cell, "neg_dhs_unmarked", g, manifest)
        manifest[f"check/{cell}_dnase_frac_h3k27ac"] = round(len(marked) / max(len(dnase), 1), 4)


HIDRA = REF / "hidra" / "GSE104001_HiDRA_counts_per_fragmentgroup.txt.gz"
HIDRA_MIN_P, HIDRA_INACTIVE_LFC, HIDRA_ACTIVE_LFC = 30, 0.0, 0.5


def _hidra_groups() -> pd.DataFrame:
    d = pd.read_csv(HIDRA, sep="\t")
    loc = d.iloc[:, 0].str.extract(r"^(chr[^:]+):(\d+)-(\d+)_(\d+)$")
    P = d[[f"P{i}" for i in range(1, 6)]].sum(axis=1).to_numpy()
    R = d[[f"R{i}" for i in range(1, 6)]].sum(axis=1).to_numpy()
    sf = R.sum() / P.sum()
    return pd.DataFrame({"chrom": loc[0], "start": loc[1].astype(np.int64), "end": loc[2].astype(np.int64),
                         "P": P, "R": R, "lfc": np.log2((R + 1) / (P * sf + 1))})


def build_pints_gm12878(manifest: dict) -> None:
    """PINTS-paper-style labels for GM12878.
    enhancers     Supp. Table 2 GM12878 rows (HiDRA STARR-seq, 3,544, hg38).
    nonenhancers  HiDRA-tested, inactive (Supp. Table 3 is K562-only): fragment groups with >= 30
                  plasmid reads and library-normalised log2(RNA/plasmid) <= 0, merged; minus any group
                  with log2FC >= 0.5 (>= 10 plasmid reads), minus positives +/-250 bp; >= 250 bp; minus
                  TSS +/-500 of all transcripts (PINTS's non-enhancer rule); restricted to regions touching a
                  GM12878 Duke/UW DHS, so these are accessible-but-inactive: hard negatives, like
                  dreg_aux dhs_unmarked.
    promoters     same recipe as K562, with GM12878 CAGE and RNA-seq."""
    t2 = _supp_table("SupplementaryTable2")
    enh38 = t2[t2["Cell line"] == "GM12878"][BED3].reset_index(drop=True)
    tx = _gtf_transcripts()
    tx_all = tx[tx.feature == "transcript"]
    prom_any = _tss_windows(tx_all)
    prom = _promoter_pool("GM12878", tx, tx_all)

    h = _hidra_groups()
    inact = h[(h.P >= HIDRA_MIN_P) & (h.lfc <= HIDRA_INACTIVE_LFC)]
    act = h[(h.P >= 10) & (h.lfc >= HIDRA_ACTIVE_LFC)]
    enh19 = clean(liftover(enh38, "hg38", "hg19"), "hg19", "none")
    pad = enh19.assign(start=lambda d: (d.start - 250).clip(lower=0), end=lambda d: d.end + 250)
    ne = merge(bf.sort_bedframe(inact[BED3]))
    ne = bf.subtract(ne, merge(bf.sort_bedframe(pd.concat([act[BED3], pad[BED3]]))))[BED3]
    ne = ne[(ne.end - ne.start) >= 250]
    # HiDRA's library spans much of the genome; keep tested-inactive regions at GM12878 DHSs (Duke or UW),
    # the analogue of PINTS's tested candidate elements.
    dhs = merge(bf.sort_bedframe(pd.concat([read_bed(REF / "dhs" / f"GM12878_{x}.hg19.narrowPeak.gz")[BED3]
                                            for x in ("Duke", "UW")])))
    manifest["check/GM12878_hidra_inactive_regions_before_dhs"] = int(len(ne))
    ne = ne[overlaps_any(ne, dhs)]
    prom_any19 = liftover(prom_any[BED3], "hg38", "hg19")
    ne19 = ne[~overlaps_any(ne, prom_any19)].reset_index(drop=True)
    manifest["check/GM12878_hidra_groups_inactive"] = int(len(inact))
    manifest["check/GM12878_hidra_groups_active_lfc0.5"] = int(len(act))
    manifest["check/GM12878_enh_frac_in_hidra_active"] = round(float(overlaps_any(enh19, act).mean()), 4)

    for genome in ("hg38", "hg19"):
        def L(df, src="hg38"):
            return clean(liftover(df, src, genome), genome)
        _save(L(enh38), "pints", "GM12878", "enhancers_hidra", genome, manifest)
        _save(L(ne19, "hg19"), "pints", "GM12878", "nonenhancers", genome, manifest)
        _save(L(prom), "pints", "GM12878", "promoter_pool", genome, manifest)


def build_jurkat(manifest: dict) -> None:
    """Jurkat (hg19). dreg_aux: UW DNase rep1 & rep2 (in place of Duke & UW), H3K27ac peaks from either
    ChIP-Atlas replicate (rep2 is sparse); no ChromHMM, so 'quiet' excludes any-rep DNase and H3K27ac.
    tss_pints: PINTS on the FANTOM5 Jurkat CAGE library (no Jurkat PRO-cap exists)."""
    g, cell = "hg19", "Jurkat"
    tss = liftover(_tss19_38()[BED3], "hg38", g)
    reps = [REF / "dhs" / f"Jurkat_UW_rep{r}.hg19.narrowPeak.gz" for r in (1, 2)]
    dnase_both = _dhs_intersect(reps)
    dnase_any = merge(bf.sort_bedframe(pd.concat([read_bed(f)[BED3] for f in reps])))
    k1, k2 = (merge(bf.sort_bedframe(read_bed(REF / "aux" / f"Jurkat_H3K27ac_rep{r}.hg19.bed")[BED3])) for r in (1, 2))
    k_rep = merge(bf.sort_bedframe(pd.concat([k1[overlaps_any(k1, k2)], k2[overlaps_any(k2, k1)]])))
    k_any = merge(bf.sort_bedframe(pd.concat([k1, k2])))
    marked = dnase_both[overlaps_any(dnase_both, k_any)]  # either rep: rep2 has half rep1's peaks
    manifest["check/Jurkat_dhs_both_frac_h3k27ac_reproduced"] = round(float(overlaps_any(dnase_both, k_rep).mean()), 4)
    _save(clean(marked, g), "dreg_aux", cell, "pos_dhs_h3k27ac", g, manifest)
    _save(clean(marked[~overlaps_any(marked, tss)], g), "dreg_aux", cell, "pos_dhs_h3k27ac_distal", g, manifest)
    _save(clean(bf.subtract(genome_view(g), merge(bf.sort_bedframe(pd.concat([dnase_any, k_any]))))[BED3], g,
                "subtract"), "dreg_aux", cell, "neg_quiet", g, manifest)
    _save(clean(dnase_both[~overlaps_any(dnase_both, k_any)], g), "dreg_aux", cell, "neg_dhs_unmarked", g, manifest)
    manifest["check/Jurkat_dhs_both_frac_h3k27ac"] = round(len(marked) / max(len(dnase_both), 1), 4)
    _tss_pints_cell(cell, TSS_LIB[cell], g, dnase_any, tss, manifest)
    build_dhs_tss_silent(manifest, cells=["Jurkat"])


# ------------------------------------------------------------------ TSS-silent DHS negatives
# Accessible but untranscribed: DHSs with no read in an independent TSS assay within +/-250 bp. This replaces
# dreg_aux dhs_unmarked (lack of H3K27ac does not imply lack of transcription) as the hard negative. Used
# only on 3' run-on libraries (on a TSS library it would be circular).
SILENT_SRC = {  # cell: (genome, DHS files intersected, TSS-assay bigWig prefix, extra TSS calls or None)
    "K562": ("hg19", [REF / "dhs" / "K562_Duke.hg19.narrowPeak.gz", REF / "dhs" / "K562_UW.hg19.narrowPeak.gz"],
             DATA / "libs" / "K562_GROcap", REF / "K562_grocap_hmm.hg19.bed"),
    "GM12878": ("hg19", [REF / "dhs" / "GM12878_Duke.hg19.narrowPeak.gz", REF / "dhs" / "GM12878_UW.hg19.narrowPeak.gz"],
                REF / "GM12878_GROcap", None),
    "HCT116": ("hg38", [REF / "dhs" / "HCT116_DNase.hg38.bed.gz"], DATA / "libs" / "HCT116_PROcap", None),
    "A549": ("hg38", [REF / "dhs" / "A549_DNase.hg38.bed.gz"], DATA / "libs" / "A549_PROcap", None),
    "Jurkat": ("hg19", [REF / "dhs" / f"Jurkat_UW_rep{r}.hg19.narrowPeak.gz" for r in (1, 2)],
               DATA / "libs" / "Jurkat_CAGE", None),
}
SILENT_PAD = 250


def build_dhs_tss_silent(manifest: dict, cells=None) -> None:
    for cell, (g, dfiles, bwp, extra) in SILENT_SRC.items():
        if cells and cell not in cells:
            continue
        bws = [Path(f"{bwp}.{s}.bw") for s in ("pl", "mn")]
        if not all(b.exists() for b in bws):
            print(f"skip silent/{cell}: missing {bwp}")
            continue
        dhs = _dhs_intersect(dfiles) if len(dfiles) == 2 else merge(bf.sort_bedframe(read_bed(dfiles[0])[BED3]))
        sig = [_nonzero_intervals(b, SILENT_PAD) for b in bws]
        if extra is not None and Path(extra).exists():
            e = read_bed(extra)[BED3]
            sig.append(e.assign(start=(e.start - SILENT_PAD).clip(lower=0), end=e.end + SILENT_PAD))
        sig = merge(bf.sort_bedframe(pd.concat(sig)))
        silent = dhs[~overlaps_any(dhs, sig)]
        _save(clean(silent, g), "silent", cell, "dhs_tss_silent", g, manifest)
        manifest[f"check/{cell}_dhs_frac_tss_silent"] = round(len(silent) / max(len(dhs), 1), 4)
        manifest[f"check/{cell}_dhs_tss_silent_source"] = Path(bwp).name


def build_phase2() -> dict:
    """Adds phase-2 label sets to the existing manifest (needs PINTS calls on the PRO-cap libraries)."""
    mf = BUILT / "manifest.json"
    manifest = json.loads(mf.read_text()) if mf.exists() else {}
    build_dreg_aux_hg38(manifest)
    build_tss_pints(manifest)
    build_pints_gm12878(manifest)
    build_jurkat(manifest)
    mf.write_text(json.dumps(manifest, indent=2))
    return manifest


def build_all() -> dict:
    manifest: dict = {}
    build_pints_k562(manifest)
    build_dreg(manifest)
    build_dreg_aux(manifest)
    BUILT.mkdir(parents=True, exist_ok=True)
    (BUILT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest
