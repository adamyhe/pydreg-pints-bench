# pydreg vs PINTS: assay-matched benchmark

See `REVIEW.md` for how the PINTS paper compared the two tools, and the problems with that comparison.
This file covers our own benchmark. Phase 1 re-runs **both papers' benchmarks on run-on data**,
in the cell lines those papers covered, with one symmetric protocol.

## Software
* uv project in this directory, Python 3.13.
* `pydreg[gpu]==1.0.0` (CuPy backend on the cbsubscb P100s).
* `pyPINTS[fast,plots]==1.3.0` (PyPI upstream). `[plots]` is required because `--debug` writes the per-strand
  candidate GTF.
* Interval arithmetic: `bioframe`. bigWig reads: `figwig`. Our code uses no pybedtools or pyBigWig. pyPINTS still pulls
  in pyBigWig internally.
* Site programs: `/programs/kentUtils/bin/{liftOver,bigWigToBedGraph}`.

## Phase 1 libraries (`config/libraries.tsv`, all hg19 bigWigs)

| lib_id | assay | cell | dREG role | source |
|---|---|---|---|---|
| K562_PROseq_G1 | PRO-seq (3′) | K562 | **training** | GSE60456 / GSM1480327 |
| K562_GROseq_G2 | GRO-seq | K562 | **training** | GSM1480325 |
| K562_PROseq_G7 | PRO-seq (3′) | K562 | holdout | Danko hub `G7.{pl,mn}.bw` |
| GM12878_GROseq | GRO-seq | GM12878 | holdout (cross cell type) | GSM1480326 |
| K562_GROcap | GRO-cap (5′ TSS) | K562 | out of domain | GSM1480321 (wTAP) |

* G1 and G2 were in dREG's training data, so dREG's K562 numbers on them are optimistic. The clean K562 test is **G7**.
* K562_GROcap reproduces PINTS's best-case TSS setting (Fig. 4g: dREG 0.776, PINTS 0.858) under a fair protocol.

## Labels (`src/ppbench/labels.py`)

**PINTS-paper labels (`pints_paper`, K562)**
* Positives: Supp. Table 2 CRISPR/CRISPRi enhancers (803), plus an equal number of expressed, CAGE-supported
  GENCODE v24 protein-coding promoters (TSS ± 500 bp). The promoters are sampled per library over 5 seeds.
* Negatives: Supp. Table 3 non-enhancers, filtered to ≥ 250 bp, merged, with every GENCODE TSS ± 500 bp removed.
* Variants:
  * `assay_filtered`: elements must have reads on both strands *in the library under test*. This is PINTS's procedure.
  * `unfiltered`: the same elements without that filter.
* Deviations from the PINTS code:
  * CAGE filter: both strands within TSS ± 1 kb. The code anchors its windows at TSS − 500, which looks like a bug.
  * OverGeneDB is approximated by GENCODE gene-body overlaps.
  * RNA-seq: ENCODE RSEM FPKM, mean of 2 replicates. Total RNA-seq is ENCSR000AEN (PINTS's own libraries); polyA is
    ENCSR000CPH.
  * CAGE: GSE118075 K562 nAnT-iCAGE CTSS (PINTS's SRR7641288/9), lifted to hg38.

**dREG-paper labels (`dreg_paper`)**
* K562:
  * Positives: Duke ∩ UW DHS ∩ GRO-cap HMM.
  * Negatives: outside all three.
  * Both come from the Danko hub copies (`k562.{positive,negative}.bed.gz`).
* GM12878:
  * Positives: from the hub.
  * Negatives: rebuilt with the same recipe. Nonzero GRO-cap signal (± 250 bp) stands in for the unavailable GM12878
    GRO-cap HMM calls.
  * The recipe is validated on K562: `manifest.json` reports the Jaccard index of the rebuilt vs. original negatives, for
    both the HMM and the signal proxy.
* Negatives are scored as 500 bp tiles that hold ≥ 1 dREG informative position in the library. This mirrors the dREG
  paper's informative-position universe. Tiles are capped at 200k, sampled with seed 0.

## Protocol (`src/ppbench/metrics.py`)

**Primary protocol (symmetric)**
* Each labelled element gets the max score of any call overlapping it by ≥ 1 bp, or −∞ if no call touches it.
* ROC and PR are computed over the labelled elements.
* Both tools are swept over their **full** candidate sets:
  * pydreg: `dREG.raw.peak.bed.gz`, ranked by max SVR score. Ranking by −log10 p is also reported.
  * PINTS: per-strand candidates from `--debug` `_peaks.gtf`, re-paired with PINTS's rule (opposite strand within 300 bp,
    lowest-p partner, merged), ranked by −log10 anchor p. This is what the PINTS paper did for PINTS.
* Reported:
  * AUROC, AUPRC, partial AUROC (FPR ≤ 0.1), TPR at FPR 1/5/10 %;
  * each tool's default-call operating point;
  * precision and recall at matched call budgets (top-N calls);
  * call widths.

**PINTS-paper replica (`pints_paper_auc`, PINTS labels only)**
* Their exact rule:
  * TP needs ≥ 50 % of the call inside a positive, from a call that touches no negative;
  * FP is any 1 bp overlap with a negative;
  * ROC is closed with (0,0) and (1,1).
* Two flavours:
  * `paper_default`: pydreg gets only its default calls with a linspace SVR sweep, i.e. the paper as run.
  * `paper_full`: pydreg gets its full candidate set.
* Comparing the replica with the primary protocol shows how much of the published gap the protocol itself creates.

**Check**
* `default_covered_by_full`: the fraction of PINTS's own bidirectional calls covered by the rebuilt candidate pairs.
  This validates the reconstruction.

## Run order (cbsubscb, account `danko_0001`)
```bash
cd pydreg-pints-bench
bash scripts/00_fetch_data.sh 2>&1 | tee logs/00_fetch_data.log   # login node, ~1.5 GB downloads
bash slurm/sub.sh slurm/00_gpu_smoke.sbatch                               # check CuPy on the P100 first
bash slurm/submit_phase1.sh                                       # labels + pydreg(GPU) + PINTS -> evaluate
.venv/bin/python scripts/04_summarize.py                          # tables in results/summary/
```
* Logs go to `$PPB_SLURM_LOGS` (default `logs/slurm/`; see `config/site.env.example`).
* Calls go to `results/calls/<lib>/{pydreg,pints}/`.
* Metrics go to `results/eval/<lib>/`.
* `uv run pytest tests` runs the unit tests on synthetic data and on the released label tables.

## Phase 2 (not yet implemented)
* **More cell types:** HCT116 (ENCODE PRO-seq ENCSR889SZI + PRO-cap ENCSR165QKR, CTCF-AID without auxin) and A549
  (PRO-seq ENCSR244HDV + PRO-cap ENCSR219IHS), hg38.
* **TSS-derived labels:** PINTS calls on matched-cell TSS data used as positives for 3′ calls (K562 GRO-cap/PRO-cap,
  GM12878 GRO-cap, HCT116/A549 PRO-cap).
* **dREG-style labels for HCT116/A549:** ENCODE DNase ∩ TSS calls.
* **GM12878 PINTS-style labels:**
  * enhancers: HiDRA STARR from Supp. Table 2;
  * promoters: GM12878 CAGE/RNA;
  * negatives: inactive HiDRA fragments (GSE104001) under PINTS's filters.
* **Extras:** ChRO-seq (Jurkat, GSE117832; needs a label source), depth downsampling, K562 mNET-seq.

## Phase 1b additions (2026-11-26)

**dREG auxiliary labels (`dreg_aux`, hg19, K562 and GM12878; `labels.build_dreg_aux`)**
* Positives:
  * `dhs_h3k27ac`: Duke ∩ UW DHS overlapping ENCODE Broad H3K27ac peaks (Wang 2019 Fig 1B).
  * `dhs_h3k27ac_distal`: the same, > 1 kb from any GENCODE TSS.
  * `chromhmm_active`: ChromHMM states 1, 4 and 5 (Danko 2015 Fig 1b).
* Negatives:
  * `quiet`: genome − (DHS ∪ H3K27ac ∪ ChromHMM states 1–8), as informative 500 bp tiles.
  * `dhs_unmarked`: Duke ∩ UW DHSs without H3K27ac and outside ChromHMM states 1–7.
* Both dREG papers report these as recall only. The negatives and the ROC framing are ours.
* No nascent or TSS data is used, so these labels are valid on GRO-cap too.
* Source files: UCSC `wgEncodeBroadHistone*H3k27acStdPk` and `wgEncodeBroadHmm*HMM`, in `data/ref/aux/`.

**Bootstrap CIs (`src/ppbench/bootstrap.py`; `slurm/40_bootstrap.sbatch`)**
* 1,000 paired replicates per label set, under two schemes:
  * `element`: stratified by label.
  * `block5Mb`: 5 Mb genomic blocks; this is the headline interval.
* PINTS-paper replicates also redraw the promoter seed.
* The PINTS-paper rule is bootstrapped through per-element indicator matrices. Tests check that unit weights reproduce
  the point metrics exactly.
* Verbose report: `scripts/06_report.py` writes `results/summary/report.md`.
* Run with `bash slurm/submit_ci_aux.sh`.

## Phase 2 labels (2026-11-26)
* **`tss_pints`**: PINTS-on-TSS labels, for the 3′ libraries of the same cell.
  * Positives: PINTS default bidirectional calls on the matched ENCODE PRO-cap library (HCT116, A549, K562). Variants:
    all calls, calls ∩ DNase, and calls > 1 kb from any TSS.
  * Negatives: genome − (DNase ∪ PRO-cap signal ± 250 ∪ positives), as informative 500 bp tiles.
  * K562 positives are lifted to hg19. The K562 hg19 negatives are the dREG hub negatives minus the lifted
    positives ± 250 bp.
  * Never scored on the PRO-cap library the labels came from.
  * The positives favour PINTS's element definition.
* **`dreg_aux` on hg38 (HCT116, A549)**: ENCODE DNase (one set, not Duke ∩ UW) with ENCODE H3K27ac. There is no
  ChromHMM, so the `quiet` negatives exclude only DNase and H3K27ac.
* **`pints_paper` for GM12878**:
  * Enhancers: Supp. Table 2 HiDRA STARR-seq rows.
  * Promoters: the K562 recipe, with GM12878 CAGE (GSM3318207/8) and ENCODE RNA-seq (ENCSR000AEE total,
    ENCSR000COQ polyA).
  * Non-enhancers: ours, because Supp. Table 3 is K562 only. These are HiDRA-tested inactive regions (≥ 30 plasmid
    reads, log2 RNA/plasmid ≤ 0) with these filters:
    * no overlap with any group at log2FC ≥ 0.5;
    * ≥ 250 bp;
    * touching a Duke/UW DHS;
    * clear of TSS ± 500.
  * 75,311 negatives in hg19, 15× the K562 set and all accessible, so this is a harder specificity test than K562's.

## Phase 3: Jurkat, hg19 (2026-11-26)
* **Libraries** (GSE117832, Chu 2018; Danko-lab pipeline bigWigs, cleaned of explicit zeros):
  * nuclei PRO-seq (~25M reads);
  * ChRO-seq batch 1 (~25M) and batch 2 (~31M);
  * leChRO-seq (~31M);
  * FANTOM5 CAGE CNhs11253 (~5M; CTSS converted to bigWigs).
  * PRO-seq and ChRO-seq from the same cells separate the assay effect from the cell effect.
* **`dreg_aux`:**
  * DNase: UW rep 1 ∩ rep 2.
  * H3K27ac: ChIP-Atlas MACS2 peaks for GSM4551183/4, from either replicate (rep 2 has 16k peaks vs. rep 1's 29k).
  * Counts: 15,613 marked DHSs, of which 2,600 are distal.
  * The H3K27ac set is promoter-heavy, so distal results rest on a small set.
* **`tss_pints`:** PINTS calls on the Jurkat CAGE library. There is no Jurkat PRO-cap or GRO-cap. CAGE also captures
  stable-RNA TSSs, and this library is shallow.


**Note (2026-11-26):** the `dhs_unmarked` negatives (DHS without H3K27ac) are demoted to diagnostic. Lacking H3K27ac does not mean untranscribed, so results there cannot separate specificity from sensitivity. Conclusions about PINTS being "more specific at accessible-but-unmarked sites" are withdrawn pending a transcription-verified hard negative: DHS with no independent TSS-assay signal.

## TSS-silent DHS negatives (`dhs_tss_silent`, replaces `dhs_unmarked` as the hard negative)
* Definition: DHSs (Duke ∩ UW, ENCODE DNase, or UW rep1 ∩ rep2) with no read in an independent TSS assay within
  ± 250 bp.
* TSS assay per cell:
  * K562: GRO-cap reads plus the GRO-cap HMM calls.
  * GM12878: GRO-cap.
  * HCT116 and A549: PRO-cap.
  * Jurkat: CAGE.
* Silent fraction of DHSs: K562 23 %, GM12878 18 %, HCT116 17 %, A549 44 %, Jurkat 63 %. A549 and Jurkat have shallow
  TSS assays, so their silent sets are contaminated by under-sampling.
* Paired with three kinds of positives, each minus any overlap with the negatives:
  * dREG-paper positives;
  * PINTS-on-TSS positives;
  * H3K27ac/DHS positives.
* Scored on 3′ run-on libraries only.

## GC matching and CpG-island split (bootstrap conditions)
* **Diagnostic** (`scripts/07_gc_diagnostic.py` writes `results/summary/gc_diagnostic.tsv`): GC alone separates
  positives from negatives with AUROC 0.73–0.91, because positives are promoter- and CpG-island-rich.
* **`|gcw`:** negatives are reweighted, in every replicate, so their GC histogram matches the positives' (20 pooled
  quantile bins).
* **`|cgi_gcw`, `|noncgi_gcw`:** positives restricted to CpG-island-like elements (GC > 0.5 and CpG o/e > 0.6) or the
  rest, with negatives GC-matched to that subset. The non-CpG-island rows isolate enhancer-like elements.
* Conditions are applied to pydreg (SVR) and PINTS. The paired difference gets CIs.

## Reporting metrics (from 2026-11-26)
* **Headline:** TPR @ 1 % FPR, TPR @ 5 % FPR, and pAUROC (FPR ≤ 0.1). They focus on the low-FPR range that matters
  genome-wide and do not depend on the class ratio.
* **AUPRC:** always shown with its baseline (the positive fraction). Negatives per positive range from 0.3 to 595 by
  construction (200k cap, read filter), so AUPRC compares tools within a label set only.
* **AUROC:** supplementary. Conclusions are not drawn from AUROC ties.
* **Matched call budgets:** recall and labelled precision for the top 10k / 20k / 40k merged calls, shown per label set.
