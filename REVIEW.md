# How PINTS (Yao et al. 2022) benchmarked dREG: a review

Paper: Yao L. et al., *Nat Biotechnol* 2022, "A comparison of experimental assays and analytical methods for genome-wide identification of active enhancers" — https://doi.org/10.1038/s41587-022-01211-7 (PMC9288987; full text from PubMed Central).
Code: https://github.com/hyulab/PINTS_analysis (cloned to `external/PINTS_analysis`). The key files are
`Analysis/tools_sensitivity_specificity.py` (label construction and driver) and `Analysis/utils_roc.py` (per-tool scoring and ROC).

## 1. What data dREG was tested on

* The tool ROC (Fig. 4g, "PINTS mean AUC 0.780 vs. dREG 0.652, second best") was run **only on TSS assays**:
  GRO-cap, CoPRO, csRNA-seq, NET-CAGE, RAMPAGE, CAGE and STRIPE-seq. It used **K562 only, hg38 only**
  (`plot_assays=(...)`, and `main()` filters to `Cell line == "K562"`). The paper says this directly: "evaluated each
  tool's performance for all TSS-assay datasets".
* The signal tracks were made with `pints_visualizer -e <EXP>` (`Protocols/Unified alignment pipeline.json`). For TSS
  assays these are single-base **5′/TSS-end counts**. dREG's SVR was trained on 3′-end Pol II run-on density
  (PRO-seq/GRO-seq). So dREG was scored on inputs it was never trained on, using the stock
  `asvm.gdm.6.6M.20170828.rdata` model.
* **No PRO-seq, GRO-seq or ChRO-seq library is in the tool ROC.** Those data appear only in the assay-level comparisons
  (Figs. 1–3), not in the head-to-head between callers. The paper's summary that "PINTS outperforms dREG" therefore has no
  direct evidence behind it for run-on data.

## 2. What the positive and negative labels were

From Methods and `_enhancer_set` / `_promoter_set`:

| Set | Source | Filters |
|---|---|---|
| Positive: enhancers | 938 CRISPR/CRISPRi-validated K562 enhancers (Fulco 2016 & 2019, Gasperini 2019, Klann, Thakore and others) | merged; **must have reads on both strands in the assay being evaluated** |
| Positive: promoters | GENCODE v24 protein-coding TSS ± 500 bp, randomly sampled to equal the enhancer count | reads on both strands in the assay; **divergent K562 CAGE signal within ±1 kb**; FPKM ≥ 1 in both total and polyA RNA-seq; overlapping genes (OverGeneDB) removed |
| Negative: non-enhancers | STARR-seq (20,941) and MPRA (17,462) segments that were not active | ≥ 250 bp; **no overlap with any ENCODE cCRE PLS/ELS**; not within ±500 bp of a GENCODE TSS; **must have reads on both strands in the assay** |

Concerns with these labels:

1. **The labels depend on the test data.** Each element is kept only if the library being evaluated has reads on both
   strands over it (`direction=2`). As a result, the labels differ per assay and are tied to the same data the callers see.
   Any negative that keeps divergent signal is still labeled negative. These are likely transcribed elements that
   happened not to score in STARR-seq or MPRA. A caller that is sensitive to bidirectional transcription gets "FPs"
   precisely where it is arguably right.
2. **The negatives are easy.** Removing every cCRE PLS/ELS overlap leaves negatives that are mostly closed chromatin.
   Telling these apart is easy, so the comparison mostly measures how cleanly each caller excludes low-level noise
   transcription.
3. **The promoter positives are selected by a TSS assay.** They must show divergent CAGE signal, and CAGE is one of the
   benchmarked assays. This favors callers built around sharp, divergent TSS peaks.
4. The positive set is small (938 enhancers before filtering) and comes from one cell line. The non-enhancer set's
   validity also rests on STARR-seq/MPRA negatives, and the paper's Discussion admits these contain false negatives.
5. Promoter sampling is **not reproducible**: `seed = int(datetime.now().strftime("%m%d%I%M"))`.

## 3. How each tool was scored. This is the biggest issue.

The ROC sweep is not run the same way for the two tools:

* **PINTS** (`ProfilerPINTS`): starts from the per-strand peak files that contain **all** candidate peaks, including
  non-significant ones. It sweeps the p-value threshold from the smallest observed value up to **1.0**
  (`_unbalanced_range_for_pqs`) and **re-pairs** divergent peaks at each threshold with `merge_opposite_peaks`. Its curve
  therefore reaches the high-sensitivity end.
* **dREG** (`ProfilerdREG`): takes only the **final, already FDR-thresholded** `dREG.peak.full.bed.gz` and sweeps 101
  values of its SVR score over that set's [min, max]. dREG never gets to show calls below its default cutoff. Its ROC
  stops early, and `_roc_formatter` / `interpolate_roc` then fill the gap to (1, 1) with a straight line, which is the
  random-classifier segment. dREG writes `dREG.raw.peak.bed.gz` with all candidate peaks, but the benchmark did not use
  it. On its own, this would lower dREG's AUC even if the two tools ranked elements equally well.
* **The overlap rules are asymmetric** (`evaluation_core_element`):
  * TP: a positive counts only if a call overlaps it with **≥ 50 % of the call's length** (`-F 0.5`), and that call
    touches no negative.
  * FP: a negative counts if any call overlaps it by **≥ 1 bp**.

  Wider calls are penalized twice: they are less likely to fall 50 % inside a positive, and more likely to graze a
  negative. Per the paper, dREG peaks average 381–460 bp and dREG.HD peaks 548–751 bp, versus 185–300 bp for PINTS.
* dREG.HD, which is what dREG users typically report, got no ROC at all, only one operating point.

## 4. Smaller code issues

* The 15 M-read downsampling builds its labels from the **19 M-read tracks**: `ds_cov_tracks_15m` is filled from
  `dspl_covs` / `dsmn_covs` instead of the `_15m` dictionaries.
* The code can't be re-run as is. Every input path points to `/fs/cbsuhy0{1,2}/storage/ly349/...` and a BioQueue
  workspace, and `config.txt` is not in the repo. pints.yulab.org returned HTTP 502 on 2026-11-26.
* **Correction:** the label sets *are* released, in the Nature Biotech supplement (`41587_2022_1211_MOESM3_ESM.xlsx`):
  * Supp. Table 2 (known enhancers):

    | Cell line | Support | Count |
    |---|---|---|
    | K562 | CRISPR/CRISPRi only | 786 |
    | K562 | CRISPR/CRISPRi + MPRA/STARR | 17 |
    | K562 | MPRA/STARR only | 533 |
    | GM12878 | STARR-seq (HiDRA) | 3,544 |

  * Supp. Table 3: 7,097 non-enhancers. There is no cell-line column; the non-enhancer methods describe K562 sources.
  * Promoter positives are not tabulated. They must be rebuilt from GENCODE v24, K562 CAGE and RNA-seq, following the
    code.
* Supp. Table 1 lists the run-on libraries in their compendium: K562 PRO-seq (SRR1554311/2), K562 GRO-seq
  (SRR1552484) and GM12878 GRO-seq (SRR1552485), all from GSE60456, plus K562 mNET-seq. None of them entered the tool ROC.
* The per-assay AUCs behind Fig. 4g are in `MOESM6` (sheet Fig_4_g). dREG and PINTS scored:

  | Assay | dREG | PINTS |
  |---|---|---|
  | GRO-cap | 0.776 | 0.858 |
  | CoPRO | 0.782 | 0.879 |
  | csRNA-seq | 0.761 | 0.866 |
  | NET-CAGE | 0.619 | 0.781 |
  | RAMPAGE | 0.557 | 0.702 |
  | CAGE | 0.553 | 0.695 |
  | STRIPE-seq | 0.515 | 0.677 |

  dREG does worst on the assays furthest from run-on signal.
* No plotting or AUC code is in the repo, so how those AUCs were integrated can't be checked. Because of the
  straight-line fill at the end of dREG's ROC, they likely understate dREG.

## Bottom line

The PINTS paper's dREG comparison:

* covers TSS assays only, with dREG fed out-of-distribution 5′-end signal;
* uses labels that are filtered by the data under test and lean toward TSS-assay behavior;
* scores PINTS across its full p-value range but dREG only above its default cutoff;
* uses overlap rules that penalize wider calls.

It says nothing about PRO-seq, GRO-seq or ChRO-seq. That case needs its own benchmark.

---

# Proposed benchmark (prepared; nothing has been run)

**Environment:** a uv project in this directory, Python 3.13, with `pydreg[gpu]==1.0.0`, `pyPINTS[fast]==1.3.0`
(upstream PyPI, 2026-11-22) and `bioframe`. Our code uses no pybedtools or pyBigWig. pyBigWig is still installed because
pyPINTS requires it. bigWig I/O uses `figwig`, which ships with pydreg. pydreg must run on a GPU node; this login node has
none.

**Data:** public bigWigs, with no realignment needed:

| Data | Accession | Genome |
|---|---|---|
| K562 PRO-seq, GRO-seq, and GRO-cap from the same cells (Core et al. 2014) | GEO GSE60456 | hg19 |
| Jurkat ChRO-seq (Chu et al. 2018) | GEO GSE117832 | — |
| K562 PRO-cap and its PINTS pairs | ENCODE ENCSR220XSM | hg38 |
| K562 ATAC peaks and nonpeaks | ENCODE ENCSR868FGK | hg38 |
| ENCODE cCREs | ENCODE SCREEN registry | hg38 |

**Truth sets, all independent of the run-on library under test:**

* **T1 initiation:** divergent TSS pairs from GRO-cap in the same study (Core 2014), and separately from ENCODE PRO-cap.
  Caveat: if these pairs are themselves PINTS calls, they lean toward PINTS's element definition. A tool-agnostic version
  can be built directly from GRO-cap read counts.
* **T2 chromatin:** K562 ATAC peaks (positives) vs. ATAC nonpeaks (negatives). These labels are cell-matched and fully
  orthogonal. Not every open region is transcribed, so recall is capped, but the cap is the same for both tools.
* **T3 functional:** CRISPRi-tested elements (e.g. the ENCODE-rE2G K562 CRISPR benchmark) plus PINTS's STARR/MPRA-style
  negatives. Each will be scored **with and without** PINTS's in-assay signal filter, to measure how much that filter
  matters.

**Metrics** (the same rules for both tools):

* Element-level ROC and PR. Each labeled element is scored by the max score of any overlapping call, and the same
  overlap rule is used for positives and negatives.
* Score sweep over the full candidate sets: pydreg `raw.peak.bed.gz`, and the PINTS per-strand peaks re-paired across
  p-values, exactly as PINTS's own ROC did.
* Precision and recall at **matched call budgets** (top-N by score), which removes the effect of each tool's default
  threshold.
* Resolution: distance from each call's center to the nearest TSS-pair center.
* Depth robustness: downsampling the K562 PRO-seq library.
* An optional "PINTS-paper metric" mode that reproduces their asymmetric rule, so we can show how much the metric alone
  changes the ranking.
