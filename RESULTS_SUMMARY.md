# pydreg 1.0.0 vs PINTS 1.3.0: consolidated results (2026-11-26)

These supersede the AUROC-based conclusions in RESULTS_PHASE1.md and RESULTS_PHASE2.md. Full per-library tables:
`results/summary/report.md`.

**Data:** 15 libraries in 5 cell types.
* 3′ run-on: K562 PRO-seq ×2, K562 GRO-seq, GM12878 GRO-seq, HCT116 PRO-seq, A549 PRO-seq, Jurkat PRO-seq,
  Jurkat ChRO-seq ×2, Jurkat leChRO-seq.
* TSS assays: K562 GRO-cap, K562 PRO-cap, HCT116 PRO-cap, A549 PRO-cap, Jurkat CAGE.

**Metrics:**
* Headline: TPR @ 5 % FPR and pAUROC (FPR ≤ 0.1). TPR @ 1 % FPR is a point estimate only.
* AUPRC is shown with its baseline. AUROC is supplementary.
* CIs: paired 5 Mb block bootstrap, 1,000 replicates.

**Sensitivity conditions:**
* GC-matched negatives.
* CpG-island and non-CpG-island positives, each scored against GC-matched negatives. The non-CpG-island rows are the
  enhancer-like elements.

## 1. Transcribed regulatory elements vs non-regulatory background, on 3′ run-on data: pydreg is better
* Label sets in this category:
  * dREG-paper negatives;
  * H3K27ac/DHS vs quiet genome;
  * PINTS-on-TSS vs TSS-silent genome.
* pydreg's TPR @ 5 % FPR lead:
  * as evaluated: +0.08 to +0.25 on essentially every run-on library and cell type, with CIs excluding 0;
  * GC-matched: the lead grows;
  * non-CpG-island (enhancer-like) positives: it grows further, typically +0.15 to +0.40.
* GC contrast had been flattering PINTS, not pydreg.
* This holds even with PINTS-on-PRO-cap positives (K562 G7: +0.23 raw, +0.40 non-CpG-island), which favour PINTS by
  construction.
* Exceptions:
  * Jurkat ChRO-seq batch 2 is the weakest case: about 0 after GC matching.
  * On Jurkat, the PINTS-on-CAGE labels favour PINTS on AUPRC.

## 2. Accessible but untranscribed DHSs (TSS-silent): close, with a small PINTS edge on enhancer-like elements
* The question here: among DHSs, can a tool tell transcribed sites from silent ones? This is the specificity test.
* Gaps are small: TPR @ 5 % within ±0.06 and AUPRC within ±0.03 in K562, GM12878 and HCT116, the cells with reliable
  TSS assays.
* On non-CpG-island positives, PINTS is modestly ahead:
  * TPR @ 5 %: K562 G7 −0.056 (dREG positives) and −0.121 (H3K27ac/DHS positives); HCT116 −0.068.
  * Exceptions: the deep K562 G1 library (pydreg ahead) and GM12878 (a tie).
* Jurkat and A549 silent sets come from shallow TSS assays (CAGE 5M, PRO-cap 9M reads) and are unreliable. Jurkat
  shows large PINTS leads there.

## 3. TSS assays: PINTS is better
* PRO-cap (K562, HCT116, A549) and GRO-cap: PINTS leads on PINTS-paper, H3K27ac/DHS and PINTS-on-TSS labels.
* Example: HCT116 PRO-cap, H3K27ac/DHS vs quiet: TPR @ 5 % −0.09 raw, −0.13 to −0.33 under GC matching and the CpG
  split.
* Exception: the shallow Jurkat CAGE library, where pydreg leads.

## 4. PINTS-paper labels (K562 CRISPR; GM12878 HiDRA)
* Without the read filter:
  * pydreg leads TPR @ 5 % on K562 run-on (+0.05 to +0.09; +0.12 to +0.19 for non-CpG-island positives);
  * PINTS leads on GRO-cap.
* GM12878 HiDRA:
  * a tie at 5 % (+0.026 raw, not significant after GC matching);
  * PINTS ahead at 1 % FPR (point estimate);
  * AUPRC tied against a baseline of 0.14.
* The PINTS paper's reported dREG deficit is largely protocol: giving dREG only its default calls, plus the
  asymmetric overlap rule (RESULTS_PHASE1.md).

## Run time
PINTS is faster overall: about 2.9 h against 4.6 h for all 15 libraries. On the deepest run-on libraries the two are
about even. On TSS assays PINTS is 4–9× faster. pydreg uses far less memory: PINTS grows with depth, up to 58 GB.

| Library | pydreg wall | PINTS wall | pydreg RAM | PINTS RAM |
|---|---|---|---|---|
| K562 PRO-seq G1 (deepest) | 44 min | 46 min | 8 GB | 58 GB |
| K562 PRO-seq G7 | 16 | 18 | 5 | 26 |
| GM12878 GRO-seq | 36 | 22 | 7 | 32 |
| Jurkat leChRO-seq | 17 | 19 | 5 | 34 |
| Jurkat PRO-seq / ChRO-seq 1 / 2 | 15 / 15 / 17 | 13 / 7 / 12 | 5 | 17–25 |
| K562 GRO-seq G2 | 15 | 10 | 5 | 16 |
| HCT116 PRO-seq | 18 | 7 | 5 | 15 |
| A549 PRO-seq | 10 | 4 | 4 | 9 |
| TSS assays (PRO-cap ×3, GRO-cap, CAGE) | 6–24 | 0.6–4 | 3–4 | 2–8 |
| **Total, 15 libraries** | **4.6 h** | **2.9 h** | | |

* Measurement: wall time and peak resident memory from `/usr/bin/time -v` in the caller jobs (Slurm 4144326/7,
  4151159/60, 4151191/2), one run each on shared cbsubscb nodes.
* Hardware: pydreg (`--backend cupy`) on one Tesla P100 with 16 CPUs; PINTS on 16 CPU threads.
* Scope: pydreg's time covers scoring every candidate site plus peak calling. CPU-only pydreg was not timed.

## Caveats
* One library per cell and assay, except Jurkat ChRO-seq ×2. Batch effects are comparable to the assay effects
  (Jurkat ChRO-seq 1 vs 2).
* Several label sets derive from dREG's authors (dREG-paper positives) or from PINTS calls (PINTS-on-TSS). The
  H3K27ac/DHS labels are tool-independent but histone-defined.
* TPR @ 1 % FPR has no CIs yet. It was added after these bootstraps started.
* The dhs_unmarked (DHS without H3K27ac) results are diagnostic only.
