# Phase 2 results (2026-11-26): HCT116, A549, K562 PRO-cap, GM12878 HiDRA, PINTS-on-TSS labels

Numbers are the gap pydreg − PINTS: estimate [95 % CI]. CIs are paired 5 Mb block bootstraps with 1,000 replicates.
Full tables: `results/summary/report.md`.

## PINTS-on-TSS labels on 3′ run-on libraries
* Positives are PINTS's own calls on the matched PRO-cap library, so these labels favour PINTS.
* pydreg nevertheless leads on every K562 run-on library. On G7, with PINTS PRO-cap calls ∩ DNase, the gaps are:
  * AUROC +0.023 [+0.021, +0.026];
  * AUPRC +0.159;
  * TPR@5 % FPR +0.229.
  * Distal elements show larger AUPRC gaps, +0.32 to +0.41.
* HCT116 and A549 PRO-seq, cell lines dREG never saw, give the same picture with pints_tss ∩ DNase:
  * AUROC +0.019 (HCT116) and +0.026 (A549);
  * TPR@5 % +0.05 and +0.12.
* Exception: on HCT116 with all PINTS PRO-cap calls (no DNase filter), PINTS is slightly ahead, AUROC −0.011.
* On a TSS assay (K562 GRO-cap), PINTS leads on these labels, AUROC −0.007 to −0.011.

## H3K27ac/DHS labels in new cell types
* *Against quiet negatives:* pydreg leads AUPRC (+0.08 to +0.15) and TPR@5 % (+0.14 to +0.18) on HCT116 and A549
  PRO-seq. AUROC +0.02 to +0.06.
* *Against DHS-unmarked negatives:* PINTS leads AUROC on HCT116 and A549 PRO-seq (−0.017 to −0.035). Every cell type
  shows this pattern: PINTS is more specific against accessible-but-unmarked sites.
* On PRO-cap (TSS assays), PINTS leads in HCT116 (AUPRC −0.05 to −0.10). A549 PRO-cap (9M reads) is near chance for
  both tools on distal elements.

## PINTS-paper labels
* K562 PRO-cap (hg38, native labels): PINTS leads, AUROC −0.018 [−0.024, −0.011] unfiltered.
* GM12878 GRO-seq with HiDRA labels (our negatives: 75k accessible tested-inactive regions):
  * AUROC: pydreg +0.015 [+0.011, +0.019];
  * AUPRC: −0.021, not significant.

## Overall pattern so far
* **3′ run-on data, any cell type and any label family:** pydreg has better early recall and AUPRC; AUROC is similar.
* **Accessible-but-unmarked negatives:** PINTS is more specific.
* **TSS assays:** PINTS is modestly better.


**Note (2026-11-26):** the `dhs_unmarked` negatives (DHS without H3K27ac) are demoted to diagnostic. Lacking H3K27ac does not mean untranscribed, so results there cannot separate specificity from sensitivity. Conclusions about PINTS being "more specific at accessible-but-unmarked sites" are withdrawn pending a transcription-verified hard negative: DHS with no independent TSS-assay signal.
