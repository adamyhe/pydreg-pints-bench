# Phase 1 results (2026-11-26)

pydreg 1.0.0 (cupy, P100) vs pyPINTS 1.3.0. Both tools are scored over their full candidate sets with the symmetric
element-level protocol described in `BENCHMARK.md`.

* Raw per-library metrics: `results/eval/<lib>/metrics.tsv`. Seed means: `results/summary/metrics_by_seed_mean.tsv`.
* PINTS-paper labels have 5 promoter-sampling seeds; the means are shown here.
* No confidence intervals yet (bootstrap is TODO).

## Validation
* The rebuilt PINTS candidate pairs cover 100 % of PINTS's own default bidirectional calls in all 5 libraries.
* At PINTS's default cutoff, the rebuilt set reproduces PINTS's output on K562 GRO-cap: 43,433 vs 43,782 elements,
  median width 180 vs 179 bp.
* The dREG negative-rebuild recipe reproduces the original K562 negatives: Jaccard 0.96 with the GRO-cap HMM calls and
  0.88 with the GRO-cap-signal proxy.

## Holdouts for dREG: K562 PRO-seq G7 and GM12878 GRO-seq (cross cell type)

| Labels | Library | pydreg AUROC / AUPRC / TPR@5 %FPR | PINTS AUROC / AUPRC / TPR@5 %FPR |
|---|---|---|---|
| dREG paper (original negatives) | K562 G7 | **0.856 / 0.734 / 0.729** | 0.821 / 0.541 / 0.459 |
| dREG paper (rebuilt negatives) | K562 G7 | **0.858 / 0.750 / 0.732** | 0.841 / 0.622 / 0.609 |
| dREG paper (rebuilt negatives) | GM12878 | **0.977 / 0.956 / 0.958** | 0.971 / 0.892 / 0.870 |
| PINTS paper, assay-filtered (467 neg) | K562 G7 | **0.903 / 0.950 / 0.643** | 0.797 / 0.898 / 0.424 |
| PINTS paper, unfiltered | K562 G7 | 0.885 / **0.809 / 0.787** | **0.897** / 0.798 / 0.733 |

The dREG training libraries (K562 G1 PRO-seq, G2 GRO-seq) show the same direction; see the tables in
`results/summary/`. pydreg's numbers on those two are optimistic.

## TSS assay (K562 GRO-cap; PINTS paper Fig. 4g reported dREG 0.776 vs PINTS 0.858)

| Labels | pydreg AUROC | PINTS AUROC |
|---|---|---|
| PINTS paper, assay-filtered | 0.868 | **0.887** |
| PINTS paper, unfiltered | 0.935 | **0.944** |
| PINTS-paper scoring rule, as run (dREG default calls only) | 0.850 | **0.886** |
| PINTS-paper scoring rule, dREG full candidates | 0.865 | **0.886** |

dREG-paper labels on GRO-cap are **circular**: the positives use GRO-cap HMM calls from this same library, and the
rebuilt negatives use its signal. So they are not interpreted.

## PINTS-paper scoring rule on run-on data (AUC)

| Library | Variant | PINTS | pydreg as run in the paper | pydreg, full candidates |
|---|---|---|---|---|
| K562 G7 | assay-filtered | 0.792 | 0.870 | **0.900** |
| K562 G7 | unfiltered | **0.896** | 0.839 | 0.881 |

* Giving pydreg its full candidate set instead of only its default calls is worth +0.015 to +0.05 AUC on every
  library.
* The rule's "≥ 50 % of the call inside a positive" requirement still penalizes pydreg's ~430 bp calls against PINTS's
  180–340 bp calls.

## Default call sets

| Library | pydreg calls (median width) | PINTS bidirectional calls (median width) |
|---|---|---|
| K562 PRO-seq G1 (deep) | 75,875 (430 bp) | 123,837 (274 bp) |
| K562 PRO-seq G7 | 34,989 (450 bp) | **7,348** (342 bp) |
| K562 GRO-seq G2 | 26,272 (430 bp) | 17,936 (194 bp) |
| GM12878 GRO-seq | 71,448 (430 bp) | 38,393 (289 bp) |
| K562 GRO-cap | 46,665 (360 bp) | 43,782 (179 bp) |

## Caveats
* Each cell type has one library per assay.
* The assay-filtered PINTS-paper negative sets are small (386–1,274) and have no CIs.
* GM12878 lacks PINTS-paper labels until phase 2 (HiDRA negatives).
* The dREG-paper labels were made by dREG's authors. They are tool-independent (DHS ∩ GRO-cap), but they come from the
  same group as the method.

## Phase 1b: bootstrap CIs and dREG auxiliary labels (2026-11-26)
The full tables are in `results/summary/report.md`. CIs come from a 5 Mb genomic block bootstrap, paired between the
tools, with 1,000 replicates.
* **PINTS-paper gaps:** most gaps are significant. On the unfiltered PINTS-paper labels, the AUROC ties on G1 and G2
  hold: the CIs span 0.
* **PINTS-paper rule, unfiltered:**
  * Giving pydreg its full candidate set closes most of the gap: it moves from −0.04…−0.06 to −0.002…−0.015.
  * The rest is significant on G7 and GRO-cap.
* **dREG auxiliary labels:** these are tool- and assay-independent.
  * *Against quiet negatives:* pydreg leads AUPRC and TPR@5%FPR on every run-on library, by +0.05…+0.29 and
    +0.09…+0.33. AUROC is within ±0.025.
  * *Against DHS-unmarked negatives:* these are accessible sites without H3K27ac, mostly CTCF. PINTS leads AUROC
    (−0.013…−0.043) and AUPRC. pydreg calls more of these sites; some may be genuinely transcribed.
* **K562 GRO-cap, non-circular auxiliary labels:** PINTS leads AUROC by 0.004–0.026. TPR@5% is a tie or a small pydreg
  edge.


**Note (2026-11-26):** the `dhs_unmarked` negatives (DHS without H3K27ac) are demoted to diagnostic. Lacking H3K27ac does not mean untranscribed, so results there cannot separate specificity from sensitivity. Conclusions about PINTS being "more specific at accessible-but-unmarked sites" are withdrawn pending a transcription-verified hard negative: DHS with no independent TSS-assay signal.
