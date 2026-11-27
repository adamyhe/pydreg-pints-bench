# Phase 1 detailed report

pydreg 1.0.0 vs pyPINTS 1.3.0. Each cell gives the estimate [95 % CI], and Δ rows add a two-sided paired-bootstrap p.

* CIs are paired bootstrap intervals (1,000 replicates; `src/ppbench/bootstrap.py`). **block5Mb** resamples 5 Mb genomic blocks and is the headline interval. **element** is stratified by label and is in `results/bootstrap/*/ci.tsv`.
* For PINTS-paper labels, each replicate also redraws the promoter sample (5 seeds), and the estimate is the mean over seeds.
* **Headline metrics:** TPR at 1 % and 5 % FPR, and pAUROC (FPR ≤ 0.1, normalised). These sit in the low-FPR range where genome-wide calling operates, and do not depend on the class ratio.
* **AUPRC** is shown with its baseline, the positive fraction, which varies widely between label sets because negative counts are capped or filtered by construction. Compare the tools within a set; do not compare AUPRC across sets.
* **AUROC** is supplementary only. It weights every FPR equally, and genome-wide only low FPRs matter.
* TPR @ 1 % FPR cells without an interval are as-evaluated point estimates: bootstraps run before that metric was added.
* All metrics are element-level and symmetric: a labelled element takes the max score of any call touching it.
* GC-matched rows reweight negatives, in every replicate, so their GC histogram (20 bins) matches the positives'. CpG-island-like = GC > 0.5 and CpG o/e > 0.6; the non-CpG-island rows isolate enhancer-like positives.
* `pints_paper_auc` uses the PINTS_analysis rule: TP needs ≥ 50 % of a call inside a positive, and FP is any touch.

## Libraries

| lib | cell | assay | signal | role for dREG | notes |
|---|---|---|---|---|---|
| K562_PROseq_G1 | K562 | PRO-seq | 3p | train | Core 2014 (GSE60456); also in PINTS compendium (SRR1554311/2); dREG training set G1 |
| K562_GROseq_G2 | K562 | GRO-seq | run-on | train | Core 2014; PINTS compendium (SRR1552484); dREG training set G2 |
| K562_PROseq_G7 | K562 | PRO-seq | 3p | holdout | dREG paper K562 holdout G7 (Danko lab hub copy; TROGDOR maps it to GSM2545325) |
| GM12878_GROseq | GM12878 | GRO-seq | run-on | holdout | dREG paper GM12878 holdout; PINTS compendium (SRR1552485) |
| K562_GROcap | K562 | GRO-cap | 5p | ood | TSS assay used in PINTS Fig 4g (their own alignment of SRR1552480); dREG out of domain |
| HCT116_PROseq | HCT116 | PRO-seq | 3p | new_cell | ENCODE ENCSR889SZI (HCT116 CTCF-AID, no auxin); 2 bio reps summed |
| HCT116_PROcap | HCT116 | PRO-cap | 5p | ood | ENCODE ENCSR165QKR; 2 bio reps summed; TSS assay + source of PINTS-on-TSS labels |
| A549_PROseq | A549 | PRO-seq | 3p | new_cell | ENCODE ENCSR244HDV (pooled reps 1+2) |
| A549_PROcap | A549 | PRO-cap | 5p | ood | ENCODE ENCSR219IHS; 2 bio reps summed |
| K562_PROcap | K562 | PRO-cap | 5p | ood | ENCODE ENCSR220XSM; TSS assay + PINTS-on-TSS labels for K562 run-on libraries |
| Jurkat_PROseq | Jurkat | PRO-seq | 3p | new_cell | GSE117832 Jurkat nuclei PRO-seq (Chu 2018); dREG 2015 applied dREG to Jurkat GRO-seq but did not train on Jurkat |
| Jurkat_ChROseq1 | Jurkat | ChRO-seq | 3p | new_cell | GSE117832 ChRO-seq batch 1 |
| Jurkat_ChROseq2 | Jurkat | ChRO-seq | 3p | new_cell | GSE117832 ChRO-seq batch 2 (no RNase) |
| Jurkat_leChROseq | Jurkat | leChRO-seq | 3p | new_cell | GSE117832 length-extension ChRO-seq (RNase) |
| Jurkat_CAGE | Jurkat | CAGE | 5p | ood | FANTOM5 hCAGE, 1 library; TSS assay + source of Jurkat PINTS-on-TSS labels |

dREG roles: train = in dREG's training data (optimistic for pydreg), holdout = held out by the dREG paper, new_cell = cell line dREG never saw, ood = TSS assay outside dREG's design.

## PINTS paper (Yao 2022) labels, with the paper's read filter

* **Positives:** **K562:** CRISPR-validated enhancers (Supp. Table 2, CRISPRi/CRISPRa rows; 801 after hg38→hg19 liftOver and blacklist) **plus** an equal number of promoters drawn from a pool (5 draws = 5 seeds). Pool: GENCODE v24 protein-coding TSS ± 500 bp with K562 CAGE reads on both strands within TSS ± 1 kb, FPKM ≥ 1 in both ENCODE total and polyA RNA-seq, excluding genes overlapping another gene (7,504 in hg19). **GM12878:** the 3,544 HiDRA STARR-seq enhancers of Supp. Table 2 (3,497 in hg19), plus promoters from the same recipe with GM12878 CAGE and ENCODE RNA-seq (pool 7,149).
* **Negatives:** **K562:** Supp. Table 3 non-enhancers (CRISPR/MPRA/STARR-tested, inactive): ≥ 250 bp, merged, minus TSS ± 500 of all transcripts (5,952 in hg19). **GM12878** (ours; Supp. Table 3 is K562 only): HiDRA (GSE104001) fragment groups with ≥ 30 plasmid reads and library-normalised log2(RNA/plasmid) ≤ 0, merged, minus any group with log2FC ≥ 0.5 and the positives ± 250 bp. Kept if ≥ 250 bp, touching a Duke/UW DHS and clear of TSS ± 500 (75,311 in hg19). These are accessible-but-inactive hard negatives.
* **Notes:** Positives and negatives are both restricted to elements with ≥ 1 read on **each** strand in the library under test (as in PINTS_analysis). This leaves few negatives on run-on data, and it favours whichever tool calls what the assay sees.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 1,444, n_neg = 1,274, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.531**; CpG-island rows 0.339 (n_pos 652), non-CpG-island rows 0.383 (n_pos 792).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.232 | 0.659 [0.586, 0.706] | 0.572 [0.517, 0.629] | 0.908 [0.891, 0.923] | 0.904 [0.889, 0.918] |  |
| pydreg, −log10 p over all candidates | 0.095 | 0.474 [0.385, 0.610] | 0.441 [0.379, 0.510] | 0.875 [0.856, 0.894] | 0.889 [0.874, 0.905] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.163 | 0.507 [0.428, 0.558] | 0.447 [0.401, 0.493] | 0.864 [0.848, 0.881] | 0.848 [0.833, 0.863] |  |
| Δ pydreg(SVR) − PINTS | 0.069 | 0.152 [0.072, 0.232] p<1e-03 | 0.126 [0.063, 0.188] p<1e-03 | 0.044 [0.023, 0.063] p<1e-03 | 0.056 [0.039, 0.073] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | -0.068 | -0.033 [-0.130, 0.118] p=0.674 | -0.005 [-0.074, 0.068] p=0.914 | 0.011 [-0.011, 0.033] p=0.294 | 0.041 [0.023, 0.060] p<1e-03 |  |
| GC-matched: pydreg (SVR) |  | 0.611 [0.439, 0.693] | 0.523 [0.434, 0.604] | 0.896 [0.871, 0.916] | 0.892 [0.870, 0.909] |  |
| GC-matched: PINTS |  | 0.509 [0.441, 0.575] | 0.466 [0.417, 0.522] | 0.873 [0.854, 0.890] | 0.851 [0.832, 0.869] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.103 [-0.088, 0.199] p=0.35 | 0.056 [-0.038, 0.140] p=0.222 | 0.023 [-0.003, 0.045] p=0.094 | 0.041 [0.014, 0.063] p=0.004 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.060 [-0.278, 0.159] p=0.858 | -0.081 [-0.273, 0.067] p=0.344 | -0.034 [-0.096, 0.016] p=0.23 | -0.005 [-0.029, 0.016] p=0.71 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.188 [0.076, 0.281] p=0.002 | 0.158 [0.086, 0.229] p<1e-03 | 0.113 [0.064, 0.160] p<1e-03 | 0.083 [0.053, 0.113] p<1e-03 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.881 [0.864, 0.896] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.899 [0.883, 0.913] |
| PINTS rule: PINTS |  |  |  |  |  | 0.841 [0.826, 0.857] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | 0.039 [0.022, 0.057] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | 0.057 [0.039, 0.075] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.399 | 0.960 | 0.324 | 0.947 |
| 20,000 | 0.623 | 0.942 | 0.456 | 0.924 |
| 40,000 | 0.787 | 0.890 | 0.578 | 0.887 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 1,088, n_neg = 386, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.738**; CpG-island rows 0.565 (n_pos 501), non-CpG-island rows 0.603 (n_pos 587).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.222 | 0.513 [0.382, 0.625] | 0.468 [0.384, 0.564] | 0.945 [0.932, 0.957] | 0.880 [0.857, 0.902] |  |
| pydreg, −log10 p over all candidates | 0.083 | 0.433 [0.263, 0.584] | 0.384 [0.277, 0.502] | 0.932 [0.914, 0.949] | 0.872 [0.848, 0.895] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.092 | 0.242 [0.185, 0.307] | 0.231 [0.175, 0.299] | 0.874 [0.850, 0.895] | 0.733 [0.704, 0.762] |  |
| Δ pydreg(SVR) − PINTS | 0.129 | 0.270 [0.130, 0.374] p<1e-03 | 0.237 [0.143, 0.325] p<1e-03 | 0.072 [0.053, 0.093] p<1e-03 | 0.147 [0.121, 0.175] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | -0.009 | 0.191 [0.036, 0.336] p=0.028 | 0.153 [0.037, 0.271] p=0.008 | 0.059 [0.036, 0.082] p<1e-03 | 0.139 [0.110, 0.169] p<1e-03 |  |
| GC-matched: pydreg (SVR) |  | 0.300 [0.247, 0.677] | 0.375 [0.255, 0.613] | 0.940 [0.917, 0.962] | 0.878 [0.838, 0.917] |  |
| GC-matched: PINTS |  | 0.286 [0.220, 0.447] | 0.303 [0.228, 0.384] | 0.894 [0.869, 0.914] | 0.769 [0.728, 0.804] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.014 [-0.146, 0.379] p=0.598 | 0.072 [-0.066, 0.301] p=0.406 | 0.047 [0.025, 0.073] p<1e-03 | 0.109 [0.076, 0.150] p<1e-03 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.187 [-0.322, 0.406] p=0.81 | -0.087 [-0.228, 0.343] p=0.856 | 0.028 [-0.018, 0.082] p=0.256 | 0.049 [0.000, 0.104] p=0.048 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.121 [0.016, 0.271] p=0.02 | 0.126 [0.046, 0.229] p=0.006 | 0.119 [0.079, 0.163] p<1e-03 | 0.160 [0.121, 0.203] p<1e-03 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.840 [0.814, 0.864] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.877 [0.854, 0.899] |
| PINTS rule: PINTS |  |  |  |  |  | 0.729 [0.699, 0.758] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | 0.111 [0.083, 0.140] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | 0.148 [0.121, 0.177] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.523 | 0.968 | 0.298 | 0.918 |
| 20,000 | 0.726 | 0.944 | 0.474 | 0.890 |
| 40,000 | 0.876 | 0.907 | 0.650 | 0.873 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 1,100, n_neg = 467, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.702**; CpG-island rows 0.512 (n_pos 489), non-CpG-island rows 0.567 (n_pos 611).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.452 | 0.643 [0.589, 0.704] | 0.603 [0.544, 0.666] | 0.950 [0.939, 0.960] | 0.903 [0.888, 0.919] |  |
| pydreg, −log10 p over all candidates | 0.223 | 0.649 [0.556, 0.698] | 0.566 [0.494, 0.642] | 0.944 [0.931, 0.955] | 0.900 [0.883, 0.915] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.198 | 0.456 [0.381, 0.512] | 0.409 [0.345, 0.482] | 0.899 [0.883, 0.914] | 0.798 [0.774, 0.822] |  |
| Δ pydreg(SVR) − PINTS | 0.254 | 0.187 [0.122, 0.273] p<1e-03 | 0.194 [0.112, 0.270] p<1e-03 | 0.052 [0.038, 0.065] p<1e-03 | 0.106 [0.082, 0.128] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | 0.024 | 0.193 [0.093, 0.274] p=0.002 | 0.158 [0.059, 0.243] p=0.004 | 0.045 [0.031, 0.060] p<1e-03 | 0.102 [0.078, 0.126] p<1e-03 |  |
| GC-matched: pydreg (SVR) |  | 0.577 [0.448, 0.677] | 0.532 [0.426, 0.648] | 0.944 [0.927, 0.958] | 0.890 [0.864, 0.914] |  |
| GC-matched: PINTS |  | 0.475 [0.386, 0.582] | 0.447 [0.346, 0.532] | 0.905 [0.885, 0.921] | 0.811 [0.785, 0.837] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.102 [-0.067, 0.224] p=0.234 | 0.085 [-0.046, 0.238] p=0.206 | 0.039 [0.019, 0.061] p=0.002 | 0.079 [0.044, 0.114] p<1e-03 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.014 [-0.413, 0.282] p=0.97 | 0.033 [-0.234, 0.351] p=0.728 | 0.019 [-0.028, 0.071] p=0.386 | 0.011 [-0.026, 0.048] p=0.496 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.149 [0.041, 0.255] p<1e-03 | 0.149 [0.062, 0.242] p<1e-03 | 0.127 [0.087, 0.169] p<1e-03 | 0.153 [0.111, 0.194] p<1e-03 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.870 [0.852, 0.887] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.900 [0.884, 0.916] |
| PINTS rule: PINTS |  |  |  |  |  | 0.793 [0.769, 0.817] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | 0.077 [0.054, 0.100] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | 0.107 [0.083, 0.131] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.453 | 0.988 | 0.369 | 0.971 |
| 20,000 | 0.654 | 0.960 | 0.475 | 0.939 |
| 40,000 | 0.835 | 0.933 | 0.658 | 0.895 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 5,268, n_neg = 32,784, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.138**; CpG-island rows 0.078 (n_pos 2,778), non-CpG-island rows 0.071 (n_pos 2,490).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.098 | 0.468 [0.444, 0.492] | 0.410 [0.391, 0.429] | 0.504 [0.483, 0.526] | 0.858 [0.852, 0.865] |  |
| pydreg, −log10 p over all candidates | 0.050 | 0.252 [0.237, 0.264] | 0.252 [0.238, 0.265] | 0.398 [0.381, 0.413] | 0.835 [0.828, 0.842] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.169 | 0.472 [0.443, 0.497] | 0.420 [0.395, 0.443] | 0.507 [0.473, 0.537] | 0.826 [0.817, 0.835] |  |
| Δ pydreg(SVR) − PINTS | -0.072 | -0.004 [-0.028, 0.026] p=0.934 | -0.010 [-0.033, 0.014] p=0.38 | -0.003 [-0.032, 0.027] p=0.79 | 0.032 [0.025, 0.040] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | -0.119 | -0.220 [-0.243, -0.196] p<1e-03 | -0.168 [-0.190, -0.146] p<1e-03 | -0.109 [-0.136, -0.082] p<1e-03 | 0.009 [0.002, 0.017] p=0.006 |  |
| GC-matched: pydreg (SVR) |  | 0.457 [0.417, 0.498] | 0.404 [0.370, 0.434] | 0.508 [0.476, 0.541] | 0.858 [0.850, 0.867] |  |
| GC-matched: PINTS |  | 0.443 [0.411, 0.480] | 0.401 [0.371, 0.432] | 0.491 [0.456, 0.525] | 0.819 [0.809, 0.829] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.014 [-0.037, 0.056] p=0.57 | 0.002 [-0.037, 0.036] p=0.9 | 0.018 [-0.021, 0.051] p=0.342 | 0.039 [0.030, 0.049] p<1e-03 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.004 [-0.113, 0.077] p=0.934 | -0.012 [-0.097, 0.055] p=0.768 | -0.002 [-0.084, 0.070] p=0.906 | 0.011 [0.000, 0.021] p=0.05 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.028 [-0.002, 0.057] p=0.08 | 0.027 [0.003, 0.049] p=0.026 | 0.030 [0.013, 0.046] p<1e-03 | 0.064 [0.052, 0.077] p<1e-03 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.826 [0.818, 0.835] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.833 [0.825, 0.841] |
| PINTS rule: PINTS |  |  |  |  |  | 0.810 [0.800, 0.819] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | 0.017 [0.009, 0.024] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | 0.023 [0.015, 0.031] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.473 | 0.600 | 0.415 | 0.645 |
| 20,000 | 0.684 | 0.460 | 0.565 | 0.516 |
| 40,000 | 0.836 | 0.324 | 0.685 | 0.391 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 1,414, n_neg = 639, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.689**; CpG-island rows 0.505 (n_pos 652), non-CpG-island rows 0.544 (n_pos 762).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.200 | 0.417 [0.317, 0.562] | 0.421 [0.352, 0.501] | 0.929 [0.915, 0.943] | 0.868 [0.850, 0.887] |  |
| pydreg, −log10 p over all candidates | 0.023 | 0.206 [0.078, 0.328] | 0.216 [0.146, 0.307] | 0.882 [0.857, 0.904] | 0.839 [0.819, 0.861] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.291 | 0.540 [0.473, 0.619] | 0.503 [0.437, 0.579] | 0.941 [0.926, 0.953] | 0.887 [0.870, 0.903] |  |
| Δ pydreg(SVR) − PINTS | -0.090 | -0.122 [-0.245, 0.027] p=0.132 | -0.083 [-0.168, 0.008] p=0.084 | -0.012 [-0.024, 0.002] p=0.098 | -0.019 [-0.035, -0.002] p=0.034 |  |
| Δ pydreg(−log p) − PINTS | -0.268 | -0.334 [-0.503, -0.209] p<1e-03 | -0.287 [-0.382, -0.182] p<1e-03 | -0.059 [-0.080, -0.039] p<1e-03 | -0.048 [-0.066, -0.027] p<1e-03 |  |
| GC-matched: pydreg (SVR) |  | 0.280 [0.217, 0.536] | 0.325 [0.215, 0.463] | 0.908 [0.881, 0.932] | 0.836 [0.804, 0.866] |  |
| GC-matched: PINTS |  | 0.521 [0.280, 0.667] | 0.476 [0.317, 0.611] | 0.937 [0.915, 0.956] | 0.882 [0.853, 0.906] |  |
| GC-matched: Δ pydreg − PINTS |  | -0.241 [-0.383, 0.134] p=0.21 | -0.151 [-0.320, 0.057] p=0.146 | -0.029 [-0.056, -0.003] p=0.03 | -0.046 [-0.079, -0.015] p=0.002 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.371 [-0.625, 0.194] p=0.154 | -0.313 [-0.554, 0.078] p=0.114 | -0.077 [-0.157, 0.005] p=0.064 | -0.066 [-0.128, -0.009] p=0.02 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.025 [-0.139, 0.113] p=0.886 | -0.008 [-0.098, 0.083] p=0.832 | -0.014 [-0.046, 0.019] p=0.392 | -0.020 [-0.043, 0.003] p=0.08 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.850 [0.829, 0.869] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.865 [0.846, 0.884] |
| PINTS rule: PINTS |  |  |  |  |  | 0.886 [0.868, 0.902] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.036 [-0.053, -0.018] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | -0.020 [-0.036, -0.004] p=0.022 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.414 | 0.951 | 0.431 | 0.976 |
| 20,000 | 0.678 | 0.933 | 0.645 | 0.947 |
| 40,000 | 0.851 | 0.879 | 0.837 | 0.907 |

### K562_PROcap (K562 PRO-cap, ood)

n_pos = 1,322, n_neg = 444, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.749**; CpG-island rows 0.580 (n_pos 613), non-CpG-island rows 0.615 (n_pos 709).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.299 | 0.639 [0.574, 0.699] | 0.581 [0.521, 0.648] | 0.957 [0.948, 0.965] | 0.890 [0.873, 0.909] |  |
| pydreg, −log10 p over all candidates | 0.032 | 0.553 [0.124, 0.669] | 0.411 [0.290, 0.521] | 0.931 [0.911, 0.948] | 0.872 [0.850, 0.893] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.509 | 0.677 [0.599, 0.734] | 0.632 [0.578, 0.688] | 0.965 [0.957, 0.972] | 0.903 [0.887, 0.920] |  |
| Δ pydreg(SVR) − PINTS | -0.210 | -0.038 [-0.107, 0.049] p=0.41 | -0.050 [-0.121, 0.020] p=0.186 | -0.008 [-0.014, -0.001] p=0.018 | -0.013 [-0.027, -0.001] p=0.046 |  |
| Δ pydreg(−log p) − PINTS | -0.477 | -0.124 [-0.538, -0.001] p=0.046 | -0.220 [-0.340, -0.102] p=0.002 | -0.034 [-0.051, -0.019] p<1e-03 | -0.032 [-0.050, -0.016] p=0.002 |  |
| GC-matched: pydreg (SVR) |  | 0.459 [0.283, 0.642] | 0.431 [0.283, 0.578] | 0.943 [0.921, 0.959] | 0.869 [0.841, 0.895] |  |
| GC-matched: PINTS |  | 0.584 [0.499, 0.700] | 0.573 [0.499, 0.651] | 0.960 [0.950, 0.968] | 0.891 [0.870, 0.911] |  |
| GC-matched: Δ pydreg − PINTS |  | -0.125 [-0.308, 0.053] p=0.256 | -0.141 [-0.303, 0.010] p=0.062 | -0.017 [-0.037, -0.003] p=0.022 | -0.022 [-0.046, -0.000] p=0.05 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.399 [-0.831, -0.069] p=0.006 | -0.438 [-0.725, -0.145] p<1e-03 | -0.064 [-0.134, -0.014] p<1e-03 | -0.054 [-0.103, -0.014] p=0.002 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.015 [-0.116, 0.175] p=0.932 | 0.026 [-0.070, 0.130] p=0.618 | -0.000 [-0.023, 0.025] p=0.98 | -0.011 [-0.033, 0.011] p=0.354 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.864 [0.846, 0.882] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.888 [0.871, 0.907] |
| PINTS rule: PINTS |  |  |  |  |  | 0.902 [0.886, 0.919] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.039 [-0.052, -0.024] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | -0.015 [-0.029, -0.002] p=0.03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.408 | 0.989 | 0.439 | 0.995 |
| 20,000 | 0.657 | 0.969 | 0.632 | 0.979 |
| 40,000 | 0.821 | 0.926 | 0.801 | 0.947 |


## PINTS paper labels, no read filter

* **Positives:** As above.
* **Negatives:** As above.
* **Notes:** No read filter. This is closer to asking whether each tool recovers the orthogonally validated elements.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 1,602, n_neg = 5,952, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.212**; CpG-island rows 0.104 (n_pos 693), non-CpG-island rows 0.132 (n_pos 909).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.579 | 0.818 [0.791, 0.850] | 0.765 [0.734, 0.796] | 0.857 [0.834, 0.878] | 0.935 [0.924, 0.946] |  |
| pydreg, −log10 p over all candidates | 0.362 | 0.818 [0.787, 0.851] | 0.730 [0.699, 0.762] | 0.824 [0.801, 0.845] | 0.932 [0.921, 0.943] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.427 | 0.725 [0.682, 0.757] | 0.663 [0.630, 0.695] | 0.827 [0.806, 0.848] | 0.938 [0.928, 0.947] |  |
| Δ pydreg(SVR) − PINTS | 0.152 | 0.093 [0.064, 0.133] p<1e-03 | 0.102 [0.072, 0.132] p<1e-03 | 0.029 [0.007, 0.051] p=0.008 | -0.003 [-0.011, 0.005] p=0.536 |  |
| Δ pydreg(−log p) − PINTS | -0.064 | 0.093 [0.062, 0.132] p<1e-03 | 0.067 [0.036, 0.100] p<1e-03 | -0.004 [-0.028, 0.017] p=0.676 | -0.006 [-0.014, 0.001] p=0.11 |  |
| GC-matched: pydreg (SVR) |  | 0.765 [0.716, 0.814] | 0.703 [0.654, 0.750] | 0.826 [0.795, 0.856] | 0.926 [0.913, 0.939] |  |
| GC-matched: PINTS |  | 0.659 [0.606, 0.720] | 0.625 [0.578, 0.672] | 0.811 [0.780, 0.836] | 0.927 [0.915, 0.938] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.106 [0.040, 0.167] p=0.002 | 0.078 [0.025, 0.128] p=0.004 | 0.016 [-0.018, 0.048] p=0.386 | -0.001 [-0.012, 0.009] p=0.882 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.031 [-0.021, 0.065] p=0.316 | -0.014 [-0.093, 0.061] p=0.708 | -0.045 [-0.124, 0.030] p=0.286 | -0.003 [-0.011, 0.006] p=0.482 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.194 [0.114, 0.248] p<1e-03 | 0.162 [0.116, 0.209] p<1e-03 | 0.104 [0.056, 0.149] p<1e-03 | -0.002 [-0.017, 0.011] p=0.782 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.898 [0.884, 0.912] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.931 [0.918, 0.942] |
| PINTS rule: PINTS |  |  |  |  |  | 0.935 [0.925, 0.944] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.037 [-0.048, -0.027] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | -0.005 [-0.013, 0.003] p=0.266 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.384 | 0.956 | 0.301 | 0.940 |
| 20,000 | 0.591 | 0.937 | 0.426 | 0.912 |
| 40,000 | 0.746 | 0.874 | 0.551 | 0.878 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 1,602, n_neg = 5,952, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.212**; CpG-island rows 0.104 (n_pos 693), non-CpG-island rows 0.132 (n_pos 909).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.647 | 0.791 [0.766, 0.815] | 0.750 [0.724, 0.777] | 0.806 [0.785, 0.827] | 0.886 [0.873, 0.899] |  |
| pydreg, −log10 p over all candidates | 0.651 | 0.791 [0.766, 0.815] | 0.744 [0.718, 0.772] | 0.795 [0.773, 0.817] | 0.885 [0.873, 0.898] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.390 | 0.730 [0.695, 0.770] | 0.670 [0.643, 0.699] | 0.766 [0.741, 0.790] | 0.886 [0.874, 0.899] |  |
| Δ pydreg(SVR) − PINTS | 0.257 | 0.061 [0.023, 0.094] p=0.002 | 0.080 [0.058, 0.102] p<1e-03 | 0.041 [0.022, 0.059] p<1e-03 | -0.001 [-0.010, 0.007] p=0.836 |  |
| Δ pydreg(−log p) − PINTS | 0.261 | 0.061 [0.023, 0.094] p=0.002 | 0.074 [0.050, 0.096] p<1e-03 | 0.030 [0.009, 0.052] p=0.006 | -0.001 [-0.010, 0.007] p=0.746 |  |
| GC-matched: pydreg (SVR) |  | 0.789 [0.765, 0.814] | 0.737 [0.701, 0.769] | 0.799 [0.771, 0.823] | 0.884 [0.871, 0.897] |  |
| GC-matched: PINTS |  | 0.696 [0.647, 0.746] | 0.656 [0.618, 0.694] | 0.766 [0.739, 0.791] | 0.883 [0.870, 0.896] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.093 [0.048, 0.139] p<1e-03 | 0.082 [0.053, 0.109] p<1e-03 | 0.033 [0.013, 0.053] p=0.002 | 0.001 [-0.009, 0.009] p=0.912 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.083 [0.035, 0.135] p<1e-03 | 0.073 [0.036, 0.113] p<1e-03 | 0.055 [0.013, 0.089] p=0.004 | 0.014 [0.005, 0.023] p=0.004 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.094 [0.048, 0.139] p<1e-03 | 0.078 [0.048, 0.108] p<1e-03 | 0.067 [0.031, 0.100] p<1e-03 | -0.011 [-0.025, 0.000] p=0.056 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.828 [0.810, 0.845] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.882 [0.869, 0.895] |
| PINTS rule: PINTS |  |  |  |  |  | 0.884 [0.872, 0.897] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.057 [-0.072, -0.044] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | -0.002 [-0.012, 0.006] p=0.598 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.436 | 0.972 | 0.236 | 0.929 |
| 20,000 | 0.600 | 0.950 | 0.396 | 0.904 |
| 40,000 | 0.737 | 0.906 | 0.540 | 0.888 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 1,602, n_neg = 5,952, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.212**; CpG-island rows 0.104 (n_pos 693), non-CpG-island rows 0.132 (n_pos 909).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.642 | 0.787 [0.760, 0.813] | 0.753 [0.725, 0.782] | 0.809 [0.787, 0.830] | 0.885 [0.870, 0.898] |  |
| pydreg, −log10 p over all candidates | 0.649 | 0.787 [0.760, 0.813] | 0.750 [0.721, 0.778] | 0.804 [0.781, 0.826] | 0.884 [0.870, 0.898] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.549 | 0.734 [0.700, 0.768] | 0.694 [0.666, 0.724] | 0.799 [0.777, 0.819] | 0.897 [0.884, 0.910] |  |
| Δ pydreg(SVR) − PINTS | 0.093 | 0.053 [0.026, 0.082] p<1e-03 | 0.059 [0.038, 0.081] p<1e-03 | 0.010 [-0.004, 0.024] p=0.174 | -0.013 [-0.020, -0.005] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | 0.100 | 0.053 [0.026, 0.082] p<1e-03 | 0.056 [0.035, 0.078] p<1e-03 | 0.005 [-0.010, 0.021] p=0.518 | -0.013 [-0.021, -0.005] p=0.002 |  |
| GC-matched: pydreg (SVR) |  | 0.781 [0.740, 0.808] | 0.715 [0.675, 0.752] | 0.791 [0.764, 0.817] | 0.879 [0.864, 0.893] |  |
| GC-matched: PINTS |  | 0.676 [0.630, 0.723] | 0.646 [0.608, 0.685] | 0.772 [0.744, 0.802] | 0.889 [0.875, 0.903] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.105 [0.055, 0.145] p<1e-03 | 0.069 [0.029, 0.106] p<1e-03 | 0.019 [-0.006, 0.045] p=0.158 | -0.010 [-0.019, -0.001] p=0.038 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.027 [-0.002, 0.047] p=0.07 | 0.031 [-0.021, 0.082] p=0.232 | 0.046 [-0.016, 0.114] p=0.164 | 0.004 [-0.003, 0.011] p=0.21 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.117 [0.065, 0.173] p<1e-03 | 0.100 [0.061, 0.136] p<1e-03 | 0.068 [0.029, 0.110] p=0.002 | -0.023 [-0.037, -0.009] p=0.002 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.839 [0.823, 0.855] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.881 [0.866, 0.895] |
| PINTS rule: PINTS |  |  |  |  |  | 0.896 [0.882, 0.909] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.057 [-0.069, -0.046] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | -0.015 [-0.023, -0.007] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.420 | 0.988 | 0.337 | 0.978 |
| 20,000 | 0.574 | 0.961 | 0.433 | 0.952 |
| 40,000 | 0.712 | 0.925 | 0.582 | 0.910 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 6,994, n_neg = 75,311, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.085**; CpG-island rows 0.043 (n_pos 3,417), non-CpG-island rows 0.045 (n_pos 3,577).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.180 | 0.576 [0.558, 0.593] | 0.505 [0.490, 0.521] | 0.468 [0.450, 0.486] | 0.863 [0.856, 0.869] |  |
| pydreg, −log10 p over all candidates | 0.093 | 0.471 [0.447, 0.495] | 0.414 [0.399, 0.428] | 0.377 [0.363, 0.392] | 0.853 [0.846, 0.859] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.289 | 0.550 [0.530, 0.569] | 0.502 [0.485, 0.519] | 0.489 [0.463, 0.517] | 0.847 [0.840, 0.855] |  |
| Δ pydreg(SVR) − PINTS | -0.109 | 0.026 [0.009, 0.043] p=0.002 | 0.003 [-0.011, 0.017] p=0.752 | -0.022 [-0.048, 0.003] p=0.084 | 0.015 [0.011, 0.019] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | -0.196 | -0.079 [-0.102, -0.053] p<1e-03 | -0.089 [-0.104, -0.074] p<1e-03 | -0.112 [-0.137, -0.091] p<1e-03 | 0.005 [0.001, 0.010] p=0.018 |  |
| GC-matched: pydreg (SVR) |  | 0.579 [0.551, 0.600] | 0.510 [0.488, 0.530] | 0.481 [0.457, 0.506] | 0.867 [0.860, 0.874] |  |
| GC-matched: PINTS |  | 0.558 [0.532, 0.580] | 0.505 [0.483, 0.525] | 0.488 [0.460, 0.519] | 0.851 [0.842, 0.859] |  |
| GC-matched: Δ pydreg − PINTS |  | 0.021 [-0.006, 0.045] p=0.152 | 0.005 [-0.016, 0.025] p=0.646 | -0.007 [-0.038, 0.018] p=0.562 | 0.017 [0.012, 0.022] p<1e-03 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.005 [-0.037, 0.043] p=0.892 | -0.012 [-0.051, 0.027] p=0.554 | -0.026 [-0.103, 0.042] p=0.398 | 0.003 [-0.002, 0.007] p=0.274 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.037 [0.014, 0.057] p=0.002 | 0.026 [0.010, 0.042] p=0.002 | 0.011 [-0.002, 0.022] p=0.102 | 0.029 [0.021, 0.035] p<1e-03 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.827 [0.818, 0.835] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.840 [0.832, 0.848] |
| PINTS rule: PINTS |  |  |  |  |  | 0.838 [0.830, 0.846] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.011 [-0.016, -0.007] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | 0.002 [-0.003, 0.006] p=0.442 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.423 | 0.607 | 0.379 | 0.677 |
| 20,000 | 0.612 | 0.470 | 0.512 | 0.550 |
| 40,000 | 0.744 | 0.331 | 0.615 | 0.421 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 1,602, n_neg = 5,952, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.212**; CpG-island rows 0.104 (n_pos 693), non-CpG-island rows 0.132 (n_pos 909).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.582 | 0.846 [0.817, 0.872] | 0.785 [0.756, 0.815] | 0.863 [0.843, 0.883] | 0.935 [0.923, 0.946] |  |
| pydreg, −log10 p over all candidates | 0.404 | 0.842 [0.815, 0.868] | 0.752 [0.720, 0.783] | 0.816 [0.791, 0.843] | 0.931 [0.920, 0.943] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.631 | 0.866 [0.839, 0.891] | 0.806 [0.780, 0.832] | 0.882 [0.861, 0.901] | 0.944 [0.933, 0.955] |  |
| Δ pydreg(SVR) − PINTS | -0.049 | -0.020 [-0.038, -0.005] p=0.006 | -0.021 [-0.037, -0.005] p=0.014 | -0.019 [-0.032, -0.005] p=0.01 | -0.009 [-0.014, -0.004] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | -0.226 | -0.025 [-0.040, -0.008] p=0.002 | -0.054 [-0.075, -0.035] p<1e-03 | -0.065 [-0.086, -0.045] p<1e-03 | -0.012 [-0.018, -0.007] p<1e-03 |  |
| GC-matched: pydreg (SVR) |  | 0.797 [0.743, 0.839] | 0.713 [0.659, 0.763] | 0.823 [0.785, 0.858] | 0.926 [0.913, 0.938] |  |
| GC-matched: PINTS |  | 0.847 [0.806, 0.872] | 0.779 [0.739, 0.815] | 0.870 [0.843, 0.893] | 0.939 [0.927, 0.951] |  |
| GC-matched: Δ pydreg − PINTS |  | -0.050 [-0.081, -0.016] p<1e-03 | -0.066 [-0.103, -0.031] p<1e-03 | -0.047 [-0.077, -0.018] p=0.002 | -0.013 [-0.019, -0.007] p<1e-03 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.030 [-0.144, 0.000] p=0.068 | -0.110 [-0.196, -0.042] p<1e-03 | -0.121 [-0.213, -0.024] p=0.004 | -0.010 [-0.020, -0.002] p=0.018 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.042 [-0.081, -0.017] p<1e-03 | -0.029 [-0.054, -0.005] p=0.022 | -0.025 [-0.058, 0.009] p=0.122 | -0.016 [-0.024, -0.008] p<1e-03 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.889 [0.874, 0.903] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.932 [0.921, 0.944] |
| PINTS rule: PINTS |  |  |  |  |  | 0.943 [0.932, 0.954] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.054 [-0.065, -0.045] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | -0.011 [-0.016, -0.006] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.369 | 0.950 | 0.391 | 0.974 |
| 20,000 | 0.619 | 0.933 | 0.599 | 0.948 |
| 40,000 | 0.778 | 0.875 | 0.768 | 0.905 |

### K562_PROcap (K562 PRO-cap, ood)

n_pos = 1,602, n_neg = 5,997, promoter seeds = 5. AUPRC baseline (positive fraction) = **0.211**; CpG-island rows 0.103 (n_pos 686), non-CpG-island rows 0.133 (n_pos 916).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) | PINTS-rule AUC (replication) |
|---|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.675 | 0.858 [0.839, 0.877] | 0.806 [0.785, 0.828] | 0.861 [0.845, 0.878] | 0.920 [0.910, 0.931] |  |
| pydreg, −log10 p over all candidates | 0.666 | 0.858 [0.839, 0.877] | 0.790 [0.767, 0.815] | 0.834 [0.814, 0.856] | 0.918 [0.909, 0.929] |  |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.726 | 0.863 [0.843, 0.882] | 0.833 [0.813, 0.853] | 0.889 [0.875, 0.904] | 0.938 [0.929, 0.946] |  |
| Δ pydreg(SVR) − PINTS | -0.051 | -0.005 [-0.018, 0.008] p=0.454 | -0.027 [-0.039, -0.015] p<1e-03 | -0.029 [-0.039, -0.018] p<1e-03 | -0.018 [-0.024, -0.011] p<1e-03 |  |
| Δ pydreg(−log p) − PINTS | -0.060 | -0.005 [-0.018, 0.008] p=0.426 | -0.043 [-0.059, -0.027] p<1e-03 | -0.056 [-0.073, -0.039] p<1e-03 | -0.020 [-0.026, -0.013] p<1e-03 |  |
| GC-matched: pydreg (SVR) |  | 0.800 [0.749, 0.843] | 0.730 [0.681, 0.778] | 0.816 [0.779, 0.851] | 0.911 [0.901, 0.922] |  |
| GC-matched: PINTS |  | 0.814 [0.778, 0.855] | 0.775 [0.739, 0.810] | 0.866 [0.847, 0.885] | 0.930 [0.920, 0.939] |  |
| GC-matched: Δ pydreg − PINTS |  | -0.015 [-0.057, 0.017] p=0.348 | -0.045 [-0.077, -0.015] p<1e-03 | -0.050 [-0.081, -0.025] p<1e-03 | -0.019 [-0.026, -0.012] p<1e-03 |  |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.010 [-0.038, 0.002] p=0.134 | -0.099 [-0.180, -0.032] p<1e-03 | -0.119 [-0.213, -0.032] p<1e-03 | -0.009 [-0.018, -0.002] p=0.006 |  |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.005 [-0.035, 0.019] p=0.722 | -0.036 [-0.057, -0.015] p=0.004 | -0.047 [-0.071, -0.022] p=0.004 | -0.031 [-0.041, -0.019] p<1e-03 |  |
| PINTS rule: pydreg default calls only (as in the paper) |  |  |  |  |  | 0.869 [0.854, 0.883] |
| PINTS rule: pydreg all candidates |  |  |  |  |  | 0.918 [0.908, 0.929] |
| PINTS rule: PINTS |  |  |  |  |  | 0.937 [0.929, 0.946] |
| Δ PINTS rule, pydreg default − PINTS |  |  |  |  |  | -0.069 [-0.080, -0.058] p<1e-03 |
| Δ PINTS rule, pydreg full − PINTS |  |  |  |  |  | -0.019 [-0.026, -0.013] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.366 | 0.990 | 0.416 | 0.996 |
| 20,000 | 0.607 | 0.968 | 0.572 | 0.981 |
| 40,000 | 0.744 | 0.920 | 0.722 | 0.951 |


## dREG paper (Wang 2019) labels, original negatives

* **Positives:** Danko-lab hub `k562.positive.bed.gz`: Duke ∩ UW DHS ∩ GRO-cap HMM transcribed sites (27,798 merged).
* **Negatives:** Hub `k562.negative.bed.gz`: genome outside Duke DHS, UW DHS and GRO-cap HMM. Tiled into 500 bp tiles; only tiles holding ≥ 1 dREG informative position in this library are kept (≤ 200k, seed 0).
* **Notes:** These are the dREG authors' own labels. The positives require GRO-cap, so they are circular on K562 GRO-cap.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 27,798, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.122**; CpG-island rows 0.050 (n_pos 10,493), non-CpG-island rows 0.080 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.804 | 0.923 [0.920, 0.927] | 0.883 [0.879, 0.886] | 0.891 [0.886, 0.895] | 0.956 [0.954, 0.957] |
| pydreg, −log10 p over all candidates | 0.800 | 0.923 [0.920, 0.927] | 0.875 [0.871, 0.879] | 0.875 [0.870, 0.879] | 0.955 [0.953, 0.957] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.334 | 0.571 [0.562, 0.582] | 0.541 [0.531, 0.550] | 0.665 [0.657, 0.673] | 0.914 [0.912, 0.917] |
| Δ pydreg(SVR) − PINTS | 0.470 | 0.352 [0.342, 0.363] p<1e-03 | 0.342 [0.333, 0.350] p<1e-03 | 0.226 [0.218, 0.233] p<1e-03 | 0.041 [0.039, 0.044] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.466 | 0.352 [0.342, 0.363] p<1e-03 | 0.334 [0.325, 0.342] p<1e-03 | 0.209 [0.202, 0.216] p<1e-03 | 0.040 [0.038, 0.043] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.911 [0.905, 0.916] | 0.859 [0.851, 0.866] | 0.873 [0.861, 0.882] | 0.953 [0.951, 0.954] |
| GC-matched: PINTS |  | 0.505 [0.486, 0.522] | 0.473 [0.459, 0.490] | 0.598 [0.581, 0.615] | 0.881 [0.876, 0.887] |
| GC-matched: Δ pydreg − PINTS |  | 0.406 [0.389, 0.423] p<1e-03 | 0.386 [0.369, 0.402] p<1e-03 | 0.275 [0.256, 0.292] p<1e-03 | 0.071 [0.065, 0.077] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.208 [0.182, 0.231] p<1e-03 | 0.219 [0.193, 0.243] p<1e-03 | 0.226 [0.162, 0.285] p<1e-03 | 0.042 [0.037, 0.048] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.523 [0.508, 0.538] p<1e-03 | 0.483 [0.469, 0.496] p<1e-03 | 0.427 [0.411, 0.441] p<1e-03 | 0.080 [0.075, 0.086] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.294 | 0.980 | 0.170 | 0.931 |
| 20,000 | 0.494 | 0.972 | 0.259 | 0.876 |
| 40,000 | 0.689 | 0.950 | 0.359 | 0.799 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 27,798, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.122**; CpG-island rows 0.050 (n_pos 10,493), non-CpG-island rows 0.080 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.618 | 0.697 [0.690, 0.704] | 0.677 [0.671, 0.683] | 0.701 [0.695, 0.708] | 0.839 [0.835, 0.842] |
| pydreg, −log10 p over all candidates | 0.617 | 0.697 [0.690, 0.704] | 0.675 [0.668, 0.681] | 0.692 [0.685, 0.698] | 0.838 [0.835, 0.842] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.219 | 0.461 [0.453, 0.469] | 0.437 [0.428, 0.446] | 0.495 [0.487, 0.503] | 0.804 [0.800, 0.807] |
| Δ pydreg(SVR) − PINTS | 0.399 | 0.236 [0.227, 0.244] p<1e-03 | 0.240 [0.230, 0.249] p<1e-03 | 0.206 [0.199, 0.213] p<1e-03 | 0.035 [0.032, 0.038] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.399 | 0.236 [0.227, 0.244] p<1e-03 | 0.238 [0.228, 0.247] p<1e-03 | 0.196 [0.189, 0.203] p<1e-03 | 0.034 [0.032, 0.037] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.694 [0.687, 0.701] | 0.669 [0.663, 0.676] | 0.698 [0.691, 0.704] | 0.836 [0.833, 0.840] |
| GC-matched: PINTS |  | 0.396 [0.386, 0.408] | 0.363 [0.352, 0.375] | 0.435 [0.422, 0.447] | 0.776 [0.771, 0.781] |
| GC-matched: Δ pydreg − PINTS |  | 0.298 [0.286, 0.307] p<1e-03 | 0.306 [0.294, 0.316] p<1e-03 | 0.263 [0.251, 0.275] p<1e-03 | 0.061 [0.056, 0.065] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.292 [0.277, 0.305] p<1e-03 | 0.337 [0.316, 0.356] p<1e-03 | 0.408 [0.381, 0.433] p<1e-03 | 0.083 [0.077, 0.088] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.297 [0.286, 0.307] p<1e-03 | 0.284 [0.273, 0.293] p<1e-03 | 0.289 [0.276, 0.299] p<1e-03 | 0.043 [0.038, 0.047] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.288 | 0.988 | 0.134 | 0.811 |
| 20,000 | 0.458 | 0.970 | 0.232 | 0.725 |
| 40,000 | 0.615 | 0.898 | 0.356 | 0.687 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 27,798, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.122**; CpG-island rows 0.050 (n_pos 10,493), non-CpG-island rows 0.080 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.686 | 0.729 [0.717, 0.738] | 0.714 [0.702, 0.723] | 0.734 [0.723, 0.742] | 0.856 [0.850, 0.861] |
| pydreg, −log10 p over all candidates | 0.689 | 0.729 [0.717, 0.738] | 0.712 [0.701, 0.721] | 0.728 [0.718, 0.737] | 0.856 [0.850, 0.861] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.311 | 0.460 [0.450, 0.469] | 0.453 [0.443, 0.463] | 0.541 [0.532, 0.549] | 0.821 [0.815, 0.826] |
| Δ pydreg(SVR) − PINTS | 0.374 | 0.270 [0.261, 0.277] p<1e-03 | 0.261 [0.251, 0.271] p<1e-03 | 0.193 [0.186, 0.199] p<1e-03 | 0.035 [0.032, 0.038] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.377 | 0.270 [0.261, 0.277] p<1e-03 | 0.259 [0.249, 0.269] p<1e-03 | 0.187 [0.180, 0.194] p<1e-03 | 0.035 [0.032, 0.038] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.726 [0.714, 0.735] | 0.704 [0.693, 0.714] | 0.726 [0.715, 0.737] | 0.854 [0.847, 0.858] |
| GC-matched: PINTS |  | 0.406 [0.396, 0.415] | 0.384 [0.374, 0.395] | 0.480 [0.470, 0.491] | 0.788 [0.781, 0.794] |
| GC-matched: Δ pydreg − PINTS |  | 0.320 [0.311, 0.329] p<1e-03 | 0.320 [0.309, 0.329] p<1e-03 | 0.247 [0.236, 0.257] p<1e-03 | 0.066 [0.060, 0.070] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.224 [0.211, 0.238] p<1e-03 | 0.254 [0.238, 0.269] p<1e-03 | 0.275 [0.246, 0.301] p<1e-03 | 0.064 [0.059, 0.069] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.371 [0.360, 0.380] p<1e-03 | 0.353 [0.342, 0.362] p<1e-03 | 0.327 [0.315, 0.339] p<1e-03 | 0.058 [0.053, 0.063] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.284 | 0.989 | 0.186 | 0.941 |
| 20,000 | 0.453 | 0.976 | 0.251 | 0.841 |
| 40,000 | 0.628 | 0.940 | 0.387 | 0.710 |

### K562_GROcap (K562 GRO-cap, ood) — **circular, do not interpret**

n_pos = 27,798, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.122**; CpG-island rows 0.050 (n_pos 10,493), non-CpG-island rows 0.080 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.736 | 0.955 [0.951, 0.958] | 0.893 [0.889, 0.896] | 0.906 [0.901, 0.911] | 0.978 [0.977, 0.979] |
| pydreg, −log10 p over all candidates | 0.691 | 0.947 [0.942, 0.951] | 0.868 [0.863, 0.873] | 0.847 [0.838, 0.854] | 0.975 [0.974, 0.976] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.726 | 0.899 [0.895, 0.903] | 0.850 [0.845, 0.855] | 0.871 [0.866, 0.876] | 0.947 [0.945, 0.949] |
| Δ pydreg(SVR) − PINTS | 0.010 | 0.056 [0.052, 0.060] p<1e-03 | 0.043 [0.039, 0.047] p<1e-03 | 0.035 [0.031, 0.040] p<1e-03 | 0.031 [0.029, 0.033] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.036 | 0.048 [0.043, 0.052] p<1e-03 | 0.018 [0.013, 0.023] p<1e-03 | -0.024 [-0.030, -0.018] p<1e-03 | 0.029 [0.027, 0.030] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.963 [0.958, 0.967] | 0.906 [0.899, 0.912] | 0.916 [0.907, 0.923] | 0.979 [0.978, 0.980] |
| GC-matched: PINTS |  | 0.901 [0.897, 0.904] | 0.859 [0.852, 0.865] | 0.881 [0.873, 0.887] | 0.947 [0.945, 0.949] |
| GC-matched: Δ pydreg − PINTS |  | 0.062 [0.057, 0.067] p<1e-03 | 0.046 [0.040, 0.053] p<1e-03 | 0.035 [0.027, 0.042] p<1e-03 | 0.032 [0.031, 0.034] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.035 [0.031, 0.040] p<1e-03 | 0.024 [0.014, 0.032] p<1e-03 | -0.005 [-0.036, 0.020] p=0.762 | 0.018 [0.016, 0.021] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.072 [0.067, 0.080] p<1e-03 | 0.058 [0.052, 0.064] p<1e-03 | 0.048 [0.039, 0.057] p<1e-03 | 0.040 [0.038, 0.043] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.311 | 0.967 | 0.251 | 0.980 |
| 20,000 | 0.547 | 0.952 | 0.441 | 0.966 |
| 40,000 | 0.800 | 0.882 | 0.686 | 0.918 |


## dREG paper labels, negatives rebuilt with one recipe for K562 and GM12878

* **Positives:** K562: hub positives as above. GM12878: hub `GM12878.positive.bed.gz` (36,627).
* **Negatives:** Genome minus (Duke ∪ UW DHS ∪ nonzero GRO-cap signal ± 250 bp ∪ positives), then tiled as above. On K562 this recipe matches the original negatives (Jaccard 0.88 with the signal proxy used here, 0.96 with the HMM).
* **Notes:** The same recipe in both cells, so the K562 → GM12878 cross-cell comparison is like for like.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 27,798, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.122**; CpG-island rows 0.050 (n_pos 10,493), non-CpG-island rows 0.080 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.880 | 0.925 [0.921, 0.928] | 0.909 [0.906, 0.913] | 0.919 [0.915, 0.922] | 0.959 [0.957, 0.961] |
| pydreg, −log10 p over all candidates | 0.880 | 0.925 [0.921, 0.928] | 0.904 [0.901, 0.908] | 0.908 [0.904, 0.913] | 0.958 [0.957, 0.960] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.474 | 0.719 [0.708, 0.729] | 0.681 [0.672, 0.689] | 0.775 [0.767, 0.783] | 0.939 [0.937, 0.941] |
| Δ pydreg(SVR) − PINTS | 0.406 | 0.206 [0.195, 0.216] p<1e-03 | 0.229 [0.220, 0.237] p<1e-03 | 0.143 [0.137, 0.151] p<1e-03 | 0.020 [0.018, 0.022] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.406 | 0.206 [0.195, 0.216] p<1e-03 | 0.224 [0.215, 0.232] p<1e-03 | 0.133 [0.126, 0.141] p<1e-03 | 0.019 [0.017, 0.021] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.923 [0.920, 0.927] | 0.890 [0.883, 0.897] | 0.902 [0.895, 0.910] | 0.956 [0.954, 0.958] |
| GC-matched: PINTS |  | 0.659 [0.648, 0.673] | 0.624 [0.611, 0.637] | 0.727 [0.711, 0.743] | 0.921 [0.918, 0.924] |
| GC-matched: Δ pydreg − PINTS |  | 0.264 [0.249, 0.275] p<1e-03 | 0.267 [0.252, 0.281] p<1e-03 | 0.175 [0.159, 0.193] p<1e-03 | 0.035 [0.032, 0.039] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.109 [0.098, 0.120] p<1e-03 | 0.125 [0.107, 0.141] p<1e-03 | 0.145 [0.093, 0.199] p<1e-03 | 0.022 [0.019, 0.025] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.355 [0.336, 0.367] p<1e-03 | 0.351 [0.337, 0.365] p<1e-03 | 0.321 [0.304, 0.338] p<1e-03 | 0.040 [0.036, 0.044] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.294 | 0.997 | 0.170 | 0.982 |
| 20,000 | 0.494 | 0.995 | 0.259 | 0.968 |
| 40,000 | 0.689 | 0.988 | 0.359 | 0.933 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 27,798, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.122**; CpG-island rows 0.050 (n_pos 10,493), non-CpG-island rows 0.080 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.688 | 0.701 [0.694, 0.708] | 0.696 [0.690, 0.703] | 0.721 [0.714, 0.727] | 0.842 [0.839, 0.846] |
| pydreg, −log10 p over all candidates | 0.688 | 0.701 [0.694, 0.708] | 0.695 [0.689, 0.702] | 0.718 [0.711, 0.725] | 0.842 [0.839, 0.846] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.393 | 0.619 [0.600, 0.633] | 0.567 [0.559, 0.575] | 0.600 [0.590, 0.609] | 0.825 [0.821, 0.828] |
| Δ pydreg(SVR) − PINTS | 0.296 | 0.082 [0.069, 0.103] p<1e-03 | 0.129 [0.122, 0.136] p<1e-03 | 0.121 [0.114, 0.127] p<1e-03 | 0.018 [0.015, 0.020] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.296 | 0.082 [0.069, 0.103] p<1e-03 | 0.128 [0.121, 0.136] p<1e-03 | 0.118 [0.112, 0.125] p<1e-03 | 0.018 [0.015, 0.020] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.699 [0.693, 0.706] | 0.692 [0.686, 0.699] | 0.718 [0.712, 0.725] | 0.841 [0.837, 0.845] |
| GC-matched: PINTS |  | 0.504 [0.493, 0.516] | 0.495 [0.482, 0.507] | 0.554 [0.540, 0.566] | 0.807 [0.803, 0.811] |
| GC-matched: Δ pydreg − PINTS |  | 0.195 [0.183, 0.207] p<1e-03 | 0.198 [0.184, 0.211] p<1e-03 | 0.165 [0.154, 0.176] p<1e-03 | 0.034 [0.030, 0.038] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.212 [0.199, 0.223] p<1e-03 | 0.231 [0.213, 0.248] p<1e-03 | 0.270 [0.245, 0.294] p<1e-03 | 0.058 [0.053, 0.062] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.171 [0.158, 0.183] p<1e-03 | 0.169 [0.157, 0.180] p<1e-03 | 0.208 [0.197, 0.219] p<1e-03 | 0.016 [0.012, 0.020] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.288 | 1.000 | 0.134 | 0.957 |
| 20,000 | 0.458 | 0.997 | 0.232 | 0.916 |
| 40,000 | 0.615 | 0.971 | 0.356 | 0.897 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 27,798, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.122**; CpG-island rows 0.050 (n_pos 10,493), non-CpG-island rows 0.080 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.720 | 0.732 [0.721, 0.741] | 0.728 [0.718, 0.737] | 0.750 [0.740, 0.758] | 0.858 [0.853, 0.863] |
| pydreg, −log10 p over all candidates | 0.720 | 0.732 [0.721, 0.741] | 0.728 [0.718, 0.737] | 0.749 [0.739, 0.757] | 0.858 [0.853, 0.863] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.407 | 0.609 [0.593, 0.626] | 0.569 [0.558, 0.578] | 0.622 [0.612, 0.632] | 0.841 [0.835, 0.846] |
| Δ pydreg(SVR) − PINTS | 0.313 | 0.122 [0.105, 0.139] p<1e-03 | 0.159 [0.152, 0.167] p<1e-03 | 0.127 [0.121, 0.134] p<1e-03 | 0.018 [0.015, 0.020] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.313 | 0.122 [0.105, 0.139] p<1e-03 | 0.159 [0.152, 0.167] p<1e-03 | 0.127 [0.120, 0.133] p<1e-03 | 0.018 [0.015, 0.020] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.729 [0.719, 0.738] | 0.722 [0.711, 0.731] | 0.745 [0.735, 0.754] | 0.857 [0.851, 0.862] |
| GC-matched: PINTS |  | 0.489 [0.478, 0.502] | 0.493 [0.480, 0.507] | 0.573 [0.561, 0.585] | 0.821 [0.815, 0.826] |
| GC-matched: Δ pydreg − PINTS |  | 0.240 [0.227, 0.251] p<1e-03 | 0.228 [0.214, 0.241] p<1e-03 | 0.172 [0.161, 0.182] p<1e-03 | 0.036 [0.032, 0.039] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.181 [0.167, 0.192] p<1e-03 | 0.185 [0.168, 0.200] p<1e-03 | 0.196 [0.174, 0.217] p<1e-03 | 0.042 [0.038, 0.046] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.261 [0.247, 0.274] p<1e-03 | 0.243 [0.230, 0.255] p<1e-03 | 0.269 [0.257, 0.281] p<1e-03 | 0.027 [0.023, 0.031] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.284 | 0.999 | 0.186 | 0.988 |
| 20,000 | 0.453 | 0.998 | 0.251 | 0.963 |
| 40,000 | 0.628 | 0.988 | 0.387 | 0.896 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 36,627, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.155**; CpG-island rows 0.074 (n_pos 16,084), non-CpG-island rows 0.093 (n_pos 20,543).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.928 | 0.958 [0.955, 0.961] | 0.946 [0.943, 0.949] | 0.956 [0.953, 0.958] | 0.977 [0.975, 0.978] |
| pydreg, −log10 p over all candidates | 0.925 | 0.958 [0.955, 0.961] | 0.941 [0.937, 0.944] | 0.948 [0.945, 0.951] | 0.976 [0.975, 0.978] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.595 | 0.871 [0.865, 0.876] | 0.819 [0.813, 0.826] | 0.892 [0.888, 0.897] | 0.971 [0.970, 0.972] |
| Δ pydreg(SVR) − PINTS | 0.332 | 0.087 [0.083, 0.093] p<1e-03 | 0.127 [0.121, 0.132] p<1e-03 | 0.063 [0.059, 0.068] p<1e-03 | 0.006 [0.004, 0.007] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.329 | 0.087 [0.083, 0.093] p<1e-03 | 0.122 [0.116, 0.127] p<1e-03 | 0.055 [0.051, 0.060] p<1e-03 | 0.005 [0.004, 0.007] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.957 [0.955, 0.960] | 0.942 [0.937, 0.946] | 0.954 [0.951, 0.957] | 0.976 [0.975, 0.978] |
| GC-matched: PINTS |  | 0.830 [0.810, 0.844] | 0.779 [0.764, 0.794] | 0.875 [0.866, 0.883] | 0.964 [0.962, 0.967] |
| GC-matched: Δ pydreg − PINTS |  | 0.128 [0.113, 0.148] p<1e-03 | 0.163 [0.148, 0.177] p<1e-03 | 0.080 [0.071, 0.088] p<1e-03 | 0.012 [0.009, 0.014] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.090 [0.067, 0.113] p<1e-03 | 0.110 [0.091, 0.126] p<1e-03 | 0.090 [0.075, 0.104] p<1e-03 | 0.010 [0.007, 0.013] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.157 [0.147, 0.167] p<1e-03 | 0.198 [0.188, 0.208] p<1e-03 | 0.169 [0.157, 0.181] p<1e-03 | 0.012 [0.009, 0.014] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.283 | 0.998 | 0.210 | 0.992 |
| 20,000 | 0.492 | 0.997 | 0.335 | 0.983 |
| 40,000 | 0.721 | 0.994 | 0.473 | 0.959 |

### K562_GROcap (K562 GRO-cap, ood) — **circular, do not interpret**

n_pos = 27,798, n_neg = 20,963, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.570**; CpG-island rows 0.334 (n_pos 10,493), non-CpG-island rows 0.452 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.974 | 0.975 [0.973, 0.977] | 0.975 [0.973, 0.977] | 0.989 [0.988, 0.990] | 0.987 [0.986, 0.988] |
| pydreg, −log10 p over all candidates | 0.974 | 0.975 [0.973, 0.977] | 0.975 [0.973, 0.977] | 0.989 [0.988, 0.990] | 0.987 [0.986, 0.988] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.925 | 0.928 [0.925, 0.931] | 0.928 [0.925, 0.931] | 0.967 [0.966, 0.969] | 0.962 [0.960, 0.964] |
| Δ pydreg(SVR) − PINTS | 0.049 | 0.047 [0.044, 0.050] p<1e-03 | 0.047 [0.044, 0.050] p<1e-03 | 0.021 [0.020, 0.023] p<1e-03 | 0.025 [0.023, 0.027] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.049 | 0.047 [0.044, 0.050] p<1e-03 | 0.047 [0.044, 0.050] p<1e-03 | 0.021 [0.020, 0.023] p<1e-03 | 0.025 [0.023, 0.027] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.975 [0.973, 0.977] | 0.975 [0.973, 0.977] | 0.989 [0.988, 0.990] | 0.987 [0.986, 0.988] |
| GC-matched: PINTS |  | 0.928 [0.925, 0.931] | 0.928 [0.925, 0.931] | 0.967 [0.966, 0.969] | 0.962 [0.960, 0.964] |
| GC-matched: Δ pydreg − PINTS |  | 0.047 [0.044, 0.050] p<1e-03 | 0.047 [0.044, 0.050] p<1e-03 | 0.021 [0.020, 0.023] p<1e-03 | 0.025 [0.023, 0.027] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.028 [0.024, 0.031] p<1e-03 | 0.028 [0.024, 0.031] p<1e-03 | 0.020 [0.017, 0.022] p<1e-03 | 0.015 [0.013, 0.016] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.059 [0.054, 0.063] p<1e-03 | 0.059 [0.054, 0.063] p<1e-03 | 0.034 [0.031, 0.037] p<1e-03 | 0.031 [0.029, 0.033] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.311 | 1.000 | 0.251 | 1.000 |
| 20,000 | 0.547 | 1.000 | 0.441 | 1.000 |
| 40,000 | 0.800 | 1.000 | 0.686 | 1.000 |


## dREG auxiliary labels: dhs_h3k27ac vs quiet

* **Positives:** Duke ∩ UW DHS overlapping an ENCODE Broad H3K27ac peak, the "DHSs marked by H3K27ac" of Wang 2019 Fig 1B (K562 34,763; GM12878 55,925).
* **Negatives:** Genome minus (any Duke/UW DHS ∪ H3K27ac ∪ ChromHMM states 1–8), in 500 bp tiles holding ≥ 1 dREG informative position (≤ 200k), like the Danko 2015 tuning negatives.
* **Notes:** Not in either dREG paper as a ROC: the papers report recall only, so the negatives are ours. No nascent or TSS data is used, so these labels are not circular for any library, GRO-cap included. H3K27ac peaks and ChromHMM segments are broad, which favours wider calls under the overlap rule. On hg38 cells (HCT116, A549) the DHS is one ENCODE DNase peak set, not Duke ∩ UW, and there is no ChromHMM: quiet excludes only DNase and H3K27ac. Jurkat (hg19) uses UW DNase rep1 ∩ rep2 with H3K27ac from either ChIP-Atlas replicate (GSM4551183/4), and has no ChromHMM.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 34,763, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.148**; CpG-island rows 0.048 (n_pos 9,989), non-CpG-island rows 0.110 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.717 | 0.789 [0.780, 0.797] | 0.769 [0.761, 0.777] | 0.799 [0.792, 0.807] | 0.887 [0.882, 0.891] |
| pydreg, −log10 p over all candidates | 0.717 | 0.789 [0.780, 0.797] | 0.766 [0.758, 0.774] | 0.792 [0.785, 0.800] | 0.886 [0.882, 0.891] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.303 | 0.537 [0.525, 0.549] | 0.511 [0.500, 0.522] | 0.670 [0.660, 0.678] | 0.891 [0.888, 0.895] |
| Δ pydreg(SVR) − PINTS | 0.414 | 0.251 [0.240, 0.263] p<1e-03 | 0.258 [0.248, 0.268] p<1e-03 | 0.130 [0.121, 0.140] p<1e-03 | -0.005 [-0.008, -0.002] p=0.004 |
| Δ pydreg(−log p) − PINTS | 0.414 | 0.251 [0.240, 0.263] p<1e-03 | 0.255 [0.245, 0.265] p<1e-03 | 0.123 [0.114, 0.132] p<1e-03 | -0.005 [-0.009, -0.002] p=0.002 |
| GC-matched: pydreg (SVR) |  | 0.785 [0.777, 0.793] | 0.751 [0.743, 0.760] | 0.788 [0.780, 0.796] | 0.883 [0.879, 0.888] |
| GC-matched: PINTS |  | 0.467 [0.452, 0.481] | 0.443 [0.429, 0.457] | 0.606 [0.591, 0.620] | 0.861 [0.855, 0.866] |
| GC-matched: Δ pydreg − PINTS |  | 0.318 [0.301, 0.336] p<1e-03 | 0.308 [0.292, 0.326] p<1e-03 | 0.182 [0.167, 0.199] p<1e-03 | 0.023 [0.017, 0.029] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.150 [0.129, 0.165] p<1e-03 | 0.172 [0.151, 0.196] p<1e-03 | 0.219 [0.170, 0.266] p<1e-03 | 0.024 [0.019, 0.028] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.384 [0.369, 0.400] p<1e-03 | 0.363 [0.347, 0.378] p<1e-03 | 0.290 [0.272, 0.306] p<1e-03 | 0.016 [0.010, 0.022] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.227 | 0.998 | 0.137 | 0.955 |
| 20,000 | 0.379 | 0.995 | 0.210 | 0.913 |
| 40,000 | 0.534 | 0.984 | 0.296 | 0.848 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 34,763, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.148**; CpG-island rows 0.048 (n_pos 9,989), non-CpG-island rows 0.110 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.515 | 0.552 [0.542, 0.562] | 0.545 [0.535, 0.554] | 0.595 [0.587, 0.604] | 0.764 [0.758, 0.769] |
| pydreg, −log10 p over all candidates | 0.515 | 0.552 [0.542, 0.562] | 0.544 [0.535, 0.553] | 0.591 [0.583, 0.600] | 0.764 [0.758, 0.769] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.200 | 0.416 [0.407, 0.426] | 0.401 [0.391, 0.411] | 0.479 [0.469, 0.489] | 0.758 [0.752, 0.762] |
| Δ pydreg(SVR) − PINTS | 0.314 | 0.136 [0.127, 0.145] p<1e-03 | 0.144 [0.134, 0.152] p<1e-03 | 0.116 [0.108, 0.123] p<1e-03 | 0.006 [0.003, 0.009] p=0.002 |
| Δ pydreg(−log p) − PINTS | 0.315 | 0.136 [0.127, 0.145] p<1e-03 | 0.143 [0.133, 0.151] p<1e-03 | 0.112 [0.104, 0.119] p<1e-03 | 0.006 [0.003, 0.009] p=0.002 |
| GC-matched: pydreg (SVR) |  | 0.548 [0.538, 0.558] | 0.536 [0.527, 0.545] | 0.590 [0.582, 0.598] | 0.761 [0.756, 0.766] |
| GC-matched: PINTS |  | 0.356 [0.342, 0.367] | 0.328 [0.319, 0.340] | 0.427 [0.416, 0.439] | 0.733 [0.728, 0.739] |
| GC-matched: Δ pydreg − PINTS |  | 0.192 [0.182, 0.209] p<1e-03 | 0.208 [0.198, 0.219] p<1e-03 | 0.163 [0.152, 0.173] p<1e-03 | 0.028 [0.024, 0.032] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.257 [0.237, 0.274] p<1e-03 | 0.301 [0.279, 0.321] p<1e-03 | 0.384 [0.351, 0.412] p<1e-03 | 0.064 [0.058, 0.070] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.168 [0.159, 0.177] p<1e-03 | 0.170 [0.160, 0.178] p<1e-03 | 0.172 [0.163, 0.181] p<1e-03 | 0.009 [0.005, 0.013] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.224 | 0.998 | 0.109 | 0.865 |
| 20,000 | 0.351 | 0.990 | 0.189 | 0.781 |
| 40,000 | 0.473 | 0.946 | 0.292 | 0.746 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 34,763, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.148**; CpG-island rows 0.048 (n_pos 9,989), non-CpG-island rows 0.110 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.581 | 0.598 [0.586, 0.609] | 0.594 [0.582, 0.605] | 0.638 [0.626, 0.648] | 0.788 [0.782, 0.794] |
| pydreg, −log10 p over all candidates | 0.581 | 0.598 [0.586, 0.609] | 0.594 [0.582, 0.605] | 0.637 [0.625, 0.647] | 0.788 [0.782, 0.794] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.298 | 0.420 [0.410, 0.430] | 0.423 [0.413, 0.433] | 0.529 [0.519, 0.539] | 0.788 [0.781, 0.794] |
| Δ pydreg(SVR) − PINTS | 0.283 | 0.178 [0.168, 0.187] p<1e-03 | 0.171 [0.161, 0.181] p<1e-03 | 0.108 [0.101, 0.115] p<1e-03 | 0.000 [-0.003, 0.004] p=0.756 |
| Δ pydreg(−log p) − PINTS | 0.283 | 0.178 [0.168, 0.187] p<1e-03 | 0.171 [0.161, 0.180] p<1e-03 | 0.108 [0.101, 0.114] p<1e-03 | 0.000 [-0.003, 0.004] p=0.772 |
| GC-matched: pydreg (SVR) |  | 0.595 [0.583, 0.607] | 0.588 [0.576, 0.599] | 0.634 [0.622, 0.644] | 0.786 [0.780, 0.792] |
| GC-matched: PINTS |  | 0.358 [0.349, 0.367] | 0.344 [0.335, 0.353] | 0.471 [0.459, 0.482] | 0.758 [0.751, 0.764] |
| GC-matched: Δ pydreg − PINTS |  | 0.237 [0.227, 0.246] p<1e-03 | 0.244 [0.233, 0.254] p<1e-03 | 0.163 [0.151, 0.172] p<1e-03 | 0.029 [0.023, 0.033] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.188 [0.175, 0.200] p<1e-03 | 0.222 [0.204, 0.239] p<1e-03 | 0.258 [0.227, 0.289] p<1e-03 | 0.045 [0.039, 0.049] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.248 [0.237, 0.257] p<1e-03 | 0.246 [0.234, 0.255] p<1e-03 | 0.221 [0.210, 0.231] p<1e-03 | 0.014 [0.010, 0.019] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.221 | 0.999 | 0.149 | 0.960 |
| 20,000 | 0.351 | 0.997 | 0.205 | 0.885 |
| 40,000 | 0.494 | 0.984 | 0.321 | 0.768 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 55,925, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.219**; CpG-island rows 0.081 (n_pos 17,611), non-CpG-island rows 0.161 (n_pos 38,314).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.756 | 0.835 [0.829, 0.840] | 0.811 [0.806, 0.816] | 0.857 [0.853, 0.861] | 0.911 [0.908, 0.913] |
| pydreg, −log10 p over all candidates | 0.753 | 0.835 [0.829, 0.840] | 0.805 [0.800, 0.811] | 0.849 [0.845, 0.854] | 0.910 [0.907, 0.913] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.374 | 0.657 [0.646, 0.667] | 0.611 [0.601, 0.620] | 0.800 [0.794, 0.807] | 0.916 [0.914, 0.919] |
| Δ pydreg(SVR) − PINTS | 0.382 | 0.178 [0.169, 0.188] p<1e-03 | 0.200 [0.192, 0.208] p<1e-03 | 0.057 [0.051, 0.062] p<1e-03 | -0.006 [-0.008, -0.004] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.378 | 0.178 [0.169, 0.188] p<1e-03 | 0.194 [0.186, 0.202] p<1e-03 | 0.049 [0.044, 0.055] p<1e-03 | -0.007 [-0.009, -0.004] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.825 [0.818, 0.833] | 0.796 [0.789, 0.804] | 0.852 [0.847, 0.856] | 0.907 [0.904, 0.910] |
| GC-matched: PINTS |  | 0.562 [0.548, 0.575] | 0.537 [0.525, 0.550] | 0.758 [0.749, 0.767] | 0.896 [0.892, 0.899] |
| GC-matched: Δ pydreg − PINTS |  | 0.263 [0.252, 0.276] p<1e-03 | 0.260 [0.248, 0.270] p<1e-03 | 0.093 [0.085, 0.101] p<1e-03 | 0.012 [0.009, 0.015] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.168 [0.146, 0.183] p<1e-03 | 0.183 [0.165, 0.201] p<1e-03 | 0.142 [0.126, 0.159] p<1e-03 | 0.017 [0.014, 0.021] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.307 [0.297, 0.316] p<1e-03 | 0.287 [0.278, 0.297] p<1e-03 | 0.159 [0.149, 0.169] p<1e-03 | 0.004 [0.001, 0.007] p=0.01 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.195 | 0.998 | 0.147 | 0.984 |
| 20,000 | 0.349 | 0.997 | 0.236 | 0.966 |
| 40,000 | 0.542 | 0.993 | 0.344 | 0.923 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 34,763, n_neg = 126,437, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.216**; CpG-island rows 0.073 (n_pos 9,989), non-CpG-island rows 0.164 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.649 | 0.793 [0.784, 0.802] | 0.754 [0.745, 0.764] | 0.815 [0.809, 0.822] | 0.887 [0.882, 0.892] |
| pydreg, −log10 p over all candidates | 0.651 | 0.793 [0.784, 0.802] | 0.747 [0.738, 0.757] | 0.799 [0.792, 0.806] | 0.886 [0.882, 0.891] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.707 | 0.769 [0.761, 0.778] | 0.754 [0.746, 0.763] | 0.828 [0.822, 0.834] | 0.891 [0.887, 0.896] |
| Δ pydreg(SVR) − PINTS | -0.059 | 0.023 [0.020, 0.027] p<1e-03 | 0.000 [-0.004, 0.004] p=0.86 | -0.013 [-0.016, -0.010] p<1e-03 | -0.004 [-0.006, -0.003] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.057 | 0.023 [0.020, 0.027] p<1e-03 | -0.007 [-0.011, -0.003] p=0.002 | -0.029 [-0.034, -0.025] p<1e-03 | -0.005 [-0.007, -0.003] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.792 [0.783, 0.801] | 0.748 [0.736, 0.760] | 0.812 [0.804, 0.821] | 0.886 [0.881, 0.891] |
| GC-matched: PINTS |  | 0.760 [0.752, 0.769] | 0.744 [0.734, 0.754] | 0.823 [0.816, 0.830] | 0.888 [0.883, 0.892] |
| GC-matched: Δ pydreg − PINTS |  | 0.031 [0.027, 0.036] p<1e-03 | 0.004 [-0.004, 0.011] p=0.296 | -0.010 [-0.015, -0.005] p=0.002 | -0.002 [-0.004, 0.000] p=0.114 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.018 [0.014, 0.021] p<1e-03 | 0.005 [-0.005, 0.015] p=0.36 | 0.005 [-0.013, 0.030] p=0.694 | 0.005 [0.003, 0.007] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.035 [0.030, 0.040] p<1e-03 | 0.004 [-0.002, 0.010] p=0.196 | -0.016 [-0.022, -0.010] p<1e-03 | -0.005 [-0.008, -0.003] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.234 | 0.997 | 0.197 | 0.995 |
| 20,000 | 0.406 | 0.990 | 0.340 | 0.995 |
| 40,000 | 0.594 | 0.963 | 0.525 | 0.991 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 32,659, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.140**; CpG-island rows 0.058 (n_pos 12,409), non-CpG-island rows 0.092 (n_pos 20,250).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.575 | 0.624 [0.615, 0.633] | 0.603 [0.594, 0.612] | 0.633 [0.623, 0.643] | 0.800 [0.795, 0.805] |
| pydreg, −log10 p over all candidates | 0.575 | 0.624 [0.615, 0.633] | 0.602 [0.593, 0.611] | 0.628 [0.618, 0.637] | 0.800 [0.795, 0.805] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.244 | 0.457 [0.448, 0.470] | 0.431 [0.423, 0.438] | 0.514 [0.504, 0.523] | 0.777 [0.772, 0.782] |
| Δ pydreg(SVR) − PINTS | 0.331 | 0.168 [0.151, 0.176] p<1e-03 | 0.173 [0.164, 0.181] p<1e-03 | 0.119 [0.112, 0.125] p<1e-03 | 0.023 [0.020, 0.026] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.331 | 0.168 [0.151, 0.176] p<1e-03 | 0.171 [0.163, 0.180] p<1e-03 | 0.114 [0.107, 0.120] p<1e-03 | 0.023 [0.020, 0.026] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.573 [0.534, 0.598] | 0.485 [0.461, 0.503] | 0.514 [0.491, 0.534] | 0.780 [0.774, 0.786] |
| GC-matched: PINTS |  | 0.298 [0.289, 0.310] | 0.277 [0.266, 0.288] | 0.372 [0.358, 0.385] | 0.730 [0.723, 0.735] |
| GC-matched: Δ pydreg − PINTS |  | 0.276 [0.240, 0.297] p<1e-03 | 0.208 [0.186, 0.227] p<1e-03 | 0.142 [0.124, 0.159] p<1e-03 | 0.051 [0.046, 0.056] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.293 [0.229, 0.341] p<1e-03 | 0.253 [0.210, 0.292] p<1e-03 | 0.218 [0.177, 0.257] p<1e-03 | 0.072 [0.066, 0.078] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.226 [0.214, 0.242] p<1e-03 | 0.199 [0.189, 0.209] p<1e-03 | 0.140 [0.129, 0.151] p<1e-03 | 0.032 [0.028, 0.036] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.286 | 0.952 | 0.121 | 0.890 |
| 20,000 | 0.442 | 0.937 | 0.202 | 0.837 |
| 40,000 | 0.590 | 0.892 | 0.317 | 0.745 |

### HCT116_PROcap (HCT116 PRO-cap, ood)

n_pos = 32,659, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.140**; CpG-island rows 0.058 (n_pos 12,409), non-CpG-island rows 0.092 (n_pos 20,250).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.208 | 0.695 [0.681, 0.708] | 0.609 [0.598, 0.619] | 0.666 [0.657, 0.675] | 0.909 [0.906, 0.912] |
| pydreg, −log10 p over all candidates | 0.064 | 0.683 [0.669, 0.696] | 0.575 [0.564, 0.586] | 0.604 [0.594, 0.613] | 0.906 [0.903, 0.908] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.427 | 0.784 [0.775, 0.792] | 0.709 [0.701, 0.717] | 0.764 [0.756, 0.772] | 0.917 [0.914, 0.920] |
| Δ pydreg(SVR) − PINTS | -0.219 | -0.088 [-0.098, -0.079] p<1e-03 | -0.100 [-0.107, -0.094] p<1e-03 | -0.099 [-0.105, -0.092] p<1e-03 | -0.008 [-0.010, -0.007] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.363 | -0.100 [-0.111, -0.091] p<1e-03 | -0.134 [-0.141, -0.127] p<1e-03 | -0.161 [-0.167, -0.155] p<1e-03 | -0.012 [-0.013, -0.010] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.346 [0.316, 0.374] | 0.337 [0.313, 0.363] | 0.474 [0.456, 0.493] | 0.867 [0.861, 0.873] |
| GC-matched: PINTS |  | 0.469 [0.441, 0.503] | 0.447 [0.422, 0.471] | 0.568 [0.548, 0.587] | 0.877 [0.872, 0.882] |
| GC-matched: Δ pydreg − PINTS |  | -0.124 [-0.145, -0.103] p<1e-03 | -0.109 [-0.122, -0.096] p<1e-03 | -0.094 [-0.107, -0.080] p<1e-03 | -0.010 [-0.012, -0.007] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.327 [-0.368, -0.293] p<1e-03 | -0.264 [-0.290, -0.234] p<1e-03 | -0.190 [-0.218, -0.158] p<1e-03 | -0.042 [-0.048, -0.037] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.132 [-0.146, -0.118] p<1e-03 | -0.088 [-0.097, -0.080] p<1e-03 | -0.071 [-0.081, -0.062] p<1e-03 | -0.004 [-0.007, -0.002] p=0.002 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.240 | 0.773 | 0.250 | 0.912 |
| 20,000 | 0.439 | 0.759 | 0.417 | 0.878 |
| 40,000 | 0.661 | 0.707 | 0.616 | 0.825 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 83,074, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.293**; CpG-island rows 0.104 (n_pos 23,272), non-CpG-island rows 0.230 (n_pos 59,802).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.388 | 0.427 [0.419, 0.434] | 0.420 [0.412, 0.427] | 0.575 [0.570, 0.579] | 0.698 [0.693, 0.701] |
| pydreg, −log10 p over all candidates | 0.390 | 0.427 [0.419, 0.434] | 0.419 [0.412, 0.427] | 0.573 [0.568, 0.577] | 0.698 [0.693, 0.701] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.167 | 0.267 [0.260, 0.273] | 0.261 [0.254, 0.267] | 0.498 [0.491, 0.505] | 0.657 [0.651, 0.662] |
| Δ pydreg(SVR) − PINTS | 0.221 | 0.160 [0.154, 0.165] p<1e-03 | 0.159 [0.153, 0.166] p<1e-03 | 0.077 [0.072, 0.082] p<1e-03 | 0.041 [0.038, 0.044] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.222 | 0.160 [0.154, 0.165] p<1e-03 | 0.158 [0.152, 0.165] p<1e-03 | 0.075 [0.070, 0.079] p<1e-03 | 0.041 [0.038, 0.044] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.425 [0.417, 0.432] | 0.417 [0.409, 0.424] | 0.574 [0.569, 0.579] | 0.696 [0.693, 0.700] |
| GC-matched: PINTS |  | 0.238 [0.231, 0.245] | 0.224 [0.216, 0.232] | 0.470 [0.462, 0.480] | 0.634 [0.628, 0.640] |
| GC-matched: Δ pydreg − PINTS |  | 0.187 [0.180, 0.192] p<1e-03 | 0.193 [0.184, 0.200] p<1e-03 | 0.103 [0.095, 0.110] p<1e-03 | 0.062 [0.057, 0.067] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.273 [0.257, 0.292] p<1e-03 | 0.311 [0.289, 0.333] p<1e-03 | 0.292 [0.263, 0.318] p<1e-03 | 0.076 [0.068, 0.083] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.153 [0.148, 0.158] p<1e-03 | 0.150 [0.144, 0.155] p<1e-03 | 0.114 [0.108, 0.119] p<1e-03 | 0.050 [0.047, 0.055] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.138 | 0.997 | 0.084 | 0.926 |
| 20,000 | 0.243 | 0.993 | 0.145 | 0.878 |
| 40,000 | 0.376 | 0.954 | 0.304 | 0.644 |

### A549_PROcap (A549 PRO-cap, ood)

n_pos = 83,074, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.293**; CpG-island rows 0.104 (n_pos 23,272), non-CpG-island rows 0.230 (n_pos 59,802).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.132 | 0.256 [0.246, 0.265] | 0.236 [0.227, 0.245] | 0.457 [0.449, 0.464] | 0.614 [0.607, 0.620] |
| pydreg, −log10 p over all candidates | 0.071 | 0.256 [0.245, 0.266] | 0.226 [0.217, 0.234] | 0.433 [0.425, 0.441] | 0.613 [0.606, 0.619] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.167 | 0.269 [0.257, 0.279] | 0.252 [0.242, 0.261] | 0.471 [0.462, 0.478] | 0.590 [0.583, 0.596] |
| Δ pydreg(SVR) − PINTS | -0.035 | -0.013 [-0.016, -0.010] p<1e-03 | -0.016 [-0.018, -0.013] p<1e-03 | -0.013 [-0.016, -0.011] p<1e-03 | 0.024 [0.022, 0.026] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.096 | -0.013 [-0.016, -0.009] p<1e-03 | -0.026 [-0.028, -0.023] p<1e-03 | -0.038 [-0.040, -0.035] p<1e-03 | 0.023 [0.021, 0.024] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.223 [0.212, 0.234] | 0.204 [0.194, 0.214] | 0.428 [0.419, 0.437] | 0.605 [0.598, 0.612] |
| GC-matched: PINTS |  | 0.239 [0.228, 0.249] | 0.221 [0.211, 0.231] | 0.440 [0.430, 0.450] | 0.571 [0.564, 0.577] |
| GC-matched: Δ pydreg − PINTS |  | -0.016 [-0.021, -0.012] p<1e-03 | -0.017 [-0.021, -0.014] p<1e-03 | -0.012 [-0.017, -0.006] p<1e-03 | 0.035 [0.032, 0.037] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.057 [-0.079, -0.024] p<1e-03 | -0.044 [-0.065, -0.022] p<1e-03 | -0.035 [-0.070, -0.001] p=0.05 | 0.015 [0.011, 0.019] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.010 [-0.013, -0.007] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 | -0.005 [-0.007, -0.002] p<1e-03 | 0.038 [0.036, 0.040] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.112 | 0.862 | 0.100 | 0.922 |
| 20,000 | 0.187 | 0.780 | 0.172 | 0.867 |
| 40,000 | 0.272 | 0.655 | 0.249 | 0.727 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 15,613, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.072**; CpG-island rows 0.056 (n_pos 11,892), non-CpG-island rows 0.018 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.694 | 0.875 [0.867, 0.884] | 0.822 [0.813, 0.832] | 0.765 [0.756, 0.776] | 0.929 [0.925, 0.933] |
| pydreg, −log10 p over all candidates | 0.683 | 0.875 [0.867, 0.884] | 0.818 [0.809, 0.828] | 0.751 [0.742, 0.761] | 0.928 [0.924, 0.933] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.502 | 0.719 [0.709, 0.731] | 0.668 [0.656, 0.681] | 0.655 [0.645, 0.666] | 0.916 [0.912, 0.920] |
| Δ pydreg(SVR) − PINTS | 0.193 | 0.155 [0.149, 0.162] p<1e-03 | 0.154 [0.148, 0.160] p<1e-03 | 0.110 [0.102, 0.117] p<1e-03 | 0.013 [0.010, 0.016] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.181 | 0.155 [0.149, 0.162] p<1e-03 | 0.149 [0.143, 0.156] p<1e-03 | 0.096 [0.088, 0.104] p<1e-03 | 0.013 [0.010, 0.015] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.847 [0.826, 0.864] | 0.755 [0.734, 0.775] | 0.705 [0.677, 0.732] | 0.920 [0.915, 0.925] |
| GC-matched: PINTS |  | 0.569 [0.557, 0.593] | 0.538 [0.513, 0.563] | 0.487 [0.452, 0.530] | 0.868 [0.862, 0.875] |
| GC-matched: Δ pydreg − PINTS |  | 0.278 [0.252, 0.296] p<1e-03 | 0.217 [0.196, 0.237] p<1e-03 | 0.218 [0.176, 0.249] p<1e-03 | 0.052 [0.047, 0.057] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.261 [0.235, 0.278] p<1e-03 | 0.215 [0.190, 0.237] p<1e-03 | 0.246 [0.188, 0.287] p<1e-03 | 0.050 [0.045, 0.055] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.330 [0.309, 0.363] p<1e-03 | 0.245 [0.228, 0.261] p<1e-03 | 0.123 [0.101, 0.142] p<1e-03 | 0.051 [0.044, 0.059] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.453 | 0.884 | 0.310 | 0.886 |
| 20,000 | 0.668 | 0.851 | 0.445 | 0.798 |
| 40,000 | 0.831 | 0.784 | 0.590 | 0.670 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 15,613, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.072**; CpG-island rows 0.056 (n_pos 11,892), non-CpG-island rows 0.018 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.715 | 0.892 [0.883, 0.900] | 0.838 [0.830, 0.847] | 0.783 [0.775, 0.792] | 0.938 [0.933, 0.942] |
| pydreg, −log10 p over all candidates | 0.697 | 0.892 [0.883, 0.900] | 0.833 [0.825, 0.841] | 0.765 [0.756, 0.774] | 0.937 [0.933, 0.942] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.572 | 0.790 [0.779, 0.802] | 0.734 [0.722, 0.747] | 0.719 [0.710, 0.730] | 0.931 [0.927, 0.935] |
| Δ pydreg(SVR) − PINTS | 0.143 | 0.102 [0.095, 0.108] p<1e-03 | 0.104 [0.097, 0.111] p<1e-03 | 0.064 [0.056, 0.072] p<1e-03 | 0.007 [0.004, 0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.125 | 0.102 [0.095, 0.108] p<1e-03 | 0.099 [0.092, 0.105] p<1e-03 | 0.046 [0.038, 0.054] p<1e-03 | 0.006 [0.004, 0.009] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.850 [0.835, 0.867] | 0.769 [0.748, 0.790] | 0.716 [0.686, 0.746] | 0.928 [0.924, 0.933] |
| GC-matched: PINTS |  | 0.687 [0.670, 0.718] | 0.632 [0.609, 0.659] | 0.597 [0.555, 0.642] | 0.895 [0.889, 0.902] |
| GC-matched: Δ pydreg − PINTS |  | 0.162 [0.128, 0.184] p<1e-03 | 0.137 [0.108, 0.161] p<1e-03 | 0.120 [0.075, 0.163] p<1e-03 | 0.033 [0.028, 0.038] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.148 [0.115, 0.169] p<1e-03 | 0.129 [0.097, 0.156] p<1e-03 | 0.134 [0.075, 0.189] p<1e-03 | 0.031 [0.026, 0.036] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.223 [0.203, 0.249] p<1e-03 | 0.175 [0.153, 0.193] p<1e-03 | 0.079 [0.049, 0.102] p<1e-03 | 0.033 [0.025, 0.040] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.461 | 0.888 | 0.309 | 0.899 |
| 20,000 | 0.681 | 0.856 | 0.492 | 0.857 |
| 40,000 | 0.846 | 0.788 | 0.675 | 0.726 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 15,613, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.072**; CpG-island rows 0.056 (n_pos 11,892), non-CpG-island rows 0.018 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.650 | 0.941 [0.936, 0.945] | 0.859 [0.852, 0.865] | 0.787 [0.780, 0.795] | 0.961 [0.958, 0.963] |
| pydreg, −log10 p over all candidates | 0.586 | 0.941 [0.936, 0.945] | 0.848 [0.841, 0.854] | 0.758 [0.749, 0.765] | 0.959 [0.957, 0.962] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.654 | 0.860 [0.852, 0.869] | 0.807 [0.798, 0.817] | 0.786 [0.777, 0.795] | 0.956 [0.953, 0.959] |
| Δ pydreg(SVR) − PINTS | -0.003 | 0.080 [0.074, 0.087] p<1e-03 | 0.052 [0.046, 0.057] p<1e-03 | 0.001 [-0.007, 0.009] p=0.826 | 0.005 [0.003, 0.006] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.068 | 0.080 [0.074, 0.087] p<1e-03 | 0.041 [0.035, 0.046] p<1e-03 | -0.028 [-0.037, -0.020] p<1e-03 | 0.003 [0.002, 0.005] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.812 [0.754, 0.850] | 0.723 [0.687, 0.754] | 0.645 [0.603, 0.685] | 0.945 [0.941, 0.949] |
| GC-matched: PINTS |  | 0.769 [0.757, 0.788] | 0.716 [0.694, 0.737] | 0.664 [0.629, 0.700] | 0.932 [0.928, 0.936] |
| GC-matched: Δ pydreg − PINTS |  | 0.042 [-0.016, 0.077] p=0.226 | 0.007 [-0.026, 0.037] p=0.668 | -0.019 [-0.070, 0.027] p=0.458 | 0.013 [0.009, 0.017] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.013 [-0.025, 0.061] p=0.458 | -0.006 [-0.042, 0.026] p=0.76 | -0.029 [-0.096, 0.032] p=0.382 | 0.011 [0.006, 0.015] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.152 [0.104, 0.187] p<1e-03 | 0.074 [0.049, 0.096] p<1e-03 | -0.011 [-0.041, 0.011] p=0.364 | 0.017 [0.012, 0.022] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.438 | 0.859 | 0.388 | 0.923 |
| 20,000 | 0.679 | 0.829 | 0.602 | 0.884 |
| 40,000 | 0.859 | 0.754 | 0.761 | 0.730 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 15,613, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.072**; CpG-island rows 0.056 (n_pos 11,892), non-CpG-island rows 0.018 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.728 | 0.887 [0.879, 0.895] | 0.839 [0.830, 0.848] | 0.790 [0.781, 0.799] | 0.936 [0.931, 0.940] |
| pydreg, −log10 p over all candidates | 0.721 | 0.887 [0.879, 0.895] | 0.835 [0.826, 0.844] | 0.776 [0.768, 0.785] | 0.935 [0.931, 0.939] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.345 | 0.709 [0.696, 0.723] | 0.641 [0.627, 0.657] | 0.611 [0.600, 0.624] | 0.922 [0.919, 0.926] |
| Δ pydreg(SVR) − PINTS | 0.382 | 0.178 [0.168, 0.188] p<1e-03 | 0.197 [0.187, 0.207] p<1e-03 | 0.178 [0.169, 0.189] p<1e-03 | 0.013 [0.010, 0.016] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.376 | 0.178 [0.168, 0.188] p<1e-03 | 0.193 [0.183, 0.203] p<1e-03 | 0.165 [0.155, 0.175] p<1e-03 | 0.013 [0.010, 0.016] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.846 [0.819, 0.865] | 0.747 [0.719, 0.771] | 0.662 [0.621, 0.705] | 0.924 [0.920, 0.929] |
| GC-matched: PINTS |  | 0.478 [0.401, 0.494] | 0.402 [0.379, 0.435] | 0.384 [0.357, 0.418] | 0.857 [0.848, 0.867] |
| GC-matched: Δ pydreg − PINTS |  | 0.368 [0.336, 0.446] p<1e-03 | 0.345 [0.301, 0.377] p<1e-03 | 0.278 [0.224, 0.325] p<1e-03 | 0.067 [0.057, 0.076] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.352 [0.317, 0.461] p<1e-03 | 0.342 [0.292, 0.377] p<1e-03 | 0.286 [0.219, 0.348] p<1e-03 | 0.067 [0.057, 0.076] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.470 [0.447, 0.489] p<1e-03 | 0.373 [0.344, 0.397] p<1e-03 | 0.155 [0.128, 0.179] p<1e-03 | 0.062 [0.052, 0.070] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.457 | 0.899 | 0.275 | 0.807 |
| 20,000 | 0.672 | 0.866 | 0.347 | 0.726 |
| 40,000 | 0.835 | 0.803 | 0.520 | 0.638 |

### Jurkat_CAGE (Jurkat CAGE, ood)

n_pos = 15,613, n_neg = 120,727, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.115**; CpG-island rows 0.090 (n_pos 11,892), non-CpG-island rows 0.030 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.196 | 0.462 [0.447, 0.478] | 0.410 [0.398, 0.423] | 0.420 [0.410, 0.430] | 0.734 [0.728, 0.741] |
| pydreg, −log10 p over all candidates | 0.147 | 0.439 [0.425, 0.455] | 0.385 [0.372, 0.399] | 0.370 [0.360, 0.379] | 0.732 [0.725, 0.739] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.226 | 0.294 [0.286, 0.303] | 0.270 [0.262, 0.279] | 0.307 [0.300, 0.315] | 0.580 [0.575, 0.586] |
| Δ pydreg(SVR) − PINTS | -0.029 | 0.168 [0.157, 0.180] p<1e-03 | 0.140 [0.131, 0.148] p<1e-03 | 0.113 [0.105, 0.121] p<1e-03 | 0.154 [0.150, 0.159] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.079 | 0.146 [0.135, 0.157] p<1e-03 | 0.115 [0.106, 0.124] p<1e-03 | 0.063 [0.054, 0.071] p<1e-03 | 0.152 [0.147, 0.156] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.266 [0.241, 0.300] | 0.253 [0.230, 0.277] | 0.288 [0.265, 0.315] | 0.677 [0.667, 0.687] |
| GC-matched: PINTS |  | 0.247 [0.236, 0.256] | 0.225 [0.210, 0.239] | 0.242 [0.221, 0.267] | 0.548 [0.540, 0.557] |
| GC-matched: Δ pydreg − PINTS |  | 0.018 [-0.003, 0.051] p=0.06 | 0.028 [0.012, 0.047] p<1e-03 | 0.046 [0.022, 0.071] p=0.002 | 0.129 [0.120, 0.136] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.015 [-0.023, 0.048] p=0.404 | 0.027 [0.006, 0.051] p=0.018 | 0.051 [0.018, 0.085] p=0.004 | 0.143 [0.134, 0.152] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.025 [0.014, 0.036] p<1e-03 | 0.021 [0.013, 0.029] p<1e-03 | -0.002 [-0.006, 0.002] p=0.414 | 0.087 [0.079, 0.095] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.385 | 0.616 | 0.276 | 0.494 |
| 20,000 | 0.527 | 0.439 | 0.345 | 0.153 |
| 40,000 | 0.558 | 0.352 | 0.350 | 0.150 |


## dREG auxiliary labels: dhs_h3k27ac_distal vs quiet

* **Positives:** As dhs_h3k27ac, but > 1 kb from any GENCODE v24 TSS, i.e. putative enhancers (K562 20,553; GM12878 28,215).
* **Negatives:** Genome minus (any Duke/UW DHS ∪ H3K27ac ∪ ChromHMM states 1–8), in 500 bp tiles holding ≥ 1 dREG informative position (≤ 200k), like the Danko 2015 tuning negatives.
* **Notes:** Not in either dREG paper as a ROC: the papers report recall only, so the negatives are ours. No nascent or TSS data is used, so these labels are not circular for any library, GRO-cap included. H3K27ac peaks and ChromHMM segments are broad, which favours wider calls under the overlap rule. On hg38 cells (HCT116, A549) the DHS is one ENCODE DNase peak set, not Duke ∩ UW, and there is no ChromHMM: quiet excludes only DNase and H3K27ac. Jurkat (hg19) uses UW DNase rep1 ∩ rep2 with H3K27ac from either ChIP-Atlas replicate (GSM4551183/4), and has no ChromHMM.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 20,553, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.093**; CpG-island rows 0.006 (n_pos 1,282), non-CpG-island rows 0.088 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.613 | 0.706 [0.695, 0.717] | 0.681 [0.670, 0.692] | 0.688 [0.678, 0.698] | 0.843 [0.837, 0.848] |
| pydreg, −log10 p over all candidates | 0.613 | 0.706 [0.695, 0.717] | 0.679 [0.668, 0.689] | 0.676 [0.666, 0.686] | 0.842 [0.837, 0.848] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.133 | 0.379 [0.365, 0.393] | 0.359 [0.347, 0.372] | 0.403 [0.388, 0.419] | 0.849 [0.844, 0.854] |
| Δ pydreg(SVR) − PINTS | 0.480 | 0.327 [0.312, 0.343] p<1e-03 | 0.322 [0.310, 0.335] p<1e-03 | 0.285 [0.269, 0.301] p<1e-03 | -0.006 [-0.011, -0.002] p=0.004 |
| Δ pydreg(−log p) − PINTS | 0.480 | 0.327 [0.312, 0.343] p<1e-03 | 0.320 [0.307, 0.333] p<1e-03 | 0.273 [0.256, 0.288] p<1e-03 | -0.006 [-0.011, -0.002] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.701 [0.691, 0.713] | 0.665 [0.654, 0.676] | 0.670 [0.660, 0.681] | 0.839 [0.833, 0.845] |
| GC-matched: PINTS |  | 0.304 [0.291, 0.320] | 0.290 [0.278, 0.304] | 0.339 [0.323, 0.357] | 0.818 [0.812, 0.824] |
| GC-matched: Δ pydreg − PINTS |  | 0.397 [0.381, 0.413] p<1e-03 | 0.375 [0.360, 0.389] p<1e-03 | 0.331 [0.314, 0.349] p<1e-03 | 0.021 [0.015, 0.027] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.348 [0.314, 0.385] p<1e-03 | 0.352 [0.320, 0.384] p<1e-03 | 0.470 [0.419, 0.515] p<1e-03 | 0.046 [0.034, 0.057] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.400 [0.383, 0.415] p<1e-03 | 0.376 [0.361, 0.390] p<1e-03 | 0.344 [0.327, 0.360] p<1e-03 | 0.019 [0.013, 0.025] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.115 | 0.994 | 0.034 | 0.760 |
| 20,000 | 0.225 | 0.986 | 0.069 | 0.669 |
| 40,000 | 0.389 | 0.964 | 0.126 | 0.584 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 20,553, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.093**; CpG-island rows 0.006 (n_pos 1,282), non-CpG-island rows 0.088 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.349 | 0.396 [0.387, 0.406] | 0.388 [0.379, 0.397] | 0.407 [0.398, 0.416] | 0.681 [0.676, 0.686] |
| pydreg, −log10 p over all candidates | 0.350 | 0.396 [0.387, 0.406] | 0.387 [0.378, 0.397] | 0.400 [0.391, 0.410] | 0.681 [0.676, 0.686] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.095 | 0.272 [0.262, 0.283] | 0.266 [0.256, 0.276] | 0.245 [0.234, 0.257] | 0.684 [0.679, 0.689] |
| Δ pydreg(SVR) − PINTS | 0.253 | 0.124 [0.114, 0.134] p<1e-03 | 0.122 [0.113, 0.131] p<1e-03 | 0.161 [0.153, 0.170] p<1e-03 | -0.003 [-0.007, 0.000] p=0.074 |
| Δ pydreg(−log p) − PINTS | 0.255 | 0.124 [0.114, 0.134] p<1e-03 | 0.122 [0.112, 0.131] p<1e-03 | 0.155 [0.146, 0.163] p<1e-03 | -0.003 [-0.007, 0.000] p=0.072 |
| GC-matched: pydreg (SVR) |  | 0.392 [0.382, 0.401] | 0.379 [0.370, 0.389] | 0.397 [0.388, 0.406] | 0.679 [0.674, 0.684] |
| GC-matched: PINTS |  | 0.222 [0.212, 0.234] | 0.211 [0.201, 0.222] | 0.212 [0.201, 0.225] | 0.665 [0.660, 0.671] |
| GC-matched: Δ pydreg − PINTS |  | 0.169 [0.159, 0.178] p<1e-03 | 0.169 [0.159, 0.177] p<1e-03 | 0.185 [0.176, 0.193] p<1e-03 | 0.013 [0.009, 0.017] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.284 [0.252, 0.314] p<1e-03 | 0.291 [0.264, 0.318] p<1e-03 | 0.399 [0.347, 0.446] p<1e-03 | 0.062 [0.049, 0.076] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.162 [0.152, 0.171] p<1e-03 | 0.160 [0.150, 0.168] p<1e-03 | 0.181 [0.173, 0.189] p<1e-03 | 0.010 [0.006, 0.014] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.079 | 0.992 | 0.043 | 0.600 |
| 20,000 | 0.177 | 0.969 | 0.088 | 0.494 |
| 40,000 | 0.302 | 0.868 | 0.153 | 0.476 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 20,553, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.093**; CpG-island rows 0.006 (n_pos 1,282), non-CpG-island rows 0.088 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.410 | 0.434 [0.425, 0.445] | 0.430 [0.420, 0.440] | 0.450 [0.441, 0.461] | 0.702 [0.697, 0.707] |
| pydreg, −log10 p over all candidates | 0.410 | 0.434 [0.425, 0.445] | 0.430 [0.420, 0.440] | 0.449 [0.439, 0.459] | 0.702 [0.697, 0.707] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.107 | 0.239 [0.230, 0.249] | 0.248 [0.239, 0.258] | 0.245 [0.235, 0.256] | 0.706 [0.701, 0.712] |
| Δ pydreg(SVR) − PINTS | 0.304 | 0.195 [0.184, 0.206] p<1e-03 | 0.182 [0.170, 0.194] p<1e-03 | 0.205 [0.195, 0.215] p<1e-03 | -0.005 [-0.009, -0.000] p=0.036 |
| Δ pydreg(−log p) − PINTS | 0.304 | 0.195 [0.184, 0.206] p<1e-03 | 0.182 [0.170, 0.193] p<1e-03 | 0.203 [0.193, 0.212] p<1e-03 | -0.005 [-0.009, -0.000] p=0.036 |
| GC-matched: pydreg (SVR) |  | 0.432 [0.422, 0.442] | 0.424 [0.415, 0.435] | 0.443 [0.433, 0.453] | 0.700 [0.695, 0.706] |
| GC-matched: PINTS |  | 0.186 [0.178, 0.194] | 0.182 [0.174, 0.191] | 0.205 [0.195, 0.216] | 0.682 [0.676, 0.688] |
| GC-matched: Δ pydreg − PINTS |  | 0.246 [0.235, 0.256] p<1e-03 | 0.242 [0.231, 0.253] p<1e-03 | 0.238 [0.227, 0.248] p<1e-03 | 0.018 [0.013, 0.023] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.332 [0.302, 0.362] p<1e-03 | 0.351 [0.321, 0.380] p<1e-03 | 0.503 [0.432, 0.556] p<1e-03 | 0.065 [0.053, 0.077] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.239 [0.228, 0.249] p<1e-03 | 0.234 [0.223, 0.245] p<1e-03 | 0.236 [0.226, 0.246] p<1e-03 | 0.015 [0.010, 0.019] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.062 | 0.997 | 0.027 | 0.720 |
| 20,000 | 0.150 | 0.989 | 0.057 | 0.560 |
| 40,000 | 0.301 | 0.957 | 0.130 | 0.442 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 28,215, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.124**; CpG-island rows 0.008 (n_pos 1,556), non-CpG-island rows 0.118 (n_pos 26,659).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.657 | 0.761 [0.754, 0.767] | 0.732 [0.725, 0.738] | 0.755 [0.749, 0.762] | 0.871 [0.868, 0.875] |
| pydreg, −log10 p over all candidates | 0.654 | 0.761 [0.754, 0.767] | 0.726 [0.719, 0.734] | 0.738 [0.730, 0.745] | 0.871 [0.867, 0.874] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.205 | 0.514 [0.502, 0.528] | 0.470 [0.459, 0.482] | 0.563 [0.550, 0.576] | 0.875 [0.871, 0.878] |
| Δ pydreg(SVR) − PINTS | 0.451 | 0.247 [0.234, 0.260] p<1e-03 | 0.262 [0.251, 0.273] p<1e-03 | 0.192 [0.181, 0.203] p<1e-03 | -0.003 [-0.007, -0.000] p=0.034 |
| Δ pydreg(−log p) − PINTS | 0.448 | 0.247 [0.234, 0.260] p<1e-03 | 0.257 [0.246, 0.267] p<1e-03 | 0.175 [0.164, 0.186] p<1e-03 | -0.004 [-0.007, -0.001] p=0.014 |
| GC-matched: pydreg (SVR) |  | 0.758 [0.751, 0.764] | 0.720 [0.712, 0.727] | 0.746 [0.739, 0.753] | 0.869 [0.865, 0.872] |
| GC-matched: PINTS |  | 0.417 [0.407, 0.433] | 0.409 [0.398, 0.421] | 0.516 [0.502, 0.529] | 0.859 [0.855, 0.863] |
| GC-matched: Δ pydreg − PINTS |  | 0.341 [0.325, 0.351] p<1e-03 | 0.311 [0.300, 0.321] p<1e-03 | 0.230 [0.218, 0.242] p<1e-03 | 0.010 [0.007, 0.013] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.276 [0.245, 0.305] p<1e-03 | 0.285 [0.259, 0.312] p<1e-03 | 0.388 [0.309, 0.477] p<1e-03 | 0.023 [0.015, 0.031] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.340 [0.317, 0.352] p<1e-03 | 0.310 [0.299, 0.320] p<1e-03 | 0.242 [0.230, 0.254] p<1e-03 | 0.008 [0.005, 0.012] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.092 | 0.991 | 0.038 | 0.891 |
| 20,000 | 0.211 | 0.992 | 0.090 | 0.847 |
| 40,000 | 0.407 | 0.981 | 0.174 | 0.755 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 20,553, n_neg = 126,437, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.140**; CpG-island rows 0.010 (n_pos 1,282), non-CpG-island rows 0.132 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.536 | 0.718 [0.707, 0.729] | 0.671 [0.660, 0.682] | 0.708 [0.699, 0.717] | 0.847 [0.841, 0.853] |
| pydreg, −log10 p over all candidates | 0.540 | 0.718 [0.707, 0.729] | 0.666 [0.655, 0.677] | 0.685 [0.675, 0.697] | 0.847 [0.841, 0.852] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.603 | 0.687 [0.676, 0.698] | 0.668 [0.657, 0.679] | 0.722 [0.713, 0.731] | 0.853 [0.847, 0.858] |
| Δ pydreg(SVR) − PINTS | -0.067 | 0.031 [0.026, 0.036] p<1e-03 | 0.003 [-0.002, 0.008] p=0.202 | -0.014 [-0.019, -0.009] p<1e-03 | -0.006 [-0.008, -0.003] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.063 | 0.031 [0.026, 0.036] p<1e-03 | -0.002 [-0.007, 0.003] p=0.552 | -0.037 [-0.044, -0.030] p<1e-03 | -0.006 [-0.009, -0.004] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.717 [0.706, 0.728] | 0.667 [0.656, 0.679] | 0.705 [0.695, 0.715] | 0.846 [0.840, 0.852] |
| GC-matched: PINTS |  | 0.677 [0.666, 0.688] | 0.658 [0.647, 0.669] | 0.715 [0.705, 0.725] | 0.849 [0.843, 0.854] |
| GC-matched: Δ pydreg − PINTS |  | 0.040 [0.034, 0.045] p<1e-03 | 0.009 [0.003, 0.015] p=0.002 | -0.010 [-0.017, -0.004] p<1e-03 | -0.003 [-0.006, -0.000] p=0.038 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.055 [0.040, 0.071] p<1e-03 | 0.026 [0.011, 0.040] p<1e-03 | -0.013 [-0.094, 0.063] p=0.704 | 0.017 [0.010, 0.025] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.039 [0.032, 0.044] p<1e-03 | 0.008 [0.002, 0.013] p=0.014 | -0.011 [-0.018, -0.005] p<1e-03 | -0.004 [-0.007, -0.001] p=0.002 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.140 | 0.991 | 0.062 | 0.976 |
| 20,000 | 0.268 | 0.975 | 0.169 | 0.983 |
| 40,000 | 0.469 | 0.924 | 0.369 | 0.979 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 16,069, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.074**; CpG-island rows 0.006 (n_pos 1,175), non-CpG-island rows 0.069 (n_pos 14,894).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.314 | 0.385 [0.371, 0.398] | 0.364 [0.352, 0.376] | 0.306 [0.291, 0.320] | 0.674 [0.667, 0.681] |
| pydreg, −log10 p over all candidates | 0.315 | 0.385 [0.371, 0.398] | 0.365 [0.352, 0.377] | 0.307 [0.292, 0.321] | 0.674 [0.667, 0.681] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.066 | 0.206 [0.198, 0.218] | 0.200 [0.191, 0.209] | 0.161 [0.152, 0.171] | 0.651 [0.644, 0.658] |
| Δ pydreg(SVR) − PINTS | 0.248 | 0.178 [0.165, 0.189] p<1e-03 | 0.165 [0.155, 0.175] p<1e-03 | 0.145 [0.135, 0.154] p<1e-03 | 0.023 [0.018, 0.028] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.249 | 0.178 [0.165, 0.189] p<1e-03 | 0.165 [0.155, 0.175] p<1e-03 | 0.145 [0.136, 0.155] p<1e-03 | 0.023 [0.018, 0.028] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.379 [0.366, 0.392] | 0.344 [0.331, 0.356] | 0.258 [0.242, 0.273] | 0.670 [0.663, 0.676] |
| GC-matched: PINTS |  | 0.163 [0.152, 0.171] | 0.154 [0.146, 0.162] | 0.133 [0.126, 0.141] | 0.633 [0.626, 0.640] |
| GC-matched: Δ pydreg − PINTS |  | 0.216 [0.205, 0.233] p<1e-03 | 0.190 [0.180, 0.201] p<1e-03 | 0.125 [0.113, 0.136] p<1e-03 | 0.037 [0.032, 0.042] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.371 [0.284, 0.416] p<1e-03 | 0.283 [0.244, 0.318] p<1e-03 | 0.050 [0.029, 0.071] p<1e-03 | 0.075 [0.062, 0.091] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.202 [0.191, 0.213] p<1e-03 | 0.182 [0.172, 0.192] p<1e-03 | 0.120 [0.110, 0.131] p<1e-03 | 0.033 [0.028, 0.038] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.061 | 0.675 | 0.023 | 0.434 |
| 20,000 | 0.162 | 0.729 | 0.049 | 0.380 |
| 40,000 | 0.333 | 0.696 | 0.102 | 0.315 |

### HCT116_PROcap (HCT116 PRO-cap, ood)

n_pos = 16,069, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.074**; CpG-island rows 0.006 (n_pos 1,175), non-CpG-island rows 0.069 (n_pos 14,894).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.119 | 0.551 [0.533, 0.568] | 0.486 [0.473, 0.498] | 0.392 [0.379, 0.406] | 0.861 [0.857, 0.865] |
| pydreg, −log10 p over all candidates | 0.060 | 0.540 [0.522, 0.558] | 0.468 [0.455, 0.480] | 0.362 [0.350, 0.374] | 0.859 [0.855, 0.863] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.167 | 0.633 [0.622, 0.647] | 0.550 [0.540, 0.561] | 0.441 [0.425, 0.457] | 0.861 [0.857, 0.866] |
| Δ pydreg(SVR) − PINTS | -0.048 | -0.082 [-0.096, -0.071] p<1e-03 | -0.064 [-0.073, -0.056] p<1e-03 | -0.049 [-0.060, -0.039] p<1e-03 | -0.001 [-0.003, 0.002] p=0.682 |
| Δ pydreg(−log p) − PINTS | -0.106 | -0.093 [-0.108, -0.082] p<1e-03 | -0.082 [-0.091, -0.074] p<1e-03 | -0.079 [-0.089, -0.071] p<1e-03 | -0.002 [-0.005, 0.000] p=0.096 |
| GC-matched: pydreg (SVR) |  | 0.505 [0.486, 0.524] | 0.450 [0.434, 0.465] | 0.365 [0.352, 0.379] | 0.856 [0.851, 0.860] |
| GC-matched: PINTS |  | 0.609 [0.595, 0.622] | 0.513 [0.499, 0.526] | 0.404 [0.387, 0.421] | 0.854 [0.850, 0.859] |
| GC-matched: Δ pydreg − PINTS |  | -0.103 [-0.118, -0.091] p<1e-03 | -0.063 [-0.072, -0.054] p<1e-03 | -0.039 [-0.050, -0.030] p<1e-03 | 0.001 [-0.002, 0.004] p=0.44 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.119 [-0.181, -0.069] p<1e-03 | -0.110 [-0.147, -0.071] p<1e-03 | -0.011 [-0.018, -0.003] p=0.004 | -0.017 [-0.025, -0.010] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.078 [-0.092, -0.067] p<1e-03 | -0.062 [-0.070, -0.054] p<1e-03 | -0.047 [-0.057, -0.038] p<1e-03 | 0.001 [-0.002, 0.004] p=0.482 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.139 | 0.492 | 0.045 | 0.480 |
| 20,000 | 0.285 | 0.501 | 0.155 | 0.568 |
| 40,000 | 0.511 | 0.478 | 0.395 | 0.597 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 49,656, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.199**; CpG-island rows 0.017 (n_pos 3,359), non-CpG-island rows 0.188 (n_pos 46,297).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.212 | 0.257 [0.251, 0.263] | 0.251 [0.245, 0.257] | 0.368 [0.360, 0.375] | 0.608 [0.605, 0.612] |
| pydreg, −log10 p over all candidates | 0.214 | 0.257 [0.251, 0.263] | 0.251 [0.245, 0.257] | 0.366 [0.358, 0.374] | 0.608 [0.605, 0.612] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.050 | 0.115 [0.110, 0.121] | 0.114 [0.109, 0.120] | 0.251 [0.240, 0.261] | 0.554 [0.548, 0.559] |
| Δ pydreg(SVR) − PINTS | 0.162 | 0.142 [0.137, 0.147] p<1e-03 | 0.137 [0.131, 0.142] p<1e-03 | 0.117 [0.112, 0.122] p<1e-03 | 0.055 [0.050, 0.059] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.164 | 0.142 [0.137, 0.147] p<1e-03 | 0.137 [0.131, 0.142] p<1e-03 | 0.115 [0.110, 0.120] p<1e-03 | 0.055 [0.050, 0.059] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.256 [0.250, 0.262] | 0.249 [0.243, 0.256] | 0.367 [0.359, 0.374] | 0.608 [0.604, 0.611] |
| GC-matched: PINTS |  | 0.106 [0.101, 0.111] | 0.103 [0.098, 0.109] | 0.243 [0.233, 0.253] | 0.546 [0.540, 0.551] |
| GC-matched: Δ pydreg − PINTS |  | 0.150 [0.145, 0.155] p<1e-03 | 0.146 [0.141, 0.151] p<1e-03 | 0.124 [0.118, 0.129] p<1e-03 | 0.062 [0.058, 0.066] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.265 [0.243, 0.290] p<1e-03 | 0.281 [0.260, 0.301] p<1e-03 | 0.356 [0.325, 0.384] p<1e-03 | 0.100 [0.088, 0.111] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.142 [0.137, 0.146] p<1e-03 | 0.137 [0.132, 0.142] p<1e-03 | 0.118 [0.113, 0.123] p<1e-03 | 0.059 [0.054, 0.063] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.028 | 0.973 | 0.020 | 0.642 |
| 20,000 | 0.085 | 0.965 | 0.043 | 0.559 |
| 40,000 | 0.199 | 0.869 | 0.142 | 0.336 |

### A549_PROcap (A549 PRO-cap, ood)

n_pos = 49,656, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.199**; CpG-island rows 0.017 (n_pos 3,359), non-CpG-island rows 0.188 (n_pos 46,297).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.025 | 0.087 [0.082, 0.093] | 0.081 [0.076, 0.086] | 0.217 [0.209, 0.225] | 0.508 [0.504, 0.512] |
| pydreg, −log10 p over all candidates | 0.018 | 0.086 [0.081, 0.092] | 0.079 [0.075, 0.084] | 0.213 [0.206, 0.221] | 0.508 [0.504, 0.512] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.032 | 0.090 [0.085, 0.096] | 0.084 [0.079, 0.088] | 0.217 [0.210, 0.226] | 0.471 [0.467, 0.475] |
| Δ pydreg(SVR) − PINTS | -0.007 | -0.003 [-0.006, -0.001] p=0.018 | -0.003 [-0.004, -0.001] p=0.004 | -0.001 [-0.003, 0.001] p=0.556 | 0.037 [0.035, 0.039] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.014 | -0.004 [-0.007, -0.001] p=0.004 | -0.005 [-0.007, -0.003] p<1e-03 | -0.005 [-0.007, -0.003] p<1e-03 | 0.037 [0.035, 0.039] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.080 [0.074, 0.085] | 0.075 [0.070, 0.079] | 0.212 [0.204, 0.220] | 0.504 [0.499, 0.508] |
| GC-matched: PINTS |  | 0.082 [0.076, 0.087] | 0.077 [0.072, 0.082] | 0.211 [0.204, 0.220] | 0.460 [0.456, 0.464] |
| GC-matched: Δ pydreg − PINTS |  | -0.002 [-0.005, 0.000] p=0.066 | -0.003 [-0.004, -0.001] p=0.006 | 0.000 [-0.001, 0.002] p=0.598 | 0.044 [0.042, 0.046] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.038 [-0.053, -0.020] p<1e-03 | -0.016 [-0.028, -0.006] p=0.004 | -0.005 [-0.011, -0.000] p=0.042 | 0.047 [0.040, 0.053] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.002 [-0.005, 0.000] p=0.068 | -0.002 [-0.004, -0.001] p=0.01 | 0.002 [0.000, 0.003] p=0.01 | 0.043 [0.041, 0.045] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.020 | 0.400 | 0.013 | 0.481 |
| 20,000 | 0.047 | 0.349 | 0.034 | 0.437 |
| 40,000 | 0.099 | 0.292 | 0.076 | 0.326 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 2,600, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.013**; CpG-island rows 0.006 (n_pos 1,212), non-CpG-island rows 0.007 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.327 | 0.603 [0.582, 0.622] | 0.545 [0.526, 0.562] | 0.176 [0.162, 0.193] | 0.785 [0.775, 0.795] |
| pydreg, −log10 p over all candidates | 0.317 | 0.603 [0.582, 0.622] | 0.543 [0.525, 0.561] | 0.170 [0.156, 0.184] | 0.785 [0.774, 0.795] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.183 | 0.388 [0.366, 0.409] | 0.351 [0.332, 0.369] | 0.095 [0.085, 0.108] | 0.774 [0.762, 0.784] |
| Δ pydreg(SVR) − PINTS | 0.143 | 0.214 [0.194, 0.237] p<1e-03 | 0.194 [0.176, 0.212] p<1e-03 | 0.081 [0.068, 0.095] p<1e-03 | 0.011 [0.002, 0.021] p=0.022 |
| Δ pydreg(−log p) − PINTS | 0.133 | 0.214 [0.194, 0.237] p<1e-03 | 0.192 [0.175, 0.211] p<1e-03 | 0.074 [0.062, 0.087] p<1e-03 | 0.011 [0.001, 0.021] p=0.022 |
| GC-matched: pydreg (SVR) |  | 0.593 [0.572, 0.612] | 0.510 [0.490, 0.529] | 0.149 [0.130, 0.173] | 0.778 [0.767, 0.788] |
| GC-matched: PINTS |  | 0.255 [0.230, 0.278] | 0.253 [0.235, 0.269] | 0.053 [0.045, 0.063] | 0.703 [0.690, 0.715] |
| GC-matched: Δ pydreg − PINTS |  | 0.339 [0.310, 0.366] p<1e-03 | 0.257 [0.237, 0.277] p<1e-03 | 0.097 [0.079, 0.116] p<1e-03 | 0.075 [0.064, 0.087] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.383 [0.343, 0.417] p<1e-03 | 0.295 [0.267, 0.323] p<1e-03 | 0.096 [0.071, 0.125] p<1e-03 | 0.086 [0.072, 0.102] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.301 [0.272, 0.331] p<1e-03 | 0.227 [0.207, 0.251] p<1e-03 | 0.035 [0.028, 0.045] p<1e-03 | 0.063 [0.050, 0.077] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.148 | 0.293 | 0.092 | 0.277 |
| 20,000 | 0.300 | 0.299 | 0.156 | 0.187 |
| 40,000 | 0.507 | 0.269 | 0.254 | 0.127 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 2,600, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.013**; CpG-island rows 0.006 (n_pos 1,212), non-CpG-island rows 0.007 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.359 | 0.641 [0.620, 0.661] | 0.579 [0.561, 0.597] | 0.201 [0.184, 0.218] | 0.805 [0.794, 0.815] |
| pydreg, −log10 p over all candidates | 0.352 | 0.641 [0.620, 0.661] | 0.578 [0.559, 0.595] | 0.190 [0.175, 0.205] | 0.805 [0.794, 0.815] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.229 | 0.454 [0.433, 0.475] | 0.415 [0.397, 0.434] | 0.124 [0.111, 0.139] | 0.796 [0.785, 0.806] |
| Δ pydreg(SVR) − PINTS | 0.130 | 0.187 [0.166, 0.209] p<1e-03 | 0.164 [0.147, 0.182] p<1e-03 | 0.076 [0.062, 0.090] p<1e-03 | 0.009 [-0.000, 0.018] p=0.054 |
| Δ pydreg(−log p) − PINTS | 0.123 | 0.187 [0.166, 0.209] p<1e-03 | 0.162 [0.145, 0.181] p<1e-03 | 0.065 [0.052, 0.079] p<1e-03 | 0.009 [-0.000, 0.018] p=0.056 |
| GC-matched: pydreg (SVR) |  | 0.623 [0.603, 0.644] | 0.538 [0.518, 0.556] | 0.169 [0.147, 0.193] | 0.797 [0.786, 0.807] |
| GC-matched: PINTS |  | 0.355 [0.325, 0.380] | 0.319 [0.299, 0.340] | 0.076 [0.064, 0.091] | 0.732 [0.719, 0.745] |
| GC-matched: Δ pydreg − PINTS |  | 0.268 [0.241, 0.301] p<1e-03 | 0.219 [0.197, 0.241] p<1e-03 | 0.093 [0.071, 0.115] p<1e-03 | 0.065 [0.053, 0.075] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.302 [0.256, 0.335] p<1e-03 | 0.247 [0.216, 0.276] p<1e-03 | 0.107 [0.076, 0.139] p<1e-03 | 0.076 [0.062, 0.091] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.253 [0.220, 0.284] p<1e-03 | 0.196 [0.172, 0.220] p<1e-03 | 0.031 [0.024, 0.040] p<1e-03 | 0.052 [0.039, 0.065] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.165 | 0.321 | 0.087 | 0.295 |
| 20,000 | 0.326 | 0.321 | 0.177 | 0.263 |
| 40,000 | 0.544 | 0.285 | 0.320 | 0.174 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 2,600, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.013**; CpG-island rows 0.006 (n_pos 1,212), non-CpG-island rows 0.007 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.370 | 0.784 [0.768, 0.801] | 0.681 [0.666, 0.697] | 0.246 [0.228, 0.268] | 0.876 [0.868, 0.885] |
| pydreg, −log10 p over all candidates | 0.332 | 0.784 [0.768, 0.801] | 0.674 [0.659, 0.689] | 0.216 [0.200, 0.233] | 0.875 [0.867, 0.884] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.321 | 0.606 [0.586, 0.624] | 0.546 [0.528, 0.562] | 0.193 [0.177, 0.211] | 0.863 [0.854, 0.871] |
| Δ pydreg(SVR) − PINTS | 0.049 | 0.178 [0.161, 0.196] p<1e-03 | 0.135 [0.122, 0.149] p<1e-03 | 0.053 [0.036, 0.070] p<1e-03 | 0.013 [0.006, 0.020] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.011 | 0.178 [0.161, 0.196] p<1e-03 | 0.128 [0.114, 0.143] p<1e-03 | 0.023 [0.008, 0.036] p=0.004 | 0.012 [0.005, 0.020] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.711 [0.675, 0.742] | 0.601 [0.577, 0.621] | 0.171 [0.150, 0.198] | 0.865 [0.855, 0.874] |
| GC-matched: PINTS |  | 0.505 [0.476, 0.539] | 0.470 [0.450, 0.490] | 0.137 [0.120, 0.158] | 0.821 [0.810, 0.831] |
| GC-matched: Δ pydreg − PINTS |  | 0.206 [0.162, 0.246] p<1e-03 | 0.131 [0.109, 0.154] p<1e-03 | 0.035 [0.011, 0.059] p=0.004 | 0.044 [0.035, 0.053] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.200 [0.151, 0.248] p<1e-03 | 0.125 [0.097, 0.154] p<1e-03 | 0.023 [-0.007, 0.057] p=0.154 | 0.054 [0.042, 0.066] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.216 [0.165, 0.260] p<1e-03 | 0.137 [0.113, 0.161] p<1e-03 | 0.018 [0.008, 0.027] p=0.004 | 0.033 [0.022, 0.044] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.215 | 0.332 | 0.121 | 0.382 |
| 20,000 | 0.398 | 0.322 | 0.273 | 0.365 |
| 40,000 | 0.612 | 0.266 | 0.457 | 0.213 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 2,600, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.013**; CpG-island rows 0.006 (n_pos 1,212), non-CpG-island rows 0.007 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.392 | 0.650 [0.631, 0.668] | 0.594 [0.577, 0.612] | 0.223 [0.205, 0.242] | 0.810 [0.801, 0.820] |
| pydreg, −log10 p over all candidates | 0.389 | 0.650 [0.631, 0.668] | 0.593 [0.576, 0.611] | 0.216 [0.200, 0.232] | 0.810 [0.800, 0.820] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.113 | 0.389 [0.368, 0.410] | 0.347 [0.328, 0.366] | 0.078 [0.069, 0.089] | 0.803 [0.794, 0.812] |
| Δ pydreg(SVR) − PINTS | 0.278 | 0.261 [0.239, 0.283] p<1e-03 | 0.247 [0.228, 0.267] p<1e-03 | 0.145 [0.129, 0.160] p<1e-03 | 0.007 [-0.002, 0.017] p=0.136 |
| Δ pydreg(−log p) − PINTS | 0.276 | 0.261 [0.239, 0.283] p<1e-03 | 0.246 [0.227, 0.266] p<1e-03 | 0.138 [0.123, 0.152] p<1e-03 | 0.007 [-0.002, 0.017] p=0.142 |
| GC-matched: pydreg (SVR) |  | 0.637 [0.618, 0.657] | 0.548 [0.528, 0.566] | 0.160 [0.138, 0.183] | 0.801 [0.791, 0.811] |
| GC-matched: PINTS |  | 0.174 [0.158, 0.191] | 0.162 [0.145, 0.181] | 0.035 [0.030, 0.043] | 0.717 [0.704, 0.729] |
| GC-matched: Δ pydreg − PINTS |  | 0.463 [0.441, 0.485] p<1e-03 | 0.385 [0.360, 0.407] p<1e-03 | 0.126 [0.104, 0.147] p<1e-03 | 0.084 [0.072, 0.098] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.505 [0.472, 0.533] p<1e-03 | 0.448 [0.415, 0.477] p<1e-03 | 0.113 [0.086, 0.143] p<1e-03 | 0.111 [0.096, 0.128] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.427 [0.398, 0.455] p<1e-03 | 0.333 [0.308, 0.358] p<1e-03 | 0.044 [0.034, 0.053] p<1e-03 | 0.059 [0.045, 0.076] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.170 | 0.357 | 0.078 | 0.165 |
| 20,000 | 0.333 | 0.348 | 0.114 | 0.127 |
| 40,000 | 0.546 | 0.307 | 0.208 | 0.105 |

### Jurkat_CAGE (Jurkat CAGE, ood)

n_pos = 2,600, n_neg = 120,727, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.021**; CpG-island rows 0.010 (n_pos 1,212), non-CpG-island rows 0.011 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.007 | 0.038 [0.030, 0.045] | 0.036 [0.029, 0.042] | 0.021 [0.019, 0.022] | 0.474 [0.468, 0.480] |
| pydreg, −log10 p over all candidates | 0.008 | 0.035 [0.028, 0.042] | 0.034 [0.028, 0.040] | 0.021 [0.019, 0.022] | 0.473 [0.467, 0.479] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.004 | 0.028 [0.021, 0.034] | 0.024 [0.019, 0.029] | 0.020 [0.018, 0.022] | 0.413 [0.407, 0.419] |
| Δ pydreg(SVR) − PINTS | 0.003 | 0.010 [0.003, 0.017] p=0.004 | 0.012 [0.007, 0.017] p<1e-03 | 0.001 [0.000, 0.001] p<1e-03 | 0.061 [0.055, 0.066] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.004 | 0.008 [0.002, 0.014] p=0.014 | 0.010 [0.006, 0.015] p<1e-03 | 0.001 [0.000, 0.001] p<1e-03 | 0.060 [0.055, 0.066] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.019 [0.014, 0.025] | 0.020 [0.015, 0.024] | 0.020 [0.018, 0.022] | 0.426 [0.419, 0.434] |
| GC-matched: PINTS |  | 0.014 [0.009, 0.020] | 0.015 [0.011, 0.019] | 0.020 [0.018, 0.022] | 0.374 [0.366, 0.382] |
| GC-matched: Δ pydreg − PINTS |  | 0.004 [-0.002, 0.011] p=0.142 | 0.004 [0.000, 0.008] p=0.044 | 0.000 [0.000, 0.000] p=0.004 | 0.052 [0.045, 0.060] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.008 [-0.003, 0.017] p=0.204 | 0.004 [-0.004, 0.011] p=0.23 | 0.000 [0.000, 0.000] p<1e-03 | 0.048 [0.038, 0.059] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.004 [-0.003, 0.012] p=0.24 | 0.005 [0.000, 0.009] p=0.032 | 0.000 [-0.000, 0.000] p=0.548 | 0.056 [0.049, 0.063] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.026 | 0.018 | 0.018 | 0.011 |
| 20,000 | 0.059 | 0.014 | 0.077 | 0.007 |
| 40,000 | 0.078 | 0.013 | 0.079 | 0.007 |


## dREG auxiliary labels: chromhmm_active vs quiet

* **Positives:** ENCODE Broad ChromHMM segments in states 1 (active promoter) and 4–5 (strong enhancer), as in Danko 2015 Fig 1b (K562 43,377; GM12878 41,415; median segment about 1.4 kb).
* **Negatives:** Genome minus (any Duke/UW DHS ∪ H3K27ac ∪ ChromHMM states 1–8), in 500 bp tiles holding ≥ 1 dREG informative position (≤ 200k), like the Danko 2015 tuning negatives.
* **Notes:** Not in either dREG paper as a ROC: the papers report recall only, so the negatives are ours. No nascent or TSS data is used, so these labels are not circular for any library, GRO-cap included. H3K27ac peaks and ChromHMM segments are broad, which favours wider calls under the overlap rule. On hg38 cells (HCT116, A549) the DHS is one ENCODE DNase peak set, not Duke ∩ UW, and there is no ChromHMM: quiet excludes only DNase and H3K27ac. Jurkat (hg19) uses UW DNase rep1 ∩ rep2 with H3K27ac from either ChIP-Atlas replicate (GSM4551183/4), and has no ChromHMM.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 43,377, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.178**; CpG-island rows 0.033 (n_pos 6,818), non-CpG-island rows 0.155 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.649 | 0.729 [0.722, 0.736] | 0.709 [0.702, 0.716] | 0.758 [0.751, 0.765] | 0.855 [0.852, 0.859] |
| pydreg, −log10 p over all candidates | 0.650 | 0.729 [0.722, 0.736] | 0.707 [0.700, 0.714] | 0.753 [0.746, 0.760] | 0.855 [0.851, 0.859] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.325 | 0.551 [0.541, 0.561] | 0.522 [0.513, 0.531] | 0.702 [0.694, 0.711] | 0.880 [0.877, 0.883] |
| Δ pydreg(SVR) − PINTS | 0.324 | 0.178 [0.170, 0.188] p<1e-03 | 0.187 [0.180, 0.195] p<1e-03 | 0.056 [0.050, 0.061] p<1e-03 | -0.025 [-0.027, -0.022] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.325 | 0.178 [0.170, 0.188] p<1e-03 | 0.185 [0.177, 0.193] p<1e-03 | 0.051 [0.045, 0.057] p<1e-03 | -0.025 [-0.027, -0.022] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.726 [0.719, 0.733] | 0.699 [0.691, 0.705] | 0.752 [0.746, 0.759] | 0.853 [0.849, 0.857] |
| GC-matched: PINTS |  | 0.501 [0.492, 0.511] | 0.478 [0.468, 0.487] | 0.667 [0.658, 0.675] | 0.863 [0.859, 0.866] |
| GC-matched: Δ pydreg − PINTS |  | 0.225 [0.215, 0.234] p<1e-03 | 0.221 [0.212, 0.231] p<1e-03 | 0.086 [0.079, 0.094] p<1e-03 | -0.010 [-0.013, -0.005] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.121 [0.111, 0.132] p<1e-03 | 0.154 [0.141, 0.168] p<1e-03 | 0.239 [0.211, 0.270] p<1e-03 | 0.024 [0.021, 0.027] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.238 [0.228, 0.250] p<1e-03 | 0.229 [0.220, 0.240] p<1e-03 | 0.110 [0.102, 0.120] p<1e-03 | -0.019 [-0.022, -0.015] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.221 | 0.999 | 0.150 | 0.967 |
| 20,000 | 0.362 | 0.996 | 0.229 | 0.934 |
| 40,000 | 0.491 | 0.986 | 0.318 | 0.882 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 43,377, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.178**; CpG-island rows 0.033 (n_pos 6,818), non-CpG-island rows 0.155 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.458 | 0.498 [0.489, 0.508] | 0.492 [0.483, 0.501] | 0.565 [0.556, 0.574] | 0.735 [0.731, 0.740] |
| pydreg, −log10 p over all candidates | 0.460 | 0.498 [0.489, 0.508] | 0.491 [0.483, 0.500] | 0.562 [0.554, 0.571] | 0.735 [0.731, 0.740] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.218 | 0.409 [0.401, 0.419] | 0.392 [0.384, 0.401] | 0.503 [0.493, 0.513] | 0.741 [0.736, 0.746] |
| Δ pydreg(SVR) − PINTS | 0.240 | 0.089 [0.081, 0.098] p<1e-03 | 0.100 [0.092, 0.108] p<1e-03 | 0.062 [0.056, 0.068] p<1e-03 | -0.006 [-0.008, -0.003] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.242 | 0.089 [0.081, 0.098] p<1e-03 | 0.099 [0.091, 0.108] p<1e-03 | 0.059 [0.053, 0.065] p<1e-03 | -0.006 [-0.008, -0.003] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.496 [0.487, 0.505] | 0.486 [0.478, 0.495] | 0.562 [0.554, 0.571] | 0.734 [0.729, 0.738] |
| GC-matched: PINTS |  | 0.372 [0.364, 0.381] | 0.352 [0.344, 0.361] | 0.477 [0.468, 0.487] | 0.729 [0.724, 0.733] |
| GC-matched: Δ pydreg − PINTS |  | 0.124 [0.115, 0.132] p<1e-03 | 0.134 [0.126, 0.142] p<1e-03 | 0.085 [0.078, 0.092] p<1e-03 | 0.005 [0.002, 0.008] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.222 [0.206, 0.237] p<1e-03 | 0.266 [0.250, 0.282] p<1e-03 | 0.410 [0.387, 0.432] p<1e-03 | 0.049 [0.044, 0.054] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.104 [0.096, 0.112] p<1e-03 | 0.108 [0.100, 0.116] p<1e-03 | 0.082 [0.075, 0.089] p<1e-03 | -0.005 [-0.008, -0.002] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.216 | 0.999 | 0.125 | 0.902 |
| 20,000 | 0.325 | 0.992 | 0.209 | 0.831 |
| 40,000 | 0.422 | 0.951 | 0.295 | 0.787 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 43,377, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.178**; CpG-island rows 0.033 (n_pos 6,818), non-CpG-island rows 0.155 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.512 | 0.532 [0.521, 0.543] | 0.529 [0.518, 0.540] | 0.595 [0.585, 0.606] | 0.753 [0.748, 0.759] |
| pydreg, −log10 p over all candidates | 0.512 | 0.532 [0.521, 0.543] | 0.529 [0.518, 0.540] | 0.594 [0.584, 0.605] | 0.753 [0.748, 0.759] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.299 | 0.424 [0.414, 0.434] | 0.417 [0.408, 0.426] | 0.548 [0.538, 0.559] | 0.769 [0.763, 0.774] |
| Δ pydreg(SVR) − PINTS | 0.214 | 0.108 [0.100, 0.116] p<1e-03 | 0.112 [0.104, 0.120] p<1e-03 | 0.046 [0.041, 0.052] p<1e-03 | -0.015 [-0.018, -0.012] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.214 | 0.108 [0.100, 0.116] p<1e-03 | 0.112 [0.104, 0.120] p<1e-03 | 0.046 [0.041, 0.051] p<1e-03 | -0.015 [-0.018, -0.012] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.530 [0.520, 0.541] | 0.526 [0.515, 0.536] | 0.593 [0.583, 0.604] | 0.752 [0.747, 0.758] |
| GC-matched: PINTS |  | 0.385 [0.377, 0.394] | 0.376 [0.367, 0.384] | 0.521 [0.512, 0.531] | 0.754 [0.748, 0.760] |
| GC-matched: Δ pydreg − PINTS |  | 0.145 [0.137, 0.153] p<1e-03 | 0.150 [0.141, 0.158] p<1e-03 | 0.072 [0.066, 0.079] p<1e-03 | -0.002 [-0.005, 0.001] p=0.276 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.143 [0.133, 0.152] p<1e-03 | 0.176 [0.163, 0.187] p<1e-03 | 0.253 [0.231, 0.276] p<1e-03 | 0.031 [0.028, 0.035] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.137 [0.129, 0.146] p<1e-03 | 0.139 [0.130, 0.147] p<1e-03 | 0.086 [0.079, 0.093] p<1e-03 | -0.011 [-0.015, -0.008] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.216 | 1.000 | 0.165 | 0.971 |
| 20,000 | 0.325 | 0.998 | 0.221 | 0.912 |
| 40,000 | 0.435 | 0.985 | 0.326 | 0.807 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 41,415, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.172**; CpG-island rows 0.034 (n_pos 7,034), non-CpG-island rows 0.147 (n_pos 34,381).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.661 | 0.747 [0.739, 0.755] | 0.724 [0.717, 0.732] | 0.771 [0.764, 0.778] | 0.865 [0.860, 0.868] |
| pydreg, −log10 p over all candidates | 0.659 | 0.747 [0.739, 0.755] | 0.720 [0.713, 0.727] | 0.763 [0.756, 0.769] | 0.864 [0.860, 0.868] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.391 | 0.640 [0.630, 0.649] | 0.595 [0.586, 0.603] | 0.734 [0.726, 0.741] | 0.885 [0.881, 0.888] |
| Δ pydreg(SVR) − PINTS | 0.269 | 0.107 [0.099, 0.116] p<1e-03 | 0.130 [0.123, 0.137] p<1e-03 | 0.038 [0.033, 0.042] p<1e-03 | -0.021 [-0.023, -0.018] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.268 | 0.107 [0.099, 0.116] p<1e-03 | 0.126 [0.119, 0.133] p<1e-03 | 0.029 [0.024, 0.034] p<1e-03 | -0.021 [-0.023, -0.019] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.744 [0.736, 0.751] | 0.716 [0.708, 0.723] | 0.767 [0.760, 0.773] | 0.863 [0.859, 0.866] |
| GC-matched: PINTS |  | 0.572 [0.562, 0.583] | 0.555 [0.546, 0.564] | 0.709 [0.700, 0.716] | 0.874 [0.871, 0.878] |
| GC-matched: Δ pydreg − PINTS |  | 0.173 [0.160, 0.181] p<1e-03 | 0.161 [0.153, 0.168] p<1e-03 | 0.058 [0.053, 0.064] p<1e-03 | -0.012 [-0.015, -0.009] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.084 [0.076, 0.093] p<1e-03 | 0.110 [0.101, 0.118] p<1e-03 | 0.152 [0.133, 0.170] p<1e-03 | 0.016 [0.014, 0.018] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.154 [0.138, 0.177] p<1e-03 | 0.160 [0.152, 0.168] p<1e-03 | 0.077 [0.070, 0.084] p<1e-03 | -0.021 [-0.024, -0.018] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.224 | 0.998 | 0.171 | 0.982 |
| 20,000 | 0.368 | 0.997 | 0.262 | 0.959 |
| 40,000 | 0.508 | 0.989 | 0.360 | 0.903 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 43,377, n_neg = 126,437, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.255**; CpG-island rows 0.051 (n_pos 6,818), non-CpG-island rows 0.224 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.551 | 0.674 [0.665, 0.681] | 0.643 [0.634, 0.651] | 0.738 [0.732, 0.745] | 0.825 [0.821, 0.829] |
| pydreg, −log10 p over all candidates | 0.554 | 0.674 [0.665, 0.681] | 0.637 [0.628, 0.645] | 0.728 [0.721, 0.735] | 0.825 [0.820, 0.829] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.596 | 0.668 [0.660, 0.676] | 0.654 [0.646, 0.662] | 0.768 [0.762, 0.775] | 0.842 [0.837, 0.846] |
| Δ pydreg(SVR) − PINTS | -0.045 | 0.005 [0.002, 0.009] p=0.002 | -0.011 [-0.014, -0.008] p<1e-03 | -0.030 [-0.032, -0.027] p<1e-03 | -0.016 [-0.018, -0.015] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.042 | 0.005 [0.002, 0.009] p=0.002 | -0.017 [-0.020, -0.013] p<1e-03 | -0.040 [-0.044, -0.037] p<1e-03 | -0.017 [-0.019, -0.015] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.673 [0.664, 0.680] | 0.640 [0.631, 0.649] | 0.738 [0.731, 0.744] | 0.824 [0.820, 0.829] |
| GC-matched: PINTS |  | 0.661 [0.653, 0.669] | 0.647 [0.639, 0.655] | 0.765 [0.758, 0.772] | 0.839 [0.834, 0.843] |
| GC-matched: Δ pydreg − PINTS |  | 0.011 [0.008, 0.015] p<1e-03 | -0.007 [-0.010, -0.003] p<1e-03 | -0.027 [-0.030, -0.025] p<1e-03 | -0.014 [-0.016, -0.012] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.010 [0.006, 0.015] p<1e-03 | 0.001 [-0.004, 0.007] p=0.558 | -0.003 [-0.016, 0.014] p=0.65 | -0.000 [-0.003, 0.002] p=0.78 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.010 [0.006, 0.014] p<1e-03 | -0.009 [-0.013, -0.006] p<1e-03 | -0.032 [-0.035, -0.028] p<1e-03 | -0.017 [-0.019, -0.015] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.217 | 0.997 | 0.191 | 0.996 |
| 20,000 | 0.369 | 0.991 | 0.308 | 0.996 |
| 40,000 | 0.510 | 0.966 | 0.445 | 0.992 |


## PINTS-on-TSS labels: pints_tss vs tss_silent

* **Positives:** PINTS default bidirectional calls on the matched-cell TSS library: ENCODE PRO-cap (hg38; K562 lifted to hg19 for the K562 run-on libraries), or FANTOM5 CAGE for Jurkat (hg19).
* **Negatives:** Genome minus (ENCODE DNase ∪ PRO-cap signal ± 250 bp ∪ positives), in informative 500 bp tiles (≤ 200k). For K562 hg19 the negatives are the dREG hub negatives minus the lifted positives ± 250 bp.
* **Notes:** Applied only to libraries other than the PRO-cap library the labels came from. The positives are PINTS calls, so they favour PINTS's element boundaries and its sensitivity profile.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 33,958, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.145**; CpG-island rows 0.091 (n_pos 20,129), non-CpG-island rows 0.065 (n_pos 13,829).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.890 | 0.959 [0.956, 0.962] | 0.930 [0.926, 0.935] | 0.941 [0.937, 0.943] | 0.976 [0.974, 0.977] |
| pydreg, −log10 p over all candidates | 0.889 | 0.959 [0.956, 0.962] | 0.921 [0.916, 0.925] | 0.926 [0.923, 0.930] | 0.975 [0.973, 0.976] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.471 | 0.754 [0.743, 0.765] | 0.705 [0.696, 0.715] | 0.819 [0.812, 0.824] | 0.958 [0.956, 0.960] |
| Δ pydreg(SVR) − PINTS | 0.419 | 0.205 [0.195, 0.215] p<1e-03 | 0.225 [0.217, 0.233] p<1e-03 | 0.122 [0.116, 0.129] p<1e-03 | 0.017 [0.016, 0.019] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.418 | 0.205 [0.195, 0.215] p<1e-03 | 0.215 [0.207, 0.224] p<1e-03 | 0.108 [0.102, 0.115] p<1e-03 | 0.016 [0.015, 0.018] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.952 [0.946, 0.957] | 0.906 [0.890, 0.919] | 0.927 [0.914, 0.937] | 0.973 [0.970, 0.975] |
| GC-matched: PINTS |  | 0.709 [0.685, 0.728] | 0.661 [0.640, 0.683] | 0.787 [0.772, 0.801] | 0.944 [0.939, 0.948] |
| GC-matched: Δ pydreg − PINTS |  | 0.243 [0.223, 0.267] p<1e-03 | 0.245 [0.219, 0.269] p<1e-03 | 0.140 [0.122, 0.157] p<1e-03 | 0.029 [0.024, 0.034] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.169 [0.146, 0.196] p<1e-03 | 0.175 [0.143, 0.204] p<1e-03 | 0.132 [0.101, 0.161] p<1e-03 | 0.023 [0.018, 0.028] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.361 [0.343, 0.379] p<1e-03 | 0.355 [0.338, 0.373] p<1e-03 | 0.342 [0.318, 0.364] p<1e-03 | 0.038 [0.033, 0.043] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.314 | 0.989 | 0.238 | 0.956 |
| 20,000 | 0.529 | 0.984 | 0.365 | 0.927 |
| 40,000 | 0.754 | 0.970 | 0.505 | 0.874 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 33,958, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.145**; CpG-island rows 0.091 (n_pos 20,129), non-CpG-island rows 0.065 (n_pos 13,829).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.707 | 0.773 [0.766, 0.779] | 0.754 [0.747, 0.760] | 0.782 [0.776, 0.788] | 0.878 [0.875, 0.882] |
| pydreg, −log10 p over all candidates | 0.707 | 0.773 [0.766, 0.779] | 0.752 [0.745, 0.758] | 0.776 [0.770, 0.781] | 0.878 [0.875, 0.881] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.254 | 0.543 [0.532, 0.553] | 0.508 [0.496, 0.520] | 0.600 [0.590, 0.610] | 0.848 [0.843, 0.851] |
| Δ pydreg(SVR) − PINTS | 0.453 | 0.230 [0.222, 0.238] p<1e-03 | 0.246 [0.236, 0.256] p<1e-03 | 0.181 [0.173, 0.190] p<1e-03 | 0.031 [0.028, 0.034] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.453 | 0.230 [0.222, 0.238] p<1e-03 | 0.243 [0.234, 0.254] p<1e-03 | 0.175 [0.167, 0.184] p<1e-03 | 0.031 [0.027, 0.034] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.770 [0.763, 0.777] | 0.744 [0.734, 0.753] | 0.777 [0.770, 0.784] | 0.876 [0.872, 0.880] |
| GC-matched: PINTS |  | 0.470 [0.452, 0.502] | 0.440 [0.420, 0.460] | 0.546 [0.527, 0.565] | 0.819 [0.811, 0.826] |
| GC-matched: Δ pydreg − PINTS |  | 0.300 [0.268, 0.316] p<1e-03 | 0.304 [0.285, 0.324] p<1e-03 | 0.232 [0.213, 0.250] p<1e-03 | 0.058 [0.051, 0.064] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.308 [0.270, 0.326] p<1e-03 | 0.318 [0.291, 0.344] p<1e-03 | 0.308 [0.275, 0.342] p<1e-03 | 0.078 [0.070, 0.086] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.285 [0.270, 0.300] p<1e-03 | 0.282 [0.269, 0.295] p<1e-03 | 0.310 [0.293, 0.327] p<1e-03 | 0.026 [0.021, 0.031] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.320 | 0.992 | 0.151 | 0.865 |
| 20,000 | 0.522 | 0.980 | 0.265 | 0.793 |
| 40,000 | 0.696 | 0.932 | 0.414 | 0.761 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 33,958, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.145**; CpG-island rows 0.091 (n_pos 20,129), non-CpG-island rows 0.065 (n_pos 13,829).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.875 | 0.893 [0.877, 0.904] | 0.880 [0.865, 0.892] | 0.894 [0.880, 0.904] | 0.942 [0.934, 0.948] |
| pydreg, −log10 p over all candidates | 0.873 | 0.893 [0.877, 0.904] | 0.879 [0.863, 0.890] | 0.891 [0.877, 0.900] | 0.942 [0.934, 0.948] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.438 | 0.655 [0.641, 0.667] | 0.639 [0.623, 0.653] | 0.745 [0.733, 0.755] | 0.921 [0.912, 0.928] |
| Δ pydreg(SVR) − PINTS | 0.438 | 0.237 [0.228, 0.246] p<1e-03 | 0.242 [0.230, 0.252] p<1e-03 | 0.148 [0.142, 0.156] p<1e-03 | 0.021 [0.019, 0.023] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.436 | 0.237 [0.228, 0.246] p<1e-03 | 0.240 [0.229, 0.251] p<1e-03 | 0.145 [0.138, 0.153] p<1e-03 | 0.021 [0.019, 0.023] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.891 [0.876, 0.903] | 0.870 [0.854, 0.886] | 0.888 [0.873, 0.902] | 0.941 [0.933, 0.947] |
| GC-matched: PINTS |  | 0.587 [0.572, 0.602] | 0.549 [0.533, 0.567] | 0.680 [0.665, 0.694] | 0.894 [0.884, 0.903] |
| GC-matched: Δ pydreg − PINTS |  | 0.304 [0.288, 0.319] p<1e-03 | 0.321 [0.304, 0.337] p<1e-03 | 0.207 [0.193, 0.224] p<1e-03 | 0.046 [0.041, 0.052] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.235 [0.220, 0.254] p<1e-03 | 0.262 [0.244, 0.281] p<1e-03 | 0.217 [0.194, 0.241] p<1e-03 | 0.045 [0.039, 0.050] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.394 [0.377, 0.412] p<1e-03 | 0.400 [0.381, 0.417] p<1e-03 | 0.406 [0.384, 0.430] p<1e-03 | 0.042 [0.036, 0.047] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.345 | 0.994 | 0.244 | 0.962 |
| 20,000 | 0.590 | 0.990 | 0.351 | 0.897 |
| 40,000 | 0.817 | 0.969 | 0.567 | 0.811 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 33,958, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.145**; CpG-island rows 0.091 (n_pos 20,129), non-CpG-island rows 0.065 (n_pos 13,829).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.667 | 0.893 [0.887, 0.898] | 0.834 [0.828, 0.839] | 0.866 [0.861, 0.872] | 0.947 [0.945, 0.949] |
| pydreg, −log10 p over all candidates | 0.645 | 0.887 [0.881, 0.892] | 0.815 [0.809, 0.821] | 0.824 [0.816, 0.832] | 0.945 [0.943, 0.947] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.802 | 0.912 [0.908, 0.916] | 0.877 [0.873, 0.882] | 0.908 [0.904, 0.912] | 0.958 [0.955, 0.960] |
| Δ pydreg(SVR) − PINTS | -0.136 | -0.019 [-0.024, -0.015] p<1e-03 | -0.044 [-0.048, -0.040] p<1e-03 | -0.042 [-0.046, -0.038] p<1e-03 | -0.011 [-0.012, -0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.157 | -0.025 [-0.030, -0.021] p<1e-03 | -0.062 [-0.066, -0.058] p<1e-03 | -0.084 [-0.091, -0.078] p<1e-03 | -0.013 [-0.014, -0.011] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.907 [0.895, 0.916] | 0.858 [0.846, 0.869] | 0.885 [0.875, 0.893] | 0.950 [0.947, 0.952] |
| GC-matched: PINTS |  | 0.914 [0.910, 0.919] | 0.896 [0.891, 0.901] | 0.922 [0.918, 0.926] | 0.959 [0.957, 0.962] |
| GC-matched: Δ pydreg − PINTS |  | -0.008 [-0.019, 0.001] p=0.082 | -0.038 [-0.049, -0.028] p<1e-03 | -0.037 [-0.046, -0.030] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.003 [-0.011, 0.008] p=0.688 | -0.025 [-0.040, -0.013] p<1e-03 | -0.035 [-0.051, -0.021] p<1e-03 | -0.005 [-0.008, -0.003] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.024 [-0.033, -0.017] p<1e-03 | -0.057 [-0.064, -0.051] p<1e-03 | -0.095 [-0.106, -0.084] p<1e-03 | -0.015 [-0.018, -0.012] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.259 | 0.973 | 0.274 | 0.991 |
| 20,000 | 0.467 | 0.961 | 0.488 | 0.983 |
| 40,000 | 0.709 | 0.903 | 0.742 | 0.951 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 46,284, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.188**; CpG-island rows 0.111 (n_pos 25,078), non-CpG-island rows 0.096 (n_pos 21,206).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.775 | 0.784 [0.780, 0.789] | 0.783 [0.779, 0.788] | 0.816 [0.811, 0.820] | 0.886 [0.884, 0.889] |
| pydreg, −log10 p over all candidates | 0.775 | 0.784 [0.780, 0.789] | 0.783 [0.778, 0.788] | 0.815 [0.810, 0.820] | 0.886 [0.884, 0.889] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.638 | 0.791 [0.784, 0.799] | 0.754 [0.747, 0.762] | 0.816 [0.807, 0.825] | 0.898 [0.894, 0.902] |
| Δ pydreg(SVR) − PINTS | 0.137 | -0.007 [-0.014, 0.000] p=0.054 | 0.029 [0.022, 0.036] p<1e-03 | -0.000 [-0.007, 0.007] p=0.964 | -0.011 [-0.015, -0.008] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.137 | -0.007 [-0.014, 0.000] p=0.054 | 0.029 [0.022, 0.036] p<1e-03 | -0.001 [-0.007, 0.006] p=0.898 | -0.011 [-0.015, -0.008] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.783 [0.778, 0.788] | 0.777 [0.767, 0.785] | 0.811 [0.803, 0.818] | 0.885 [0.883, 0.888] |
| GC-matched: PINTS |  | 0.739 [0.721, 0.752] | 0.704 [0.688, 0.719] | 0.792 [0.779, 0.804] | 0.889 [0.885, 0.893] |
| GC-matched: Δ pydreg − PINTS |  | 0.044 [0.032, 0.061] p<1e-03 | 0.073 [0.061, 0.087] p<1e-03 | 0.020 [0.010, 0.030] p<1e-03 | -0.004 [-0.007, 0.000] p=0.076 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.047 [0.034, 0.061] p<1e-03 | 0.071 [0.057, 0.086] p<1e-03 | 0.050 [0.037, 0.064] p<1e-03 | 0.009 [0.005, 0.012] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.038 [0.025, 0.051] p<1e-03 | 0.069 [0.057, 0.081] p<1e-03 | 0.076 [0.059, 0.094] p<1e-03 | -0.020 [-0.025, -0.014] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.317 | 1.000 | 0.151 | 0.997 |
| 20,000 | 0.538 | 1.000 | 0.260 | 0.996 |
| 40,000 | 0.752 | 0.993 | 0.419 | 0.985 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 17,548, n_neg = 162,340, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.098**; CpG-island rows 0.074 (n_pos 13,033), non-CpG-island rows 0.027 (n_pos 4,515).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.607 | 0.756 [0.742, 0.770] | 0.714 [0.700, 0.727] | 0.712 [0.698, 0.723] | 0.868 [0.860, 0.875] |
| pydreg, −log10 p over all candidates | 0.587 | 0.756 [0.742, 0.770] | 0.701 [0.688, 0.714] | 0.655 [0.641, 0.668] | 0.866 [0.859, 0.873] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.527 | 0.644 [0.629, 0.659] | 0.626 [0.611, 0.640] | 0.627 [0.609, 0.643] | 0.844 [0.836, 0.852] |
| Δ pydreg(SVR) − PINTS | 0.080 | 0.112 [0.101, 0.122] p<1e-03 | 0.088 [0.080, 0.097] p<1e-03 | 0.085 [0.074, 0.096] p<1e-03 | 0.023 [0.019, 0.028] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.060 | 0.112 [0.101, 0.122] p<1e-03 | 0.075 [0.067, 0.084] p<1e-03 | 0.028 [0.017, 0.039] p<1e-03 | 0.022 [0.018, 0.027] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.698 [0.669, 0.728] | 0.640 [0.609, 0.666] | 0.628 [0.583, 0.666] | 0.854 [0.846, 0.862] |
| GC-matched: PINTS |  | 0.585 [0.568, 0.600] | 0.544 [0.522, 0.566] | 0.560 [0.526, 0.592] | 0.816 [0.806, 0.825] |
| GC-matched: Δ pydreg − PINTS |  | 0.113 [0.091, 0.137] p<1e-03 | 0.097 [0.073, 0.117] p<1e-03 | 0.068 [0.027, 0.107] p=0.002 | 0.038 [0.033, 0.044] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.119 [0.093, 0.143] p<1e-03 | 0.103 [0.070, 0.130] p<1e-03 | 0.077 [0.013, 0.137] p=0.024 | 0.038 [0.033, 0.044] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.115 [0.100, 0.129] p<1e-03 | 0.073 [0.061, 0.085] p<1e-03 | 0.032 [0.017, 0.045] p<1e-03 | 0.028 [0.021, 0.036] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.443 | 0.929 | 0.286 | 0.913 |
| 20,000 | 0.628 | 0.856 | 0.435 | 0.865 |
| 40,000 | 0.740 | 0.722 | 0.646 | 0.579 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 4,706, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.023**; CpG-island rows 0.019 (n_pos 3,939), non-CpG-island rows 0.004 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.844 | 0.892 [0.880, 0.903] | 0.864 [0.853, 0.876] | 0.717 [0.700, 0.735] | 0.940 [0.934, 0.946] |
| pydreg, −log10 p over all candidates | 0.837 | 0.892 [0.880, 0.903] | 0.857 [0.845, 0.868] | 0.666 [0.650, 0.682] | 0.940 [0.934, 0.945] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.799 | 0.888 [0.878, 0.898] | 0.861 [0.850, 0.871] | 0.778 [0.765, 0.791] | 0.948 [0.942, 0.953] |
| Δ pydreg(SVR) − PINTS | 0.044 | 0.004 [-0.004, 0.012] p=0.31 | 0.004 [-0.004, 0.011] p=0.328 | -0.060 [-0.076, -0.045] p<1e-03 | -0.007 [-0.012, -0.003] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.038 | 0.004 [-0.004, 0.012] p=0.31 | -0.004 [-0.011, 0.004] p=0.316 | -0.112 [-0.127, -0.098] p<1e-03 | -0.008 [-0.012, -0.004] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.888 [0.876, 0.899] | 0.843 [0.830, 0.855] | 0.659 [0.617, 0.698] | 0.937 [0.931, 0.943] |
| GC-matched: PINTS |  | 0.846 [0.835, 0.866] | 0.813 [0.801, 0.826] | 0.678 [0.643, 0.712] | 0.933 [0.927, 0.939] |
| GC-matched: Δ pydreg − PINTS |  | 0.043 [0.025, 0.052] p<1e-03 | 0.030 [0.017, 0.041] p<1e-03 | -0.020 [-0.070, 0.032] p=0.472 | 0.004 [-0.001, 0.008] p=0.118 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.053 [0.040, 0.065] p<1e-03 | 0.040 [0.026, 0.053] p<1e-03 | -0.011 [-0.074, 0.056] p=0.768 | 0.006 [0.002, 0.009] p=0.002 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.034 [-0.065, -0.003] p=0.034 | -0.033 [-0.061, -0.005] p=0.016 | -0.093 [-0.134, -0.056] p<1e-03 | -0.031 [-0.050, -0.012] p=0.002 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.673 | 0.798 | 0.574 | 0.910 |
| 20,000 | 0.825 | 0.702 | 0.716 | 0.761 |
| 40,000 | 0.878 | 0.546 | 0.819 | 0.546 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 4,706, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.023**; CpG-island rows 0.019 (n_pos 3,939), non-CpG-island rows 0.004 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.850 | 0.897 [0.887, 0.909] | 0.868 [0.858, 0.880] | 0.708 [0.691, 0.725] | 0.943 [0.938, 0.949] |
| pydreg, −log10 p over all candidates | 0.841 | 0.897 [0.887, 0.909] | 0.861 [0.851, 0.872] | 0.663 [0.647, 0.679] | 0.942 [0.937, 0.948] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.839 | 0.895 [0.884, 0.906] | 0.875 [0.864, 0.886] | 0.802 [0.790, 0.817] | 0.949 [0.942, 0.955] |
| Δ pydreg(SVR) − PINTS | 0.012 | 0.003 [-0.003, 0.008] p=0.374 | -0.007 [-0.012, -0.001] p=0.022 | -0.095 [-0.110, -0.079] p<1e-03 | -0.005 [-0.009, -0.002] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.003 | 0.003 [-0.003, 0.008] p=0.374 | -0.014 [-0.020, -0.008] p<1e-03 | -0.140 [-0.153, -0.126] p<1e-03 | -0.006 [-0.010, -0.003] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.888 [0.875, 0.899] | 0.830 [0.813, 0.846] | 0.575 [0.525, 0.632] | 0.937 [0.931, 0.944] |
| GC-matched: PINTS |  | 0.875 [0.864, 0.886] | 0.837 [0.824, 0.850] | 0.701 [0.662, 0.737] | 0.935 [0.928, 0.942] |
| GC-matched: Δ pydreg − PINTS |  | 0.013 [0.005, 0.020] p=0.004 | -0.007 [-0.021, 0.004] p=0.258 | -0.126 [-0.179, -0.063] p<1e-03 | 0.002 [-0.002, 0.006] p=0.326 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.023 [0.017, 0.030] p<1e-03 | 0.000 [-0.015, 0.015] p=0.972 | -0.134 [-0.200, -0.052] p<1e-03 | 0.002 [-0.000, 0.005] p=0.116 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.041 [-0.075, -0.010] p=0.012 | -0.047 [-0.077, -0.020] p<1e-03 | -0.125 [-0.168, -0.083] p<1e-03 | -0.023 [-0.041, -0.005] p=0.008 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.669 | 0.801 | 0.556 | 0.917 |
| 20,000 | 0.828 | 0.711 | 0.747 | 0.826 |
| 40,000 | 0.883 | 0.543 | 0.854 | 0.598 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 4,706, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.023**; CpG-island rows 0.019 (n_pos 3,939), non-CpG-island rows 0.004 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.813 | 0.903 [0.893, 0.914] | 0.856 [0.846, 0.867] | 0.619 [0.603, 0.638] | 0.944 [0.939, 0.950] |
| pydreg, −log10 p over all candidates | 0.752 | 0.903 [0.893, 0.914] | 0.843 [0.833, 0.854] | 0.564 [0.548, 0.580] | 0.943 [0.938, 0.949] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.884 | 0.917 [0.907, 0.926] | 0.901 [0.892, 0.911] | 0.844 [0.834, 0.856] | 0.956 [0.951, 0.962] |
| Δ pydreg(SVR) − PINTS | -0.070 | -0.013 [-0.019, -0.009] p<1e-03 | -0.045 [-0.051, -0.040] p<1e-03 | -0.225 [-0.241, -0.209] p<1e-03 | -0.012 [-0.015, -0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.132 | -0.013 [-0.019, -0.009] p<1e-03 | -0.058 [-0.064, -0.053] p<1e-03 | -0.281 [-0.296, -0.266] p<1e-03 | -0.013 [-0.016, -0.010] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.885 [0.871, 0.897] | 0.788 [0.769, 0.805] | 0.434 [0.393, 0.476] | 0.935 [0.929, 0.941] |
| GC-matched: PINTS |  | 0.906 [0.897, 0.916] | 0.881 [0.871, 0.892] | 0.789 [0.756, 0.821] | 0.948 [0.942, 0.954] |
| GC-matched: Δ pydreg − PINTS |  | -0.021 [-0.032, -0.013] p<1e-03 | -0.093 [-0.111, -0.078] p<1e-03 | -0.355 [-0.407, -0.301] p<1e-03 | -0.013 [-0.017, -0.009] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.008 [-0.016, 0.001] p=0.08 | -0.093 [-0.113, -0.076] p<1e-03 | -0.401 [-0.464, -0.333] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.095 [-0.125, -0.065] p<1e-03 | -0.117 [-0.145, -0.091] p<1e-03 | -0.168 [-0.217, -0.125] p<1e-03 | -0.049 [-0.069, -0.032] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.593 | 0.715 | 0.671 | 0.932 |
| 20,000 | 0.830 | 0.637 | 0.834 | 0.843 |
| 40,000 | 0.890 | 0.468 | 0.889 | 0.584 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 4,706, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.023**; CpG-island rows 0.019 (n_pos 3,939), non-CpG-island rows 0.004 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.844 | 0.891 [0.881, 0.903] | 0.865 [0.854, 0.876] | 0.721 [0.705, 0.738] | 0.940 [0.935, 0.946] |
| pydreg, −log10 p over all candidates | 0.840 | 0.891 [0.881, 0.903] | 0.858 [0.847, 0.869] | 0.673 [0.659, 0.687] | 0.939 [0.934, 0.945] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.752 | 0.893 [0.883, 0.904] | 0.852 [0.841, 0.864] | 0.720 [0.702, 0.737] | 0.949 [0.944, 0.955] |
| Δ pydreg(SVR) − PINTS | 0.092 | -0.002 [-0.009, 0.006] p=0.672 | 0.012 [0.005, 0.019] p=0.002 | 0.001 [-0.017, 0.018] p=0.878 | -0.009 [-0.013, -0.005] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.089 | -0.002 [-0.009, 0.006] p=0.672 | 0.005 [-0.002, 0.013] p=0.16 | -0.046 [-0.064, -0.030] p<1e-03 | -0.010 [-0.013, -0.006] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.884 [0.872, 0.896] | 0.823 [0.806, 0.839] | 0.565 [0.508, 0.623] | 0.934 [0.928, 0.941] |
| GC-matched: PINTS |  | 0.774 [0.755, 0.820] | 0.739 [0.713, 0.764] | 0.508 [0.464, 0.556] | 0.929 [0.922, 0.936] |
| GC-matched: Δ pydreg − PINTS |  | 0.109 [0.067, 0.127] p<1e-03 | 0.084 [0.058, 0.109] p<1e-03 | 0.058 [-0.018, 0.125] p=0.122 | 0.005 [0.001, 0.010] p=0.026 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.130 [0.105, 0.144] p<1e-03 | 0.109 [0.078, 0.138] p<1e-03 | 0.078 [-0.012, 0.161] p=0.092 | 0.011 [0.007, 0.015] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.044 [-0.080, -0.010] p=0.014 | -0.036 [-0.068, -0.006] p=0.024 | -0.069 [-0.111, -0.032] p<1e-03 | -0.042 [-0.061, -0.024] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.667 | 0.809 | 0.507 | 0.813 |
| 20,000 | 0.814 | 0.718 | 0.592 | 0.705 |
| 40,000 | 0.874 | 0.566 | 0.767 | 0.557 |


## PINTS-on-TSS labels: pints_tss_dhs vs tss_silent

* **Positives:** pints_tss overlapping the ENCODE DNase peaks of that cell: the dREG-paper recipe (DHS ∩ TSS-assay transcription), with PINTS-on-PRO-cap in place of the GRO-cap HMM.
* **Negatives:** Genome minus (ENCODE DNase ∪ PRO-cap signal ± 250 bp ∪ positives), in informative 500 bp tiles (≤ 200k). For K562 hg19 the negatives are the dREG hub negatives minus the lifted positives ± 250 bp.
* **Notes:** Applied only to libraries other than the PRO-cap library the labels came from. The positives are PINTS calls, so they favour PINTS's element boundaries and its sensitivity profile.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 29,622, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.129**; CpG-island rows 0.082 (n_pos 17,778), non-CpG-island rows 0.056 (n_pos 11,844).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.901 | 0.966 [0.963, 0.970] | 0.938 [0.934, 0.943] | 0.944 [0.940, 0.947] | 0.979 [0.978, 0.981] |
| pydreg, −log10 p over all candidates | 0.901 | 0.966 [0.963, 0.970] | 0.929 [0.924, 0.933] | 0.928 [0.924, 0.932] | 0.979 [0.977, 0.980] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.486 | 0.772 [0.762, 0.784] | 0.721 [0.712, 0.731] | 0.812 [0.806, 0.819] | 0.962 [0.960, 0.964] |
| Δ pydreg(SVR) − PINTS | 0.415 | 0.194 [0.183, 0.204] p<1e-03 | 0.217 [0.208, 0.225] p<1e-03 | 0.131 [0.125, 0.139] p<1e-03 | 0.018 [0.016, 0.020] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.415 | 0.194 [0.183, 0.204] p<1e-03 | 0.207 [0.199, 0.215] p<1e-03 | 0.116 [0.109, 0.123] p<1e-03 | 0.017 [0.015, 0.019] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.959 [0.953, 0.964] | 0.918 [0.905, 0.930] | 0.932 [0.920, 0.941] | 0.977 [0.975, 0.979] |
| GC-matched: PINTS |  | 0.722 [0.693, 0.745] | 0.673 [0.647, 0.698] | 0.778 [0.759, 0.796] | 0.947 [0.942, 0.952] |
| GC-matched: Δ pydreg − PINTS |  | 0.237 [0.214, 0.265] p<1e-03 | 0.245 [0.217, 0.273] p<1e-03 | 0.153 [0.132, 0.175] p<1e-03 | 0.030 [0.025, 0.035] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.158 [0.133, 0.188] p<1e-03 | 0.172 [0.139, 0.206] p<1e-03 | 0.141 [0.107, 0.175] p<1e-03 | 0.022 [0.017, 0.028] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.362 [0.341, 0.385] p<1e-03 | 0.361 [0.342, 0.380] p<1e-03 | 0.370 [0.343, 0.396] p<1e-03 | 0.040 [0.036, 0.045] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.323 | 0.988 | 0.245 | 0.952 |
| 20,000 | 0.549 | 0.983 | 0.376 | 0.920 |
| 40,000 | 0.773 | 0.967 | 0.521 | 0.861 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 29,622, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.129**; CpG-island rows 0.082 (n_pos 17,778), non-CpG-island rows 0.056 (n_pos 11,844).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.734 | 0.799 [0.793, 0.805] | 0.780 [0.773, 0.787] | 0.799 [0.794, 0.805] | 0.892 [0.889, 0.896] |
| pydreg, −log10 p over all candidates | 0.734 | 0.799 [0.793, 0.805] | 0.777 [0.771, 0.784] | 0.792 [0.786, 0.799] | 0.892 [0.889, 0.895] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.276 | 0.576 [0.565, 0.587] | 0.539 [0.526, 0.551] | 0.609 [0.599, 0.619] | 0.864 [0.860, 0.867] |
| Δ pydreg(SVR) − PINTS | 0.458 | 0.223 [0.214, 0.232] p<1e-03 | 0.241 [0.231, 0.251] p<1e-03 | 0.191 [0.182, 0.199] p<1e-03 | 0.029 [0.026, 0.031] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.457 | 0.223 [0.214, 0.232] p<1e-03 | 0.239 [0.229, 0.248] p<1e-03 | 0.184 [0.175, 0.192] p<1e-03 | 0.028 [0.025, 0.031] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.796 [0.790, 0.803] | 0.772 [0.763, 0.781] | 0.796 [0.789, 0.803] | 0.890 [0.887, 0.894] |
| GC-matched: PINTS |  | 0.505 [0.488, 0.537] | 0.471 [0.453, 0.492] | 0.557 [0.538, 0.575] | 0.837 [0.830, 0.844] |
| GC-matched: Δ pydreg − PINTS |  | 0.291 [0.260, 0.307] p<1e-03 | 0.301 [0.282, 0.319] p<1e-03 | 0.240 [0.221, 0.257] p<1e-03 | 0.053 [0.047, 0.060] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.291 [0.254, 0.310] p<1e-03 | 0.307 [0.282, 0.331] p<1e-03 | 0.308 [0.277, 0.339] p<1e-03 | 0.070 [0.062, 0.078] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.290 [0.274, 0.304] p<1e-03 | 0.288 [0.275, 0.300] p<1e-03 | 0.330 [0.310, 0.346] p<1e-03 | 0.026 [0.021, 0.030] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.336 | 0.992 | 0.166 | 0.860 |
| 20,000 | 0.547 | 0.979 | 0.288 | 0.785 |
| 40,000 | 0.723 | 0.926 | 0.447 | 0.750 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 29,622, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.129**; CpG-island rows 0.082 (n_pos 17,778), non-CpG-island rows 0.056 (n_pos 11,844).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.888 | 0.905 [0.889, 0.917] | 0.893 [0.877, 0.905] | 0.902 [0.888, 0.912] | 0.949 [0.941, 0.955] |
| pydreg, −log10 p over all candidates | 0.887 | 0.905 [0.889, 0.917] | 0.892 [0.876, 0.903] | 0.898 [0.884, 0.909] | 0.949 [0.940, 0.955] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.462 | 0.676 [0.662, 0.688] | 0.658 [0.643, 0.672] | 0.743 [0.731, 0.753] | 0.926 [0.917, 0.932] |
| Δ pydreg(SVR) − PINTS | 0.426 | 0.229 [0.220, 0.238] p<1e-03 | 0.235 [0.223, 0.246] p<1e-03 | 0.159 [0.152, 0.167] p<1e-03 | 0.023 [0.021, 0.026] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.424 | 0.229 [0.220, 0.238] p<1e-03 | 0.233 [0.221, 0.244] p<1e-03 | 0.156 [0.148, 0.163] p<1e-03 | 0.023 [0.021, 0.025] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.904 [0.888, 0.915] | 0.885 [0.869, 0.899] | 0.897 [0.881, 0.910] | 0.948 [0.939, 0.954] |
| GC-matched: PINTS |  | 0.608 [0.590, 0.623] | 0.567 [0.550, 0.587] | 0.674 [0.657, 0.691] | 0.899 [0.889, 0.908] |
| GC-matched: Δ pydreg − PINTS |  | 0.296 [0.281, 0.313] p<1e-03 | 0.318 [0.299, 0.334] p<1e-03 | 0.223 [0.205, 0.242] p<1e-03 | 0.048 [0.042, 0.054] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.219 [0.203, 0.236] p<1e-03 | 0.252 [0.231, 0.272] p<1e-03 | 0.226 [0.197, 0.259] p<1e-03 | 0.044 [0.039, 0.050] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.402 [0.384, 0.418] p<1e-03 | 0.407 [0.388, 0.424] p<1e-03 | 0.432 [0.407, 0.456] p<1e-03 | 0.046 [0.040, 0.052] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.356 | 0.993 | 0.262 | 0.960 |
| 20,000 | 0.608 | 0.988 | 0.372 | 0.889 |
| 40,000 | 0.832 | 0.965 | 0.592 | 0.796 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 29,622, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.129**; CpG-island rows 0.082 (n_pos 17,778), non-CpG-island rows 0.056 (n_pos 11,844).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.717 | 0.920 [0.915, 0.924] | 0.864 [0.859, 0.870] | 0.885 [0.879, 0.890] | 0.960 [0.958, 0.962] |
| pydreg, −log10 p over all candidates | 0.692 | 0.916 [0.911, 0.920] | 0.845 [0.840, 0.851] | 0.835 [0.826, 0.844] | 0.958 [0.956, 0.960] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.843 | 0.932 [0.928, 0.936] | 0.902 [0.898, 0.906] | 0.922 [0.918, 0.926] | 0.967 [0.965, 0.969] |
| Δ pydreg(SVR) − PINTS | -0.126 | -0.013 [-0.017, -0.009] p<1e-03 | -0.038 [-0.042, -0.034] p<1e-03 | -0.038 [-0.042, -0.033] p<1e-03 | -0.007 [-0.008, -0.005] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.151 | -0.017 [-0.021, -0.013] p<1e-03 | -0.057 [-0.062, -0.052] p<1e-03 | -0.087 [-0.094, -0.080] p<1e-03 | -0.009 [-0.010, -0.007] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.934 [0.923, 0.940] | 0.889 [0.878, 0.897] | 0.905 [0.895, 0.913] | 0.963 [0.961, 0.965] |
| GC-matched: PINTS |  | 0.934 [0.931, 0.938] | 0.920 [0.915, 0.924] | 0.937 [0.933, 0.940] | 0.968 [0.966, 0.970] |
| GC-matched: Δ pydreg − PINTS |  | -0.001 [-0.011, 0.006] p=0.83 | -0.031 [-0.040, -0.023] p<1e-03 | -0.032 [-0.040, -0.025] p<1e-03 | -0.005 [-0.007, -0.004] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.007 [-0.004, 0.010] p=0.152 | -0.019 [-0.031, -0.010] p<1e-03 | -0.029 [-0.044, -0.018] p<1e-03 | -0.002 [-0.004, 0.000] p=0.06 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.017 [-0.024, -0.010] p<1e-03 | -0.050 [-0.057, -0.044] p<1e-03 | -0.095 [-0.106, -0.083] p<1e-03 | -0.011 [-0.013, -0.008] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.284 | 0.972 | 0.304 | 0.990 |
| 20,000 | 0.511 | 0.959 | 0.536 | 0.982 |
| 40,000 | 0.757 | 0.896 | 0.788 | 0.947 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 24,837, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.110**; CpG-island rows 0.066 (n_pos 14,094), non-CpG-island rows 0.051 (n_pos 10,743).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.828 | 0.835 [0.829, 0.841] | 0.834 [0.828, 0.840] | 0.845 [0.839, 0.851] | 0.913 [0.910, 0.916] |
| pydreg, −log10 p over all candidates | 0.828 | 0.835 [0.829, 0.841] | 0.834 [0.828, 0.840] | 0.844 [0.838, 0.851] | 0.913 [0.910, 0.916] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.656 | 0.788 [0.779, 0.796] | 0.757 [0.749, 0.765] | 0.780 [0.770, 0.790] | 0.895 [0.890, 0.899] |
| Δ pydreg(SVR) − PINTS | 0.172 | 0.047 [0.041, 0.054] p<1e-03 | 0.077 [0.071, 0.084] p<1e-03 | 0.065 [0.058, 0.071] p<1e-03 | 0.019 [0.016, 0.022] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.172 | 0.047 [0.041, 0.054] p<1e-03 | 0.077 [0.071, 0.083] p<1e-03 | 0.064 [0.057, 0.071] p<1e-03 | 0.018 [0.015, 0.022] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.834 [0.828, 0.840] | 0.829 [0.820, 0.837] | 0.839 [0.829, 0.848] | 0.912 [0.909, 0.915] |
| GC-matched: PINTS |  | 0.732 [0.720, 0.746] | 0.704 [0.690, 0.717] | 0.741 [0.725, 0.757] | 0.885 [0.880, 0.889] |
| GC-matched: Δ pydreg − PINTS |  | 0.102 [0.088, 0.113] p<1e-03 | 0.125 [0.113, 0.137] p<1e-03 | 0.098 [0.085, 0.112] p<1e-03 | 0.027 [0.024, 0.031] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.075 [0.066, 0.090] p<1e-03 | 0.099 [0.086, 0.112] p<1e-03 | 0.106 [0.085, 0.133] p<1e-03 | 0.028 [0.025, 0.031] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.113 [0.100, 0.125] p<1e-03 | 0.143 [0.132, 0.156] p<1e-03 | 0.221 [0.203, 0.238] p<1e-03 | 0.023 [0.018, 0.028] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.389 | 1.000 | 0.199 | 0.995 |
| 20,000 | 0.613 | 1.000 | 0.320 | 0.994 |
| 40,000 | 0.807 | 0.988 | 0.476 | 0.975 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 13,470, n_neg = 162,340, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.077**; CpG-island rows 0.064 (n_pos 11,145), non-CpG-island rows 0.014 (n_pos 2,325).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.716 | 0.872 [0.858, 0.884] | 0.826 [0.812, 0.839] | 0.801 [0.787, 0.814] | 0.928 [0.920, 0.934] |
| pydreg, −log10 p over all candidates | 0.695 | 0.872 [0.858, 0.884] | 0.811 [0.796, 0.823] | 0.729 [0.712, 0.744] | 0.927 [0.919, 0.933] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.636 | 0.752 [0.737, 0.767] | 0.729 [0.714, 0.743] | 0.698 [0.681, 0.713] | 0.902 [0.894, 0.910] |
| Δ pydreg(SVR) − PINTS | 0.080 | 0.120 [0.110, 0.129] p<1e-03 | 0.097 [0.089, 0.105] p<1e-03 | 0.104 [0.092, 0.116] p<1e-03 | 0.026 [0.022, 0.029] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.059 | 0.120 [0.110, 0.129] p<1e-03 | 0.082 [0.073, 0.090] p<1e-03 | 0.031 [0.018, 0.043] p<1e-03 | 0.024 [0.020, 0.027] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.831 [0.803, 0.850] | 0.761 [0.731, 0.785] | 0.713 [0.661, 0.758] | 0.919 [0.911, 0.926] |
| GC-matched: PINTS |  | 0.694 [0.678, 0.711] | 0.642 [0.613, 0.669] | 0.612 [0.563, 0.658] | 0.881 [0.870, 0.890] |
| GC-matched: Δ pydreg − PINTS |  | 0.137 [0.114, 0.152] p<1e-03 | 0.119 [0.094, 0.144] p<1e-03 | 0.102 [0.048, 0.156] p<1e-03 | 0.038 [0.033, 0.043] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.129 [0.113, 0.148] p<1e-03 | 0.120 [0.091, 0.150] p<1e-03 | 0.110 [0.040, 0.180] p=0.002 | 0.040 [0.034, 0.045] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.186 [0.164, 0.205] p<1e-03 | 0.118 [0.100, 0.137] p<1e-03 | 0.037 [0.017, 0.059] p<1e-03 | 0.023 [0.013, 0.032] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.531 | 0.923 | 0.353 | 0.909 |
| 20,000 | 0.739 | 0.843 | 0.529 | 0.857 |
| 40,000 | 0.857 | 0.697 | 0.754 | 0.552 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 4,154, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.020**; CpG-island rows 0.019 (n_pos 3,803), non-CpG-island rows 0.002 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.932 | 0.979 [0.973, 0.985] | 0.949 [0.942, 0.956] | 0.785 [0.769, 0.802] | 0.986 [0.983, 0.989] |
| pydreg, −log10 p over all candidates | 0.925 | 0.979 [0.973, 0.985] | 0.941 [0.934, 0.948] | 0.728 [0.714, 0.743] | 0.985 [0.982, 0.988] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.881 | 0.963 [0.957, 0.970] | 0.936 [0.929, 0.943] | 0.845 [0.833, 0.857] | 0.988 [0.985, 0.990] |
| Δ pydreg(SVR) − PINTS | 0.051 | 0.016 [0.008, 0.023] p<1e-03 | 0.013 [0.006, 0.020] p<1e-03 | -0.059 [-0.077, -0.042] p<1e-03 | -0.002 [-0.005, 0.001] p=0.302 |
| Δ pydreg(−log p) − PINTS | 0.044 | 0.016 [0.008, 0.023] p<1e-03 | 0.005 [-0.002, 0.012] p=0.148 | -0.117 [-0.134, -0.101] p<1e-03 | -0.002 [-0.006, 0.001] p=0.112 |
| GC-matched: pydreg (SVR) |  | 0.978 [0.970, 0.984] | 0.926 [0.915, 0.935] | 0.716 [0.667, 0.762] | 0.983 [0.980, 0.987] |
| GC-matched: PINTS |  | 0.923 [0.912, 0.939] | 0.887 [0.875, 0.900] | 0.731 [0.691, 0.770] | 0.979 [0.976, 0.982] |
| GC-matched: Δ pydreg − PINTS |  | 0.054 [0.039, 0.064] p<1e-03 | 0.039 [0.026, 0.051] p<1e-03 | -0.014 [-0.078, 0.046] p=0.662 | 0.004 [0.001, 0.008] p=0.006 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.050 [0.038, 0.061] p<1e-03 | 0.037 [0.024, 0.050] p<1e-03 | -0.013 [-0.085, 0.055] p=0.744 | 0.005 [0.002, 0.008] p=0.004 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.058 [0.013, 0.099] p=0.012 | 0.048 [0.008, 0.085] p=0.018 | -0.146 [-0.214, -0.083] p<1e-03 | -0.006 [-0.025, 0.014] p=0.56 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.747 | 0.795 | 0.636 | 0.908 |
| 20,000 | 0.912 | 0.697 | 0.791 | 0.756 |
| 40,000 | 0.967 | 0.539 | 0.900 | 0.538 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 4,154, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.020**; CpG-island rows 0.019 (n_pos 3,803), non-CpG-island rows 0.002 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.940 | 0.985 [0.981, 0.989] | 0.953 [0.948, 0.959] | 0.774 [0.760, 0.791] | 0.989 [0.987, 0.991] |
| pydreg, −log10 p over all candidates | 0.931 | 0.985 [0.981, 0.989] | 0.945 [0.940, 0.951] | 0.724 [0.710, 0.737] | 0.988 [0.986, 0.991] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.920 | 0.972 [0.966, 0.978] | 0.951 [0.945, 0.958] | 0.871 [0.858, 0.883] | 0.990 [0.988, 0.992] |
| Δ pydreg(SVR) − PINTS | 0.020 | 0.013 [0.008, 0.017] p<1e-03 | 0.002 [-0.002, 0.007] p=0.364 | -0.096 [-0.113, -0.079] p<1e-03 | -0.001 [-0.003, 0.001] p=0.216 |
| Δ pydreg(−log p) − PINTS | 0.012 | 0.013 [0.008, 0.017] p<1e-03 | -0.006 [-0.011, -0.002] p=0.006 | -0.147 [-0.162, -0.131] p<1e-03 | -0.002 [-0.004, -0.000] p=0.028 |
| GC-matched: pydreg (SVR) |  | 0.976 [0.968, 0.983] | 0.912 [0.897, 0.925] | 0.620 [0.563, 0.682] | 0.984 [0.982, 0.987] |
| GC-matched: PINTS |  | 0.954 [0.946, 0.962] | 0.913 [0.902, 0.925] | 0.756 [0.711, 0.801] | 0.984 [0.981, 0.986] |
| GC-matched: Δ pydreg − PINTS |  | 0.022 [0.015, 0.029] p<1e-03 | -0.001 [-0.015, 0.012] p=0.9 | -0.136 [-0.200, -0.069] p<1e-03 | 0.001 [-0.002, 0.003] p=0.51 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.022 [0.016, 0.029] p<1e-03 | -0.000 [-0.015, 0.013] p=0.92 | -0.138 [-0.206, -0.065] p<1e-03 | 0.001 [-0.001, 0.004] p=0.214 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.046 [0.011, 0.090] p=0.022 | 0.005 [-0.031, 0.040] p=0.712 | -0.214 [-0.292, -0.144] p<1e-03 | -0.011 [-0.028, 0.004] p=0.224 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.744 | 0.798 | 0.615 | 0.915 |
| 20,000 | 0.917 | 0.707 | 0.825 | 0.822 |
| 40,000 | 0.973 | 0.536 | 0.935 | 0.590 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 4,154, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.020**; CpG-island rows 0.019 (n_pos 3,803), non-CpG-island rows 0.002 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.903 | 0.992 [0.989, 0.995] | 0.941 [0.936, 0.946] | 0.675 [0.658, 0.692] | 0.991 [0.989, 0.992] |
| pydreg, −log10 p over all candidates | 0.834 | 0.992 [0.989, 0.995] | 0.926 [0.922, 0.932] | 0.613 [0.597, 0.628] | 0.989 [0.987, 0.991] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.968 | 0.989 [0.986, 0.992] | 0.976 [0.972, 0.979] | 0.914 [0.905, 0.924] | 0.996 [0.994, 0.997] |
| Δ pydreg(SVR) − PINTS | -0.065 | 0.003 [-0.000, 0.006] p=0.08 | -0.035 [-0.039, -0.031] p<1e-03 | -0.240 [-0.258, -0.220] p<1e-03 | -0.005 [-0.006, -0.004] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.134 | 0.003 [-0.000, 0.006] p=0.076 | -0.049 [-0.053, -0.045] p<1e-03 | -0.301 [-0.319, -0.286] p<1e-03 | -0.006 [-0.008, -0.005] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.972 [0.961, 0.981] | 0.863 [0.843, 0.881] | 0.458 [0.414, 0.513] | 0.983 [0.980, 0.985] |
| GC-matched: PINTS |  | 0.984 [0.979, 0.988] | 0.958 [0.952, 0.964] | 0.852 [0.813, 0.891] | 0.993 [0.991, 0.994] |
| GC-matched: Δ pydreg − PINTS |  | -0.011 [-0.021, -0.003] p=0.002 | -0.095 [-0.115, -0.078] p<1e-03 | -0.395 [-0.453, -0.330] p<1e-03 | -0.010 [-0.013, -0.008] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.008 [-0.015, -0.001] p=0.028 | -0.094 [-0.114, -0.077] p<1e-03 | -0.409 [-0.471, -0.339] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.022 [-0.062, 0.013] p=0.29 | -0.091 [-0.125, -0.057] p<1e-03 | -0.284 [-0.374, -0.209] p<1e-03 | -0.025 [-0.040, -0.009] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.659 | 0.711 | 0.743 | 0.931 |
| 20,000 | 0.920 | 0.632 | 0.920 | 0.840 |
| 40,000 | 0.981 | 0.461 | 0.972 | 0.575 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 4,154, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.020**; CpG-island rows 0.019 (n_pos 3,803), non-CpG-island rows 0.002 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.933 | 0.978 [0.972, 0.984] | 0.949 [0.943, 0.956] | 0.790 [0.774, 0.806] | 0.985 [0.983, 0.989] |
| pydreg, −log10 p over all candidates | 0.929 | 0.978 [0.972, 0.984] | 0.942 [0.935, 0.949] | 0.736 [0.722, 0.750] | 0.985 [0.982, 0.988] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.822 | 0.965 [0.958, 0.970] | 0.923 [0.914, 0.931] | 0.775 [0.757, 0.793] | 0.987 [0.985, 0.989] |
| Δ pydreg(SVR) − PINTS | 0.111 | 0.013 [0.007, 0.020] p<1e-03 | 0.026 [0.020, 0.033] p<1e-03 | 0.014 [-0.006, 0.034] p=0.154 | -0.001 [-0.004, 0.001] p=0.322 |
| Δ pydreg(−log p) − PINTS | 0.107 | 0.013 [0.007, 0.020] p<1e-03 | 0.019 [0.013, 0.025] p<1e-03 | -0.039 [-0.059, -0.020] p<1e-03 | -0.002 [-0.005, 0.000] p=0.11 |
| GC-matched: pydreg (SVR) |  | 0.972 [0.964, 0.978] | 0.903 [0.887, 0.917] | 0.607 [0.545, 0.672] | 0.980 [0.977, 0.984] |
| GC-matched: PINTS |  | 0.837 [0.824, 0.875] | 0.797 [0.769, 0.822] | 0.533 [0.486, 0.584] | 0.971 [0.967, 0.975] |
| GC-matched: Δ pydreg − PINTS |  | 0.135 [0.099, 0.148] p<1e-03 | 0.106 [0.077, 0.136] p<1e-03 | 0.073 [-0.008, 0.153] p=0.084 | 0.010 [0.006, 0.013] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.128 [0.103, 0.141] p<1e-03 | 0.108 [0.078, 0.139] p<1e-03 | 0.079 [-0.008, 0.167] p=0.072 | 0.010 [0.007, 0.014] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.094 [0.047, 0.156] p<1e-03 | 0.092 [0.050, 0.142] p<1e-03 | -0.099 [-0.172, -0.032] p=0.004 | -0.004 [-0.023, 0.016] p=0.738 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.742 | 0.806 | 0.559 | 0.809 |
| 20,000 | 0.902 | 0.713 | 0.649 | 0.699 |
| 40,000 | 0.963 | 0.559 | 0.838 | 0.548 |


## PINTS-on-TSS labels: pints_tss_distal vs tss_silent

* **Positives:** pints_tss > 1 kb from any GENCODE v24 TSS.
* **Negatives:** Genome minus (ENCODE DNase ∪ PRO-cap signal ± 250 bp ∪ positives), in informative 500 bp tiles (≤ 200k). For K562 hg19 the negatives are the dREG hub negatives minus the lifted positives ± 250 bp.
* **Notes:** Applied only to libraries other than the PRO-cap library the labels came from. The positives are PINTS calls, so they favour PINTS's element boundaries and its sensitivity profile.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 9,359, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.045**; CpG-island rows 0.010 (n_pos 1,971), non-CpG-island rows 0.036 (n_pos 7,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.835 | 0.931 [0.924, 0.938] | 0.896 [0.888, 0.903] | 0.838 [0.827, 0.847] | 0.960 [0.957, 0.964] |
| pydreg, −log10 p over all candidates | 0.835 | 0.931 [0.924, 0.938] | 0.887 [0.879, 0.894] | 0.790 [0.778, 0.801] | 0.959 [0.956, 0.963] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.280 | 0.623 [0.606, 0.640] | 0.574 [0.559, 0.588] | 0.453 [0.429, 0.480] | 0.937 [0.933, 0.940] |
| Δ pydreg(SVR) − PINTS | 0.555 | 0.309 [0.293, 0.325] p<1e-03 | 0.322 [0.309, 0.336] p<1e-03 | 0.384 [0.362, 0.407] p<1e-03 | 0.023 [0.019, 0.027] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.555 | 0.309 [0.293, 0.325] p<1e-03 | 0.313 [0.301, 0.326] p<1e-03 | 0.337 [0.314, 0.359] p<1e-03 | 0.022 [0.018, 0.026] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.925 [0.917, 0.932] | 0.877 [0.869, 0.886] | 0.808 [0.793, 0.824] | 0.958 [0.954, 0.962] |
| GC-matched: PINTS |  | 0.540 [0.515, 0.563] | 0.496 [0.476, 0.515] | 0.388 [0.360, 0.420] | 0.913 [0.908, 0.918] |
| GC-matched: Δ pydreg − PINTS |  | 0.385 [0.362, 0.409] p<1e-03 | 0.382 [0.362, 0.402] p<1e-03 | 0.420 [0.392, 0.449] p<1e-03 | 0.044 [0.039, 0.050] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.343 [0.300, 0.393] p<1e-03 | 0.343 [0.311, 0.377] p<1e-03 | 0.417 [0.345, 0.485] p<1e-03 | 0.049 [0.040, 0.059] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.397 [0.375, 0.417] p<1e-03 | 0.391 [0.373, 0.410] p<1e-03 | 0.455 [0.426, 0.481] p<1e-03 | 0.042 [0.037, 0.048] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.230 | 0.947 | 0.107 | 0.730 |
| 20,000 | 0.429 | 0.934 | 0.190 | 0.646 |
| 40,000 | 0.674 | 0.889 | 0.314 | 0.542 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 9,359, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.045**; CpG-island rows 0.010 (n_pos 1,971), non-CpG-island rows 0.036 (n_pos 7,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.594 | 0.676 [0.664, 0.688] | 0.654 [0.642, 0.666] | 0.602 [0.588, 0.616] | 0.827 [0.821, 0.833] |
| pydreg, −log10 p over all candidates | 0.591 | 0.676 [0.664, 0.688] | 0.652 [0.640, 0.664] | 0.581 [0.565, 0.597] | 0.827 [0.821, 0.833] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.174 | 0.445 [0.431, 0.460] | 0.419 [0.404, 0.434] | 0.282 [0.261, 0.307] | 0.809 [0.802, 0.815] |
| Δ pydreg(SVR) − PINTS | 0.420 | 0.232 [0.219, 0.244] p<1e-03 | 0.236 [0.223, 0.248] p<1e-03 | 0.321 [0.304, 0.335] p<1e-03 | 0.019 [0.013, 0.024] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.417 | 0.232 [0.219, 0.244] p<1e-03 | 0.234 [0.222, 0.246] p<1e-03 | 0.300 [0.283, 0.316] p<1e-03 | 0.019 [0.013, 0.024] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.672 [0.660, 0.684] | 0.642 [0.629, 0.654] | 0.581 [0.562, 0.600] | 0.825 [0.818, 0.831] |
| GC-matched: PINTS |  | 0.365 [0.348, 0.383] | 0.339 [0.323, 0.356] | 0.228 [0.205, 0.257] | 0.784 [0.776, 0.791] |
| GC-matched: Δ pydreg − PINTS |  | 0.308 [0.292, 0.321] p<1e-03 | 0.303 [0.289, 0.314] p<1e-03 | 0.354 [0.333, 0.371] p<1e-03 | 0.041 [0.035, 0.047] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.357 [0.332, 0.380] p<1e-03 | 0.356 [0.332, 0.379] p<1e-03 | 0.385 [0.326, 0.436] p<1e-03 | 0.083 [0.071, 0.095] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.290 [0.273, 0.304] p<1e-03 | 0.286 [0.272, 0.297] p<1e-03 | 0.349 [0.328, 0.365] p<1e-03 | 0.030 [0.024, 0.035] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.181 | 0.952 | 0.101 | 0.542 |
| 20,000 | 0.380 | 0.910 | 0.185 | 0.425 |
| 40,000 | 0.579 | 0.759 | 0.311 | 0.398 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 9,359, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.045**; CpG-island rows 0.010 (n_pos 1,971), non-CpG-island rows 0.036 (n_pos 7,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.768 | 0.800 [0.787, 0.813] | 0.783 [0.770, 0.795] | 0.741 [0.727, 0.754] | 0.893 [0.886, 0.900] |
| pydreg, −log10 p over all candidates | 0.766 | 0.800 [0.787, 0.813] | 0.782 [0.768, 0.794] | 0.729 [0.715, 0.742] | 0.893 [0.886, 0.900] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.234 | 0.461 [0.447, 0.475] | 0.456 [0.441, 0.472] | 0.334 [0.317, 0.351] | 0.864 [0.856, 0.871] |
| Δ pydreg(SVR) − PINTS | 0.534 | 0.340 [0.324, 0.354] p<1e-03 | 0.327 [0.311, 0.343] p<1e-03 | 0.407 [0.384, 0.427] p<1e-03 | 0.029 [0.024, 0.035] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.532 | 0.340 [0.324, 0.354] p<1e-03 | 0.326 [0.310, 0.342] p<1e-03 | 0.395 [0.373, 0.415] p<1e-03 | 0.029 [0.024, 0.034] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.799 [0.786, 0.811] | 0.778 [0.764, 0.790] | 0.733 [0.720, 0.747] | 0.892 [0.885, 0.898] |
| GC-matched: PINTS |  | 0.382 [0.369, 0.396] | 0.359 [0.346, 0.373] | 0.270 [0.254, 0.289] | 0.836 [0.828, 0.844] |
| GC-matched: Δ pydreg − PINTS |  | 0.417 [0.398, 0.433] p<1e-03 | 0.419 [0.401, 0.434] p<1e-03 | 0.463 [0.440, 0.484] p<1e-03 | 0.056 [0.049, 0.062] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.392 [0.363, 0.422] p<1e-03 | 0.411 [0.386, 0.436] p<1e-03 | 0.546 [0.490, 0.603] p<1e-03 | 0.075 [0.064, 0.085] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.420 [0.401, 0.438] p<1e-03 | 0.416 [0.398, 0.433] p<1e-03 | 0.475 [0.452, 0.495] p<1e-03 | 0.049 [0.042, 0.056] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.169 | 0.955 | 0.091 | 0.723 |
| 20,000 | 0.382 | 0.944 | 0.171 | 0.539 |
| 40,000 | 0.674 | 0.877 | 0.344 | 0.418 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 9,359, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.045**; CpG-island rows 0.010 (n_pos 1,971), non-CpG-island rows 0.036 (n_pos 7,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.630 | 0.855 [0.846, 0.864] | 0.798 [0.790, 0.806] | 0.702 [0.690, 0.714] | 0.929 [0.925, 0.933] |
| pydreg, −log10 p over all candidates | 0.602 | 0.850 [0.841, 0.859] | 0.780 [0.772, 0.789] | 0.601 [0.584, 0.618] | 0.927 [0.923, 0.932] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.750 | 0.870 [0.862, 0.878] | 0.832 [0.824, 0.840] | 0.776 [0.763, 0.789] | 0.936 [0.932, 0.940] |
| Δ pydreg(SVR) − PINTS | -0.120 | -0.015 [-0.023, -0.007] p<1e-03 | -0.034 [-0.041, -0.027] p<1e-03 | -0.073 [-0.085, -0.060] p<1e-03 | -0.007 [-0.010, -0.003] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.148 | -0.020 [-0.027, -0.012] p<1e-03 | -0.052 [-0.059, -0.044] p<1e-03 | -0.175 [-0.188, -0.161] p<1e-03 | -0.008 [-0.012, -0.005] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.864 [0.855, 0.873] | 0.810 [0.800, 0.818] | 0.720 [0.704, 0.735] | 0.931 [0.926, 0.935] |
| GC-matched: PINTS |  | 0.871 [0.863, 0.879] | 0.841 [0.833, 0.850] | 0.799 [0.785, 0.812] | 0.936 [0.931, 0.940] |
| GC-matched: Δ pydreg − PINTS |  | -0.007 [-0.015, 0.001] p=0.096 | -0.032 [-0.040, -0.024] p<1e-03 | -0.079 [-0.093, -0.065] p<1e-03 | -0.005 [-0.009, -0.002] p=0.008 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.023 [0.006, 0.038] p=0.008 | -0.009 [-0.025, 0.007] p=0.27 | -0.168 [-0.229, -0.110] p<1e-03 | 0.004 [-0.004, 0.012] p=0.372 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.015 [-0.023, -0.007] p<1e-03 | -0.038 [-0.045, -0.030] p<1e-03 | -0.095 [-0.110, -0.079] p<1e-03 | -0.008 [-0.012, -0.004] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.223 | 0.895 | 0.161 | 0.945 |
| 20,000 | 0.412 | 0.856 | 0.380 | 0.926 |
| 40,000 | 0.674 | 0.709 | 0.680 | 0.830 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 14,595, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.068**; CpG-island rows 0.014 (n_pos 2,806), non-CpG-island rows 0.056 (n_pos 11,789).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.598 | 0.614 [0.604, 0.625] | 0.613 [0.603, 0.623] | 0.619 [0.608, 0.630] | 0.797 [0.792, 0.802] |
| pydreg, −log10 p over all candidates | 0.598 | 0.614 [0.604, 0.625] | 0.613 [0.603, 0.623] | 0.617 [0.606, 0.628] | 0.797 [0.792, 0.802] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.381 | 0.598 [0.584, 0.613] | 0.554 [0.541, 0.568] | 0.506 [0.485, 0.526] | 0.803 [0.796, 0.811] |
| Δ pydreg(SVR) − PINTS | 0.218 | 0.016 [0.005, 0.027] p=0.002 | 0.059 [0.049, 0.069] p<1e-03 | 0.114 [0.099, 0.128] p<1e-03 | -0.006 [-0.011, -0.001] p=0.014 |
| Δ pydreg(−log p) − PINTS | 0.218 | 0.016 [0.005, 0.027] p=0.002 | 0.059 [0.049, 0.069] p<1e-03 | 0.112 [0.097, 0.126] p<1e-03 | -0.007 [-0.011, -0.001] p=0.014 |
| GC-matched: pydreg (SVR) |  | 0.613 [0.603, 0.623] | 0.609 [0.599, 0.619] | 0.611 [0.598, 0.623] | 0.796 [0.791, 0.801] |
| GC-matched: PINTS |  | 0.525 [0.510, 0.541] | 0.494 [0.479, 0.509] | 0.447 [0.424, 0.471] | 0.792 [0.784, 0.799] |
| GC-matched: Δ pydreg − PINTS |  | 0.087 [0.075, 0.100] p<1e-03 | 0.115 [0.103, 0.126] p<1e-03 | 0.164 [0.146, 0.182] p<1e-03 | 0.004 [-0.001, 0.009] p=0.138 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.120 [0.097, 0.142] p<1e-03 | 0.151 [0.129, 0.174] p<1e-03 | 0.346 [0.286, 0.404] p<1e-03 | 0.020 [0.010, 0.030] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.076 [0.063, 0.089] p<1e-03 | 0.104 [0.092, 0.115] p<1e-03 | 0.185 [0.168, 0.204] p<1e-03 | -0.000 [-0.006, 0.005] p=0.928 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.101 | 0.999 | 0.041 | 0.961 |
| 20,000 | 0.280 | 0.999 | 0.087 | 0.963 |
| 40,000 | 0.561 | 0.972 | 0.184 | 0.900 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 3,093, n_neg = 162,340, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.019**; CpG-island rows 0.005 (n_pos 797), non-CpG-island rows 0.014 (n_pos 2,296).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.192 | 0.337 [0.316, 0.359] | 0.306 [0.287, 0.326] | 0.113 [0.100, 0.129] | 0.648 [0.637, 0.659] |
| pydreg, −log10 p over all candidates | 0.177 | 0.337 [0.316, 0.359] | 0.301 [0.283, 0.321] | 0.084 [0.073, 0.096] | 0.647 [0.637, 0.659] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.167 | 0.269 [0.250, 0.290] | 0.262 [0.244, 0.282] | 0.081 [0.069, 0.096] | 0.617 [0.604, 0.632] |
| Δ pydreg(SVR) − PINTS | 0.025 | 0.068 [0.052, 0.083] p<1e-03 | 0.044 [0.031, 0.056] p<1e-03 | 0.032 [0.020, 0.044] p<1e-03 | 0.030 [0.022, 0.039] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.009 | 0.068 [0.052, 0.083] p<1e-03 | 0.039 [0.026, 0.052] p<1e-03 | 0.003 [-0.009, 0.012] p=0.592 | 0.030 [0.021, 0.038] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.330 [0.307, 0.351] | 0.290 [0.271, 0.310] | 0.090 [0.078, 0.107] | 0.643 [0.632, 0.655] |
| GC-matched: PINTS |  | 0.236 [0.219, 0.255] | 0.229 [0.213, 0.247] | 0.067 [0.057, 0.081] | 0.598 [0.585, 0.612] |
| GC-matched: Δ pydreg − PINTS |  | 0.094 [0.079, 0.109] p<1e-03 | 0.061 [0.049, 0.074] p<1e-03 | 0.023 [0.010, 0.037] p<1e-03 | 0.045 [0.036, 0.053] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.160 [0.126, 0.192] p<1e-03 | 0.124 [0.093, 0.154] p<1e-03 | 0.041 [0.005, 0.086] p=0.026 | 0.046 [0.030, 0.061] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.070 [0.054, 0.086] p<1e-03 | 0.042 [0.029, 0.055] p<1e-03 | 0.003 [-0.003, 0.008] p=0.27 | 0.042 [0.032, 0.051] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.112 | 0.367 | 0.070 | 0.313 |
| 20,000 | 0.205 | 0.255 | 0.129 | 0.251 |
| 40,000 | 0.310 | 0.161 | 0.272 | 0.093 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 336, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.002**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.002 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.074 | 0.127 [0.096, 0.160] | 0.124 [0.095, 0.154] | 0.004 [0.002, 0.007] | 0.540 [0.524, 0.557] |
| pydreg, −log10 p over all candidates | 0.080 | 0.127 [0.096, 0.160] | 0.123 [0.094, 0.155] | 0.003 [0.002, 0.005] | 0.540 [0.524, 0.557] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.115 | 0.223 [0.175, 0.272] | 0.204 [0.162, 0.249] | 0.016 [0.007, 0.030] | 0.593 [0.556, 0.629] |
| Δ pydreg(SVR) − PINTS | -0.041 | -0.096 [-0.138, -0.056] p<1e-03 | -0.081 [-0.117, -0.048] p<1e-03 | -0.012 [-0.026, -0.003] p<1e-03 | -0.053 [-0.084, -0.021] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.035 | -0.096 [-0.138, -0.056] p<1e-03 | -0.081 [-0.117, -0.048] p<1e-03 | -0.012 [-0.027, -0.004] p<1e-03 | -0.053 [-0.084, -0.021] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.124 [0.093, 0.157] | 0.121 [0.092, 0.152] | 0.004 [0.002, 0.007] | 0.539 [0.523, 0.556] |
| GC-matched: PINTS |  | 0.204 [0.163, 0.256] | 0.194 [0.153, 0.237] | 0.015 [0.006, 0.029] | 0.585 [0.549, 0.621] |
| GC-matched: Δ pydreg − PINTS |  | -0.080 [-0.126, -0.046] p<1e-03 | -0.073 [-0.108, -0.041] p<1e-03 | -0.011 [-0.025, -0.002] p<1e-03 | -0.046 [-0.078, -0.014] p=0.004 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.129 [-0.037, 0.278] p=0.176 | 0.072 [-0.034, 0.168] p=0.184 | -0.033 [-0.119, 0.055] p=0.574 | 0.046 [-0.030, 0.115] p=0.192 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.108 [-0.151, -0.062] p<1e-03 | -0.090 [-0.126, -0.054] p<1e-03 | -0.008 [-0.019, -0.002] p<1e-03 | -0.056 [-0.090, -0.022] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.060 | 0.024 | 0.039 | 0.046 |
| 20,000 | 0.071 | 0.014 | 0.092 | 0.028 |
| 40,000 | 0.101 | 0.010 | 0.137 | 0.014 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 336, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.002**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.002 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.074 | 0.134 [0.099, 0.171] | 0.128 [0.097, 0.161] | 0.004 [0.002, 0.009] | 0.543 [0.526, 0.563] |
| pydreg, −log10 p over all candidates | 0.077 | 0.134 [0.099, 0.171] | 0.128 [0.096, 0.161] | 0.003 [0.002, 0.004] | 0.543 [0.526, 0.563] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.119 | 0.219 [0.172, 0.269] | 0.201 [0.159, 0.246] | 0.015 [0.007, 0.031] | 0.584 [0.552, 0.617] |
| Δ pydreg(SVR) − PINTS | -0.045 | -0.086 [-0.126, -0.046] p<1e-03 | -0.073 [-0.108, -0.043] p<1e-03 | -0.011 [-0.025, -0.002] p=0.024 | -0.040 [-0.068, -0.014] p=0.002 |
| Δ pydreg(−log p) − PINTS | -0.042 | -0.086 [-0.126, -0.046] p<1e-03 | -0.073 [-0.108, -0.042] p<1e-03 | -0.012 [-0.027, -0.004] p<1e-03 | -0.040 [-0.068, -0.014] p=0.002 |
| GC-matched: pydreg (SVR) |  | 0.130 [0.096, 0.167] | 0.124 [0.093, 0.157] | 0.004 [0.002, 0.009] | 0.541 [0.524, 0.560] |
| GC-matched: PINTS |  | 0.210 [0.165, 0.263] | 0.192 [0.152, 0.238] | 0.014 [0.006, 0.029] | 0.576 [0.543, 0.609] |
| GC-matched: Δ pydreg − PINTS |  | -0.080 [-0.122, -0.042] p<1e-03 | -0.069 [-0.102, -0.036] p<1e-03 | -0.011 [-0.023, -0.001] p=0.028 | -0.035 [-0.063, -0.008] p=0.008 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.039 [-0.081, 0.191] p=0.684 | 0.030 [-0.080, 0.157] p=0.652 | -0.015 [-0.099, 0.009] p=0.336 | 0.051 [-0.019, 0.109] p=0.154 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.095 [-0.137, -0.056] p<1e-03 | -0.081 [-0.116, -0.049] p<1e-03 | -0.008 [-0.020, 0.000] p=0.06 | -0.044 [-0.073, -0.017] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.036 | 0.015 | 0.042 | 0.056 |
| 20,000 | 0.071 | 0.015 | 0.077 | 0.034 |
| 40,000 | 0.095 | 0.009 | 0.131 | 0.016 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 336, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.002**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.002 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.074 | 0.119 [0.082, 0.157] | 0.114 [0.082, 0.147] | 0.003 [0.002, 0.004] | 0.535 [0.516, 0.555] |
| pydreg, −log10 p over all candidates | 0.051 | 0.119 [0.082, 0.157] | 0.113 [0.081, 0.145] | 0.002 [0.002, 0.003] | 0.535 [0.516, 0.555] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.134 | 0.234 [0.188, 0.288] | 0.212 [0.170, 0.260] | 0.015 [0.008, 0.030] | 0.589 [0.557, 0.626] |
| Δ pydreg(SVR) − PINTS | -0.060 | -0.115 [-0.163, -0.074] p<1e-03 | -0.098 [-0.139, -0.064] p<1e-03 | -0.012 [-0.027, -0.005] p<1e-03 | -0.054 [-0.086, -0.025] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.083 | -0.115 [-0.163, -0.074] p<1e-03 | -0.100 [-0.140, -0.066] p<1e-03 | -0.013 [-0.028, -0.006] p<1e-03 | -0.055 [-0.087, -0.025] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.119 [0.082, 0.157] | 0.110 [0.079, 0.142] | 0.003 [0.002, 0.004] | 0.532 [0.514, 0.552] |
| GC-matched: PINTS |  | 0.230 [0.184, 0.283] | 0.206 [0.165, 0.252] | 0.014 [0.007, 0.030] | 0.583 [0.551, 0.620] |
| GC-matched: Δ pydreg − PINTS |  | -0.111 [-0.159, -0.070] p<1e-03 | -0.096 [-0.136, -0.062] p<1e-03 | -0.012 [-0.027, -0.005] p<1e-03 | -0.051 [-0.083, -0.021] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.000 [-0.154, 0.160] p=1 | -0.022 [-0.113, 0.068] p=0.602 | -0.027 [-0.121, -0.001] p=0.008 | 0.031 [-0.039, 0.093] p=0.36 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.129 [-0.178, -0.084] p<1e-03 | -0.109 [-0.151, -0.075] p<1e-03 | -0.007 [-0.018, -0.003] p<1e-03 | -0.060 [-0.095, -0.029] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.036 | 0.011 | 0.042 | 0.057 |
| 20,000 | 0.074 | 0.011 | 0.092 | 0.041 |
| 40,000 | 0.095 | 0.007 | 0.143 | 0.016 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 336, n_neg = 200,000, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.002**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.002 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.077 | 0.126 [0.093, 0.161] | 0.124 [0.092, 0.156] | 0.003 [0.002, 0.006] | 0.540 [0.522, 0.558] |
| pydreg, −log10 p over all candidates | 0.086 | 0.126 [0.093, 0.161] | 0.124 [0.092, 0.156] | 0.003 [0.002, 0.005] | 0.540 [0.522, 0.558] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.122 | 0.248 [0.201, 0.300] | 0.223 [0.182, 0.269] | 0.013 [0.006, 0.026] | 0.597 [0.564, 0.632] |
| Δ pydreg(SVR) − PINTS | -0.045 | -0.122 [-0.169, -0.079] p<1e-03 | -0.100 [-0.143, -0.063] p<1e-03 | -0.010 [-0.022, -0.003] p<1e-03 | -0.057 [-0.089, -0.028] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.036 | -0.122 [-0.169, -0.079] p<1e-03 | -0.099 [-0.143, -0.062] p<1e-03 | -0.010 [-0.023, -0.003] p<1e-03 | -0.057 [-0.089, -0.028] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.122 [0.089, 0.157] | 0.120 [0.089, 0.152] | 0.003 [0.002, 0.005] | 0.538 [0.520, 0.556] |
| GC-matched: PINTS |  | 0.229 [0.189, 0.283] | 0.209 [0.170, 0.251] | 0.011 [0.005, 0.023] | 0.589 [0.557, 0.621] |
| GC-matched: Δ pydreg − PINTS |  | -0.107 [-0.157, -0.069] p<1e-03 | -0.089 [-0.131, -0.053] p<1e-03 | -0.008 [-0.019, -0.002] p<1e-03 | -0.051 [-0.082, -0.022] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.032 [-0.087, 0.167] p=0.852 | 0.034 [-0.028, 0.093] p=0.252 | -0.014 [-0.102, 0.003] p=0.658 | 0.015 [-0.053, 0.075] p=0.642 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.124 [-0.174, -0.080] p<1e-03 | -0.101 [-0.146, -0.062] p<1e-03 | -0.007 [-0.017, -0.002] p<1e-03 | -0.058 [-0.091, -0.027] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.048 | 0.021 | 0.062 | 0.037 |
| 20,000 | 0.074 | 0.016 | 0.086 | 0.024 |
| 40,000 | 0.101 | 0.011 | 0.140 | 0.016 |


## dREG paper positives vs TSS-silent DHSs

* **Positives:** dREG hub positives (DHS ∩ GRO-cap HMM; K562, GM12878).
* **Negatives:** DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).
* **Notes:** Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces dhs_unmarked as the specificity test.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 27,798, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.613**; CpG-island rows 0.375 (n_pos 10,493), non-CpG-island rows 0.497 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.805 | 0.924 [0.920, 0.927] | 0.889 [0.883, 0.895] | 0.967 [0.965, 0.969] | 0.956 [0.954, 0.958] |
| pydreg, −log10 p over all candidates | 0.811 | 0.924 [0.920, 0.927] | 0.882 [0.874, 0.889] | 0.966 [0.964, 0.968] | 0.956 [0.953, 0.958] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.517 | 0.791 [0.779, 0.804] | 0.742 [0.728, 0.756] | 0.963 [0.961, 0.965] | 0.950 [0.947, 0.952] |
| Δ pydreg(SVR) − PINTS | 0.288 | 0.133 [0.119, 0.144] p<1e-03 | 0.147 [0.133, 0.160] p<1e-03 | 0.004 [0.002, 0.006] p<1e-03 | 0.007 [0.004, 0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.295 | 0.133 [0.119, 0.144] p<1e-03 | 0.140 [0.126, 0.154] p<1e-03 | 0.003 [0.001, 0.005] p=0.004 | 0.006 [0.004, 0.008] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.921 [0.911, 0.926] | 0.870 [0.855, 0.883] | 0.966 [0.963, 0.968] | 0.954 [0.951, 0.956] |
| GC-matched: PINTS |  | 0.789 [0.770, 0.804] | 0.739 [0.721, 0.757] | 0.963 [0.960, 0.965] | 0.948 [0.945, 0.951] |
| GC-matched: Δ pydreg − PINTS |  | 0.132 [0.116, 0.148] p<1e-03 | 0.131 [0.116, 0.146] p<1e-03 | 0.003 [0.001, 0.005] p=0.002 | 0.006 [0.003, 0.009] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.042 [0.033, 0.058] p<1e-03 | 0.046 [0.031, 0.061] p<1e-03 | 0.011 [0.006, 0.017] p<1e-03 | 0.008 [0.005, 0.010] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.185 [0.165, 0.202] p<1e-03 | 0.197 [0.179, 0.215] p<1e-03 | 0.014 [0.010, 0.018] p<1e-03 | 0.006 [0.002, 0.009] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.294 | 1.000 | 0.170 | 0.999 |
| 20,000 | 0.494 | 0.999 | 0.259 | 0.998 |
| 40,000 | 0.689 | 0.996 | 0.359 | 0.995 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 27,798, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.613**; CpG-island rows 0.375 (n_pos 10,493), non-CpG-island rows 0.497 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.690 | 0.703 [0.696, 0.709] | 0.701 [0.694, 0.708] | 0.879 [0.874, 0.884] | 0.844 [0.840, 0.847] |
| pydreg, −log10 p over all candidates | 0.690 | 0.703 [0.696, 0.709] | 0.701 [0.694, 0.708] | 0.879 [0.874, 0.884] | 0.843 [0.840, 0.847] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.638 | 0.713 [0.706, 0.720] | 0.688 [0.680, 0.695] | 0.882 [0.877, 0.887] | 0.847 [0.843, 0.850] |
| Δ pydreg(SVR) − PINTS | 0.053 | -0.010 [-0.015, -0.005] p<1e-03 | 0.013 [0.008, 0.020] p<1e-03 | -0.002 [-0.004, -0.000] p=0.028 | -0.003 [-0.006, -0.000] p=0.03 |
| Δ pydreg(−log p) − PINTS | 0.053 | -0.010 [-0.015, -0.005] p<1e-03 | 0.013 [0.007, 0.019] p<1e-03 | -0.002 [-0.005, -0.000] p=0.02 | -0.003 [-0.006, -0.000] p=0.028 |
| GC-matched: pydreg (SVR) |  | 0.703 [0.696, 0.709] | 0.702 [0.694, 0.708] | 0.879 [0.874, 0.884] | 0.844 [0.840, 0.847] |
| GC-matched: PINTS |  | 0.714 [0.707, 0.721] | 0.691 [0.683, 0.698] | 0.882 [0.877, 0.887] | 0.847 [0.843, 0.851] |
| GC-matched: Δ pydreg − PINTS |  | -0.011 [-0.016, -0.006] p<1e-03 | 0.011 [0.004, 0.017] p=0.002 | -0.003 [-0.005, -0.000] p=0.012 | -0.003 [-0.006, -0.001] p=0.014 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.047 [0.041, 0.053] p<1e-03 | 0.061 [0.052, 0.070] p<1e-03 | 0.034 [0.030, 0.038] p<1e-03 | 0.026 [0.023, 0.029] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.045 [-0.052, -0.038] p<1e-03 | -0.018 [-0.025, -0.010] p<1e-03 | -0.019 [-0.023, -0.015] p<1e-03 | -0.021 [-0.025, -0.018] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.288 | 1.000 | 0.134 | 0.998 |
| 20,000 | 0.458 | 1.000 | 0.232 | 0.998 |
| 40,000 | 0.615 | 0.999 | 0.356 | 0.998 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 27,798, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.613**; CpG-island rows 0.375 (n_pos 10,493), non-CpG-island rows 0.497 (n_pos 17,305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.721 | 0.732 [0.721, 0.741] | 0.728 [0.717, 0.737] | 0.891 [0.886, 0.896] | 0.859 [0.853, 0.863] |
| pydreg, −log10 p over all candidates | 0.721 | 0.732 [0.721, 0.741] | 0.728 [0.717, 0.737] | 0.891 [0.885, 0.896] | 0.859 [0.853, 0.863] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.577 | 0.756 [0.744, 0.766] | 0.712 [0.701, 0.722] | 0.898 [0.893, 0.903] | 0.867 [0.861, 0.872] |
| Δ pydreg(SVR) − PINTS | 0.144 | -0.024 [-0.029, -0.019] p<1e-03 | 0.016 [0.009, 0.024] p<1e-03 | -0.007 [-0.009, -0.005] p<1e-03 | -0.008 [-0.011, -0.006] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.144 | -0.024 [-0.029, -0.019] p<1e-03 | 0.016 [0.008, 0.024] p=0.002 | -0.007 [-0.009, -0.005] p<1e-03 | -0.009 [-0.011, -0.006] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.731 [0.719, 0.740] | 0.720 [0.709, 0.731] | 0.891 [0.885, 0.896] | 0.857 [0.851, 0.862] |
| GC-matched: PINTS |  | 0.754 [0.743, 0.764] | 0.702 [0.688, 0.714] | 0.897 [0.891, 0.902] | 0.865 [0.860, 0.871] |
| GC-matched: Δ pydreg − PINTS |  | -0.024 [-0.029, -0.018] p<1e-03 | 0.019 [0.009, 0.029] p<1e-03 | -0.006 [-0.009, -0.004] p<1e-03 | -0.008 [-0.011, -0.005] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.031 [0.026, 0.041] p<1e-03 | 0.057 [0.043, 0.072] p<1e-03 | 0.026 [0.021, 0.032] p<1e-03 | 0.019 [0.016, 0.022] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.056 [-0.064, -0.049] p<1e-03 | -0.003 [-0.013, 0.007] p=0.618 | -0.020 [-0.025, -0.017] p<1e-03 | -0.024 [-0.028, -0.020] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.284 | 1.000 | 0.186 | 1.000 |
| 20,000 | 0.453 | 1.000 | 0.251 | 0.999 |
| 40,000 | 0.628 | 0.998 | 0.387 | 0.998 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 36,627, n_neg = 16,074, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.695**; CpG-island rows 0.500 (n_pos 16,084), non-CpG-island rows 0.561 (n_pos 20,543).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.867 | 0.957 [0.955, 0.960] | 0.931 [0.926, 0.936] | 0.986 [0.985, 0.986] | 0.975 [0.973, 0.977] |
| pydreg, −log10 p over all candidates | 0.870 | 0.957 [0.955, 0.960] | 0.919 [0.911, 0.925] | 0.984 [0.983, 0.985] | 0.974 [0.972, 0.975] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.731 | 0.937 [0.932, 0.941] | 0.886 [0.877, 0.894] | 0.988 [0.987, 0.989] | 0.979 [0.977, 0.980] |
| Δ pydreg(SVR) − PINTS | 0.136 | 0.021 [0.016, 0.026] p<1e-03 | 0.046 [0.038, 0.054] p<1e-03 | -0.003 [-0.004, -0.002] p<1e-03 | -0.004 [-0.005, -0.002] p=0.002 |
| Δ pydreg(−log p) − PINTS | 0.139 | 0.021 [0.016, 0.026] p<1e-03 | 0.033 [0.023, 0.043] p<1e-03 | -0.005 [-0.006, -0.003] p<1e-03 | -0.005 [-0.007, -0.003] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.957 [0.954, 0.960] | 0.930 [0.923, 0.937] | 0.986 [0.985, 0.987] | 0.975 [0.973, 0.977] |
| GC-matched: PINTS |  | 0.938 [0.929, 0.950] | 0.888 [0.874, 0.900] | 0.989 [0.987, 0.989] | 0.979 [0.977, 0.980] |
| GC-matched: Δ pydreg − PINTS |  | 0.020 [0.007, 0.027] p<1e-03 | 0.043 [0.030, 0.056] p<1e-03 | -0.003 [-0.004, -0.002] p<1e-03 | -0.004 [-0.006, -0.002] p=0.002 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.010 [-0.003, 0.016] p=0.21 | 0.027 [0.014, 0.041] p<1e-03 | -0.000 [-0.002, 0.002] p=0.684 | -0.001 [-0.003, 0.001] p=0.404 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.029 [0.022, 0.038] p<1e-03 | 0.059 [0.048, 0.072] p<1e-03 | -0.004 [-0.006, -0.002] p=0.002 | -0.006 [-0.008, -0.004] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.283 | 1.000 | 0.210 | 1.000 |
| 20,000 | 0.492 | 1.000 | 0.335 | 0.999 |
| 40,000 | 0.721 | 0.999 | 0.473 | 0.998 |


## dREG auxiliary labels: dhs_h3k27ac vs TSS-silent DHSs

* **Positives:** Duke ∩ UW DHS overlapping an ENCODE Broad H3K27ac peak, the "DHSs marked by H3K27ac" of Wang 2019 Fig 1B (K562 34,763; GM12878 55,925).
* **Negatives:** DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).
* **Notes:** Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces dhs_unmarked as the specificity test.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 33,874, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.659**; CpG-island rows 0.362 (n_pos 9,941), non-CpG-island rows 0.577 (n_pos 23,933).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.656 | 0.803 [0.796, 0.810] | 0.765 [0.756, 0.774] | 0.929 [0.926, 0.932] | 0.892 [0.889, 0.896] |
| pydreg, −log10 p over all candidates | 0.663 | 0.803 [0.796, 0.810] | 0.759 [0.750, 0.768] | 0.928 [0.925, 0.931] | 0.892 [0.888, 0.896] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.450 | 0.740 [0.724, 0.756] | 0.692 [0.676, 0.706] | 0.954 [0.952, 0.957] | 0.927 [0.924, 0.930] |
| Δ pydreg(SVR) − PINTS | 0.207 | 0.063 [0.047, 0.078] p<1e-03 | 0.073 [0.058, 0.088] p<1e-03 | -0.025 [-0.028, -0.023] p<1e-03 | -0.035 [-0.038, -0.032] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.213 | 0.063 [0.047, 0.078] p<1e-03 | 0.067 [0.053, 0.082] p<1e-03 | -0.026 [-0.029, -0.024] p<1e-03 | -0.036 [-0.039, -0.032] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.801 [0.793, 0.809] | 0.750 [0.737, 0.763] | 0.928 [0.924, 0.931] | 0.890 [0.886, 0.894] |
| GC-matched: PINTS |  | 0.735 [0.715, 0.753] | 0.688 [0.672, 0.706] | 0.954 [0.951, 0.956] | 0.926 [0.922, 0.929] |
| GC-matched: Δ pydreg − PINTS |  | 0.066 [0.047, 0.084] p<1e-03 | 0.062 [0.046, 0.077] p<1e-03 | -0.026 [-0.028, -0.023] p<1e-03 | -0.036 [-0.039, -0.032] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.021 [0.011, 0.035] p<1e-03 | 0.025 [0.009, 0.043] p<1e-03 | 0.000 [-0.005, 0.008] p=0.952 | -0.002 [-0.004, 0.001] p=0.132 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.081 [0.059, 0.100] p<1e-03 | 0.088 [0.070, 0.107] p<1e-03 | -0.038 [-0.042, -0.033] p<1e-03 | -0.049 [-0.053, -0.044] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.233 | 1.000 | 0.141 | 0.999 |
| 20,000 | 0.389 | 0.999 | 0.216 | 0.998 |
| 40,000 | 0.548 | 0.996 | 0.303 | 0.995 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 33,874, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.659**; CpG-island rows 0.362 (n_pos 9,941), non-CpG-island rows 0.577 (n_pos 23,933).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.553 | 0.571 [0.562, 0.580] | 0.569 [0.561, 0.579] | 0.846 [0.840, 0.852] | 0.774 [0.770, 0.779] |
| pydreg, −log10 p over all candidates | 0.553 | 0.571 [0.562, 0.580] | 0.569 [0.561, 0.578] | 0.846 [0.840, 0.852] | 0.774 [0.770, 0.779] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.553 | 0.629 [0.620, 0.638] | 0.606 [0.597, 0.615] | 0.866 [0.860, 0.872] | 0.802 [0.798, 0.807] |
| Δ pydreg(SVR) − PINTS | -0.001 | -0.058 [-0.063, -0.052] p<1e-03 | -0.036 [-0.042, -0.030] p<1e-03 | -0.020 [-0.021, -0.018] p<1e-03 | -0.028 [-0.031, -0.025] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.001 | -0.058 [-0.063, -0.052] p<1e-03 | -0.037 [-0.042, -0.031] p<1e-03 | -0.020 [-0.022, -0.018] p<1e-03 | -0.028 [-0.031, -0.025] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.571 [0.563, 0.580] | 0.569 [0.561, 0.579] | 0.846 [0.840, 0.852] | 0.774 [0.770, 0.779] |
| GC-matched: PINTS |  | 0.629 [0.621, 0.639] | 0.607 [0.598, 0.617] | 0.866 [0.861, 0.872] | 0.803 [0.798, 0.808] |
| GC-matched: Δ pydreg − PINTS |  | -0.058 [-0.064, -0.053] p<1e-03 | -0.038 [-0.044, -0.032] p<1e-03 | -0.020 [-0.022, -0.018] p<1e-03 | -0.029 [-0.031, -0.026] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.029 [0.022, 0.036] p<1e-03 | 0.043 [0.034, 0.052] p<1e-03 | 0.023 [0.018, 0.027] p<1e-03 | 0.017 [0.013, 0.020] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.094 [-0.101, -0.087] p<1e-03 | -0.070 [-0.077, -0.063] p<1e-03 | -0.039 [-0.042, -0.036] p<1e-03 | -0.047 [-0.051, -0.044] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.230 | 1.000 | 0.112 | 0.998 |
| 20,000 | 0.360 | 1.000 | 0.194 | 0.998 |
| 40,000 | 0.485 | 0.999 | 0.299 | 0.998 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 33,874, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.659**; CpG-island rows 0.362 (n_pos 9,941), non-CpG-island rows 0.577 (n_pos 23,933).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.597 | 0.613 [0.601, 0.624] | 0.609 [0.597, 0.621] | 0.861 [0.855, 0.867] | 0.796 [0.790, 0.802] |
| pydreg, −log10 p over all candidates | 0.597 | 0.613 [0.601, 0.624] | 0.609 [0.597, 0.620] | 0.861 [0.855, 0.867] | 0.796 [0.790, 0.802] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.519 | 0.696 [0.683, 0.707] | 0.653 [0.641, 0.664] | 0.889 [0.883, 0.894] | 0.836 [0.829, 0.842] |
| Δ pydreg(SVR) − PINTS | 0.078 | -0.083 [-0.089, -0.078] p<1e-03 | -0.043 [-0.051, -0.036] p<1e-03 | -0.028 [-0.030, -0.026] p<1e-03 | -0.040 [-0.043, -0.037] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.078 | -0.083 [-0.089, -0.078] p<1e-03 | -0.044 [-0.052, -0.036] p<1e-03 | -0.028 [-0.030, -0.026] p<1e-03 | -0.040 [-0.043, -0.037] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.612 [0.599, 0.623] | 0.604 [0.591, 0.616] | 0.861 [0.855, 0.866] | 0.795 [0.788, 0.801] |
| GC-matched: PINTS |  | 0.695 [0.681, 0.706] | 0.645 [0.632, 0.657] | 0.888 [0.883, 0.893] | 0.834 [0.828, 0.840] |
| GC-matched: Δ pydreg − PINTS |  | -0.083 [-0.089, -0.077] p<1e-03 | -0.041 [-0.050, -0.032] p<1e-03 | -0.027 [-0.030, -0.025] p<1e-03 | -0.039 [-0.043, -0.036] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.010 [0.006, 0.018] p<1e-03 | 0.038 [0.025, 0.054] p<1e-03 | 0.014 [0.009, 0.021] p<1e-03 | 0.008 [0.006, 0.011] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.121 [-0.129, -0.114] p<1e-03 | -0.073 [-0.083, -0.064] p<1e-03 | -0.048 [-0.052, -0.045] p<1e-03 | -0.059 [-0.063, -0.055] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.227 | 1.000 | 0.153 | 1.000 |
| 20,000 | 0.360 | 1.000 | 0.210 | 0.999 |
| 40,000 | 0.506 | 0.998 | 0.329 | 0.998 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 54,478, n_neg = 16,074, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.772**; CpG-island rows 0.522 (n_pos 17,574), non-CpG-island rows 0.697 (n_pos 36,904).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.713 | 0.851 [0.846, 0.856] | 0.814 [0.807, 0.821] | 0.964 [0.963, 0.966] | 0.918 [0.915, 0.921] |
| pydreg, −log10 p over all candidates | 0.715 | 0.851 [0.846, 0.856] | 0.804 [0.795, 0.813] | 0.963 [0.961, 0.965] | 0.917 [0.914, 0.920] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.591 | 0.845 [0.836, 0.854] | 0.789 [0.778, 0.800] | 0.979 [0.977, 0.980] | 0.947 [0.945, 0.949] |
| Δ pydreg(SVR) − PINTS | 0.122 | 0.006 [-0.002, 0.015] p=0.124 | 0.025 [0.015, 0.036] p<1e-03 | -0.014 [-0.015, -0.013] p<1e-03 | -0.029 [-0.031, -0.027] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.124 | 0.006 [-0.002, 0.015] p=0.126 | 0.015 [0.003, 0.027] p=0.01 | -0.015 [-0.017, -0.014] p<1e-03 | -0.030 [-0.032, -0.027] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.851 [0.845, 0.856] | 0.813 [0.803, 0.822] | 0.964 [0.963, 0.966] | 0.918 [0.915, 0.921] |
| GC-matched: PINTS |  | 0.845 [0.833, 0.858] | 0.790 [0.773, 0.805] | 0.979 [0.977, 0.980] | 0.947 [0.945, 0.950] |
| GC-matched: Δ pydreg − PINTS |  | 0.006 [-0.008, 0.018] p=0.316 | 0.023 [0.009, 0.038] p<1e-03 | -0.014 [-0.015, -0.013] p<1e-03 | -0.029 [-0.032, -0.027] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.000 [-0.017, 0.010] p=0.832 | 0.018 [0.000, 0.035] p=0.048 | -0.009 [-0.011, -0.006] p<1e-03 | -0.010 [-0.013, -0.007] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.010 [-0.002, 0.021] p=0.082 | 0.029 [0.016, 0.043] p<1e-03 | -0.023 [-0.025, -0.021] p<1e-03 | -0.038 [-0.041, -0.035] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.201 | 1.000 | 0.151 | 1.000 |
| 20,000 | 0.358 | 1.000 | 0.243 | 0.999 |
| 40,000 | 0.557 | 0.999 | 0.352 | 0.998 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 32,404, n_neg = 11,691, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.735**; CpG-island rows 0.515 (n_pos 12,409), non-CpG-island rows 0.631 (n_pos 19,995).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.619 | 0.635 [0.624, 0.645] | 0.635 [0.624, 0.645] | 0.898 [0.890, 0.905] | 0.808 [0.802, 0.813] |
| pydreg, −log10 p over all candidates | 0.619 | 0.635 [0.624, 0.645] | 0.635 [0.624, 0.645] | 0.898 [0.890, 0.905] | 0.808 [0.802, 0.813] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.655 | 0.669 [0.658, 0.681] | 0.664 [0.654, 0.676] | 0.908 [0.900, 0.915] | 0.825 [0.820, 0.832] |
| Δ pydreg(SVR) − PINTS | -0.036 | -0.034 [-0.039, -0.029] p<1e-03 | -0.030 [-0.035, -0.024] p<1e-03 | -0.010 [-0.011, -0.008] p<1e-03 | -0.018 [-0.020, -0.015] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.036 | -0.034 [-0.039, -0.029] p<1e-03 | -0.030 [-0.035, -0.024] p<1e-03 | -0.010 [-0.011, -0.008] p<1e-03 | -0.018 [-0.020, -0.015] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.635 [0.624, 0.645] | 0.634 [0.624, 0.645] | 0.898 [0.890, 0.905] | 0.808 [0.802, 0.813] |
| GC-matched: PINTS |  | 0.667 [0.656, 0.679] | 0.655 [0.638, 0.671] | 0.907 [0.900, 0.914] | 0.823 [0.817, 0.830] |
| GC-matched: Δ pydreg − PINTS |  | -0.032 [-0.038, -0.026] p<1e-03 | -0.020 [-0.031, -0.003] p=0.016 | -0.009 [-0.011, -0.008] p<1e-03 | -0.016 [-0.019, -0.012] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.023 [0.017, 0.029] p<1e-03 | 0.038 [0.022, 0.062] p<1e-03 | 0.013 [0.010, 0.016] p<1e-03 | 0.014 [0.010, 0.017] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.068 [-0.075, -0.060] p<1e-03 | -0.060 [-0.068, -0.051] p<1e-03 | -0.026 [-0.029, -0.023] p<1e-03 | -0.035 [-0.039, -0.031] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.289 | 1.000 | 0.122 | 1.000 |
| 20,000 | 0.446 | 1.000 | 0.203 | 1.000 |
| 40,000 | 0.594 | 1.000 | 0.319 | 1.000 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 57,093, n_neg = 62,494, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.477**; CpG-island rows 0.260 (n_pos 21,963), non-CpG-island rows 0.360 (n_pos 35,130).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.226 | 0.414 [0.402, 0.427] | 0.389 [0.379, 0.399] | 0.708 [0.698, 0.718] | 0.721 [0.716, 0.725] |
| pydreg, −log10 p over all candidates | 0.186 | 0.411 [0.400, 0.423] | 0.376 [0.366, 0.387] | 0.693 [0.682, 0.704] | 0.719 [0.715, 0.724] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.178 | 0.427 [0.408, 0.448] | 0.393 [0.380, 0.405] | 0.731 [0.719, 0.742] | 0.754 [0.748, 0.759] |
| Δ pydreg(SVR) − PINTS | 0.048 | -0.013 [-0.027, 0.001] p=0.078 | -0.003 [-0.010, 0.003] p=0.368 | -0.023 [-0.026, -0.020] p<1e-03 | -0.033 [-0.036, -0.030] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.007 | -0.016 [-0.030, -0.002] p=0.02 | -0.016 [-0.023, -0.010] p<1e-03 | -0.038 [-0.041, -0.034] p<1e-03 | -0.034 [-0.037, -0.031] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.275 [0.249, 0.295] | 0.253 [0.237, 0.269] | 0.661 [0.645, 0.676] | 0.685 [0.679, 0.692] |
| GC-matched: PINTS |  | 0.285 [0.264, 0.303] | 0.259 [0.241, 0.276] | 0.680 [0.664, 0.695] | 0.714 [0.707, 0.720] |
| GC-matched: Δ pydreg − PINTS |  | -0.010 [-0.031, 0.006] p=0.216 | -0.006 [-0.019, 0.007] p=0.412 | -0.019 [-0.027, -0.011] p<1e-03 | -0.028 [-0.032, -0.025] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.005 [-0.080, 0.051] p=0.76 | -0.001 [-0.042, 0.038] p=0.924 | 0.018 [-0.012, 0.049] p=0.24 | 0.001 [-0.005, 0.008] p=0.748 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.027 [-0.036, -0.018] p<1e-03 | -0.025 [-0.032, -0.018] p<1e-03 | -0.044 [-0.050, -0.038] p<1e-03 | -0.049 [-0.053, -0.045] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.193 | 0.963 | 0.117 | 0.954 |
| 20,000 | 0.325 | 0.920 | 0.198 | 0.941 |
| 40,000 | 0.473 | 0.861 | 0.395 | 0.890 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 13,004, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.173**; CpG-island rows 0.146 (n_pos 10,613), non-CpG-island rows 0.037 (n_pos 2,391).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.412 | 0.771 [0.755, 0.787] | 0.700 [0.686, 0.714] | 0.798 [0.789, 0.808] | 0.936 [0.931, 0.940] |
| pydreg, −log10 p over all candidates | 0.223 | 0.755 [0.738, 0.770] | 0.652 [0.637, 0.668] | 0.748 [0.738, 0.759] | 0.931 [0.926, 0.935] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.417 | 0.735 [0.712, 0.761] | 0.670 [0.654, 0.687] | 0.787 [0.774, 0.800] | 0.939 [0.935, 0.942] |
| Δ pydreg(SVR) − PINTS | -0.005 | 0.036 [0.019, 0.051] p<1e-03 | 0.030 [0.019, 0.040] p<1e-03 | 0.011 [-0.001, 0.024] p=0.074 | -0.003 [-0.006, 0.000] p=0.062 |
| Δ pydreg(−log p) − PINTS | -0.193 | 0.020 [0.003, 0.035] p=0.026 | -0.018 [-0.029, -0.006] p<1e-03 | -0.039 [-0.051, -0.026] p<1e-03 | -0.008 [-0.011, -0.005] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.568 [0.537, 0.609] | 0.511 [0.482, 0.542] | 0.663 [0.637, 0.692] | 0.898 [0.892, 0.905] |
| GC-matched: PINTS |  | 0.555 [0.521, 0.600] | 0.517 [0.494, 0.541] | 0.684 [0.664, 0.705] | 0.899 [0.892, 0.905] |
| GC-matched: Δ pydreg − PINTS |  | 0.013 [-0.020, 0.057] p=0.404 | -0.006 [-0.035, 0.023] p=0.704 | -0.021 [-0.052, 0.007] p=0.14 | -0.000 [-0.006, 0.005] p=0.882 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.021 [-0.020, 0.060] p=0.298 | -0.007 [-0.042, 0.025] p=0.65 | -0.025 [-0.062, 0.013] p=0.182 | 0.005 [-0.000, 0.011] p=0.06 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.002 [-0.029, 0.028] p=0.902 | -0.005 [-0.026, 0.017] p=0.64 | -0.023 [-0.046, -0.002] p=0.032 | -0.025 [-0.034, -0.016] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.525 | 0.879 | 0.363 | 0.906 |
| 20,000 | 0.753 | 0.778 | 0.514 | 0.850 |
| 40,000 | 0.899 | 0.620 | 0.672 | 0.786 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 13,004, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.173**; CpG-island rows 0.146 (n_pos 10,613), non-CpG-island rows 0.037 (n_pos 2,391).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.434 | 0.784 [0.767, 0.801] | 0.714 [0.699, 0.729] | 0.814 [0.804, 0.824] | 0.941 [0.937, 0.946] |
| pydreg, −log10 p over all candidates | 0.221 | 0.761 [0.742, 0.780] | 0.659 [0.643, 0.674] | 0.754 [0.744, 0.765] | 0.936 [0.931, 0.940] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.368 | 0.755 [0.739, 0.775] | 0.678 [0.659, 0.696] | 0.795 [0.782, 0.809] | 0.944 [0.940, 0.947] |
| Δ pydreg(SVR) − PINTS | 0.066 | 0.029 [0.016, 0.040] p<1e-03 | 0.037 [0.026, 0.048] p<1e-03 | 0.019 [0.008, 0.032] p<1e-03 | -0.002 [-0.005, 0.000] p=0.076 |
| Δ pydreg(−log p) − PINTS | -0.148 | 0.006 [-0.008, 0.018] p=0.47 | -0.019 [-0.029, -0.008] p<1e-03 | -0.041 [-0.052, -0.029] p<1e-03 | -0.008 [-0.011, -0.005] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.537 [0.497, 0.578] | 0.481 [0.452, 0.512] | 0.648 [0.622, 0.676] | 0.896 [0.890, 0.903] |
| GC-matched: PINTS |  | 0.554 [0.528, 0.587] | 0.510 [0.484, 0.539] | 0.689 [0.665, 0.713] | 0.903 [0.897, 0.909] |
| GC-matched: Δ pydreg − PINTS |  | -0.017 [-0.055, 0.023] p=0.434 | -0.030 [-0.061, 0.003] p=0.08 | -0.041 [-0.071, -0.006] p=0.028 | -0.007 [-0.011, -0.002] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.019 [-0.062, 0.026] p=0.438 | -0.035 [-0.073, 0.005] p=0.088 | -0.046 [-0.083, 0.001] p=0.054 | -0.002 [-0.007, 0.004] p=0.502 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.001 [-0.027, 0.027] p=0.91 | -0.015 [-0.034, 0.003] p=0.112 | -0.031 [-0.051, -0.012] p=0.004 | -0.028 [-0.036, -0.020] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.532 | 0.881 | 0.361 | 0.889 |
| 20,000 | 0.761 | 0.781 | 0.565 | 0.838 |
| 40,000 | 0.910 | 0.617 | 0.756 | 0.758 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 13,004, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.173**; CpG-island rows 0.146 (n_pos 10,613), non-CpG-island rows 0.037 (n_pos 2,391).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.274 | 0.737 [0.722, 0.752] | 0.651 [0.636, 0.665] | 0.769 [0.758, 0.780] | 0.945 [0.942, 0.948] |
| pydreg, −log10 p over all candidates | 0.137 | 0.684 [0.663, 0.700] | 0.574 [0.560, 0.588] | 0.704 [0.692, 0.716] | 0.937 [0.934, 0.940] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.474 | 0.756 [0.742, 0.771] | 0.705 [0.689, 0.721] | 0.825 [0.810, 0.838] | 0.954 [0.951, 0.956] |
| Δ pydreg(SVR) − PINTS | -0.200 | -0.019 [-0.031, -0.008] p<1e-03 | -0.054 [-0.066, -0.041] p<1e-03 | -0.056 [-0.069, -0.040] p<1e-03 | -0.008 [-0.010, -0.006] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.336 | -0.072 [-0.088, -0.060] p<1e-03 | -0.131 [-0.144, -0.118] p<1e-03 | -0.121 [-0.135, -0.106] p<1e-03 | -0.016 [-0.018, -0.014] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.348 [0.306, 0.407] | 0.339 [0.308, 0.371] | 0.538 [0.513, 0.567] | 0.879 [0.873, 0.886] |
| GC-matched: PINTS |  | 0.595 [0.568, 0.611] | 0.542 [0.521, 0.563] | 0.718 [0.701, 0.735] | 0.912 [0.908, 0.917] |
| GC-matched: Δ pydreg − PINTS |  | -0.247 [-0.288, -0.191] p<1e-03 | -0.202 [-0.233, -0.173] p<1e-03 | -0.180 [-0.205, -0.155] p<1e-03 | -0.033 [-0.038, -0.029] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.283 [-0.316, -0.220] p<1e-03 | -0.226 [-0.261, -0.191] p<1e-03 | -0.209 [-0.238, -0.179] p<1e-03 | -0.032 [-0.038, -0.027] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.112 [-0.146, -0.079] p<1e-03 | -0.114 [-0.136, -0.093] p<1e-03 | -0.101 [-0.120, -0.085] p<1e-03 | -0.035 [-0.041, -0.029] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.498 | 0.834 | 0.455 | 0.912 |
| 20,000 | 0.754 | 0.745 | 0.687 | 0.815 |
| 40,000 | 0.916 | 0.582 | 0.838 | 0.691 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 13,004, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.173**; CpG-island rows 0.146 (n_pos 10,613), non-CpG-island rows 0.037 (n_pos 2,391).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.453 | 0.773 [0.756, 0.788] | 0.708 [0.693, 0.723] | 0.811 [0.802, 0.820] | 0.938 [0.934, 0.943] |
| pydreg, −log10 p over all candidates | 0.249 | 0.758 [0.738, 0.776] | 0.664 [0.648, 0.680] | 0.761 [0.751, 0.771] | 0.934 [0.929, 0.938] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.387 | 0.765 [0.745, 0.779] | 0.679 [0.661, 0.696] | 0.793 [0.777, 0.806] | 0.943 [0.940, 0.947] |
| Δ pydreg(SVR) − PINTS | 0.065 | 0.008 [-0.003, 0.021] p=0.178 | 0.029 [0.018, 0.041] p<1e-03 | 0.018 [0.006, 0.033] p<1e-03 | -0.005 [-0.009, -0.002] p=0.002 |
| Δ pydreg(−log p) − PINTS | -0.138 | -0.007 [-0.021, 0.007] p=0.334 | -0.015 [-0.027, -0.001] p=0.036 | -0.032 [-0.044, -0.017] p<1e-03 | -0.009 [-0.013, -0.006] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.520 [0.482, 0.553] | 0.464 [0.435, 0.494] | 0.639 [0.614, 0.664] | 0.891 [0.885, 0.898] |
| GC-matched: PINTS |  | 0.554 [0.519, 0.570] | 0.484 [0.462, 0.508] | 0.666 [0.647, 0.685] | 0.900 [0.894, 0.905] |
| GC-matched: Δ pydreg − PINTS |  | -0.034 [-0.064, 0.006] p=0.1 | -0.020 [-0.051, 0.008] p=0.168 | -0.027 [-0.051, -0.001] p=0.036 | -0.008 [-0.014, -0.003] p=0.002 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.033 [-0.075, 0.022] p=0.204 | -0.024 [-0.059, 0.011] p=0.182 | -0.032 [-0.062, 0.000] p=0.056 | -0.004 [-0.010, 0.002] p=0.192 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.014 [-0.015, 0.042] p=0.34 | 0.005 [-0.017, 0.025] p=0.662 | -0.014 [-0.032, 0.003] p=0.112 | -0.026 [-0.035, -0.017] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.527 | 0.886 | 0.322 | 0.907 |
| 20,000 | 0.750 | 0.780 | 0.403 | 0.879 |
| 40,000 | 0.899 | 0.619 | 0.595 | 0.824 |


## dREG auxiliary labels: dhs_h3k27ac_distal vs TSS-silent DHSs

* **Positives:** As dhs_h3k27ac, but > 1 kb from any GENCODE v24 TSS, i.e. putative enhancers (K562 20,553; GM12878 28,215).
* **Negatives:** DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).
* **Notes:** Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces dhs_unmarked as the specificity test.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 19,783, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.530**; CpG-island rows 0.067 (n_pos 1,252), non-CpG-island rows 0.514 (n_pos 18,531).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.533 | 0.724 [0.714, 0.733] | 0.676 [0.665, 0.688] | 0.861 [0.856, 0.867] | 0.850 [0.845, 0.855] |
| pydreg, −log10 p over all candidates | 0.542 | 0.724 [0.714, 0.733] | 0.672 [0.660, 0.683] | 0.859 [0.854, 0.865] | 0.850 [0.844, 0.854] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.276 | 0.629 [0.608, 0.650] | 0.576 [0.556, 0.595] | 0.897 [0.892, 0.902] | 0.896 [0.891, 0.900] |
| Δ pydreg(SVR) − PINTS | 0.257 | 0.096 [0.075, 0.116] p<1e-03 | 0.100 [0.083, 0.119] p<1e-03 | -0.036 [-0.040, -0.031] p<1e-03 | -0.046 [-0.050, -0.041] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.266 | 0.096 [0.075, 0.116] p<1e-03 | 0.095 [0.078, 0.115] p<1e-03 | -0.038 [-0.043, -0.033] p<1e-03 | -0.046 [-0.050, -0.041] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.724 [0.715, 0.733] | 0.677 [0.666, 0.688] | 0.861 [0.856, 0.867] | 0.850 [0.845, 0.855] |
| GC-matched: PINTS |  | 0.626 [0.607, 0.648] | 0.575 [0.555, 0.595] | 0.897 [0.892, 0.902] | 0.896 [0.891, 0.899] |
| GC-matched: Δ pydreg − PINTS |  | 0.098 [0.075, 0.119] p<1e-03 | 0.102 [0.084, 0.120] p<1e-03 | -0.035 [-0.040, -0.031] p<1e-03 | -0.046 [-0.050, -0.041] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.062 [0.026, 0.091] p<1e-03 | 0.069 [0.043, 0.098] p<1e-03 | 0.065 [0.033, 0.100] p<1e-03 | -0.014 [-0.023, -0.005] p=0.004 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.100 [0.077, 0.121] p<1e-03 | 0.104 [0.086, 0.123] p<1e-03 | -0.036 [-0.042, -0.031] p<1e-03 | -0.048 [-0.052, -0.043] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.120 | 0.999 | 0.036 | 0.990 |
| 20,000 | 0.233 | 0.997 | 0.071 | 0.990 |
| 40,000 | 0.403 | 0.991 | 0.130 | 0.979 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 19,783, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.530**; CpG-island rows 0.067 (n_pos 1,252), non-CpG-island rows 0.514 (n_pos 18,531).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.392 | 0.417 [0.407, 0.426] | 0.415 [0.406, 0.425] | 0.712 [0.704, 0.720] | 0.693 [0.688, 0.698] |
| pydreg, −log10 p over all candidates | 0.392 | 0.417 [0.407, 0.426] | 0.415 [0.406, 0.425] | 0.712 [0.704, 0.720] | 0.693 [0.688, 0.698] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.411 | 0.505 [0.494, 0.514] | 0.480 [0.469, 0.491] | 0.752 [0.744, 0.759] | 0.737 [0.731, 0.742] |
| Δ pydreg(SVR) − PINTS | -0.019 | -0.088 [-0.094, -0.082] p<1e-03 | -0.065 [-0.072, -0.058] p<1e-03 | -0.040 [-0.043, -0.036] p<1e-03 | -0.044 [-0.047, -0.041] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.019 | -0.088 [-0.094, -0.082] p<1e-03 | -0.066 [-0.072, -0.058] p<1e-03 | -0.040 [-0.043, -0.036] p<1e-03 | -0.044 [-0.047, -0.041] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.416 [0.407, 0.426] | 0.415 [0.406, 0.425] | 0.712 [0.704, 0.720] | 0.693 [0.688, 0.698] |
| GC-matched: PINTS |  | 0.504 [0.494, 0.514] | 0.479 [0.468, 0.489] | 0.751 [0.744, 0.759] | 0.737 [0.731, 0.742] |
| GC-matched: Δ pydreg − PINTS |  | -0.088 [-0.094, -0.081] p<1e-03 | -0.065 [-0.072, -0.057] p<1e-03 | -0.039 [-0.043, -0.036] p<1e-03 | -0.044 [-0.047, -0.040] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.011 [-0.035, 0.012] p=0.354 | 0.009 [-0.014, 0.032] p=0.416 | 0.044 [0.018, 0.073] p=0.002 | -0.004 [-0.016, 0.008] p=0.552 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.093 [-0.100, -0.086] p<1e-03 | -0.069 [-0.077, -0.062] p<1e-03 | -0.043 [-0.046, -0.039] p<1e-03 | -0.047 [-0.050, -0.043] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.082 | 1.000 | 0.045 | 0.993 |
| 20,000 | 0.184 | 1.000 | 0.091 | 0.994 |
| 40,000 | 0.313 | 0.998 | 0.159 | 0.993 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 19,783, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.530**; CpG-island rows 0.067 (n_pos 1,252), non-CpG-island rows 0.514 (n_pos 18,531).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.427 | 0.450 [0.441, 0.460] | 0.446 [0.436, 0.456] | 0.728 [0.721, 0.736] | 0.710 [0.705, 0.715] |
| pydreg, −log10 p over all candidates | 0.427 | 0.450 [0.441, 0.460] | 0.446 [0.436, 0.456] | 0.728 [0.720, 0.735] | 0.710 [0.705, 0.715] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.344 | 0.566 [0.555, 0.577] | 0.514 [0.503, 0.525] | 0.778 [0.771, 0.785] | 0.767 [0.761, 0.772] |
| Δ pydreg(SVR) − PINTS | 0.083 | -0.116 [-0.124, -0.108] p<1e-03 | -0.068 [-0.078, -0.059] p<1e-03 | -0.050 [-0.054, -0.046] p<1e-03 | -0.056 [-0.060, -0.052] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.083 | -0.116 [-0.124, -0.108] p<1e-03 | -0.068 [-0.078, -0.059] p<1e-03 | -0.050 [-0.054, -0.046] p<1e-03 | -0.056 [-0.061, -0.052] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.450 [0.441, 0.460] | 0.446 [0.436, 0.456] | 0.728 [0.721, 0.736] | 0.710 [0.705, 0.716] |
| GC-matched: PINTS |  | 0.567 [0.556, 0.578] | 0.515 [0.504, 0.526] | 0.778 [0.771, 0.785] | 0.767 [0.761, 0.772] |
| GC-matched: Δ pydreg − PINTS |  | -0.116 [-0.124, -0.109] p<1e-03 | -0.069 [-0.079, -0.059] p<1e-03 | -0.050 [-0.054, -0.046] p<1e-03 | -0.056 [-0.061, -0.052] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.021 [-0.042, -0.000] p=0.05 | 0.021 [-0.002, 0.044] p=0.084 | 0.071 [0.036, 0.112] p<1e-03 | -0.007 [-0.018, 0.004] p=0.212 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.123 [-0.131, -0.115] p<1e-03 | -0.075 [-0.085, -0.065] p<1e-03 | -0.054 [-0.058, -0.050] p<1e-03 | -0.060 [-0.064, -0.055] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.064 | 0.999 | 0.028 | 0.998 |
| 20,000 | 0.156 | 0.999 | 0.060 | 0.995 |
| 40,000 | 0.313 | 0.996 | 0.135 | 0.991 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 26,959, n_neg = 16,074, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.626**; CpG-island rows 0.088 (n_pos 1,544), non-CpG-island rows 0.613 (n_pos 25,415).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.608 | 0.787 [0.781, 0.794] | 0.740 [0.731, 0.747] | 0.916 [0.912, 0.919] | 0.883 [0.880, 0.886] |
| pydreg, −log10 p over all candidates | 0.609 | 0.787 [0.781, 0.794] | 0.731 [0.722, 0.740] | 0.913 [0.909, 0.916] | 0.882 [0.879, 0.886] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.429 | 0.752 [0.739, 0.767] | 0.692 [0.679, 0.707] | 0.940 [0.937, 0.943] | 0.917 [0.914, 0.921] |
| Δ pydreg(SVR) − PINTS | 0.179 | 0.035 [0.020, 0.048] p<1e-03 | 0.047 [0.033, 0.060] p<1e-03 | -0.024 [-0.027, -0.022] p<1e-03 | -0.034 [-0.038, -0.031] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.180 | 0.035 [0.020, 0.048] p<1e-03 | 0.039 [0.023, 0.052] p<1e-03 | -0.027 [-0.030, -0.024] p<1e-03 | -0.035 [-0.039, -0.032] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.787 [0.781, 0.794] | 0.740 [0.732, 0.748] | 0.916 [0.912, 0.919] | 0.883 [0.880, 0.887] |
| GC-matched: PINTS |  | 0.752 [0.739, 0.766] | 0.692 [0.679, 0.707] | 0.940 [0.937, 0.943] | 0.917 [0.914, 0.920] |
| GC-matched: Δ pydreg − PINTS |  | 0.035 [0.021, 0.049] p<1e-03 | 0.048 [0.034, 0.061] p<1e-03 | -0.024 [-0.027, -0.022] p<1e-03 | -0.034 [-0.038, -0.030] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.008 [-0.020, 0.027] p=0.622 | 0.026 [0.000, 0.048] p=0.046 | 0.023 [-0.001, 0.047] p=0.054 | -0.022 [-0.030, -0.015] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.037 [0.021, 0.051] p<1e-03 | 0.050 [0.035, 0.063] p<1e-03 | -0.025 [-0.028, -0.022] p<1e-03 | -0.035 [-0.038, -0.031] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.096 | 1.000 | 0.040 | 0.997 |
| 20,000 | 0.221 | 0.999 | 0.095 | 0.996 |
| 40,000 | 0.426 | 0.997 | 0.181 | 0.992 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 15,834, n_neg = 11,691, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.575**; CpG-island rows 0.091 (n_pos 1,175), non-CpG-island rows 0.556 (n_pos 14,659).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.374 | 0.399 [0.386, 0.411] | 0.399 [0.386, 0.410] | 0.731 [0.716, 0.745] | 0.684 [0.677, 0.690] |
| pydreg, −log10 p over all candidates | 0.374 | 0.399 [0.386, 0.411] | 0.399 [0.386, 0.410] | 0.731 [0.716, 0.745] | 0.684 [0.677, 0.690] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.438 | 0.461 [0.445, 0.474] | 0.456 [0.440, 0.469] | 0.759 [0.743, 0.773] | 0.716 [0.707, 0.723] |
| Δ pydreg(SVR) − PINTS | -0.064 | -0.062 [-0.070, -0.052] p<1e-03 | -0.057 [-0.065, -0.047] p<1e-03 | -0.028 [-0.031, -0.023] p<1e-03 | -0.032 [-0.036, -0.027] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.064 | -0.062 [-0.070, -0.052] p<1e-03 | -0.057 [-0.065, -0.047] p<1e-03 | -0.028 [-0.031, -0.023] p<1e-03 | -0.032 [-0.036, -0.027] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.399 [0.386, 0.411] | 0.399 [0.386, 0.410] | 0.731 [0.716, 0.745] | 0.684 [0.677, 0.690] |
| GC-matched: PINTS |  | 0.460 [0.444, 0.473] | 0.452 [0.437, 0.465] | 0.758 [0.743, 0.772] | 0.715 [0.707, 0.722] |
| GC-matched: Δ pydreg − PINTS |  | -0.061 [-0.069, -0.051] p<1e-03 | -0.053 [-0.062, -0.043] p<1e-03 | -0.027 [-0.031, -0.022] p<1e-03 | -0.031 [-0.035, -0.026] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.005 [-0.030, 0.021] p=0.672 | 0.012 [-0.018, 0.049] p=0.49 | 0.018 [-0.016, 0.058] p=0.386 | -0.001 [-0.014, 0.013] p=0.88 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.065 [-0.073, -0.056] p<1e-03 | -0.058 [-0.066, -0.049] p<1e-03 | -0.030 [-0.034, -0.026] p<1e-03 | -0.034 [-0.038, -0.029] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.062 | 1.000 | 0.023 | 1.000 |
| 20,000 | 0.164 | 1.000 | 0.050 | 0.999 |
| 40,000 | 0.338 | 1.000 | 0.103 | 0.999 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 27,173, n_neg = 62,494, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.303**; CpG-island rows 0.042 (n_pos 2,766), non-CpG-island rows 0.281 (n_pos 24,407).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.057 | 0.187 [0.179, 0.196] | 0.180 [0.174, 0.187] | 0.403 [0.389, 0.417] | 0.601 [0.597, 0.605] |
| pydreg, −log10 p over all candidates | 0.050 | 0.185 [0.177, 0.194] | 0.178 [0.171, 0.185] | 0.396 [0.382, 0.409] | 0.601 [0.597, 0.605] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.054 | 0.213 [0.200, 0.228] | 0.201 [0.192, 0.211] | 0.437 [0.421, 0.453] | 0.641 [0.636, 0.647] |
| Δ pydreg(SVR) − PINTS | 0.003 | -0.026 [-0.036, -0.016] p<1e-03 | -0.021 [-0.027, -0.015] p<1e-03 | -0.034 [-0.040, -0.029] p<1e-03 | -0.040 [-0.044, -0.037] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.005 | -0.028 [-0.038, -0.018] p<1e-03 | -0.023 [-0.029, -0.018] p<1e-03 | -0.042 [-0.048, -0.036] p<1e-03 | -0.041 [-0.044, -0.037] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.142 [0.134, 0.149] | 0.140 [0.133, 0.147] | 0.377 [0.363, 0.391] | 0.589 [0.585, 0.593] |
| GC-matched: PINTS |  | 0.159 [0.150, 0.169] | 0.156 [0.148, 0.165] | 0.408 [0.393, 0.424] | 0.626 [0.621, 0.631] |
| GC-matched: Δ pydreg − PINTS |  | -0.017 [-0.024, -0.010] p<1e-03 | -0.016 [-0.022, -0.011] p<1e-03 | -0.031 [-0.036, -0.026] p<1e-03 | -0.037 [-0.041, -0.034] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.012 [-0.040, 0.024] p=0.604 | 0.001 [-0.020, 0.022] p=0.898 | 0.005 [-0.007, 0.021] p=0.372 | -0.010 [-0.019, -0.001] p=0.03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.022 [-0.029, -0.014] p<1e-03 | -0.020 [-0.025, -0.014] p<1e-03 | -0.036 [-0.041, -0.031] p<1e-03 | -0.040 [-0.044, -0.037] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.043 | 0.736 | 0.029 | 0.714 |
| 20,000 | 0.114 | 0.657 | 0.061 | 0.700 |
| 40,000 | 0.246 | 0.605 | 0.190 | 0.649 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 1,175, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.019**; CpG-island rows 0.010 (n_pos 613), non-CpG-island rows 0.009 (n_pos 562).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.164 | 0.449 [0.416, 0.483] | 0.413 [0.388, 0.441] | 0.137 [0.121, 0.155] | 0.810 [0.794, 0.827] |
| pydreg, −log10 p over all candidates | 0.098 | 0.433 [0.403, 0.465] | 0.388 [0.363, 0.413] | 0.106 [0.094, 0.118] | 0.807 [0.791, 0.824] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.177 | 0.446 [0.415, 0.486] | 0.413 [0.388, 0.439] | 0.140 [0.122, 0.161] | 0.846 [0.833, 0.859] |
| Δ pydreg(SVR) − PINTS | -0.012 | 0.002 [-0.032, 0.033] p=0.99 | 0.001 [-0.022, 0.024] p=0.94 | -0.002 [-0.021, 0.013] p=0.808 | -0.036 [-0.050, -0.024] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.079 | -0.013 [-0.047, 0.019] p=0.386 | -0.025 [-0.047, -0.001] p=0.044 | -0.033 [-0.052, -0.020] p<1e-03 | -0.039 [-0.052, -0.027] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.305 [0.275, 0.340] | 0.281 [0.257, 0.307] | 0.080 [0.068, 0.093] | 0.768 [0.751, 0.785] |
| GC-matched: PINTS |  | 0.312 [0.274, 0.345] | 0.287 [0.262, 0.314] | 0.088 [0.077, 0.103] | 0.795 [0.781, 0.809] |
| GC-matched: Δ pydreg − PINTS |  | -0.008 [-0.039, 0.034] p=0.894 | -0.006 [-0.031, 0.017] p=0.6 | -0.008 [-0.023, 0.005] p=0.236 | -0.027 [-0.041, -0.014] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.012 [-0.040, 0.066] p=0.606 | 0.003 [-0.034, 0.039] p=0.838 | -0.009 [-0.031, 0.009] p=0.308 | -0.006 [-0.022, 0.012] p=0.496 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.029 [-0.066, 0.013] p=0.172 | -0.026 [-0.053, 0.004] p=0.088 | -0.002 [-0.008, 0.005] p=0.556 | -0.051 [-0.070, -0.032] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.228 | 0.221 | 0.151 | 0.265 |
| 20,000 | 0.426 | 0.152 | 0.249 | 0.199 |
| 40,000 | 0.652 | 0.097 | 0.393 | 0.163 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 1,175, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.019**; CpG-island rows 0.010 (n_pos 613), non-CpG-island rows 0.009 (n_pos 562).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.179 | 0.466 [0.437, 0.503] | 0.429 [0.404, 0.458] | 0.156 [0.137, 0.177] | 0.818 [0.803, 0.832] |
| pydreg, −log10 p over all candidates | 0.099 | 0.452 [0.419, 0.487] | 0.402 [0.379, 0.428] | 0.111 [0.099, 0.124] | 0.815 [0.799, 0.830] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.142 | 0.444 [0.417, 0.476] | 0.406 [0.384, 0.430] | 0.132 [0.117, 0.150] | 0.844 [0.831, 0.857] |
| Δ pydreg(SVR) − PINTS | 0.037 | 0.022 [-0.005, 0.053] p=0.116 | 0.023 [0.003, 0.044] p=0.012 | 0.024 [0.004, 0.040] p=0.014 | -0.026 [-0.039, -0.014] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.043 | 0.008 [-0.021, 0.036] p=0.646 | -0.004 [-0.023, 0.016] p=0.658 | -0.021 [-0.036, -0.010] p<1e-03 | -0.029 [-0.041, -0.017] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.303 [0.271, 0.334] | 0.274 [0.251, 0.299] | 0.080 [0.069, 0.097] | 0.769 [0.754, 0.784] |
| GC-matched: PINTS |  | 0.287 [0.261, 0.316] | 0.272 [0.250, 0.296] | 0.084 [0.073, 0.099] | 0.791 [0.776, 0.805] |
| GC-matched: Δ pydreg − PINTS |  | 0.016 [-0.014, 0.046] p=0.302 | 0.003 [-0.019, 0.024] p=0.802 | -0.004 [-0.019, 0.011] p=0.676 | -0.022 [-0.034, -0.010] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.018 [-0.028, 0.077] p=0.382 | 0.017 [-0.019, 0.053] p=0.37 | -0.001 [-0.024, 0.022] p=0.972 | -0.003 [-0.017, 0.012] p=0.77 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.018 [-0.054, 0.021] p=0.37 | -0.018 [-0.045, 0.007] p=0.166 | -0.004 [-0.009, 0.007] p=0.43 | -0.044 [-0.063, -0.024] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.248 | 0.238 | 0.140 | 0.219 |
| 20,000 | 0.442 | 0.158 | 0.268 | 0.182 |
| 40,000 | 0.675 | 0.098 | 0.444 | 0.143 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 1,175, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.019**; CpG-island rows 0.010 (n_pos 613), non-CpG-island rows 0.009 (n_pos 562).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.166 | 0.492 [0.460, 0.525] | 0.438 [0.412, 0.462] | 0.152 [0.135, 0.171] | 0.859 [0.848, 0.871] |
| pydreg, −log10 p over all candidates | 0.082 | 0.452 [0.419, 0.482] | 0.389 [0.366, 0.413] | 0.109 [0.097, 0.120] | 0.854 [0.843, 0.866] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.208 | 0.483 [0.450, 0.515] | 0.457 [0.431, 0.484] | 0.167 [0.148, 0.194] | 0.875 [0.864, 0.886] |
| Δ pydreg(SVR) − PINTS | -0.042 | 0.009 [-0.025, 0.041] p=0.638 | -0.019 [-0.043, 0.003] p=0.086 | -0.015 [-0.038, 0.004] p=0.136 | -0.016 [-0.025, -0.008] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.126 | -0.031 [-0.068, -0.002] p=0.038 | -0.068 [-0.091, -0.045] p<1e-03 | -0.058 [-0.081, -0.042] p<1e-03 | -0.021 [-0.030, -0.013] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.257 [0.228, 0.287] | 0.249 [0.224, 0.274] | 0.072 [0.063, 0.085] | 0.795 [0.783, 0.808] |
| GC-matched: PINTS |  | 0.343 [0.312, 0.370] | 0.316 [0.291, 0.342] | 0.101 [0.087, 0.119] | 0.822 [0.809, 0.834] |
| GC-matched: Δ pydreg − PINTS |  | -0.086 [-0.119, -0.050] p<1e-03 | -0.068 [-0.095, -0.040] p<1e-03 | -0.029 [-0.045, -0.015] p<1e-03 | -0.026 [-0.036, -0.017] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.105 [-0.160, -0.050] p<1e-03 | -0.087 [-0.128, -0.046] p<1e-03 | -0.029 [-0.053, -0.011] p<1e-03 | -0.012 [-0.024, 0.000] p=0.064 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.065 [-0.105, -0.018] p<1e-03 | -0.050 [-0.079, -0.021] p<1e-03 | -0.008 [-0.014, -0.002] p=0.012 | -0.042 [-0.057, -0.028] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.280 | 0.203 | 0.196 | 0.287 |
| 20,000 | 0.517 | 0.153 | 0.395 | 0.186 |
| 40,000 | 0.727 | 0.091 | 0.608 | 0.128 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 1,175, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.019**; CpG-island rows 0.010 (n_pos 613), non-CpG-island rows 0.009 (n_pos 562).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.204 | 0.470 [0.435, 0.503] | 0.439 [0.412, 0.466] | 0.164 [0.145, 0.189] | 0.828 [0.814, 0.842] |
| pydreg, −log10 p over all candidates | 0.121 | 0.464 [0.432, 0.499] | 0.416 [0.390, 0.443] | 0.123 [0.109, 0.138] | 0.826 [0.811, 0.840] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.161 | 0.494 [0.459, 0.527] | 0.434 [0.408, 0.460] | 0.141 [0.124, 0.164] | 0.864 [0.851, 0.876] |
| Δ pydreg(SVR) − PINTS | 0.044 | -0.024 [-0.057, 0.011] p=0.148 | 0.005 [-0.018, 0.030] p=0.644 | 0.024 [0.003, 0.042] p=0.014 | -0.036 [-0.049, -0.022] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.040 | -0.030 [-0.061, 0.006] p=0.108 | -0.018 [-0.042, 0.008] p=0.162 | -0.018 [-0.036, -0.004] p=0.014 | -0.039 [-0.052, -0.024] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.300 [0.269, 0.329] | 0.277 [0.254, 0.303] | 0.084 [0.072, 0.103] | 0.780 [0.765, 0.794] |
| GC-matched: PINTS |  | 0.272 [0.242, 0.302] | 0.264 [0.239, 0.290] | 0.081 [0.070, 0.096] | 0.810 [0.796, 0.824] |
| GC-matched: Δ pydreg − PINTS |  | 0.028 [-0.002, 0.060] p=0.076 | 0.013 [-0.011, 0.039] p=0.302 | 0.003 [-0.010, 0.018] p=0.608 | -0.031 [-0.044, -0.017] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.020 [-0.027, 0.077] p=0.394 | 0.034 [-0.003, 0.075] p=0.084 | 0.007 [-0.012, 0.035] p=0.37 | -0.006 [-0.021, 0.009] p=0.474 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.002 [-0.046, 0.047] p=0.948 | -0.016 [-0.046, 0.016] p=0.308 | -0.004 [-0.010, 0.001] p=0.104 | -0.058 [-0.078, -0.037] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.248 | 0.247 | 0.123 | 0.250 |
| 20,000 | 0.449 | 0.160 | 0.172 | 0.219 |
| 40,000 | 0.680 | 0.100 | 0.294 | 0.173 |


## dREG auxiliary labels: chromhmm_active vs TSS-silent DHSs

* **Positives:** ENCODE Broad ChromHMM segments in states 1 (active promoter) and 4–5 (strong enhancer), as in Danko 2015 Fig 1b (K562 43,377; GM12878 41,415; median segment about 1.4 kb).
* **Negatives:** DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).
* **Notes:** Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces dhs_unmarked as the specificity test.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 43,027, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.711**; CpG-island rows 0.280 (n_pos 6,802), non-CpG-island rows 0.674 (n_pos 36,225).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.583 | 0.728 [0.720, 0.735] | 0.691 [0.682, 0.700] | 0.918 [0.914, 0.921] | 0.853 [0.849, 0.857] |
| pydreg, −log10 p over all candidates | 0.589 | 0.728 [0.720, 0.735] | 0.686 [0.677, 0.695] | 0.917 [0.913, 0.920] | 0.853 [0.849, 0.856] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.462 | 0.714 [0.701, 0.727] | 0.671 [0.657, 0.685] | 0.952 [0.950, 0.955] | 0.908 [0.904, 0.911] |
| Δ pydreg(SVR) − PINTS | 0.121 | 0.014 [0.001, 0.025] p=0.034 | 0.019 [0.007, 0.031] p=0.004 | -0.034 [-0.037, -0.033] p<1e-03 | -0.055 [-0.057, -0.052] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.127 | 0.014 [0.001, 0.025] p=0.034 | 0.015 [0.003, 0.027] p=0.016 | -0.035 [-0.037, -0.033] p<1e-03 | -0.055 [-0.058, -0.052] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.728 [0.721, 0.736] | 0.694 [0.686, 0.702] | 0.918 [0.914, 0.921] | 0.854 [0.850, 0.857] |
| GC-matched: PINTS |  | 0.717 [0.704, 0.729] | 0.675 [0.661, 0.688] | 0.952 [0.950, 0.955] | 0.909 [0.905, 0.913] |
| GC-matched: Δ pydreg − PINTS |  | 0.012 [-0.001, 0.023] p=0.064 | 0.019 [0.007, 0.031] p=0.004 | -0.034 [-0.037, -0.033] p<1e-03 | -0.055 [-0.058, -0.052] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.027 [0.022, 0.035] p<1e-03 | 0.036 [0.028, 0.044] p<1e-03 | 0.005 [0.001, 0.009] p=0.002 | 0.002 [-0.000, 0.004] p=0.068 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.008 [-0.008, 0.022] p=0.354 | 0.016 [0.002, 0.031] p=0.034 | -0.045 [-0.048, -0.042] p<1e-03 | -0.066 [-0.069, -0.063] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.223 | 1.000 | 0.151 | 0.999 |
| 20,000 | 0.364 | 0.999 | 0.230 | 0.999 |
| 40,000 | 0.494 | 0.997 | 0.319 | 0.996 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 43,027, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.711**; CpG-island rows 0.280 (n_pos 6,802), non-CpG-island rows 0.674 (n_pos 36,225).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.486 | 0.507 [0.498, 0.516] | 0.506 [0.497, 0.515] | 0.850 [0.844, 0.856] | 0.741 [0.736, 0.745] |
| pydreg, −log10 p over all candidates | 0.486 | 0.507 [0.498, 0.516] | 0.506 [0.496, 0.515] | 0.850 [0.844, 0.856] | 0.740 [0.736, 0.745] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.510 | 0.585 [0.575, 0.594] | 0.565 [0.556, 0.574] | 0.873 [0.868, 0.878] | 0.779 [0.774, 0.784] |
| Δ pydreg(SVR) − PINTS | -0.024 | -0.077 [-0.081, -0.073] p<1e-03 | -0.059 [-0.063, -0.054] p<1e-03 | -0.023 [-0.024, -0.022] p<1e-03 | -0.039 [-0.041, -0.037] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.024 | -0.077 [-0.081, -0.073] p<1e-03 | -0.059 [-0.064, -0.054] p<1e-03 | -0.023 [-0.024, -0.022] p<1e-03 | -0.039 [-0.041, -0.037] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.507 [0.498, 0.516] | 0.506 [0.497, 0.515] | 0.850 [0.844, 0.856] | 0.740 [0.736, 0.745] |
| GC-matched: PINTS |  | 0.585 [0.575, 0.594] | 0.565 [0.555, 0.575] | 0.873 [0.868, 0.878] | 0.779 [0.774, 0.784] |
| GC-matched: Δ pydreg − PINTS |  | -0.077 [-0.081, -0.073] p<1e-03 | -0.059 [-0.064, -0.054] p<1e-03 | -0.023 [-0.024, -0.022] p<1e-03 | -0.039 [-0.041, -0.037] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.008 [0.001, 0.015] p=0.024 | 0.024 [0.016, 0.032] p<1e-03 | 0.012 [0.006, 0.018] p<1e-03 | 0.006 [0.002, 0.010] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.093 [-0.097, -0.089] p<1e-03 | -0.075 [-0.080, -0.069] p<1e-03 | -0.031 [-0.033, -0.030] p<1e-03 | -0.047 [-0.049, -0.045] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.218 | 1.000 | 0.126 | 0.999 |
| 20,000 | 0.327 | 1.000 | 0.210 | 0.999 |
| 40,000 | 0.425 | 0.999 | 0.297 | 0.998 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 43,027, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.711**; CpG-island rows 0.280 (n_pos 6,802), non-CpG-island rows 0.674 (n_pos 36,225).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.516 | 0.536 [0.525, 0.548] | 0.533 [0.522, 0.545] | 0.859 [0.853, 0.864] | 0.755 [0.750, 0.762] |
| pydreg, −log10 p over all candidates | 0.516 | 0.536 [0.525, 0.548] | 0.533 [0.521, 0.545] | 0.859 [0.853, 0.864] | 0.755 [0.750, 0.762] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.493 | 0.645 [0.633, 0.657] | 0.608 [0.597, 0.619] | 0.892 [0.887, 0.896] | 0.809 [0.803, 0.815] |
| Δ pydreg(SVR) − PINTS | 0.023 | -0.109 [-0.113, -0.105] p<1e-03 | -0.075 [-0.081, -0.069] p<1e-03 | -0.033 [-0.034, -0.031] p<1e-03 | -0.054 [-0.056, -0.052] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.023 | -0.109 [-0.113, -0.105] p<1e-03 | -0.075 [-0.082, -0.069] p<1e-03 | -0.033 [-0.034, -0.031] p<1e-03 | -0.054 [-0.056, -0.052] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.536 [0.525, 0.548] | 0.534 [0.523, 0.546] | 0.859 [0.853, 0.864] | 0.756 [0.750, 0.762] |
| GC-matched: PINTS |  | 0.646 [0.634, 0.658] | 0.611 [0.600, 0.622] | 0.892 [0.887, 0.896] | 0.810 [0.804, 0.816] |
| GC-matched: Δ pydreg − PINTS |  | -0.110 [-0.114, -0.105] p<1e-03 | -0.077 [-0.084, -0.071] p<1e-03 | -0.033 [-0.034, -0.031] p<1e-03 | -0.055 [-0.057, -0.052] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.000 [-0.005, 0.005] p=0.952 | 0.019 [0.012, 0.027] p<1e-03 | 0.006 [0.002, 0.011] p=0.004 | 0.002 [-0.001, 0.004] p=0.14 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.131 [-0.136, -0.125] p<1e-03 | -0.096 [-0.103, -0.089] p<1e-03 | -0.044 [-0.046, -0.041] p<1e-03 | -0.065 [-0.068, -0.063] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.217 | 1.000 | 0.166 | 1.000 |
| 20,000 | 0.327 | 1.000 | 0.222 | 0.999 |
| 40,000 | 0.438 | 0.999 | 0.328 | 0.998 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 40,690, n_neg = 16,074, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.717**; CpG-island rows 0.303 (n_pos 6,994), non-CpG-island rows 0.677 (n_pos 33,696).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.621 | 0.751 [0.744, 0.758] | 0.717 [0.708, 0.725] | 0.927 [0.924, 0.930] | 0.865 [0.862, 0.869] |
| pydreg, −log10 p over all candidates | 0.621 | 0.751 [0.744, 0.758] | 0.708 [0.699, 0.718] | 0.926 [0.923, 0.929] | 0.865 [0.861, 0.869] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.568 | 0.777 [0.767, 0.788] | 0.736 [0.726, 0.747] | 0.953 [0.951, 0.956] | 0.909 [0.905, 0.912] |
| Δ pydreg(SVR) − PINTS | 0.053 | -0.026 [-0.036, -0.017] p<1e-03 | -0.019 [-0.029, -0.011] p<1e-03 | -0.026 [-0.028, -0.025] p<1e-03 | -0.043 [-0.046, -0.041] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.053 | -0.026 [-0.036, -0.017] p<1e-03 | -0.027 [-0.038, -0.018] p<1e-03 | -0.028 [-0.029, -0.026] p<1e-03 | -0.044 [-0.047, -0.042] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.751 [0.744, 0.758] | 0.717 [0.708, 0.725] | 0.927 [0.924, 0.930] | 0.865 [0.862, 0.870] |
| GC-matched: PINTS |  | 0.777 [0.767, 0.788] | 0.736 [0.726, 0.747] | 0.953 [0.951, 0.956] | 0.909 [0.905, 0.912] |
| GC-matched: Δ pydreg − PINTS |  | -0.026 [-0.037, -0.017] p<1e-03 | -0.020 [-0.029, -0.011] p<1e-03 | -0.026 [-0.028, -0.025] p<1e-03 | -0.044 [-0.046, -0.041] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.013 [0.007, 0.018] p<1e-03 | 0.023 [0.017, 0.031] p<1e-03 | 0.004 [0.001, 0.006] p=0.014 | 0.001 [-0.001, 0.003] p=0.2 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.035 [-0.046, -0.025] p<1e-03 | -0.029 [-0.040, -0.020] p<1e-03 | -0.036 [-0.038, -0.034] p<1e-03 | -0.053 [-0.056, -0.050] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.227 | 1.000 | 0.173 | 1.000 |
| 20,000 | 0.373 | 1.000 | 0.265 | 0.999 |
| 40,000 | 0.514 | 0.998 | 0.365 | 0.997 |


## PINTS-on-TSS labels: pints_tss vs TSS-silent DHSs

* **Positives:** PINTS default bidirectional calls on the matched-cell TSS library: ENCODE PRO-cap (hg38; K562 lifted to hg19 for the K562 run-on libraries), or FANTOM5 CAGE for Jurkat (hg19).
* **Negatives:** DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).
* **Notes:** Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces dhs_unmarked as the specificity test. The negatives come from the same TSS library as the positives, so the two sets are disjoint by construction.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 33,933, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.659**; CpG-island rows 0.534 (n_pos 20,118), non-CpG-island rows 0.441 (n_pos 13,815).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.876 | 0.959 [0.956, 0.963] | 0.932 [0.926, 0.938] | 0.984 [0.983, 0.985] | 0.976 [0.974, 0.978] |
| pydreg, −log10 p over all candidates | 0.882 | 0.959 [0.956, 0.963] | 0.924 [0.916, 0.931] | 0.983 [0.982, 0.984] | 0.975 [0.973, 0.977] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.695 | 0.925 [0.916, 0.934] | 0.873 [0.862, 0.884] | 0.988 [0.986, 0.989] | 0.980 [0.978, 0.982] |
| Δ pydreg(SVR) − PINTS | 0.181 | 0.035 [0.026, 0.044] p<1e-03 | 0.059 [0.049, 0.068] p<1e-03 | -0.003 [-0.005, -0.002] p<1e-03 | -0.004 [-0.006, -0.002] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.186 | 0.035 [0.026, 0.044] p<1e-03 | 0.051 [0.040, 0.061] p<1e-03 | -0.005 [-0.006, -0.003] p<1e-03 | -0.005 [-0.007, -0.003] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.952 [0.942, 0.960] | 0.903 [0.883, 0.921] | 0.982 [0.980, 0.984] | 0.973 [0.970, 0.975] |
| GC-matched: PINTS |  | 0.924 [0.905, 0.938] | 0.874 [0.856, 0.891] | 0.988 [0.986, 0.989] | 0.980 [0.977, 0.982] |
| GC-matched: Δ pydreg − PINTS |  | 0.029 [0.016, 0.045] p<1e-03 | 0.029 [0.015, 0.042] p=0.002 | -0.005 [-0.007, -0.003] p<1e-03 | -0.007 [-0.009, -0.005] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.016 [0.004, 0.033] p=0.016 | 0.007 [-0.012, 0.025] p=0.496 | -0.004 [-0.007, -0.000] p=0.05 | -0.004 [-0.006, -0.001] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.052 [0.036, 0.064] p<1e-03 | 0.078 [0.064, 0.093] p<1e-03 | -0.005 [-0.009, -0.001] p=0.018 | -0.010 [-0.013, -0.006] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.314 | 1.000 | 0.238 | 0.999 |
| 20,000 | 0.529 | 0.999 | 0.366 | 0.999 |
| 40,000 | 0.754 | 0.997 | 0.505 | 0.997 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 33,933, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.659**; CpG-island rows 0.534 (n_pos 20,118), non-CpG-island rows 0.441 (n_pos 13,815).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.768 | 0.778 [0.771, 0.784] | 0.776 [0.770, 0.783] | 0.920 [0.917, 0.924] | 0.883 [0.879, 0.886] |
| pydreg, −log10 p over all candidates | 0.768 | 0.778 [0.771, 0.784] | 0.776 [0.769, 0.782] | 0.920 [0.916, 0.924] | 0.883 [0.879, 0.886] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.726 | 0.789 [0.782, 0.797] | 0.764 [0.756, 0.772] | 0.923 [0.919, 0.928] | 0.887 [0.883, 0.890] |
| Δ pydreg(SVR) − PINTS | 0.042 | -0.012 [-0.018, -0.005] p<1e-03 | 0.012 [0.006, 0.020] p<1e-03 | -0.003 [-0.005, -0.000] p=0.016 | -0.004 [-0.007, -0.001] p=0.024 |
| Δ pydreg(−log p) − PINTS | 0.042 | -0.012 [-0.018, -0.005] p<1e-03 | 0.012 [0.005, 0.019] p<1e-03 | -0.003 [-0.005, -0.001] p=0.01 | -0.004 [-0.007, -0.001] p=0.024 |
| GC-matched: pydreg (SVR) |  | 0.778 [0.771, 0.784] | 0.776 [0.770, 0.783] | 0.920 [0.917, 0.924] | 0.883 [0.879, 0.886] |
| GC-matched: PINTS |  | 0.790 [0.783, 0.797] | 0.769 [0.760, 0.778] | 0.924 [0.919, 0.928] | 0.887 [0.884, 0.891] |
| GC-matched: Δ pydreg − PINTS |  | -0.013 [-0.018, -0.006] p<1e-03 | 0.008 [0.000, 0.016] p=0.03 | -0.003 [-0.005, -0.001] p=0.01 | -0.005 [-0.008, -0.001] p=0.012 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.024 [0.018, 0.031] p<1e-03 | 0.040 [0.032, 0.049] p<1e-03 | 0.014 [0.011, 0.017] p<1e-03 | 0.014 [0.011, 0.018] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.066 [-0.074, -0.057] p<1e-03 | -0.038 [-0.046, -0.029] p<1e-03 | -0.033 [-0.037, -0.028] p<1e-03 | -0.032 [-0.036, -0.027] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.320 | 1.000 | 0.151 | 0.999 |
| 20,000 | 0.522 | 1.000 | 0.265 | 0.999 |
| 40,000 | 0.696 | 0.999 | 0.415 | 0.998 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 33,933, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.659**; CpG-island rows 0.534 (n_pos 20,118), non-CpG-island rows 0.441 (n_pos 13,815).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.890 | 0.894 [0.878, 0.905] | 0.891 [0.875, 0.902] | 0.962 [0.956, 0.966] | 0.944 [0.935, 0.950] |
| pydreg, −log10 p over all candidates | 0.890 | 0.894 [0.878, 0.905] | 0.890 [0.874, 0.901] | 0.962 [0.956, 0.966] | 0.944 [0.935, 0.950] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.794 | 0.917 [0.899, 0.929] | 0.883 [0.866, 0.895] | 0.968 [0.963, 0.972] | 0.953 [0.944, 0.959] |
| Δ pydreg(SVR) − PINTS | 0.096 | -0.023 [-0.027, -0.019] p<1e-03 | 0.008 [0.002, 0.014] p=0.006 | -0.006 [-0.008, -0.005] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.096 | -0.023 [-0.027, -0.019] p<1e-03 | 0.007 [0.001, 0.014] p=0.024 | -0.007 [-0.008, -0.005] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.893 [0.877, 0.904] | 0.879 [0.860, 0.892] | 0.961 [0.956, 0.965] | 0.942 [0.933, 0.948] |
| GC-matched: PINTS |  | 0.916 [0.899, 0.928] | 0.868 [0.846, 0.886] | 0.967 [0.961, 0.971] | 0.951 [0.942, 0.957] |
| GC-matched: Δ pydreg − PINTS |  | -0.023 [-0.027, -0.019] p<1e-03 | 0.011 [-0.003, 0.026] p=0.102 | -0.006 [-0.008, -0.003] p<1e-03 | -0.009 [-0.011, -0.006] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.001 [-0.005, 0.006] p=0.732 | 0.026 [0.009, 0.047] p<1e-03 | 0.003 [0.000, 0.008] p=0.044 | 0.002 [-0.001, 0.005] p=0.1 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.055 [-0.063, -0.047] p<1e-03 | -0.007 [-0.017, 0.004] p=0.214 | -0.023 [-0.027, -0.018] p<1e-03 | -0.024 [-0.028, -0.020] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.345 | 1.000 | 0.244 | 1.000 |
| 20,000 | 0.590 | 1.000 | 0.351 | 0.999 |
| 40,000 | 0.817 | 0.999 | 0.567 | 0.999 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 46,284, n_neg = 11,691, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.798**; CpG-island rows 0.682 (n_pos 25,078), non-CpG-island rows 0.645 (n_pos 21,206).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.776 | 0.785 [0.780, 0.789] | 0.785 [0.780, 0.789] | 0.954 [0.951, 0.957] | 0.887 [0.884, 0.889] |
| pydreg, −log10 p over all candidates | 0.776 | 0.785 [0.780, 0.789] | 0.785 [0.780, 0.789] | 0.954 [0.951, 0.957] | 0.887 [0.884, 0.889] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.822 | 0.829 [0.822, 0.835] | 0.824 [0.817, 0.831] | 0.964 [0.960, 0.967] | 0.910 [0.906, 0.913] |
| Δ pydreg(SVR) − PINTS | -0.046 | -0.044 [-0.050, -0.038] p<1e-03 | -0.039 [-0.046, -0.033] p<1e-03 | -0.009 [-0.010, -0.008] p<1e-03 | -0.023 [-0.026, -0.019] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.046 | -0.044 [-0.050, -0.038] p<1e-03 | -0.039 [-0.046, -0.033] p<1e-03 | -0.009 [-0.010, -0.008] p<1e-03 | -0.023 [-0.026, -0.019] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.785 [0.780, 0.789] | 0.785 [0.780, 0.789] | 0.954 [0.951, 0.957] | 0.887 [0.884, 0.889] |
| GC-matched: PINTS |  | 0.828 [0.821, 0.835] | 0.814 [0.791, 0.827] | 0.963 [0.960, 0.967] | 0.908 [0.904, 0.912] |
| GC-matched: Δ pydreg − PINTS |  | -0.043 [-0.049, -0.036] p<1e-03 | -0.029 [-0.042, -0.007] p=0.014 | -0.009 [-0.010, -0.008] p<1e-03 | -0.021 [-0.025, -0.017] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.006 [-0.011, -0.000] p=0.04 | 0.005 [-0.007, 0.029] p=0.65 | -0.002 [-0.003, 0.000] p=0.108 | -0.002 [-0.005, 0.001] p=0.264 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.087 [-0.096, -0.077] p<1e-03 | -0.072 [-0.084, -0.056] p<1e-03 | -0.032 [-0.035, -0.028] p<1e-03 | -0.044 [-0.049, -0.039] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.317 | 1.000 | 0.151 | 1.000 |
| 20,000 | 0.538 | 1.000 | 0.260 | 1.000 |
| 40,000 | 0.752 | 1.000 | 0.419 | 1.000 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 17,548, n_neg = 62,494, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.219**; CpG-island rows 0.173 (n_pos 13,033), non-CpG-island rows 0.067 (n_pos 4,515).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.501 | 0.704 [0.686, 0.720] | 0.659 [0.643, 0.674] | 0.755 [0.744, 0.766] | 0.857 [0.849, 0.865] |
| pydreg, −log10 p over all candidates | 0.404 | 0.699 [0.682, 0.716] | 0.629 [0.612, 0.644] | 0.700 [0.686, 0.713] | 0.854 [0.846, 0.862] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.396 | 0.680 [0.656, 0.703] | 0.625 [0.610, 0.639] | 0.732 [0.719, 0.744] | 0.858 [0.849, 0.865] |
| Δ pydreg(SVR) − PINTS | 0.105 | 0.024 [0.006, 0.042] p=0.004 | 0.034 [0.024, 0.044] p<1e-03 | 0.023 [0.014, 0.032] p<1e-03 | -0.000 [-0.004, 0.004] p=0.866 |
| Δ pydreg(−log p) − PINTS | 0.009 | 0.019 [0.002, 0.039] p=0.024 | 0.004 [-0.006, 0.014] p=0.488 | -0.032 [-0.041, -0.022] p<1e-03 | -0.003 [-0.007, 0.001] p=0.118 |
| GC-matched: pydreg (SVR) |  | 0.435 [0.364, 0.495] | 0.408 [0.369, 0.452] | 0.607 [0.575, 0.638] | 0.806 [0.795, 0.816] |
| GC-matched: PINTS |  | 0.458 [0.394, 0.525] | 0.412 [0.374, 0.451] | 0.588 [0.552, 0.624] | 0.804 [0.793, 0.814] |
| GC-matched: Δ pydreg − PINTS |  | -0.023 [-0.088, 0.045] p=0.46 | -0.004 [-0.041, 0.037] p=0.824 | 0.019 [-0.016, 0.056] p=0.294 | 0.002 [-0.004, 0.009] p=0.518 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.029 [-0.115, 0.049] p=0.496 | -0.010 [-0.064, 0.051] p=0.716 | 0.036 [-0.018, 0.089] p=0.222 | 0.012 [0.004, 0.020] p=0.006 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.002 [-0.019, 0.012] p=0.706 | -0.010 [-0.023, 0.003] p=0.144 | -0.008 [-0.021, 0.005] p=0.21 | -0.028 [-0.036, -0.021] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.443 | 0.949 | 0.286 | 0.940 |
| 20,000 | 0.628 | 0.872 | 0.435 | 0.915 |
| 40,000 | 0.740 | 0.748 | 0.646 | 0.803 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 4,706, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.070**; CpG-island rows 0.060 (n_pos 3,939), non-CpG-island rows 0.012 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.563 | 0.832 [0.819, 0.847] | 0.765 [0.751, 0.779] | 0.702 [0.683, 0.717] | 0.923 [0.917, 0.930] |
| pydreg, −log10 p over all candidates | 0.296 | 0.822 [0.807, 0.837] | 0.711 [0.697, 0.727] | 0.598 [0.582, 0.614] | 0.918 [0.912, 0.925] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.633 | 0.852 [0.836, 0.867] | 0.792 [0.778, 0.806] | 0.748 [0.722, 0.770] | 0.940 [0.934, 0.945] |
| Δ pydreg(SVR) − PINTS | -0.070 | -0.020 [-0.031, -0.007] p<1e-03 | -0.027 [-0.038, -0.015] p<1e-03 | -0.046 [-0.070, -0.017] p=0.006 | -0.016 [-0.021, -0.012] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.337 | -0.030 [-0.043, -0.018] p<1e-03 | -0.081 [-0.093, -0.069] p<1e-03 | -0.150 [-0.173, -0.122] p<1e-03 | -0.022 [-0.026, -0.017] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.741 [0.710, 0.762] | 0.651 [0.628, 0.675] | 0.566 [0.530, 0.602] | 0.902 [0.895, 0.909] |
| GC-matched: PINTS |  | 0.775 [0.747, 0.796] | 0.705 [0.686, 0.725] | 0.670 [0.639, 0.700] | 0.920 [0.914, 0.927] |
| GC-matched: Δ pydreg − PINTS |  | -0.033 [-0.060, -0.007] p=0.01 | -0.054 [-0.077, -0.032] p<1e-03 | -0.104 [-0.147, -0.060] p<1e-03 | -0.019 [-0.024, -0.013] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.031 [-0.065, 0.008] p=0.142 | -0.055 [-0.084, -0.027] p<1e-03 | -0.111 [-0.165, -0.056] p<1e-03 | -0.009 [-0.013, -0.004] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.088 [-0.119, -0.053] p<1e-03 | -0.087 [-0.113, -0.061] p<1e-03 | -0.069 [-0.098, -0.046] p<1e-03 | -0.078 [-0.096, -0.058] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.673 | 0.770 | 0.574 | 0.846 |
| 20,000 | 0.825 | 0.581 | 0.716 | 0.741 |
| 40,000 | 0.878 | 0.366 | 0.819 | 0.619 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 4,706, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.070**; CpG-island rows 0.060 (n_pos 3,939), non-CpG-island rows 0.012 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.568 | 0.838 [0.823, 0.852] | 0.768 [0.753, 0.783] | 0.709 [0.690, 0.726] | 0.925 [0.919, 0.932] |
| pydreg, −log10 p over all candidates | 0.288 | 0.822 [0.806, 0.839] | 0.710 [0.695, 0.724] | 0.595 [0.579, 0.610] | 0.920 [0.914, 0.926] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.564 | 0.854 [0.842, 0.867] | 0.781 [0.767, 0.796] | 0.731 [0.709, 0.752] | 0.936 [0.930, 0.942] |
| Δ pydreg(SVR) − PINTS | 0.004 | -0.016 [-0.027, -0.007] p=0.004 | -0.013 [-0.024, -0.001] p=0.036 | -0.022 [-0.045, 0.004] p=0.096 | -0.011 [-0.014, -0.007] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.277 | -0.032 [-0.043, -0.020] p<1e-03 | -0.071 [-0.082, -0.060] p<1e-03 | -0.136 [-0.155, -0.114] p<1e-03 | -0.016 [-0.020, -0.013] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.705 [0.674, 0.731] | 0.615 [0.587, 0.642] | 0.528 [0.491, 0.562] | 0.898 [0.891, 0.906] |
| GC-matched: PINTS |  | 0.758 [0.732, 0.775] | 0.690 [0.669, 0.712] | 0.656 [0.622, 0.689] | 0.916 [0.909, 0.922] |
| GC-matched: Δ pydreg − PINTS |  | -0.052 [-0.082, -0.023] p<1e-03 | -0.075 [-0.106, -0.047] p<1e-03 | -0.128 [-0.173, -0.080] p<1e-03 | -0.018 [-0.023, -0.013] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.053 [-0.092, -0.019] p=0.002 | -0.082 [-0.120, -0.046] p<1e-03 | -0.142 [-0.200, -0.080] p<1e-03 | -0.010 [-0.015, -0.005] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.075 [-0.109, -0.043] p<1e-03 | -0.079 [-0.107, -0.051] p<1e-03 | -0.071 [-0.103, -0.042] p<1e-03 | -0.061 [-0.078, -0.042] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.669 | 0.771 | 0.556 | 0.817 |
| 20,000 | 0.828 | 0.584 | 0.747 | 0.713 |
| 40,000 | 0.883 | 0.361 | 0.854 | 0.562 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 4,706, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.070**; CpG-island rows 0.060 (n_pos 3,939), non-CpG-island rows 0.012 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.319 | 0.821 [0.806, 0.835] | 0.712 [0.696, 0.726] | 0.601 [0.579, 0.620] | 0.918 [0.911, 0.924] |
| pydreg, −log10 p over all candidates | 0.159 | 0.778 [0.756, 0.795] | 0.628 [0.612, 0.642] | 0.494 [0.476, 0.510] | 0.909 [0.903, 0.915] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.690 | 0.860 [0.847, 0.872] | 0.814 [0.800, 0.826] | 0.779 [0.749, 0.802] | 0.941 [0.935, 0.947] |
| Δ pydreg(SVR) − PINTS | -0.371 | -0.038 [-0.050, -0.029] p<1e-03 | -0.102 [-0.114, -0.090] p<1e-03 | -0.178 [-0.206, -0.148] p<1e-03 | -0.023 [-0.027, -0.020] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.531 | -0.081 [-0.101, -0.067] p<1e-03 | -0.186 [-0.200, -0.171] p<1e-03 | -0.285 [-0.313, -0.254] p<1e-03 | -0.032 [-0.035, -0.028] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.478 [0.407, 0.520] | 0.441 [0.406, 0.474] | 0.362 [0.335, 0.392] | 0.875 [0.867, 0.883] |
| GC-matched: PINTS |  | 0.791 [0.777, 0.805] | 0.731 [0.714, 0.748] | 0.704 [0.677, 0.728] | 0.922 [0.915, 0.928] |
| GC-matched: Δ pydreg − PINTS |  | -0.314 [-0.384, -0.272] p<1e-03 | -0.289 [-0.324, -0.257] p<1e-03 | -0.342 [-0.374, -0.305] p<1e-03 | -0.047 [-0.052, -0.042] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.405 [-0.465, -0.328] p<1e-03 | -0.341 [-0.386, -0.302] p<1e-03 | -0.389 [-0.429, -0.345] p<1e-03 | -0.041 [-0.047, -0.036] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.082 [-0.131, -0.050] p<1e-03 | -0.116 [-0.149, -0.085] p<1e-03 | -0.107 [-0.139, -0.080] p<1e-03 | -0.088 [-0.106, -0.069] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.593 | 0.684 | 0.671 | 0.847 |
| 20,000 | 0.830 | 0.538 | 0.834 | 0.660 |
| 40,000 | 0.890 | 0.328 | 0.889 | 0.462 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 4,706, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.070**; CpG-island rows 0.060 (n_pos 3,939), non-CpG-island rows 0.012 (n_pos 767).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.597 | 0.826 [0.811, 0.841] | 0.766 [0.750, 0.781] | 0.716 [0.699, 0.733] | 0.922 [0.916, 0.929] |
| pydreg, −log10 p over all candidates | 0.317 | 0.817 [0.800, 0.832] | 0.715 [0.699, 0.730] | 0.607 [0.590, 0.622] | 0.917 [0.911, 0.924] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.579 | 0.871 [0.858, 0.883] | 0.795 [0.781, 0.809] | 0.746 [0.719, 0.769] | 0.940 [0.934, 0.946] |
| Δ pydreg(SVR) − PINTS | 0.018 | -0.045 [-0.057, -0.032] p<1e-03 | -0.030 [-0.042, -0.017] p<1e-03 | -0.030 [-0.054, -0.001] p=0.048 | -0.018 [-0.022, -0.013] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.262 | -0.055 [-0.068, -0.041] p<1e-03 | -0.080 [-0.093, -0.066] p<1e-03 | -0.139 [-0.163, -0.112] p<1e-03 | -0.023 [-0.027, -0.019] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.688 [0.658, 0.717] | 0.606 [0.580, 0.632] | 0.535 [0.502, 0.571] | 0.894 [0.887, 0.901] |
| GC-matched: PINTS |  | 0.749 [0.732, 0.765] | 0.668 [0.648, 0.690] | 0.633 [0.609, 0.662] | 0.918 [0.911, 0.924] |
| GC-matched: Δ pydreg − PINTS |  | -0.061 [-0.089, -0.033] p<1e-03 | -0.062 [-0.089, -0.036] p<1e-03 | -0.098 [-0.133, -0.061] p<1e-03 | -0.024 [-0.029, -0.018] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.061 [-0.102, -0.023] p<1e-03 | -0.057 [-0.091, -0.027] p<1e-03 | -0.098 [-0.143, -0.053] p<1e-03 | -0.013 [-0.018, -0.008] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.116 [-0.155, -0.078] p<1e-03 | -0.115 [-0.147, -0.084] p<1e-03 | -0.079 [-0.112, -0.050] p<1e-03 | -0.082 [-0.101, -0.062] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.667 | 0.780 | 0.507 | 0.847 |
| 20,000 | 0.814 | 0.581 | 0.592 | 0.795 |
| 40,000 | 0.874 | 0.363 | 0.767 | 0.686 |


## PINTS-on-TSS labels: pints_tss_dhs vs TSS-silent DHSs

* **Positives:** pints_tss overlapping the ENCODE DNase peaks of that cell: the dREG-paper recipe (DHS ∩ TSS-assay transcription), with PINTS-on-PRO-cap in place of the GRO-cap HMM.
* **Negatives:** DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).
* **Notes:** Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces dhs_unmarked as the specificity test. The negatives come from the same TSS library as the positives, so the two sets are disjoint by construction.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 29,601, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.628**; CpG-island rows 0.503 (n_pos 17,769), non-CpG-island rows 0.403 (n_pos 11,832).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.888 | 0.967 [0.963, 0.970] | 0.940 [0.934, 0.946] | 0.986 [0.984, 0.987] | 0.980 [0.978, 0.982] |
| pydreg, −log10 p over all candidates | 0.894 | 0.967 [0.963, 0.970] | 0.932 [0.925, 0.939] | 0.984 [0.983, 0.986] | 0.979 [0.977, 0.981] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.713 | 0.935 [0.927, 0.943] | 0.884 [0.874, 0.894] | 0.988 [0.986, 0.989] | 0.982 [0.981, 0.984] |
| Δ pydreg(SVR) − PINTS | 0.175 | 0.032 [0.024, 0.039] p<1e-03 | 0.057 [0.047, 0.065] p<1e-03 | -0.002 [-0.004, -0.001] p<1e-03 | -0.003 [-0.004, -0.001] p=0.002 |
| Δ pydreg(−log p) − PINTS | 0.181 | 0.032 [0.024, 0.039] p<1e-03 | 0.049 [0.039, 0.058] p<1e-03 | -0.003 [-0.005, -0.002] p<1e-03 | -0.003 [-0.005, -0.002] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.960 [0.951, 0.967] | 0.912 [0.893, 0.931] | 0.983 [0.982, 0.985] | 0.977 [0.974, 0.979] |
| GC-matched: PINTS |  | 0.934 [0.918, 0.946] | 0.884 [0.866, 0.901] | 0.988 [0.985, 0.989] | 0.982 [0.980, 0.984] |
| GC-matched: Δ pydreg − PINTS |  | 0.026 [0.014, 0.040] p<1e-03 | 0.028 [0.012, 0.044] p<1e-03 | -0.004 [-0.006, -0.002] p<1e-03 | -0.005 [-0.007, -0.003] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.011 [0.001, 0.027] p=0.034 | 0.005 [-0.014, 0.025] p=0.598 | -0.003 [-0.007, 0.001] p=0.15 | -0.003 [-0.006, -0.001] p=0.006 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.050 [0.037, 0.062] p<1e-03 | 0.080 [0.065, 0.094] p<1e-03 | 0.000 [-0.004, 0.004] p=0.998 | -0.006 [-0.009, -0.004] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.324 | 1.000 | 0.245 | 0.999 |
| 20,000 | 0.550 | 0.999 | 0.376 | 0.999 |
| 40,000 | 0.774 | 0.997 | 0.521 | 0.997 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 29,601, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.628**; CpG-island rows 0.503 (n_pos 17,769), non-CpG-island rows 0.403 (n_pos 11,832).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.795 | 0.803 [0.797, 0.809] | 0.801 [0.795, 0.808] | 0.923 [0.919, 0.927] | 0.896 [0.893, 0.899] |
| pydreg, −log10 p over all candidates | 0.795 | 0.803 [0.797, 0.809] | 0.801 [0.795, 0.807] | 0.923 [0.919, 0.927] | 0.896 [0.893, 0.899] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.757 | 0.816 [0.809, 0.822] | 0.791 [0.783, 0.799] | 0.926 [0.922, 0.931] | 0.901 [0.897, 0.904] |
| Δ pydreg(SVR) − PINTS | 0.038 | -0.013 [-0.018, -0.007] p<1e-03 | 0.011 [0.004, 0.017] p<1e-03 | -0.003 [-0.005, -0.001] p=0.002 | -0.004 [-0.007, -0.002] p=0.002 |
| Δ pydreg(−log p) − PINTS | 0.038 | -0.013 [-0.018, -0.007] p<1e-03 | 0.010 [0.004, 0.016] p=0.002 | -0.003 [-0.005, -0.001] p<1e-03 | -0.004 [-0.007, -0.002] p=0.002 |
| GC-matched: pydreg (SVR) |  | 0.803 [0.797, 0.809] | 0.802 [0.796, 0.808] | 0.923 [0.919, 0.927] | 0.896 [0.893, 0.899] |
| GC-matched: PINTS |  | 0.817 [0.810, 0.823] | 0.797 [0.789, 0.805] | 0.927 [0.922, 0.931] | 0.901 [0.898, 0.905] |
| GC-matched: Δ pydreg − PINTS |  | -0.014 [-0.019, -0.008] p<1e-03 | 0.005 [-0.002, 0.013] p=0.152 | -0.004 [-0.006, -0.002] p<1e-03 | -0.005 [-0.008, -0.002] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.020 [0.014, 0.026] p<1e-03 | 0.034 [0.026, 0.042] p<1e-03 | 0.012 [0.009, 0.016] p<1e-03 | 0.012 [0.009, 0.015] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.063 [-0.072, -0.055] p<1e-03 | -0.036 [-0.045, -0.028] p<1e-03 | -0.033 [-0.038, -0.028] p<1e-03 | -0.031 [-0.035, -0.026] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.336 | 1.000 | 0.166 | 0.999 |
| 20,000 | 0.548 | 1.000 | 0.289 | 0.999 |
| 40,000 | 0.723 | 0.999 | 0.448 | 0.998 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 29,601, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.628**; CpG-island rows 0.503 (n_pos 17,769), non-CpG-island rows 0.403 (n_pos 11,832).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.903 | 0.907 [0.892, 0.918] | 0.903 [0.888, 0.914] | 0.963 [0.958, 0.968] | 0.951 [0.943, 0.956] |
| pydreg, −log10 p over all candidates | 0.903 | 0.907 [0.892, 0.918] | 0.902 [0.888, 0.913] | 0.963 [0.957, 0.967] | 0.950 [0.943, 0.956] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.807 | 0.922 [0.906, 0.933] | 0.890 [0.874, 0.902] | 0.967 [0.962, 0.972] | 0.956 [0.948, 0.962] |
| Δ pydreg(SVR) − PINTS | 0.095 | -0.016 [-0.020, -0.012] p<1e-03 | 0.014 [0.008, 0.020] p<1e-03 | -0.004 [-0.006, -0.003] p<1e-03 | -0.005 [-0.007, -0.003] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.095 | -0.016 [-0.020, -0.012] p<1e-03 | 0.013 [0.007, 0.019] p<1e-03 | -0.004 [-0.006, -0.003] p<1e-03 | -0.005 [-0.007, -0.003] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.906 [0.891, 0.917] | 0.891 [0.875, 0.906] | 0.963 [0.957, 0.967] | 0.949 [0.941, 0.955] |
| GC-matched: PINTS |  | 0.921 [0.905, 0.932] | 0.875 [0.855, 0.894] | 0.966 [0.960, 0.971] | 0.954 [0.946, 0.960] |
| GC-matched: Δ pydreg − PINTS |  | -0.016 [-0.020, -0.011] p<1e-03 | 0.016 [0.001, 0.030] p=0.026 | -0.003 [-0.005, -0.001] p<1e-03 | -0.005 [-0.008, -0.003] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.004 [0.000, 0.010] p=0.038 | 0.027 [0.010, 0.044] p=0.002 | 0.006 [0.002, 0.010] p<1e-03 | 0.004 [0.002, 0.007] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.044 [-0.053, -0.036] p<1e-03 | 0.003 [-0.007, 0.013] p=0.668 | -0.017 [-0.022, -0.012] p<1e-03 | -0.019 [-0.023, -0.014] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.356 | 1.000 | 0.263 | 1.000 |
| 20,000 | 0.609 | 1.000 | 0.372 | 0.999 |
| 40,000 | 0.832 | 0.999 | 0.592 | 0.999 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 24,837, n_neg = 11,691, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.680**; CpG-island rows 0.547 (n_pos 14,094), non-CpG-island rows 0.479 (n_pos 10,743).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.829 | 0.836 [0.830, 0.841] | 0.835 [0.829, 0.841] | 0.945 [0.940, 0.949] | 0.913 [0.910, 0.917] |
| pydreg, −log10 p over all candidates | 0.829 | 0.836 [0.830, 0.841] | 0.835 [0.829, 0.841] | 0.945 [0.940, 0.949] | 0.913 [0.910, 0.917] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.814 | 0.822 [0.814, 0.829] | 0.817 [0.810, 0.825] | 0.940 [0.935, 0.945] | 0.906 [0.902, 0.910] |
| Δ pydreg(SVR) − PINTS | 0.015 | 0.014 [0.008, 0.020] p<1e-03 | 0.018 [0.012, 0.024] p<1e-03 | 0.005 [0.003, 0.007] p<1e-03 | 0.008 [0.005, 0.011] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.015 | 0.014 [0.008, 0.020] p<1e-03 | 0.018 [0.012, 0.024] p<1e-03 | 0.005 [0.003, 0.007] p<1e-03 | 0.008 [0.005, 0.011] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.836 [0.830, 0.841] | 0.835 [0.829, 0.841] | 0.945 [0.940, 0.949] | 0.913 [0.910, 0.917] |
| GC-matched: PINTS |  | 0.820 [0.813, 0.828] | 0.808 [0.788, 0.821] | 0.939 [0.934, 0.944] | 0.904 [0.900, 0.908] |
| GC-matched: Δ pydreg − PINTS |  | 0.015 [0.010, 0.021] p<1e-03 | 0.027 [0.015, 0.047] p<1e-03 | 0.005 [0.003, 0.008] p<1e-03 | 0.009 [0.006, 0.013] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.032 [0.027, 0.038] p<1e-03 | 0.041 [0.030, 0.058] p<1e-03 | 0.016 [0.013, 0.019] p<1e-03 | 0.018 [0.015, 0.021] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.008 [-0.018, 0.003] p=0.14 | 0.006 [-0.007, 0.020] p=0.402 | -0.003 [-0.008, 0.003] p=0.322 | -0.003 [-0.008, 0.003] p=0.324 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.389 | 1.000 | 0.199 | 1.000 |
| 20,000 | 0.613 | 1.000 | 0.320 | 1.000 |
| 40,000 | 0.807 | 1.000 | 0.476 | 1.000 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 13,470, n_neg = 62,494, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.177**; CpG-island rows 0.151 (n_pos 11,145), non-CpG-island rows 0.036 (n_pos 2,325).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.597 | 0.821 [0.804, 0.835] | 0.769 [0.753, 0.783] | 0.831 [0.818, 0.843] | 0.920 [0.912, 0.927] |
| pydreg, −log10 p over all candidates | 0.485 | 0.816 [0.800, 0.831] | 0.733 [0.716, 0.749] | 0.760 [0.744, 0.775] | 0.916 [0.909, 0.923] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.483 | 0.787 [0.762, 0.810] | 0.724 [0.708, 0.738] | 0.790 [0.774, 0.803] | 0.910 [0.902, 0.917] |
| Δ pydreg(SVR) − PINTS | 0.114 | 0.033 [0.017, 0.051] p<1e-03 | 0.045 [0.036, 0.054] p<1e-03 | 0.041 [0.032, 0.050] p<1e-03 | 0.010 [0.006, 0.013] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.002 | 0.028 [0.013, 0.047] p<1e-03 | 0.009 [-0.002, 0.019] p=0.092 | -0.029 [-0.040, -0.019] p<1e-03 | 0.006 [0.002, 0.009] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.501 [0.418, 0.575] | 0.468 [0.419, 0.521] | 0.639 [0.599, 0.682] | 0.870 [0.859, 0.882] |
| GC-matched: PINTS |  | 0.545 [0.480, 0.629] | 0.491 [0.445, 0.533] | 0.621 [0.583, 0.664] | 0.862 [0.851, 0.872] |
| GC-matched: Δ pydreg − PINTS |  | -0.043 [-0.118, 0.037] p=0.236 | -0.022 [-0.066, 0.025] p=0.354 | 0.018 [-0.026, 0.060] p=0.426 | 0.008 [0.002, 0.015] p=0.014 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.039 [-0.131, 0.042] p=0.302 | -0.025 [-0.080, 0.036] p=0.408 | 0.027 [-0.030, 0.082] p=0.346 | 0.014 [0.007, 0.022] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.012 [-0.036, 0.016] p=0.432 | -0.012 [-0.031, 0.008] p=0.224 | -0.012 [-0.031, 0.005] p=0.198 | -0.024 [-0.033, -0.014] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.531 | 0.944 | 0.353 | 0.937 |
| 20,000 | 0.739 | 0.861 | 0.529 | 0.910 |
| 40,000 | 0.857 | 0.726 | 0.754 | 0.785 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 4,154, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.063**; CpG-island rows 0.058 (n_pos 3,803), non-CpG-island rows 0.006 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.625 | 0.920 [0.908, 0.933] | 0.846 [0.832, 0.860] | 0.764 [0.747, 0.783] | 0.974 [0.970, 0.978] |
| pydreg, −log10 p over all candidates | 0.330 | 0.909 [0.894, 0.922] | 0.787 [0.772, 0.801] | 0.648 [0.632, 0.664] | 0.968 [0.964, 0.972] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.702 | 0.933 [0.919, 0.945] | 0.867 [0.854, 0.879] | 0.802 [0.776, 0.827] | 0.980 [0.978, 0.983] |
| Δ pydreg(SVR) − PINTS | -0.077 | -0.012 [-0.023, -0.000] p=0.046 | -0.021 [-0.032, -0.010] p=0.002 | -0.039 [-0.064, -0.011] p=0.008 | -0.006 [-0.010, -0.003] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.372 | -0.024 [-0.036, -0.012] p<1e-03 | -0.080 [-0.092, -0.068] p<1e-03 | -0.155 [-0.178, -0.128] p<1e-03 | -0.012 [-0.016, -0.009] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.803 [0.768, 0.836] | 0.710 [0.679, 0.738] | 0.600 [0.561, 0.649] | 0.955 [0.949, 0.960] |
| GC-matched: PINTS |  | 0.834 [0.805, 0.872] | 0.768 [0.746, 0.791] | 0.714 [0.679, 0.749] | 0.965 [0.962, 0.969] |
| GC-matched: Δ pydreg − PINTS |  | -0.032 [-0.069, 0.003] p=0.078 | -0.059 [-0.085, -0.030] p<1e-03 | -0.113 [-0.161, -0.056] p<1e-03 | -0.011 [-0.015, -0.006] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.032 [-0.066, 0.008] p=0.116 | -0.058 [-0.087, -0.028] p<1e-03 | -0.116 [-0.171, -0.052] p<1e-03 | -0.009 [-0.014, -0.004] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.066 [-0.113, 0.004] p=0.072 | -0.068 [-0.109, -0.023] p=0.004 | -0.101 [-0.145, -0.055] p<1e-03 | -0.031 [-0.051, -0.011] p=0.004 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.747 | 0.767 | 0.636 | 0.844 |
| 20,000 | 0.912 | 0.576 | 0.791 | 0.736 |
| 40,000 | 0.967 | 0.360 | 0.900 | 0.611 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 4,154, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.063**; CpG-island rows 0.058 (n_pos 3,803), non-CpG-island rows 0.006 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.631 | 0.928 [0.916, 0.939] | 0.851 [0.837, 0.864] | 0.773 [0.753, 0.791] | 0.977 [0.974, 0.980] |
| pydreg, −log10 p over all candidates | 0.320 | 0.911 [0.896, 0.926] | 0.786 [0.771, 0.800] | 0.645 [0.629, 0.661] | 0.971 [0.967, 0.974] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.625 | 0.935 [0.925, 0.946] | 0.856 [0.840, 0.869] | 0.786 [0.762, 0.808] | 0.980 [0.977, 0.982] |
| Δ pydreg(SVR) − PINTS | 0.007 | -0.006 [-0.017, 0.003] p=0.172 | -0.005 [-0.016, 0.006] p=0.402 | -0.013 [-0.038, 0.011] p=0.318 | -0.003 [-0.005, -0.000] p=0.034 |
| Δ pydreg(−log p) − PINTS | -0.305 | -0.024 [-0.036, -0.012] p<1e-03 | -0.069 [-0.080, -0.058] p<1e-03 | -0.141 [-0.162, -0.118] p<1e-03 | -0.009 [-0.012, -0.007] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.762 [0.726, 0.797] | 0.666 [0.634, 0.695] | 0.555 [0.519, 0.598] | 0.952 [0.946, 0.957] |
| GC-matched: PINTS |  | 0.823 [0.796, 0.847] | 0.751 [0.722, 0.776] | 0.698 [0.654, 0.739] | 0.964 [0.960, 0.968] |
| GC-matched: Δ pydreg − PINTS |  | -0.061 [-0.097, -0.026] p<1e-03 | -0.086 [-0.118, -0.052] p<1e-03 | -0.143 [-0.196, -0.085] p<1e-03 | -0.012 [-0.017, -0.007] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.055 [-0.097, -0.020] p=0.004 | -0.084 [-0.119, -0.046] p<1e-03 | -0.146 [-0.204, -0.080] p<1e-03 | -0.010 [-0.015, -0.005] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.046 [-0.099, 0.009] p=0.122 | -0.082 [-0.127, -0.038] p=0.002 | -0.122 [-0.179, -0.066] p<1e-03 | -0.031 [-0.050, -0.013] p=0.004 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.744 | 0.768 | 0.615 | 0.813 |
| 20,000 | 0.917 | 0.579 | 0.825 | 0.708 |
| 40,000 | 0.973 | 0.355 | 0.935 | 0.553 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 4,154, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.063**; CpG-island rows 0.058 (n_pos 3,803), non-CpG-island rows 0.006 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.354 | 0.911 [0.898, 0.923] | 0.789 [0.776, 0.804] | 0.651 [0.631, 0.673] | 0.973 [0.971, 0.976] |
| pydreg, −log10 p over all candidates | 0.176 | 0.864 [0.842, 0.883] | 0.696 [0.681, 0.712] | 0.532 [0.514, 0.549] | 0.964 [0.961, 0.967] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.764 | 0.945 [0.937, 0.953] | 0.894 [0.880, 0.905] | 0.840 [0.807, 0.866] | 0.986 [0.984, 0.988] |
| Δ pydreg(SVR) − PINTS | -0.411 | -0.034 [-0.045, -0.024] p<1e-03 | -0.104 [-0.118, -0.090] p<1e-03 | -0.189 [-0.220, -0.152] p<1e-03 | -0.013 [-0.015, -0.010] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.588 | -0.081 [-0.101, -0.064] p<1e-03 | -0.198 [-0.213, -0.182] p<1e-03 | -0.308 [-0.337, -0.275] p<1e-03 | -0.022 [-0.025, -0.020] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.487 [0.409, 0.550] | 0.462 [0.421, 0.501] | 0.369 [0.341, 0.404] | 0.930 [0.924, 0.936] |
| GC-matched: PINTS |  | 0.865 [0.852, 0.881] | 0.798 [0.778, 0.817] | 0.752 [0.722, 0.783] | 0.972 [0.969, 0.975] |
| GC-matched: Δ pydreg − PINTS |  | -0.378 [-0.457, -0.319] p<1e-03 | -0.336 [-0.377, -0.296] p<1e-03 | -0.384 [-0.425, -0.342] p<1e-03 | -0.042 [-0.048, -0.036] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.417 [-0.470, -0.338] p<1e-03 | -0.349 [-0.393, -0.305] p<1e-03 | -0.397 [-0.439, -0.351] p<1e-03 | -0.041 [-0.047, -0.035] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.191 [-0.262, -0.097] p<1e-03 | -0.201 [-0.250, -0.148] p<1e-03 | -0.181 [-0.233, -0.138] p<1e-03 | -0.048 [-0.069, -0.028] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.659 | 0.680 | 0.743 | 0.844 |
| 20,000 | 0.920 | 0.533 | 0.920 | 0.654 |
| 40,000 | 0.981 | 0.322 | 0.972 | 0.453 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 4,154, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.063**; CpG-island rows 0.058 (n_pos 3,803), non-CpG-island rows 0.006 (n_pos 351).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.665 | 0.914 [0.900, 0.928] | 0.848 [0.834, 0.861] | 0.782 [0.765, 0.799] | 0.973 [0.969, 0.977] |
| pydreg, −log10 p over all candidates | 0.353 | 0.905 [0.888, 0.918] | 0.792 [0.776, 0.806] | 0.659 [0.641, 0.674] | 0.968 [0.963, 0.971] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.636 | 0.948 [0.938, 0.956] | 0.864 [0.850, 0.877] | 0.797 [0.766, 0.823] | 0.980 [0.977, 0.983] |
| Δ pydreg(SVR) − PINTS | 0.030 | -0.034 [-0.046, -0.021] p<1e-03 | -0.017 [-0.029, -0.003] p=0.012 | -0.015 [-0.043, 0.017] p=0.296 | -0.007 [-0.010, -0.003] p=0.002 |
| Δ pydreg(−log p) − PINTS | -0.283 | -0.044 [-0.058, -0.030] p<1e-03 | -0.073 [-0.086, -0.059] p<1e-03 | -0.138 [-0.164, -0.109] p<1e-03 | -0.013 [-0.016, -0.009] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.747 [0.711, 0.782] | 0.656 [0.625, 0.686] | 0.565 [0.528, 0.603] | 0.946 [0.941, 0.952] |
| GC-matched: PINTS |  | 0.816 [0.788, 0.829] | 0.718 [0.693, 0.740] | 0.665 [0.631, 0.695] | 0.961 [0.957, 0.964] |
| GC-matched: Δ pydreg − PINTS |  | -0.069 [-0.097, -0.030] p<1e-03 | -0.061 [-0.092, -0.031] p<1e-03 | -0.100 [-0.141, -0.058] p<1e-03 | -0.014 [-0.019, -0.009] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.058 [-0.102, -0.022] p=0.004 | -0.057 [-0.090, -0.023] p<1e-03 | -0.099 [-0.144, -0.052] p<1e-03 | -0.012 [-0.017, -0.007] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.078 [-0.132, -0.012] p=0.02 | -0.080 [-0.125, -0.031] p<1e-03 | -0.096 [-0.143, -0.052] p<1e-03 | -0.032 [-0.055, -0.011] p=0.004 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.742 | 0.777 | 0.559 | 0.843 |
| 20,000 | 0.902 | 0.576 | 0.649 | 0.789 |
| 40,000 | 0.963 | 0.357 | 0.838 | 0.678 |


## PINTS-on-TSS labels: pints_tss_distal vs TSS-silent DHSs

* **Positives:** pints_tss > 1 kb from any GENCODE v24 TSS.
* **Negatives:** DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).
* **Notes:** Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces dhs_unmarked as the specificity test. The negatives come from the same TSS library as the positives, so the two sets are disjoint by construction.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 9,343, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.348**; CpG-island rows 0.101 (n_pos 1,967), non-CpG-island rows 0.296 (n_pos 7,376).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.818 | 0.931 [0.924, 0.939] | 0.897 [0.887, 0.906] | 0.946 [0.941, 0.952] | 0.960 [0.957, 0.964] |
| pydreg, −log10 p over all candidates | 0.826 | 0.931 [0.924, 0.939] | 0.889 [0.879, 0.899] | 0.942 [0.936, 0.948] | 0.960 [0.956, 0.964] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.543 | 0.870 [0.857, 0.886] | 0.803 [0.786, 0.819] | 0.943 [0.936, 0.949] | 0.969 [0.966, 0.971] |
| Δ pydreg(SVR) − PINTS | 0.275 | 0.061 [0.045, 0.076] p<1e-03 | 0.094 [0.078, 0.109] p<1e-03 | 0.004 [-0.002, 0.010] p=0.234 | -0.008 [-0.012, -0.005] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.284 | 0.061 [0.045, 0.076] p<1e-03 | 0.086 [0.071, 0.101] p<1e-03 | -0.001 [-0.007, 0.005] p=0.772 | -0.009 [-0.013, -0.005] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.931 [0.924, 0.938] | 0.889 [0.877, 0.901] | 0.944 [0.939, 0.950] | 0.960 [0.956, 0.964] |
| GC-matched: PINTS |  | 0.867 [0.853, 0.885] | 0.800 [0.782, 0.818] | 0.942 [0.935, 0.949] | 0.968 [0.965, 0.971] |
| GC-matched: Δ pydreg − PINTS |  | 0.064 [0.047, 0.079] p<1e-03 | 0.089 [0.072, 0.105] p<1e-03 | 0.002 [-0.004, 0.009] p=0.47 | -0.009 [-0.013, -0.005] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.044 [0.021, 0.074] p<1e-03 | 0.045 [0.019, 0.070] p=0.002 | 0.028 [-0.003, 0.065] p=0.084 | -0.007 [-0.014, 0.000] p=0.072 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.069 [0.049, 0.082] p<1e-03 | 0.102 [0.084, 0.119] p<1e-03 | 0.011 [0.003, 0.019] p=0.004 | -0.009 [-0.014, -0.005] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.230 | 0.999 | 0.107 | 0.993 |
| 20,000 | 0.430 | 0.997 | 0.190 | 0.992 |
| 40,000 | 0.674 | 0.989 | 0.314 | 0.982 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 9,343, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.348**; CpG-island rows 0.101 (n_pos 1,967), non-CpG-island rows 0.296 (n_pos 7,376).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.670 | 0.683 [0.672, 0.696] | 0.681 [0.670, 0.694] | 0.783 [0.772, 0.794] | 0.833 [0.827, 0.840] |
| pydreg, −log10 p over all candidates | 0.670 | 0.683 [0.672, 0.696] | 0.681 [0.669, 0.694] | 0.782 [0.771, 0.794] | 0.833 [0.827, 0.840] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.650 | 0.728 [0.716, 0.740] | 0.699 [0.687, 0.713] | 0.804 [0.793, 0.816] | 0.854 [0.848, 0.860] |
| Δ pydreg(SVR) − PINTS | 0.020 | -0.045 [-0.055, -0.035] p<1e-03 | -0.018 [-0.028, -0.008] p=0.002 | -0.022 [-0.029, -0.015] p<1e-03 | -0.021 [-0.026, -0.016] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.020 | -0.045 [-0.055, -0.035] p<1e-03 | -0.019 [-0.029, -0.009] p<1e-03 | -0.022 [-0.029, -0.015] p<1e-03 | -0.021 [-0.026, -0.016] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.683 [0.672, 0.696] | 0.682 [0.670, 0.694] | 0.783 [0.772, 0.794] | 0.833 [0.827, 0.840] |
| GC-matched: PINTS |  | 0.728 [0.716, 0.740] | 0.701 [0.688, 0.714] | 0.805 [0.794, 0.816] | 0.854 [0.848, 0.860] |
| GC-matched: Δ pydreg − PINTS |  | -0.045 [-0.055, -0.035] p<1e-03 | -0.019 [-0.030, -0.009] p=0.002 | -0.022 [-0.029, -0.016] p<1e-03 | -0.021 [-0.026, -0.016] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.000 [-0.021, 0.019] p=0.994 | 0.021 [-0.000, 0.041] p=0.052 | 0.036 [0.014, 0.059] p=0.002 | 0.002 [-0.009, 0.012] p=0.718 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.057 [-0.068, -0.046] p<1e-03 | -0.030 [-0.041, -0.019] p<1e-03 | -0.030 [-0.039, -0.023] p<1e-03 | -0.027 [-0.033, -0.022] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.181 | 1.000 | 0.101 | 0.994 |
| 20,000 | 0.381 | 1.000 | 0.186 | 0.994 |
| 40,000 | 0.580 | 0.998 | 0.312 | 0.993 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 9,343, n_neg = 17,523, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.348**; CpG-island rows 0.101 (n_pos 1,967), non-CpG-island rows 0.296 (n_pos 7,376).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.795 | 0.803 [0.789, 0.815] | 0.798 [0.784, 0.810] | 0.864 [0.855, 0.872] | 0.896 [0.889, 0.902] |
| pydreg, −log10 p over all candidates | 0.795 | 0.803 [0.789, 0.815] | 0.797 [0.784, 0.810] | 0.863 [0.854, 0.871] | 0.896 [0.889, 0.902] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.645 | 0.842 [0.827, 0.854] | 0.790 [0.774, 0.804] | 0.876 [0.867, 0.885] | 0.911 [0.904, 0.918] |
| Δ pydreg(SVR) − PINTS | 0.150 | -0.039 [-0.047, -0.030] p<1e-03 | 0.008 [-0.003, 0.020] p=0.14 | -0.013 [-0.019, -0.006] p<1e-03 | -0.016 [-0.020, -0.011] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.150 | -0.039 [-0.047, -0.030] p<1e-03 | 0.008 [-0.003, 0.020] p=0.186 | -0.013 [-0.019, -0.006] p<1e-03 | -0.016 [-0.020, -0.011] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.803 [0.789, 0.815] | 0.794 [0.780, 0.807] | 0.863 [0.854, 0.871] | 0.895 [0.888, 0.902] |
| GC-matched: PINTS |  | 0.841 [0.826, 0.853] | 0.783 [0.767, 0.799] | 0.874 [0.864, 0.883] | 0.911 [0.903, 0.917] |
| GC-matched: Δ pydreg − PINTS |  | -0.039 [-0.047, -0.029] p<1e-03 | 0.011 [-0.000, 0.024] p=0.056 | -0.011 [-0.018, -0.003] p=0.006 | -0.015 [-0.020, -0.010] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.003 [-0.019, 0.017] p=0.798 | 0.039 [0.015, 0.067] p<1e-03 | 0.063 [0.026, 0.107] p<1e-03 | 0.003 [-0.006, 0.012] p=0.554 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.048 [-0.059, -0.038] p<1e-03 | 0.003 [-0.009, 0.016] p=0.6 | -0.015 [-0.023, -0.007] p<1e-03 | -0.020 [-0.026, -0.015] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.170 | 0.999 | 0.091 | 0.999 |
| 20,000 | 0.383 | 0.999 | 0.171 | 0.996 |
| 40,000 | 0.675 | 0.996 | 0.345 | 0.992 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 14,595, n_neg = 11,691, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.555**; CpG-island rows 0.194 (n_pos 2,806), non-CpG-island rows 0.502 (n_pos 11,789).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.599 | 0.616 [0.605, 0.625] | 0.615 [0.604, 0.625] | 0.820 [0.808, 0.831] | 0.798 [0.792, 0.803] |
| pydreg, −log10 p over all candidates | 0.599 | 0.616 [0.605, 0.625] | 0.615 [0.604, 0.625] | 0.820 [0.808, 0.831] | 0.798 [0.792, 0.803] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.651 | 0.666 [0.651, 0.679] | 0.659 [0.643, 0.672] | 0.843 [0.830, 0.855] | 0.823 [0.815, 0.830] |
| Δ pydreg(SVR) − PINTS | -0.052 | -0.050 [-0.060, -0.039] p<1e-03 | -0.043 [-0.053, -0.033] p<1e-03 | -0.023 [-0.027, -0.018] p<1e-03 | -0.026 [-0.031, -0.020] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.052 | -0.050 [-0.060, -0.039] p<1e-03 | -0.043 [-0.053, -0.033] p<1e-03 | -0.023 [-0.027, -0.018] p<1e-03 | -0.026 [-0.031, -0.020] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.615 [0.605, 0.625] | 0.615 [0.604, 0.625] | 0.820 [0.808, 0.831] | 0.798 [0.792, 0.803] |
| GC-matched: PINTS |  | 0.664 [0.649, 0.677] | 0.650 [0.633, 0.665] | 0.842 [0.829, 0.853] | 0.822 [0.814, 0.829] |
| GC-matched: Δ pydreg − PINTS |  | -0.049 [-0.058, -0.038] p<1e-03 | -0.035 [-0.048, -0.020] p<1e-03 | -0.022 [-0.026, -0.017] p<1e-03 | -0.024 [-0.030, -0.018] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.015 [-0.034, 0.004] p=0.132 | 0.001 [-0.025, 0.031] p=0.998 | -0.004 [-0.024, 0.017] p=0.67 | -0.006 [-0.016, 0.004] p=0.23 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.057 [-0.066, -0.046] p<1e-03 | -0.044 [-0.056, -0.031] p<1e-03 | -0.028 [-0.033, -0.023] p<1e-03 | -0.029 [-0.034, -0.023] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.101 | 1.000 | 0.041 | 1.000 |
| 20,000 | 0.280 | 1.000 | 0.087 | 0.999 |
| 40,000 | 0.561 | 1.000 | 0.184 | 0.999 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 3,093, n_neg = 62,494, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.047**; CpG-island rows 0.013 (n_pos 797), non-CpG-island rows 0.035 (n_pos 2,296).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.134 | 0.267 [0.245, 0.285] | 0.250 [0.231, 0.267] | 0.152 [0.135, 0.171] | 0.630 [0.618, 0.640] |
| pydreg, −log10 p over all candidates | 0.098 | 0.260 [0.239, 0.279] | 0.239 [0.221, 0.256] | 0.113 [0.100, 0.128] | 0.629 [0.617, 0.639] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.116 | 0.303 [0.278, 0.327] | 0.269 [0.250, 0.289] | 0.146 [0.127, 0.166] | 0.651 [0.638, 0.663] |
| Δ pydreg(SVR) − PINTS | 0.018 | -0.037 [-0.055, -0.016] p=0.002 | -0.020 [-0.032, -0.008] p=0.002 | 0.006 [-0.007, 0.021] p=0.394 | -0.021 [-0.029, -0.014] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.018 | -0.044 [-0.061, -0.021] p<1e-03 | -0.030 [-0.042, -0.018] p<1e-03 | -0.033 [-0.045, -0.021] p<1e-03 | -0.022 [-0.030, -0.015] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.195 [0.178, 0.213] | 0.179 [0.164, 0.195] | 0.095 [0.084, 0.110] | 0.605 [0.594, 0.616] |
| GC-matched: PINTS |  | 0.208 [0.190, 0.227] | 0.192 [0.177, 0.210] | 0.099 [0.087, 0.115] | 0.621 [0.609, 0.633] |
| GC-matched: Δ pydreg − PINTS |  | -0.013 [-0.028, 0.001] p=0.062 | -0.013 [-0.026, -0.001] p=0.024 | -0.004 [-0.014, 0.007] p=0.448 | -0.016 [-0.024, -0.009] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.001 [-0.058, 0.052] p=0.96 | 0.007 [-0.028, 0.039] p=0.672 | 0.005 [-0.010, 0.029] p=0.47 | 0.000 [-0.014, 0.016] p=0.99 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.024 [-0.037, -0.008] p=0.002 | -0.025 [-0.037, -0.014] p<1e-03 | -0.008 [-0.014, -0.003] p<1e-03 | -0.023 [-0.032, -0.015] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.112 | 0.450 | 0.070 | 0.404 |
| 20,000 | 0.205 | 0.282 | 0.129 | 0.360 |
| 40,000 | 0.310 | 0.180 | 0.272 | 0.231 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 336, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.005**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.005 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.033 | 0.071 [0.044, 0.103] | 0.071 [0.045, 0.100] | 0.007 [0.006, 0.010] | 0.480 [0.462, 0.499] |
| pydreg, −log10 p over all candidates | 0.022 | 0.080 [0.050, 0.114] | 0.070 [0.044, 0.099] | 0.006 [0.005, 0.007] | 0.480 [0.461, 0.499] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.054 | 0.155 [0.115, 0.200] | 0.151 [0.117, 0.187] | 0.012 [0.008, 0.021] | 0.582 [0.550, 0.618] |
| Δ pydreg(SVR) − PINTS | -0.021 | -0.083 [-0.122, -0.050] p<1e-03 | -0.079 [-0.108, -0.053] p<1e-03 | -0.005 [-0.013, -0.001] p=0.01 | -0.103 [-0.133, -0.075] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.031 | -0.074 [-0.111, -0.042] p<1e-03 | -0.080 [-0.108, -0.054] p<1e-03 | -0.006 [-0.014, -0.003] p<1e-03 | -0.103 [-0.133, -0.075] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.071 [0.045, 0.103] | 0.074 [0.048, 0.104] | 0.007 [0.006, 0.011] | 0.488 [0.470, 0.507] |
| GC-matched: PINTS |  | 0.173 [0.129, 0.218] | 0.163 [0.128, 0.202] | 0.013 [0.009, 0.023] | 0.594 [0.561, 0.629] |
| GC-matched: Δ pydreg − PINTS |  | -0.101 [-0.139, -0.066] p<1e-03 | -0.089 [-0.121, -0.061] p<1e-03 | -0.006 [-0.014, -0.002] p=0.008 | -0.106 [-0.136, -0.079] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.065 [-0.100, 0.200] p=0.602 | 0.020 [-0.087, 0.127] p=0.7 | -0.004 [-0.060, 0.007] p=0.53 | -0.028 [-0.097, 0.031] p=0.378 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.112 [-0.151, -0.075] p<1e-03 | -0.103 [-0.136, -0.074] p<1e-03 | -0.005 [-0.012, -0.002] p<1e-03 | -0.114 [-0.146, -0.084] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.060 | 0.021 | 0.039 | 0.026 |
| 20,000 | 0.071 | 0.009 | 0.092 | 0.026 |
| 40,000 | 0.101 | 0.005 | 0.137 | 0.019 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 336, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.005**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.005 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.027 | 0.071 [0.043, 0.104] | 0.064 [0.040, 0.092] | 0.008 [0.005, 0.018] | 0.476 [0.458, 0.497] |
| pydreg, −log10 p over all candidates | 0.018 | 0.077 [0.047, 0.110] | 0.063 [0.039, 0.090] | 0.006 [0.005, 0.007] | 0.476 [0.458, 0.497] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.042 | 0.130 [0.094, 0.169] | 0.125 [0.095, 0.159] | 0.011 [0.007, 0.022] | 0.554 [0.525, 0.588] |
| Δ pydreg(SVR) − PINTS | -0.015 | -0.059 [-0.087, -0.032] p<1e-03 | -0.061 [-0.084, -0.040] p<1e-03 | -0.003 [-0.014, 0.007] p=0.462 | -0.078 [-0.106, -0.054] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.023 | -0.053 [-0.085, -0.024] p<1e-03 | -0.062 [-0.086, -0.039] p<1e-03 | -0.005 [-0.016, -0.002] p<1e-03 | -0.078 [-0.106, -0.053] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.074 [0.045, 0.105] | 0.068 [0.043, 0.095] | 0.008 [0.006, 0.018] | 0.486 [0.468, 0.506] |
| GC-matched: PINTS |  | 0.146 [0.105, 0.184] | 0.138 [0.105, 0.173] | 0.012 [0.008, 0.024] | 0.567 [0.538, 0.600] |
| GC-matched: Δ pydreg − PINTS |  | -0.071 [-0.102, -0.043] p<1e-03 | -0.070 [-0.095, -0.047] p<1e-03 | -0.004 [-0.015, 0.006] p=0.4 | -0.081 [-0.109, -0.057] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.000 [-0.083, 0.111] p=1 | -0.015 [-0.100, 0.072] p=0.784 | -0.001 [-0.031, 0.025] p=0.722 | -0.012 [-0.077, 0.040] p=0.718 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.079 [-0.116, -0.049] p<1e-03 | -0.081 [-0.108, -0.058] p<1e-03 | -0.003 [-0.015, 0.006] p=0.42 | -0.088 [-0.117, -0.063] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.036 | 0.013 | 0.042 | 0.023 |
| 20,000 | 0.071 | 0.009 | 0.077 | 0.018 |
| 40,000 | 0.095 | 0.004 | 0.131 | 0.014 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 336, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.005**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.005 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.021 | 0.074 [0.048, 0.104] | 0.061 [0.040, 0.086] | 0.006 [0.005, 0.007] | 0.430 [0.411, 0.451] |
| pydreg, −log10 p over all candidates | 0.010 | 0.060 [0.035, 0.088] | 0.052 [0.033, 0.074] | 0.005 [0.005, 0.006] | 0.429 [0.410, 0.450] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.042 | 0.098 [0.067, 0.135] | 0.107 [0.077, 0.140] | 0.012 [0.007, 0.023] | 0.525 [0.493, 0.561] |
| Δ pydreg(SVR) − PINTS | -0.021 | -0.024 [-0.057, 0.006] p=0.12 | -0.046 [-0.071, -0.024] p<1e-03 | -0.007 [-0.017, -0.002] p<1e-03 | -0.095 [-0.126, -0.067] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.032 | -0.039 [-0.071, -0.009] p=0.014 | -0.055 [-0.077, -0.034] p<1e-03 | -0.007 [-0.018, -0.002] p<1e-03 | -0.096 [-0.127, -0.068] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.074 [0.048, 0.104] | 0.066 [0.044, 0.091] | 0.006 [0.005, 0.007] | 0.445 [0.427, 0.466] |
| GC-matched: PINTS |  | 0.119 [0.081, 0.159] | 0.118 [0.087, 0.151] | 0.013 [0.008, 0.025] | 0.542 [0.510, 0.577] |
| GC-matched: Δ pydreg − PINTS |  | -0.045 [-0.077, -0.012] p=0.006 | -0.052 [-0.077, -0.029] p<1e-03 | -0.007 [-0.019, -0.002] p<1e-03 | -0.097 [-0.128, -0.069] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.065 [-0.167, 0.000] p=0.482 | -0.055 [-0.137, 0.004] p=0.076 | -0.013 [-0.088, -0.000] p=0.004 | -0.029 [-0.092, 0.019] p=0.268 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.052 [-0.087, -0.019] p=0.004 | -0.058 [-0.087, -0.033] p<1e-03 | -0.004 [-0.014, -0.001] p<1e-03 | -0.105 [-0.137, -0.075] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.036 | 0.009 | 0.042 | 0.024 |
| 20,000 | 0.074 | 0.007 | 0.092 | 0.015 |
| 40,000 | 0.095 | 0.004 | 0.143 | 0.010 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 336, n_neg = 62,150, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.005**; CpG-island rows 0.000 (n_pos 31), non-CpG-island rows 0.005 (n_pos 305).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.042 | 0.077 [0.048, 0.109] | 0.069 [0.045, 0.096] | 0.007 [0.005, 0.009] | 0.474 [0.456, 0.493] |
| pydreg, −log10 p over all candidates | 0.025 | 0.080 [0.052, 0.113] | 0.071 [0.046, 0.097] | 0.006 [0.005, 0.007] | 0.474 [0.455, 0.493] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.078 | 0.198 [0.153, 0.247] | 0.182 [0.143, 0.224] | 0.013 [0.010, 0.021] | 0.569 [0.535, 0.605] |
| Δ pydreg(SVR) − PINTS | -0.036 | -0.121 [-0.162, -0.080] p<1e-03 | -0.112 [-0.150, -0.078] p<1e-03 | -0.007 [-0.013, -0.003] p<1e-03 | -0.095 [-0.126, -0.064] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.053 | -0.118 [-0.164, -0.075] p<1e-03 | -0.110 [-0.149, -0.075] p<1e-03 | -0.007 [-0.015, -0.004] p<1e-03 | -0.095 [-0.126, -0.064] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.077 [0.049, 0.109] | 0.073 [0.048, 0.100] | 0.007 [0.006, 0.010] | 0.484 [0.466, 0.503] |
| GC-matched: PINTS |  | 0.211 [0.169, 0.258] | 0.192 [0.154, 0.233] | 0.015 [0.011, 0.023] | 0.580 [0.546, 0.615] |
| GC-matched: Δ pydreg − PINTS |  | -0.134 [-0.176, -0.094] p<1e-03 | -0.119 [-0.157, -0.083] p<1e-03 | -0.008 [-0.015, -0.004] p<1e-03 | -0.097 [-0.128, -0.067] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.000 [-0.128, 0.121] p=1 | -0.036 [-0.125, 0.042] p=0.376 | -0.001 [-0.009, 0.006] p=0.652 | -0.047 [-0.124, 0.014] p=0.152 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.150 [-0.194, -0.103] p<1e-03 | -0.129 [-0.168, -0.089] p<1e-03 | -0.007 [-0.015, -0.003] p<1e-03 | -0.102 [-0.133, -0.069] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.048 | 0.018 | 0.062 | 0.046 |
| 20,000 | 0.074 | 0.009 | 0.086 | 0.039 |
| 40,000 | 0.101 | 0.005 | 0.140 | 0.028 |


## dREG auxiliary labels: dhs_h3k27ac vs dhs_unmarked

* **Positives:** Duke ∩ UW DHS overlapping an ENCODE Broad H3K27ac peak, the "DHSs marked by H3K27ac" of Wang 2019 Fig 1B (K562 34,763; GM12878 55,925).
* **Negatives:** **Diagnostic only, not a headline contrast.** Duke ∩ UW DHSs with no H3K27ac and outside ChromHMM states 1–7 (K562 23,291, GM12878 15,556; CTCF/insulator-rich). Lacking H3K27ac does not mean untranscribed, so a call here may be a true weak element. Results on this set cannot distinguish specificity from sensitivity.
* **Notes:** Not in either dREG paper as a ROC: the papers report recall only, so the negatives are ours. No nascent or TSS data is used, so these labels are not circular for any library, GRO-cap included. H3K27ac peaks and ChromHMM segments are broad, which favours wider calls under the overlap rule. On hg38 cells (HCT116, A549) the DHS is one ENCODE DNase peak set, not Duke ∩ UW, and there is no ChromHMM: quiet excludes only DNase and H3K27ac. Jurkat (hg19) uses UW DNase rep1 ∩ rep2 with H3K27ac from either ChIP-Atlas replicate (GSM4551183/4), and has no ChromHMM.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 34,763, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.599**; CpG-island rows 0.300 (n_pos 9,989), non-CpG-island rows 0.515 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.532 | 0.735 [0.724, 0.747] | 0.692 [0.680, 0.703] | 0.904 [0.900, 0.908] | 0.874 [0.870, 0.879] |
| pydreg, −log10 p over all candidates | 0.491 | 0.731 [0.721, 0.743] | 0.675 [0.662, 0.687] | 0.899 [0.895, 0.903] | 0.872 [0.868, 0.877] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.245 | 0.505 [0.482, 0.529] | 0.477 [0.458, 0.497] | 0.909 [0.904, 0.913] | 0.887 [0.883, 0.892] |
| Δ pydreg(SVR) − PINTS | 0.287 | 0.231 [0.208, 0.253] p<1e-03 | 0.215 [0.196, 0.233] p<1e-03 | -0.005 [-0.009, -0.000] p=0.032 | -0.013 [-0.017, -0.008] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.245 | 0.226 [0.204, 0.249] p<1e-03 | 0.198 [0.179, 0.217] p<1e-03 | -0.009 [-0.014, -0.005] p<1e-03 | -0.015 [-0.019, -0.010] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.692 [0.675, 0.710] | 0.650 [0.634, 0.666] | 0.899 [0.895, 0.903] | 0.868 [0.863, 0.873] |
| GC-matched: PINTS |  | 0.478 [0.454, 0.504] | 0.454 [0.433, 0.477] | 0.904 [0.898, 0.909] | 0.880 [0.875, 0.886] |
| GC-matched: Δ pydreg − PINTS |  | 0.214 [0.188, 0.238] p<1e-03 | 0.196 [0.173, 0.219] p<1e-03 | -0.004 [-0.009, 0.001] p=0.12 | -0.012 [-0.018, -0.007] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.095 [0.067, 0.125] p<1e-03 | 0.114 [0.082, 0.152] p<1e-03 | 0.031 [0.019, 0.045] p<1e-03 | 0.011 [0.006, 0.016] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.280 [0.254, 0.306] p<1e-03 | 0.250 [0.229, 0.272] p<1e-03 | 0.007 [-0.001, 0.015] p=0.096 | -0.023 [-0.029, -0.018] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.227 | 0.998 | 0.137 | 0.985 |
| 20,000 | 0.379 | 0.996 | 0.210 | 0.977 |
| 40,000 | 0.534 | 0.987 | 0.296 | 0.967 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 34,763, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.599**; CpG-island rows 0.300 (n_pos 9,989), non-CpG-island rows 0.515 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.493 | 0.552 [0.543, 0.561] | 0.539 [0.530, 0.548] | 0.812 [0.806, 0.818] | 0.763 [0.758, 0.768] |
| pydreg, −log10 p over all candidates | 0.492 | 0.552 [0.543, 0.561] | 0.537 [0.528, 0.547] | 0.811 [0.805, 0.817] | 0.763 [0.758, 0.768] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.299 | 0.581 [0.563, 0.601] | 0.506 [0.495, 0.519] | 0.823 [0.817, 0.829] | 0.780 [0.775, 0.785] |
| Δ pydreg(SVR) − PINTS | 0.194 | -0.029 [-0.049, -0.011] p<1e-03 | 0.032 [0.021, 0.043] p<1e-03 | -0.012 [-0.014, -0.009] p<1e-03 | -0.017 [-0.020, -0.014] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.193 | -0.029 [-0.049, -0.011] p<1e-03 | 0.031 [0.020, 0.042] p<1e-03 | -0.013 [-0.015, -0.010] p<1e-03 | -0.017 [-0.020, -0.014] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.550 [0.540, 0.559] | 0.531 [0.522, 0.541] | 0.811 [0.805, 0.817] | 0.761 [0.756, 0.766] |
| GC-matched: PINTS |  | 0.561 [0.543, 0.583] | 0.499 [0.487, 0.512] | 0.822 [0.816, 0.828] | 0.778 [0.773, 0.783] |
| GC-matched: Δ pydreg − PINTS |  | -0.012 [-0.033, 0.008] p=0.228 | 0.032 [0.019, 0.044] p<1e-03 | -0.011 [-0.014, -0.009] p<1e-03 | -0.017 [-0.020, -0.014] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.081 [0.058, 0.124] p<1e-03 | 0.113 [0.091, 0.134] p<1e-03 | 0.052 [0.043, 0.060] p<1e-03 | 0.026 [0.022, 0.030] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.067 [-0.084, -0.045] p<1e-03 | 0.001 [-0.010, 0.011] p=0.972 | -0.024 [-0.028, -0.020] p<1e-03 | -0.034 [-0.038, -0.031] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.224 | 0.999 | 0.109 | 0.992 |
| 20,000 | 0.351 | 0.996 | 0.189 | 0.985 |
| 40,000 | 0.473 | 0.989 | 0.292 | 0.979 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 34,763, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.599**; CpG-island rows 0.300 (n_pos 9,989), non-CpG-island rows 0.515 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.539 | 0.595 [0.582, 0.606] | 0.580 [0.567, 0.590] | 0.829 [0.823, 0.835] | 0.785 [0.779, 0.791] |
| pydreg, −log10 p over all candidates | 0.536 | 0.595 [0.582, 0.606] | 0.578 [0.566, 0.589] | 0.829 [0.822, 0.834] | 0.785 [0.779, 0.791] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.317 | 0.592 [0.573, 0.610] | 0.538 [0.525, 0.552] | 0.850 [0.844, 0.855] | 0.812 [0.806, 0.818] |
| Δ pydreg(SVR) − PINTS | 0.222 | 0.003 [-0.018, 0.023] p=0.81 | 0.042 [0.027, 0.056] p<1e-03 | -0.020 [-0.023, -0.018] p<1e-03 | -0.027 [-0.030, -0.024] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.219 | 0.003 [-0.018, 0.023] p=0.81 | 0.040 [0.025, 0.055] p<1e-03 | -0.021 [-0.024, -0.018] p<1e-03 | -0.027 [-0.030, -0.024] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.590 [0.578, 0.601] | 0.556 [0.543, 0.568] | 0.826 [0.819, 0.832] | 0.781 [0.774, 0.787] |
| GC-matched: PINTS |  | 0.545 [0.519, 0.578] | 0.512 [0.495, 0.528] | 0.846 [0.840, 0.851] | 0.808 [0.802, 0.813] |
| GC-matched: Δ pydreg − PINTS |  | 0.046 [0.010, 0.074] p<1e-03 | 0.044 [0.026, 0.060] p<1e-03 | -0.020 [-0.023, -0.017] p<1e-03 | -0.027 [-0.030, -0.023] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.109 [0.074, 0.126] p<1e-03 | 0.071 [0.046, 0.093] p<1e-03 | 0.023 [0.013, 0.032] p<1e-03 | 0.013 [0.009, 0.017] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.027 [-0.056, -0.003] p=0.02 | 0.022 [0.006, 0.037] p=0.004 | -0.029 [-0.033, -0.025] p<1e-03 | -0.044 [-0.048, -0.040] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.221 | 0.999 | 0.149 | 0.996 |
| 20,000 | 0.351 | 0.997 | 0.205 | 0.991 |
| 40,000 | 0.494 | 0.991 | 0.321 | 0.979 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 55,925, n_neg = 15,556, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.782**; CpG-island rows 0.531 (n_pos 17,611), non-CpG-island rows 0.711 (n_pos 38,314).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.534 | 0.788 [0.774, 0.799] | 0.728 [0.716, 0.739] | 0.959 [0.957, 0.960] | 0.899 [0.896, 0.902] |
| pydreg, −log10 p over all candidates | 0.465 | 0.785 [0.771, 0.796] | 0.709 [0.694, 0.723] | 0.957 [0.955, 0.958] | 0.897 [0.894, 0.900] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.351 | 0.669 [0.650, 0.694] | 0.621 [0.602, 0.639] | 0.968 [0.967, 0.970] | 0.919 [0.916, 0.922] |
| Δ pydreg(SVR) − PINTS | 0.183 | 0.119 [0.095, 0.136] p<1e-03 | 0.108 [0.092, 0.122] p<1e-03 | -0.010 [-0.011, -0.008] p<1e-03 | -0.020 [-0.022, -0.017] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.114 | 0.116 [0.092, 0.132] p<1e-03 | 0.088 [0.071, 0.104] p<1e-03 | -0.012 [-0.013, -0.011] p<1e-03 | -0.022 [-0.024, -0.019] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.734 [0.712, 0.759] | 0.674 [0.647, 0.698] | 0.956 [0.954, 0.958] | 0.892 [0.888, 0.896] |
| GC-matched: PINTS |  | 0.614 [0.579, 0.647] | 0.573 [0.548, 0.599] | 0.966 [0.963, 0.967] | 0.911 [0.907, 0.915] |
| GC-matched: Δ pydreg − PINTS |  | 0.121 [0.089, 0.153] p<1e-03 | 0.101 [0.073, 0.126] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 | -0.019 [-0.023, -0.015] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.085 [0.016, 0.133] p=0.016 | 0.057 [-0.006, 0.112] p=0.074 | -0.001 [-0.011, 0.009] p=0.894 | -0.003 [-0.010, 0.004] p=0.402 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.133 [0.107, 0.154] p<1e-03 | 0.120 [0.102, 0.137] p<1e-03 | -0.015 [-0.017, -0.012] p<1e-03 | -0.028 [-0.031, -0.024] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.195 | 0.999 | 0.147 | 0.997 |
| 20,000 | 0.349 | 0.997 | 0.236 | 0.996 |
| 40,000 | 0.542 | 0.994 | 0.344 | 0.993 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 34,763, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.599**; CpG-island rows 0.300 (n_pos 9,989), non-CpG-island rows 0.515 (n_pos 24,774).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.571 | 0.782 [0.766, 0.794] | 0.724 [0.712, 0.735] | 0.909 [0.905, 0.913] | 0.883 [0.878, 0.887] |
| pydreg, −log10 p over all candidates | 0.572 | 0.779 [0.766, 0.793] | 0.710 [0.698, 0.723] | 0.902 [0.898, 0.907] | 0.881 [0.877, 0.886] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.611 | 0.774 [0.763, 0.786] | 0.740 [0.728, 0.751] | 0.919 [0.916, 0.923] | 0.894 [0.890, 0.898] |
| Δ pydreg(SVR) − PINTS | -0.040 | 0.008 [-0.002, 0.015] p=0.118 | -0.016 [-0.022, -0.011] p<1e-03 | -0.010 [-0.012, -0.009] p<1e-03 | -0.011 [-0.013, -0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.039 | 0.005 [-0.003, 0.014] p=0.19 | -0.030 [-0.036, -0.023] p<1e-03 | -0.017 [-0.019, -0.015] p<1e-03 | -0.013 [-0.015, -0.011] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.743 [0.726, 0.758] | 0.694 [0.678, 0.710] | 0.906 [0.901, 0.910] | 0.879 [0.873, 0.883] |
| GC-matched: PINTS |  | 0.754 [0.741, 0.766] | 0.713 [0.697, 0.726] | 0.916 [0.912, 0.920] | 0.890 [0.886, 0.894] |
| GC-matched: Δ pydreg − PINTS |  | -0.011 [-0.022, -0.002] p=0.018 | -0.018 [-0.026, -0.010] p<1e-03 | -0.010 [-0.012, -0.008] p<1e-03 | -0.012 [-0.014, -0.009] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.013 [-0.024, -0.003] p=0.008 | -0.020 [-0.043, 0.003] p=0.108 | -0.002 [-0.011, 0.007] p=0.676 | 0.001 [-0.002, 0.004] p=0.518 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.010 [-0.001, 0.017] p=0.1 | -0.021 [-0.028, -0.015] p<1e-03 | -0.017 [-0.020, -0.015] p<1e-03 | -0.017 [-0.019, -0.014] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.234 | 0.998 | 0.197 | 0.999 |
| 20,000 | 0.406 | 0.996 | 0.340 | 0.997 |
| 40,000 | 0.594 | 0.986 | 0.525 | 0.994 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 32,659, n_neg = 34,256, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.488**; CpG-island rows 0.266 (n_pos 12,409), non-CpG-island rows 0.372 (n_pos 20,250).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.171 | 0.467 [0.453, 0.484] | 0.429 [0.418, 0.440] | 0.754 [0.743, 0.765] | 0.771 [0.766, 0.776] |
| pydreg, −log10 p over all candidates | 0.151 | 0.469 [0.454, 0.485] | 0.425 [0.414, 0.436] | 0.750 [0.739, 0.761] | 0.770 [0.765, 0.776] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.173 | 0.495 [0.480, 0.509] | 0.447 [0.435, 0.459] | 0.770 [0.758, 0.781] | 0.787 [0.781, 0.793] |
| Δ pydreg(SVR) − PINTS | -0.002 | -0.027 [-0.040, -0.015] p<1e-03 | -0.019 [-0.027, -0.011] p<1e-03 | -0.016 [-0.019, -0.012] p<1e-03 | -0.017 [-0.020, -0.013] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.022 | -0.026 [-0.039, -0.013] p<1e-03 | -0.022 [-0.032, -0.014] p<1e-03 | -0.020 [-0.023, -0.016] p<1e-03 | -0.017 [-0.020, -0.014] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.256 [0.242, 0.269] | 0.237 [0.226, 0.248] | 0.690 [0.677, 0.705] | 0.720 [0.714, 0.725] |
| GC-matched: PINTS |  | 0.278 [0.262, 0.292] | 0.258 [0.244, 0.271] | 0.708 [0.693, 0.723] | 0.742 [0.736, 0.748] |
| GC-matched: Δ pydreg − PINTS |  | -0.022 [-0.038, -0.006] p=0.008 | -0.021 [-0.033, -0.008] p=0.002 | -0.018 [-0.025, -0.011] p<1e-03 | -0.022 [-0.026, -0.018] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.025 [-0.027, 0.051] p=0.408 | 0.008 [-0.022, 0.040] p=0.64 | 0.014 [-0.010, 0.037] p=0.26 | 0.004 [-0.003, 0.012] p=0.282 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.032 [-0.046, -0.018] p<1e-03 | -0.038 [-0.047, -0.030] p<1e-03 | -0.041 [-0.047, -0.035] p<1e-03 | -0.032 [-0.036, -0.028] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.286 | 0.928 | 0.121 | 0.952 |
| 20,000 | 0.442 | 0.903 | 0.202 | 0.943 |
| 40,000 | 0.590 | 0.877 | 0.317 | 0.928 |

### HCT116_PROcap (HCT116 PRO-cap, ood)

n_pos = 32,659, n_neg = 34,256, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.488**; CpG-island rows 0.266 (n_pos 12,409), non-CpG-island rows 0.372 (n_pos 20,250).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.057 | 0.279 [0.264, 0.294] | 0.274 [0.263, 0.286] | 0.787 [0.778, 0.796] | 0.840 [0.836, 0.844] |
| pydreg, −log10 p over all candidates | 0.024 | 0.234 [0.217, 0.252] | 0.237 [0.224, 0.251] | 0.766 [0.755, 0.776] | 0.837 [0.833, 0.840] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.142 | 0.429 [0.417, 0.443] | 0.397 [0.385, 0.410] | 0.843 [0.835, 0.851] | 0.872 [0.868, 0.875] |
| Δ pydreg(SVR) − PINTS | -0.085 | -0.150 [-0.167, -0.133] p<1e-03 | -0.123 [-0.135, -0.110] p<1e-03 | -0.056 [-0.061, -0.051] p<1e-03 | -0.032 [-0.034, -0.029] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.118 | -0.195 [-0.212, -0.177] p<1e-03 | -0.160 [-0.173, -0.147] p<1e-03 | -0.077 [-0.083, -0.071] p<1e-03 | -0.035 [-0.038, -0.033] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.112 [0.101, 0.121] | 0.114 [0.106, 0.123] | 0.669 [0.656, 0.682] | 0.752 [0.747, 0.758] |
| GC-matched: PINTS |  | 0.231 [0.216, 0.244] | 0.215 [0.204, 0.226] | 0.745 [0.732, 0.758] | 0.796 [0.791, 0.800] |
| GC-matched: Δ pydreg − PINTS |  | -0.119 [-0.136, -0.101] p<1e-03 | -0.101 [-0.115, -0.088] p<1e-03 | -0.076 [-0.085, -0.066] p<1e-03 | -0.043 [-0.048, -0.039] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.205 [-0.238, -0.167] p<1e-03 | -0.184 [-0.211, -0.158] p<1e-03 | -0.200 [-0.224, -0.177] p<1e-03 | -0.115 [-0.125, -0.105] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.072 [-0.088, -0.056] p<1e-03 | -0.060 [-0.070, -0.048] p<1e-03 | -0.047 [-0.055, -0.039] p<1e-03 | -0.032 [-0.035, -0.029] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.240 | 0.843 | 0.250 | 0.917 |
| 20,000 | 0.439 | 0.838 | 0.417 | 0.894 |
| 40,000 | 0.661 | 0.814 | 0.616 | 0.857 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 83,074, n_neg = 57,867, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.589**; CpG-island rows 0.287 (n_pos 23,272), non-CpG-island rows 0.508 (n_pos 59,802).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.401 | 0.429 [0.421, 0.436] | 0.420 [0.412, 0.426] | 0.753 [0.746, 0.761] | 0.699 [0.695, 0.702] |
| pydreg, −log10 p over all candidates | 0.398 | 0.429 [0.421, 0.436] | 0.419 [0.411, 0.426] | 0.752 [0.744, 0.760] | 0.698 [0.694, 0.702] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.291 | 0.500 [0.490, 0.509] | 0.452 [0.444, 0.459] | 0.779 [0.771, 0.787] | 0.732 [0.727, 0.736] |
| Δ pydreg(SVR) − PINTS | 0.109 | -0.071 [-0.075, -0.066] p<1e-03 | -0.033 [-0.038, -0.027] p<1e-03 | -0.026 [-0.028, -0.024] p<1e-03 | -0.033 [-0.035, -0.031] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.107 | -0.071 [-0.075, -0.066] p<1e-03 | -0.034 [-0.039, -0.028] p<1e-03 | -0.027 [-0.029, -0.025] p<1e-03 | -0.033 [-0.036, -0.031] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.423 [0.416, 0.430] | 0.402 [0.393, 0.410] | 0.750 [0.742, 0.757] | 0.694 [0.690, 0.698] |
| GC-matched: PINTS |  | 0.464 [0.443, 0.488] | 0.412 [0.399, 0.423] | 0.772 [0.764, 0.780] | 0.723 [0.718, 0.728] |
| GC-matched: Δ pydreg − PINTS |  | -0.041 [-0.067, -0.018] p<1e-03 | -0.010 [-0.023, 0.003] p=0.112 | -0.022 [-0.025, -0.020] p<1e-03 | -0.029 [-0.032, -0.026] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.070 [0.036, 0.134] p<1e-03 | 0.063 [0.029, 0.093] p<1e-03 | 0.014 [-0.001, 0.027] p=0.074 | 0.001 [-0.004, 0.006] p=0.842 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.089 [-0.094, -0.083] p<1e-03 | -0.045 [-0.051, -0.038] p<1e-03 | -0.034 [-0.037, -0.031] p<1e-03 | -0.042 [-0.045, -0.039] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.138 | 0.997 | 0.084 | 0.991 |
| 20,000 | 0.243 | 0.995 | 0.145 | 0.990 |
| 40,000 | 0.376 | 0.986 | 0.304 | 0.976 |

### A549_PROcap (A549 PRO-cap, ood)

n_pos = 83,074, n_neg = 57,867, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.589**; CpG-island rows 0.287 (n_pos 23,272), non-CpG-island rows 0.508 (n_pos 59,802).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.180 | 0.339 [0.329, 0.348] | 0.302 [0.293, 0.310] | 0.717 [0.707, 0.726] | 0.651 [0.645, 0.656] |
| pydreg, −log10 p over all candidates | 0.169 | 0.338 [0.328, 0.347] | 0.298 [0.289, 0.306] | 0.712 [0.702, 0.721] | 0.650 [0.645, 0.655] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.198 | 0.341 [0.331, 0.350] | 0.315 [0.306, 0.324] | 0.727 [0.717, 0.736] | 0.661 [0.656, 0.667] |
| Δ pydreg(SVR) − PINTS | -0.018 | -0.001 [-0.006, 0.002] p=0.368 | -0.013 [-0.016, -0.011] p<1e-03 | -0.010 [-0.011, -0.009] p<1e-03 | -0.011 [-0.012, -0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.029 | -0.003 [-0.007, 0.001] p=0.14 | -0.017 [-0.020, -0.015] p<1e-03 | -0.014 [-0.016, -0.013] p<1e-03 | -0.011 [-0.012, -0.010] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.253 [0.236, 0.265] | 0.234 [0.223, 0.246] | 0.699 [0.688, 0.708] | 0.633 [0.627, 0.638] |
| GC-matched: PINTS |  | 0.266 [0.251, 0.280] | 0.241 [0.228, 0.254] | 0.705 [0.694, 0.716] | 0.642 [0.636, 0.647] |
| GC-matched: Δ pydreg − PINTS |  | -0.014 [-0.024, -0.005] p=0.002 | -0.007 [-0.013, 0.000] p=0.052 | -0.007 [-0.011, -0.003] p<1e-03 | -0.009 [-0.011, -0.007] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.009 [-0.078, 0.086] p=0.808 | 0.002 [-0.036, 0.047] p=0.906 | 0.005 [-0.024, 0.037] p=0.732 | -0.005 [-0.009, 0.001] p=0.08 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.007 [-0.010, -0.002] p=0.002 | -0.012 [-0.015, -0.010] p<1e-03 | -0.010 [-0.012, -0.008] p<1e-03 | -0.012 [-0.013, -0.010] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.112 | 0.978 | 0.100 | 0.978 |
| 20,000 | 0.187 | 0.961 | 0.172 | 0.971 |
| 40,000 | 0.272 | 0.933 | 0.249 | 0.948 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 15,613, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.159**; CpG-island rows 0.126 (n_pos 11,892), non-CpG-island rows 0.043 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.155 | 0.521 [0.502, 0.543] | 0.470 [0.453, 0.487] | 0.597 [0.586, 0.610] | 0.875 [0.869, 0.882] |
| pydreg, −log10 p over all candidates | 0.101 | 0.477 [0.454, 0.501] | 0.430 [0.412, 0.448] | 0.560 [0.548, 0.572] | 0.871 [0.865, 0.877] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.112 | 0.445 [0.420, 0.475] | 0.411 [0.390, 0.432] | 0.549 [0.532, 0.566] | 0.872 [0.866, 0.878] |
| Δ pydreg(SVR) − PINTS | 0.042 | 0.076 [0.059, 0.092] p<1e-03 | 0.059 [0.048, 0.071] p<1e-03 | 0.048 [0.036, 0.060] p<1e-03 | 0.003 [0.000, 0.007] p=0.038 |
| Δ pydreg(−log p) − PINTS | -0.011 | 0.032 [0.015, 0.047] p<1e-03 | 0.019 [0.007, 0.031] p<1e-03 | 0.010 [-0.002, 0.023] p=0.096 | -0.001 [-0.004, 0.002] p=0.608 |
| GC-matched: pydreg (SVR) |  | 0.249 [0.234, 0.266] | 0.234 [0.220, 0.248] | 0.391 [0.372, 0.409] | 0.788 [0.782, 0.794] |
| GC-matched: PINTS |  | 0.155 [0.138, 0.173] | 0.161 [0.147, 0.176] | 0.330 [0.316, 0.345] | 0.774 [0.767, 0.781] |
| GC-matched: Δ pydreg − PINTS |  | 0.094 [0.070, 0.119] p<1e-03 | 0.072 [0.055, 0.090] p<1e-03 | 0.060 [0.045, 0.077] p<1e-03 | 0.014 [0.009, 0.020] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.121 [0.093, 0.144] p<1e-03 | 0.093 [0.071, 0.114] p<1e-03 | 0.074 [0.056, 0.093] p<1e-03 | 0.023 [0.016, 0.030] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.024 [0.008, 0.044] p=0.004 | 0.023 [0.011, 0.036] p=0.002 | 0.006 [0.002, 0.011] p=0.012 | -0.012 [-0.019, -0.005] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.453 | 0.689 | 0.310 | 0.676 |
| 20,000 | 0.668 | 0.598 | 0.445 | 0.628 |
| 40,000 | 0.831 | 0.487 | 0.590 | 0.576 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 15,613, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.159**; CpG-island rows 0.126 (n_pos 11,892), non-CpG-island rows 0.043 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.178 | 0.543 [0.522, 0.564] | 0.493 [0.476, 0.510] | 0.626 [0.615, 0.639] | 0.886 [0.880, 0.891] |
| pydreg, −log10 p over all candidates | 0.106 | 0.495 [0.470, 0.521] | 0.444 [0.426, 0.462] | 0.576 [0.563, 0.588] | 0.880 [0.875, 0.887] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.106 | 0.465 [0.437, 0.492] | 0.423 [0.401, 0.445] | 0.566 [0.547, 0.582] | 0.881 [0.875, 0.886] |
| Δ pydreg(SVR) − PINTS | 0.072 | 0.078 [0.062, 0.094] p<1e-03 | 0.070 [0.058, 0.083] p<1e-03 | 0.061 [0.048, 0.075] p<1e-03 | 0.005 [0.002, 0.008] p=0.004 |
| Δ pydreg(−log p) − PINTS | 0.000 | 0.031 [0.013, 0.047] p=0.002 | 0.020 [0.009, 0.033] p=0.002 | 0.010 [-0.002, 0.023] p=0.098 | -0.000 [-0.003, 0.003] p=0.944 |
| GC-matched: pydreg (SVR) |  | 0.255 [0.236, 0.268] | 0.233 [0.219, 0.246] | 0.397 [0.378, 0.416] | 0.791 [0.785, 0.797] |
| GC-matched: PINTS |  | 0.183 [0.161, 0.202] | 0.180 [0.167, 0.193] | 0.350 [0.335, 0.365] | 0.787 [0.781, 0.793] |
| GC-matched: Δ pydreg − PINTS |  | 0.072 [0.049, 0.097] p<1e-03 | 0.053 [0.037, 0.069] p<1e-03 | 0.047 [0.033, 0.062] p<1e-03 | 0.004 [-0.000, 0.009] p=0.054 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.093 [0.057, 0.122] p<1e-03 | 0.068 [0.048, 0.088] p<1e-03 | 0.060 [0.043, 0.079] p<1e-03 | 0.012 [0.006, 0.017] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.008 [-0.007, 0.024] p=0.298 | 0.010 [-0.002, 0.022] p=0.118 | 0.003 [-0.001, 0.008] p=0.182 | -0.016 [-0.022, -0.009] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.461 | 0.701 | 0.309 | 0.671 |
| 20,000 | 0.681 | 0.612 | 0.492 | 0.633 |
| 40,000 | 0.846 | 0.491 | 0.675 | 0.571 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 15,613, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.159**; CpG-island rows 0.126 (n_pos 11,892), non-CpG-island rows 0.043 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.140 | 0.512 [0.495, 0.531] | 0.459 [0.445, 0.475] | 0.604 [0.593, 0.617] | 0.896 [0.892, 0.901] |
| pydreg, −log10 p over all candidates | 0.080 | 0.421 [0.399, 0.444] | 0.389 [0.375, 0.404] | 0.548 [0.536, 0.560] | 0.888 [0.884, 0.893] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.122 | 0.503 [0.473, 0.525] | 0.444 [0.423, 0.465] | 0.586 [0.566, 0.604] | 0.892 [0.887, 0.897] |
| Δ pydreg(SVR) − PINTS | 0.018 | 0.009 [-0.007, 0.030] p=0.298 | 0.014 [0.001, 0.030] p=0.038 | 0.018 [0.003, 0.035] p=0.012 | 0.004 [0.001, 0.007] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.042 | -0.082 [-0.101, -0.060] p<1e-03 | -0.055 [-0.070, -0.040] p<1e-03 | -0.039 [-0.053, -0.022] p<1e-03 | -0.004 [-0.007, -0.001] p=0.006 |
| GC-matched: pydreg (SVR) |  | 0.185 [0.168, 0.203] | 0.178 [0.163, 0.192] | 0.350 [0.331, 0.367] | 0.785 [0.779, 0.791] |
| GC-matched: PINTS |  | 0.173 [0.154, 0.190] | 0.166 [0.152, 0.181] | 0.341 [0.326, 0.357] | 0.788 [0.781, 0.795] |
| GC-matched: Δ pydreg − PINTS |  | 0.012 [-0.011, 0.037] p=0.31 | 0.011 [-0.009, 0.032] p=0.278 | 0.009 [-0.007, 0.025] p=0.282 | -0.004 [-0.009, 0.002] p=0.238 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.013 [-0.018, 0.044] p=0.394 | 0.011 [-0.012, 0.034] p=0.384 | 0.009 [-0.007, 0.027] p=0.284 | -0.001 [-0.008, 0.006] p=0.722 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.015 [-0.005, 0.036] p=0.164 | 0.011 [-0.002, 0.026] p=0.116 | 0.003 [-0.002, 0.008] p=0.232 | -0.009 [-0.015, -0.004] p=0.002 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.438 | 0.680 | 0.388 | 0.690 |
| 20,000 | 0.679 | 0.592 | 0.602 | 0.610 |
| 40,000 | 0.859 | 0.472 | 0.761 | 0.520 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 15,613, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.159**; CpG-island rows 0.126 (n_pos 11,892), non-CpG-island rows 0.043 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.190 | 0.543 [0.522, 0.568] | 0.492 [0.475, 0.510] | 0.626 [0.614, 0.638] | 0.883 [0.877, 0.888] |
| pydreg, −log10 p over all candidates | 0.107 | 0.509 [0.486, 0.533] | 0.448 [0.431, 0.466] | 0.576 [0.564, 0.587] | 0.878 [0.872, 0.884] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.110 | 0.474 [0.439, 0.493] | 0.414 [0.390, 0.436] | 0.552 [0.533, 0.569] | 0.877 [0.872, 0.882] |
| Δ pydreg(SVR) − PINTS | 0.080 | 0.069 [0.058, 0.092] p<1e-03 | 0.078 [0.066, 0.094] p<1e-03 | 0.074 [0.061, 0.090] p<1e-03 | 0.006 [0.002, 0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.003 | 0.035 [0.023, 0.057] p<1e-03 | 0.034 [0.023, 0.050] p<1e-03 | 0.024 [0.012, 0.038] p<1e-03 | 0.001 [-0.002, 0.005] p=0.458 |
| GC-matched: pydreg (SVR) |  | 0.218 [0.196, 0.234] | 0.204 [0.191, 0.219] | 0.375 [0.358, 0.393] | 0.786 [0.779, 0.792] |
| GC-matched: PINTS |  | 0.180 [0.163, 0.199] | 0.175 [0.162, 0.188] | 0.339 [0.324, 0.353] | 0.779 [0.772, 0.786] |
| GC-matched: Δ pydreg − PINTS |  | 0.038 [0.009, 0.060] p=0.01 | 0.029 [0.012, 0.046] p<1e-03 | 0.037 [0.022, 0.052] p<1e-03 | 0.007 [0.001, 0.013] p=0.018 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.039 [0.007, 0.065] p=0.02 | 0.033 [0.012, 0.054] p=0.002 | 0.042 [0.024, 0.061] p<1e-03 | 0.013 [0.006, 0.020] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.020 [0.002, 0.037] p=0.024 | 0.022 [0.009, 0.034] p<1e-03 | 0.007 [0.003, 0.011] p=0.002 | -0.011 [-0.018, -0.003] p=0.006 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.457 | 0.706 | 0.275 | 0.676 |
| 20,000 | 0.672 | 0.611 | 0.347 | 0.648 |
| 40,000 | 0.835 | 0.489 | 0.520 | 0.609 |

### Jurkat_CAGE (Jurkat CAGE, ood)

n_pos = 15,613, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.159**; CpG-island rows 0.126 (n_pos 11,892), non-CpG-island rows 0.043 (n_pos 3,721).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.185 | 0.559 [0.547, 0.570] | 0.457 [0.445, 0.471] | 0.494 [0.483, 0.505] | 0.758 [0.752, 0.764] |
| pydreg, −log10 p over all candidates | 0.169 | 0.559 [0.547, 0.570] | 0.453 [0.440, 0.466] | 0.479 [0.468, 0.491] | 0.757 [0.751, 0.763] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.190 | 0.361 [0.351, 0.369] | 0.325 [0.316, 0.334] | 0.365 [0.356, 0.375] | 0.660 [0.655, 0.664] |
| Δ pydreg(SVR) − PINTS | -0.005 | 0.198 [0.190, 0.207] p<1e-03 | 0.132 [0.123, 0.140] p<1e-03 | 0.129 [0.121, 0.137] p<1e-03 | 0.098 [0.094, 0.102] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.021 | 0.198 [0.190, 0.207] p<1e-03 | 0.128 [0.119, 0.136] p<1e-03 | 0.115 [0.106, 0.124] p<1e-03 | 0.097 [0.093, 0.102] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.177 [0.162, 0.188] | 0.163 [0.153, 0.173] | 0.279 [0.265, 0.292] | 0.674 [0.668, 0.680] |
| GC-matched: PINTS |  | 0.196 [0.180, 0.211] | 0.183 [0.172, 0.194] | 0.252 [0.241, 0.264] | 0.622 [0.617, 0.627] |
| GC-matched: Δ pydreg − PINTS |  | -0.019 [-0.032, -0.007] p=0.004 | -0.020 [-0.027, -0.011] p<1e-03 | 0.026 [0.016, 0.037] p<1e-03 | 0.052 [0.048, 0.057] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.022 [-0.040, -0.007] p=0.006 | -0.022 [-0.032, -0.010] p<1e-03 | 0.027 [0.015, 0.040] p<1e-03 | 0.064 [0.059, 0.070] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.025 [-0.038, -0.013] p<1e-03 | -0.009 [-0.017, -0.000] p=0.05 | 0.001 [-0.002, 0.003] p=0.536 | 0.018 [0.012, 0.025] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.385 | 0.732 | 0.276 | 0.753 |
| 20,000 | 0.527 | 0.701 | 0.345 | 0.667 |
| 40,000 | 0.558 | 0.691 | 0.350 | 0.663 |


## dREG auxiliary labels: dhs_h3k27ac_distal vs dhs_unmarked

* **Positives:** As dhs_h3k27ac, but > 1 kb from any GENCODE v24 TSS, i.e. putative enhancers (K562 20,553; GM12878 28,215).
* **Negatives:** **Diagnostic only, not a headline contrast.** Duke ∩ UW DHSs with no H3K27ac and outside ChromHMM states 1–7 (K562 23,291, GM12878 15,556; CTCF/insulator-rich). Lacking H3K27ac does not mean untranscribed, so a call here may be a true weak element. Results on this set cannot distinguish specificity from sensitivity.
* **Notes:** Not in either dREG paper as a ROC: the papers report recall only, so the negatives are ours. No nascent or TSS data is used, so these labels are not circular for any library, GRO-cap included. H3K27ac peaks and ChromHMM segments are broad, which favours wider calls under the overlap rule. On hg38 cells (HCT116, A549) the DHS is one ENCODE DNase peak set, not Duke ∩ UW, and there is no ChromHMM: quiet excludes only DNase and H3K27ac. Jurkat (hg19) uses UW DNase rep1 ∩ rep2 with H3K27ac from either ChIP-Atlas replicate (GSM4551183/4), and has no ChromHMM.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 20,553, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.469**; CpG-island rows 0.052 (n_pos 1,282), non-CpG-island rows 0.453 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.386 | 0.636 [0.619, 0.652] | 0.585 [0.571, 0.600] | 0.817 [0.809, 0.825] | 0.826 [0.820, 0.832] |
| pydreg, −log10 p over all candidates | 0.341 | 0.631 [0.617, 0.647] | 0.572 [0.556, 0.587] | 0.809 [0.800, 0.817] | 0.825 [0.819, 0.831] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.091 | 0.341 [0.312, 0.371] | 0.325 [0.303, 0.348] | 0.793 [0.781, 0.806] | 0.844 [0.837, 0.850] |
| Δ pydreg(SVR) − PINTS | 0.295 | 0.295 [0.267, 0.322] p<1e-03 | 0.260 [0.238, 0.283] p<1e-03 | 0.024 [0.013, 0.034] p<1e-03 | -0.018 [-0.024, -0.011] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.250 | 0.290 [0.263, 0.318] p<1e-03 | 0.247 [0.225, 0.269] p<1e-03 | 0.016 [0.006, 0.026] p<1e-03 | -0.019 [-0.025, -0.012] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.644 [0.629, 0.660] | 0.593 [0.580, 0.608] | 0.818 [0.810, 0.826] | 0.828 [0.821, 0.834] |
| GC-matched: PINTS |  | 0.351 [0.324, 0.379] | 0.334 [0.312, 0.356] | 0.797 [0.785, 0.809] | 0.847 [0.841, 0.853] |
| GC-matched: Δ pydreg − PINTS |  | 0.292 [0.266, 0.320] p<1e-03 | 0.260 [0.238, 0.282] p<1e-03 | 0.022 [0.012, 0.031] p<1e-03 | -0.019 [-0.025, -0.014] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.225 [0.182, 0.279] p<1e-03 | 0.212 [0.173, 0.251] p<1e-03 | 0.215 [0.168, 0.261] p<1e-03 | 0.011 [-0.001, 0.022] p=0.072 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.299 [0.269, 0.325] p<1e-03 | 0.263 [0.242, 0.286] p<1e-03 | 0.026 [0.016, 0.037] p<1e-03 | -0.022 [-0.027, -0.016] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.115 | 0.994 | 0.034 | 0.909 |
| 20,000 | 0.225 | 0.987 | 0.069 | 0.893 |
| 40,000 | 0.389 | 0.971 | 0.126 | 0.879 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 20,553, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.469**; CpG-island rows 0.052 (n_pos 1,282), non-CpG-island rows 0.453 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.324 | 0.395 [0.385, 0.405] | 0.381 [0.371, 0.392] | 0.661 [0.652, 0.671] | 0.680 [0.675, 0.686] |
| pydreg, −log10 p over all candidates | 0.323 | 0.395 [0.385, 0.405] | 0.381 [0.371, 0.391] | 0.660 [0.650, 0.669] | 0.680 [0.675, 0.685] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.157 | 0.447 [0.421, 0.470] | 0.375 [0.362, 0.388] | 0.680 [0.669, 0.691] | 0.711 [0.705, 0.716] |
| Δ pydreg(SVR) − PINTS | 0.166 | -0.051 [-0.074, -0.027] p<1e-03 | 0.007 [-0.005, 0.018] p=0.242 | -0.019 [-0.024, -0.014] p<1e-03 | -0.030 [-0.034, -0.027] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.166 | -0.051 [-0.074, -0.027] p<1e-03 | 0.006 [-0.005, 0.018] p=0.3 | -0.021 [-0.025, -0.016] p<1e-03 | -0.031 [-0.034, -0.027] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.396 [0.386, 0.406] | 0.383 [0.373, 0.393] | 0.662 [0.652, 0.671] | 0.681 [0.675, 0.686] |
| GC-matched: PINTS |  | 0.456 [0.432, 0.475] | 0.379 [0.365, 0.391] | 0.681 [0.670, 0.691] | 0.712 [0.706, 0.717] |
| GC-matched: Δ pydreg − PINTS |  | -0.060 [-0.077, -0.036] p<1e-03 | 0.004 [-0.006, 0.016] p=0.41 | -0.020 [-0.024, -0.015] p<1e-03 | -0.031 [-0.035, -0.027] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.052 [0.023, 0.092] p<1e-03 | 0.087 [0.059, 0.112] p<1e-03 | 0.159 [0.123, 0.193] p<1e-03 | 0.009 [-0.003, 0.021] p=0.126 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.071 [-0.084, -0.047] p<1e-03 | -0.001 [-0.012, 0.010] p=0.812 | -0.021 [-0.026, -0.017] p<1e-03 | -0.034 [-0.037, -0.030] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.079 | 0.995 | 0.043 | 0.965 |
| 20,000 | 0.177 | 0.988 | 0.088 | 0.949 |
| 40,000 | 0.302 | 0.972 | 0.153 | 0.934 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 20,553, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.469**; CpG-island rows 0.052 (n_pos 1,282), non-CpG-island rows 0.453 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.358 | 0.430 [0.420, 0.440] | 0.413 [0.403, 0.422] | 0.679 [0.671, 0.687] | 0.698 [0.693, 0.703] |
| pydreg, −log10 p over all candidates | 0.354 | 0.430 [0.420, 0.440] | 0.412 [0.403, 0.422] | 0.678 [0.669, 0.685] | 0.698 [0.693, 0.703] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.126 | 0.432 [0.404, 0.460] | 0.382 [0.366, 0.398] | 0.701 [0.692, 0.711] | 0.738 [0.733, 0.743] |
| Δ pydreg(SVR) − PINTS | 0.232 | -0.002 [-0.031, 0.027] p=0.888 | 0.031 [0.014, 0.048] p<1e-03 | -0.023 [-0.028, -0.017] p<1e-03 | -0.040 [-0.044, -0.036] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.228 | -0.002 [-0.031, 0.027] p=0.888 | 0.030 [0.014, 0.047] p<1e-03 | -0.024 [-0.029, -0.019] p<1e-03 | -0.040 [-0.044, -0.036] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.431 [0.421, 0.440] | 0.416 [0.406, 0.426] | 0.680 [0.672, 0.687] | 0.699 [0.694, 0.704] |
| GC-matched: PINTS |  | 0.453 [0.428, 0.482] | 0.393 [0.378, 0.409] | 0.704 [0.696, 0.713] | 0.741 [0.735, 0.746] |
| GC-matched: Δ pydreg − PINTS |  | -0.022 [-0.052, 0.007] p=0.122 | 0.023 [0.007, 0.039] p=0.004 | -0.025 [-0.030, -0.019] p<1e-03 | -0.042 [-0.045, -0.038] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.144 [0.105, 0.189] p<1e-03 | 0.119 [0.088, 0.150] p<1e-03 | 0.135 [0.092, 0.177] p<1e-03 | 0.007 [-0.004, 0.017] p=0.244 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.039 [-0.071, -0.010] p=0.004 | 0.015 [-0.000, 0.031] p=0.06 | -0.026 [-0.031, -0.020] p<1e-03 | -0.045 [-0.049, -0.041] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.062 | 0.994 | 0.027 | 0.963 |
| 20,000 | 0.150 | 0.988 | 0.057 | 0.946 |
| 40,000 | 0.301 | 0.976 | 0.130 | 0.917 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 28,215, n_neg = 15,556, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.645**; CpG-island rows 0.091 (n_pos 1,556), non-CpG-island rows 0.632 (n_pos 26,659).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.398 | 0.697 [0.679, 0.712] | 0.632 [0.618, 0.645] | 0.901 [0.896, 0.905] | 0.857 [0.853, 0.861] |
| pydreg, −log10 p over all candidates | 0.327 | 0.696 [0.676, 0.709] | 0.616 [0.600, 0.632] | 0.896 [0.891, 0.901] | 0.855 [0.851, 0.859] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.181 | 0.529 [0.505, 0.562] | 0.483 [0.462, 0.505] | 0.911 [0.906, 0.915] | 0.878 [0.874, 0.882] |
| Δ pydreg(SVR) − PINTS | 0.217 | 0.168 [0.139, 0.188] p<1e-03 | 0.148 [0.131, 0.166] p<1e-03 | -0.010 [-0.013, -0.007] p<1e-03 | -0.021 [-0.025, -0.018] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.146 | 0.166 [0.137, 0.185] p<1e-03 | 0.133 [0.115, 0.151] p<1e-03 | -0.015 [-0.018, -0.011] p<1e-03 | -0.023 [-0.027, -0.019] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.708 [0.693, 0.721] | 0.641 [0.628, 0.654] | 0.902 [0.897, 0.906] | 0.858 [0.854, 0.862] |
| GC-matched: PINTS |  | 0.558 [0.525, 0.580] | 0.501 [0.481, 0.523] | 0.914 [0.909, 0.918] | 0.882 [0.877, 0.886] |
| GC-matched: Δ pydreg − PINTS |  | 0.150 [0.129, 0.178] p<1e-03 | 0.140 [0.122, 0.157] p<1e-03 | -0.012 [-0.015, -0.009] p<1e-03 | -0.023 [-0.027, -0.020] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.130 [0.088, 0.173] p<1e-03 | 0.107 [0.068, 0.148] p<1e-03 | 0.073 [0.026, 0.121] p=0.002 | -0.010 [-0.019, -0.001] p=0.018 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.152 [0.128, 0.176] p<1e-03 | 0.140 [0.120, 0.158] p<1e-03 | -0.012 [-0.015, -0.008] p<1e-03 | -0.025 [-0.028, -0.021] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.092 | 0.995 | 0.038 | 0.978 |
| 20,000 | 0.211 | 0.990 | 0.090 | 0.981 |
| 40,000 | 0.407 | 0.985 | 0.174 | 0.971 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 20,553, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.469**; CpG-island rows 0.052 (n_pos 1,282), non-CpG-island rows 0.453 (n_pos 19,271).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.442 | 0.704 [0.683, 0.722] | 0.635 [0.620, 0.649] | 0.832 [0.824, 0.839] | 0.842 [0.835, 0.848] |
| pydreg, −log10 p over all candidates | 0.446 | 0.701 [0.682, 0.721] | 0.626 [0.611, 0.642] | 0.823 [0.814, 0.830] | 0.841 [0.834, 0.847] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.479 | 0.694 [0.678, 0.710] | 0.651 [0.636, 0.666] | 0.848 [0.841, 0.855] | 0.856 [0.851, 0.862] |
| Δ pydreg(SVR) − PINTS | -0.037 | 0.010 [-0.002, 0.021] p=0.11 | -0.016 [-0.023, -0.009] p<1e-03 | -0.016 [-0.019, -0.013] p<1e-03 | -0.015 [-0.017, -0.012] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.033 | 0.007 [-0.003, 0.020] p=0.172 | -0.025 [-0.032, -0.017] p<1e-03 | -0.026 [-0.030, -0.022] p<1e-03 | -0.015 [-0.018, -0.013] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.712 [0.692, 0.725] | 0.640 [0.626, 0.654] | 0.833 [0.825, 0.840] | 0.843 [0.836, 0.848] |
| GC-matched: PINTS |  | 0.700 [0.684, 0.716] | 0.656 [0.641, 0.670] | 0.849 [0.842, 0.856] | 0.857 [0.851, 0.863] |
| GC-matched: Δ pydreg − PINTS |  | 0.012 [0.000, 0.021] p=0.05 | -0.016 [-0.022, -0.009] p<1e-03 | -0.016 [-0.019, -0.013] p<1e-03 | -0.015 [-0.017, -0.012] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.000 [-0.022, 0.021] p=0.894 | -0.011 [-0.031, 0.010] p=0.342 | -0.000 [-0.039, 0.041] p=0.99 | 0.007 [-0.001, 0.016] p=0.086 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.009 [-0.002, 0.018] p=0.134 | -0.018 [-0.024, -0.011] p<1e-03 | -0.018 [-0.021, -0.015] p<1e-03 | -0.016 [-0.019, -0.013] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.140 | 0.996 | 0.062 | 0.993 |
| 20,000 | 0.268 | 0.990 | 0.169 | 0.990 |
| 40,000 | 0.469 | 0.971 | 0.369 | 0.985 |

### HCT116_PROseq (HCT116 PRO-seq, new_cell)

n_pos = 16,069, n_neg = 34,256, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.319**; CpG-island rows 0.033 (n_pos 1,175), non-CpG-island rows 0.303 (n_pos 14,894).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.028 | 0.186 [0.171, 0.201] | 0.192 [0.181, 0.202] | 0.431 [0.414, 0.446] | 0.637 [0.631, 0.643] |
| pydreg, −log10 p over all candidates | 0.037 | 0.188 [0.175, 0.205] | 0.196 [0.185, 0.206] | 0.435 [0.419, 0.450] | 0.638 [0.631, 0.643] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.039 | 0.242 [0.226, 0.258] | 0.232 [0.220, 0.243] | 0.475 [0.457, 0.492] | 0.668 [0.660, 0.675] |
| Δ pydreg(SVR) − PINTS | -0.012 | -0.056 [-0.070, -0.043] p<1e-03 | -0.040 [-0.048, -0.032] p<1e-03 | -0.044 [-0.049, -0.037] p<1e-03 | -0.031 [-0.035, -0.026] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.002 | -0.053 [-0.066, -0.039] p<1e-03 | -0.036 [-0.044, -0.028] p<1e-03 | -0.039 [-0.046, -0.033] p<1e-03 | -0.030 [-0.034, -0.026] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.247 [0.227, 0.269] | 0.228 [0.217, 0.239] | 0.451 [0.436, 0.465] | 0.645 [0.639, 0.650] |
| GC-matched: PINTS |  | 0.283 [0.268, 0.301] | 0.264 [0.253, 0.275] | 0.493 [0.477, 0.509] | 0.674 [0.667, 0.680] |
| GC-matched: Δ pydreg − PINTS |  | -0.035 [-0.051, -0.020] p<1e-03 | -0.036 [-0.044, -0.028] p<1e-03 | -0.042 [-0.048, -0.035] p<1e-03 | -0.029 [-0.033, -0.025] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.026 [-0.067, 0.016] p=0.258 | -0.028 [-0.057, 0.002] p=0.07 | -0.011 [-0.026, 0.005] p=0.178 | -0.027 [-0.043, -0.011] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.003 [-0.023, 0.007] p=0.386 | -0.035 [-0.042, -0.027] p<1e-03 | -0.042 [-0.048, -0.035] p<1e-03 | -0.030 [-0.034, -0.025] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.061 | 0.571 | 0.023 | 0.650 |
| 20,000 | 0.162 | 0.628 | 0.049 | 0.662 |
| 40,000 | 0.333 | 0.665 | 0.102 | 0.671 |

### HCT116_PROcap (HCT116 PRO-cap, ood)

n_pos = 16,069, n_neg = 34,256, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.319**; CpG-island rows 0.033 (n_pos 1,175), non-CpG-island rows 0.303 (n_pos 14,894).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.027 | 0.165 [0.154, 0.176] | 0.170 [0.160, 0.180] | 0.553 [0.539, 0.566] | 0.778 [0.772, 0.783] |
| pydreg, −log10 p over all candidates | 0.026 | 0.161 [0.150, 0.172] | 0.166 [0.156, 0.176] | 0.552 [0.537, 0.565] | 0.778 [0.772, 0.784] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.017 | 0.169 [0.155, 0.181] | 0.180 [0.168, 0.191] | 0.570 [0.554, 0.584] | 0.801 [0.796, 0.806] |
| Δ pydreg(SVR) − PINTS | 0.009 | -0.003 [-0.018, 0.010] p=0.588 | -0.009 [-0.019, 0.001] p=0.08 | -0.017 [-0.025, -0.010] p<1e-03 | -0.023 [-0.026, -0.020] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.009 | -0.007 [-0.022, 0.006] p=0.298 | -0.014 [-0.024, -0.002] p=0.004 | -0.018 [-0.026, -0.010] p<1e-03 | -0.023 [-0.026, -0.020] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.199 [0.184, 0.220] | 0.201 [0.187, 0.216] | 0.575 [0.562, 0.590] | 0.786 [0.780, 0.793] |
| GC-matched: PINTS |  | 0.226 [0.202, 0.250] | 0.221 [0.204, 0.238] | 0.597 [0.582, 0.612] | 0.810 [0.805, 0.816] |
| GC-matched: Δ pydreg − PINTS |  | -0.027 [-0.041, -0.009] p<1e-03 | -0.020 [-0.030, -0.009] p<1e-03 | -0.021 [-0.029, -0.013] p<1e-03 | -0.024 [-0.027, -0.021] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.021 [-0.056, 0.009] p=0.17 | -0.030 [-0.057, -0.001] p=0.042 | -0.014 [-0.022, -0.005] p=0.002 | -0.058 [-0.073, -0.043] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.038 [-0.055, -0.023] p<1e-03 | -0.027 [-0.038, -0.016] p<1e-03 | -0.022 [-0.030, -0.014] p<1e-03 | -0.024 [-0.028, -0.021] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.139 | 0.604 | 0.045 | 0.497 |
| 20,000 | 0.285 | 0.623 | 0.155 | 0.607 |
| 40,000 | 0.511 | 0.625 | 0.395 | 0.655 |

### A549_PROseq (A549 PRO-seq, new_cell)

n_pos = 49,656, n_neg = 57,867, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.462**; CpG-island rows 0.055 (n_pos 3,359), non-CpG-island rows 0.444 (n_pos 46,297).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.224 | 0.260 [0.254, 0.266] | 0.252 [0.246, 0.258] | 0.578 [0.568, 0.588] | 0.610 [0.607, 0.613] |
| pydreg, −log10 p over all candidates | 0.222 | 0.260 [0.254, 0.266] | 0.252 [0.246, 0.258] | 0.577 [0.567, 0.587] | 0.610 [0.607, 0.613] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.133 | 0.332 [0.324, 0.341] | 0.290 [0.282, 0.297] | 0.606 [0.594, 0.617] | 0.644 [0.640, 0.649] |
| Δ pydreg(SVR) − PINTS | 0.091 | -0.072 [-0.078, -0.067] p<1e-03 | -0.038 [-0.043, -0.032] p<1e-03 | -0.028 [-0.031, -0.025] p<1e-03 | -0.035 [-0.037, -0.032] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.089 | -0.072 [-0.078, -0.067] p<1e-03 | -0.038 [-0.044, -0.032] p<1e-03 | -0.029 [-0.032, -0.025] p<1e-03 | -0.035 [-0.037, -0.032] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.259 [0.253, 0.265] | 0.249 [0.243, 0.255] | 0.576 [0.566, 0.586] | 0.609 [0.606, 0.612] |
| GC-matched: PINTS |  | 0.329 [0.321, 0.338] | 0.280 [0.273, 0.288] | 0.602 [0.590, 0.612] | 0.642 [0.638, 0.646] |
| GC-matched: Δ pydreg − PINTS |  | -0.070 [-0.075, -0.065] p<1e-03 | -0.031 [-0.037, -0.025] p<1e-03 | -0.026 [-0.029, -0.022] p<1e-03 | -0.033 [-0.036, -0.030] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.075 [0.044, 0.122] p<1e-03 | 0.066 [0.041, 0.093] p<1e-03 | 0.087 [0.050, 0.128] p<1e-03 | -0.008 [-0.016, 0.000] p=0.054 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.073 [-0.079, -0.068] p<1e-03 | -0.038 [-0.044, -0.032] p<1e-03 | -0.028 [-0.032, -0.024] p<1e-03 | -0.035 [-0.038, -0.032] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.028 | 0.975 | 0.020 | 0.939 |
| 20,000 | 0.085 | 0.976 | 0.043 | 0.944 |
| 40,000 | 0.199 | 0.957 | 0.142 | 0.918 |

### A549_PROcap (A549 PRO-cap, ood)

n_pos = 49,656, n_neg = 57,867, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.462**; CpG-island rows 0.055 (n_pos 3,359), non-CpG-island rows 0.444 (n_pos 46,297).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.044 | 0.148 [0.142, 0.154] | 0.130 [0.125, 0.135] | 0.510 [0.499, 0.521] | 0.551 [0.549, 0.554] |
| pydreg, −log10 p over all candidates | 0.041 | 0.148 [0.142, 0.154] | 0.129 [0.125, 0.134] | 0.507 [0.496, 0.518] | 0.551 [0.548, 0.554] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.047 | 0.150 [0.144, 0.156] | 0.137 [0.132, 0.143] | 0.515 [0.503, 0.526] | 0.559 [0.556, 0.563] |
| Δ pydreg(SVR) − PINTS | -0.003 | -0.002 [-0.006, 0.002] p=0.324 | -0.007 [-0.010, -0.005] p<1e-03 | -0.005 [-0.007, -0.004] p<1e-03 | -0.008 [-0.009, -0.006] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.006 | -0.002 [-0.006, 0.002] p=0.268 | -0.008 [-0.010, -0.005] p<1e-03 | -0.008 [-0.010, -0.006] p<1e-03 | -0.008 [-0.009, -0.006] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.136 [0.131, 0.142] | 0.122 [0.118, 0.126] | 0.504 [0.493, 0.516] | 0.548 [0.546, 0.551] |
| GC-matched: PINTS |  | 0.139 [0.133, 0.144] | 0.128 [0.123, 0.133] | 0.509 [0.497, 0.520] | 0.555 [0.552, 0.559] |
| GC-matched: Δ pydreg − PINTS |  | -0.002 [-0.006, 0.001] p=0.186 | -0.006 [-0.009, -0.003] p<1e-03 | -0.004 [-0.006, -0.003] p<1e-03 | -0.007 [-0.009, -0.006] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.025 [-0.047, -0.003] p=0.032 | -0.001 [-0.015, 0.013] p=0.896 | 0.005 [-0.007, 0.019] p=0.386 | -0.002 [-0.010, 0.004] p=0.554 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.004 [-0.007, 0.000] p=0.062 | -0.007 [-0.010, -0.005] p<1e-03 | -0.005 [-0.006, -0.003] p<1e-03 | -0.008 [-0.009, -0.006] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.020 | 0.823 | 0.013 | 0.779 |
| 20,000 | 0.047 | 0.786 | 0.034 | 0.800 |
| 40,000 | 0.099 | 0.752 | 0.076 | 0.769 |

### Jurkat_PROseq (Jurkat PRO-seq, new_cell)

n_pos = 2,600, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.031**; CpG-island rows 0.014 (n_pos 1,212), non-CpG-island rows 0.017 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.043 | 0.185 [0.167, 0.205] | 0.183 [0.169, 0.200] | 0.073 [0.067, 0.081] | 0.701 [0.689, 0.711] |
| pydreg, −log10 p over all candidates | 0.034 | 0.173 [0.156, 0.192] | 0.173 [0.159, 0.188] | 0.069 [0.063, 0.075] | 0.699 [0.688, 0.710] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.027 | 0.156 [0.138, 0.174] | 0.156 [0.142, 0.171] | 0.067 [0.062, 0.074] | 0.719 [0.708, 0.729] |
| Δ pydreg(SVR) − PINTS | 0.015 | 0.029 [0.014, 0.046] p<1e-03 | 0.027 [0.015, 0.040] p<1e-03 | 0.006 [0.001, 0.011] p=0.006 | -0.018 [-0.028, -0.009] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.007 | 0.017 [-0.000, 0.034] p=0.056 | 0.017 [0.005, 0.029] p=0.012 | 0.001 [-0.003, 0.005] p=0.496 | -0.019 [-0.029, -0.010] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.076 [0.066, 0.089] | 0.076 [0.067, 0.086] | 0.041 [0.038, 0.047] | 0.612 [0.602, 0.623] |
| GC-matched: PINTS |  | 0.050 [0.038, 0.062] | 0.054 [0.045, 0.064] | 0.039 [0.036, 0.043] | 0.619 [0.608, 0.629] |
| GC-matched: Δ pydreg − PINTS |  | 0.026 [0.011, 0.043] p<1e-03 | 0.022 [0.010, 0.033] p<1e-03 | 0.002 [-0.000, 0.005] p=0.052 | -0.007 [-0.017, 0.002] p=0.152 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.043 [0.018, 0.064] p<1e-03 | 0.034 [0.014, 0.052] p<1e-03 | 0.003 [-0.000, 0.006] p=0.062 | 0.008 [-0.006, 0.022] p=0.26 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.009 [-0.006, 0.025] p=0.244 | 0.013 [0.002, 0.025] p=0.018 | 0.000 [-0.001, 0.003] p=0.734 | -0.020 [-0.034, -0.008] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.148 | 0.107 | 0.092 | 0.093 |
| 20,000 | 0.300 | 0.100 | 0.156 | 0.090 |
| 40,000 | 0.507 | 0.088 | 0.254 | 0.089 |

### Jurkat_ChROseq1 (Jurkat ChRO-seq, new_cell)

n_pos = 2,600, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.031**; CpG-island rows 0.014 (n_pos 1,212), non-CpG-island rows 0.017 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.055 | 0.208 [0.192, 0.230] | 0.208 [0.193, 0.224] | 0.086 [0.078, 0.095] | 0.720 [0.708, 0.731] |
| pydreg, −log10 p over all candidates | 0.039 | 0.198 [0.179, 0.217] | 0.195 [0.181, 0.209] | 0.077 [0.071, 0.084] | 0.718 [0.707, 0.730] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.027 | 0.160 [0.141, 0.181] | 0.164 [0.150, 0.179] | 0.071 [0.065, 0.078] | 0.727 [0.717, 0.738] |
| Δ pydreg(SVR) − PINTS | 0.028 | 0.048 [0.031, 0.069] p<1e-03 | 0.044 [0.032, 0.055] p<1e-03 | 0.015 [0.010, 0.021] p<1e-03 | -0.007 [-0.017, 0.002] p=0.122 |
| Δ pydreg(−log p) − PINTS | 0.012 | 0.038 [0.017, 0.057] p<1e-03 | 0.031 [0.018, 0.043] p<1e-03 | 0.006 [0.002, 0.010] p=0.008 | -0.009 [-0.018, 0.000] p=0.056 |
| GC-matched: pydreg (SVR) |  | 0.086 [0.075, 0.098] | 0.083 [0.074, 0.094] | 0.045 [0.040, 0.051] | 0.624 [0.614, 0.635] |
| GC-matched: PINTS |  | 0.062 [0.049, 0.071] | 0.059 [0.050, 0.067] | 0.041 [0.037, 0.045] | 0.628 [0.617, 0.638] |
| GC-matched: Δ pydreg − PINTS |  | 0.025 [0.011, 0.041] p<1e-03 | 0.025 [0.013, 0.037] p<1e-03 | 0.005 [0.002, 0.009] p=0.004 | -0.003 [-0.013, 0.005] p=0.492 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.046 [0.015, 0.074] p=0.004 | 0.039 [0.017, 0.061] p<1e-03 | 0.006 [0.002, 0.012] p<1e-03 | 0.012 [-0.002, 0.025] p=0.088 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.007 [-0.007, 0.023] p=0.266 | 0.009 [-0.002, 0.021] p=0.102 | 0.000 [-0.001, 0.004] p=0.558 | -0.017 [-0.028, -0.006] p=0.004 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.165 | 0.123 | 0.087 | 0.087 |
| 20,000 | 0.326 | 0.112 | 0.177 | 0.093 |
| 40,000 | 0.544 | 0.093 | 0.320 | 0.095 |

### Jurkat_ChROseq2 (Jurkat ChRO-seq, new_cell)

n_pos = 2,600, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.031**; CpG-island rows 0.014 (n_pos 1,212), non-CpG-island rows 0.017 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.080 | 0.258 [0.239, 0.278] | 0.250 [0.234, 0.266] | 0.109 [0.100, 0.120] | 0.771 [0.760, 0.781] |
| pydreg, −log10 p over all candidates | 0.042 | 0.216 [0.197, 0.233] | 0.214 [0.201, 0.227] | 0.091 [0.084, 0.098] | 0.767 [0.757, 0.777] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.026 | 0.193 [0.170, 0.215] | 0.183 [0.167, 0.200] | 0.080 [0.074, 0.087] | 0.761 [0.751, 0.771] |
| Δ pydreg(SVR) − PINTS | 0.053 | 0.065 [0.042, 0.089] p<1e-03 | 0.067 [0.052, 0.082] p<1e-03 | 0.029 [0.022, 0.037] p<1e-03 | 0.010 [0.003, 0.017] p=0.014 |
| Δ pydreg(−log p) − PINTS | 0.016 | 0.023 [-0.002, 0.045] p=0.068 | 0.031 [0.016, 0.044] p<1e-03 | 0.011 [0.006, 0.015] p<1e-03 | 0.006 [-0.001, 0.013] p=0.086 |
| GC-matched: pydreg (SVR) |  | 0.112 [0.098, 0.126] | 0.108 [0.097, 0.120] | 0.053 [0.047, 0.059] | 0.665 [0.653, 0.676] |
| GC-matched: PINTS |  | 0.042 [0.033, 0.053] | 0.050 [0.041, 0.059] | 0.043 [0.039, 0.046] | 0.652 [0.641, 0.663] |
| GC-matched: Δ pydreg − PINTS |  | 0.070 [0.052, 0.089] p<1e-03 | 0.059 [0.044, 0.074] p<1e-03 | 0.010 [0.007, 0.015] p<1e-03 | 0.013 [0.005, 0.020] p=0.002 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.082 [0.053, 0.109] p<1e-03 | 0.074 [0.053, 0.097] p<1e-03 | 0.009 [0.005, 0.017] p<1e-03 | 0.034 [0.023, 0.046] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.040 [0.023, 0.059] p<1e-03 | 0.036 [0.022, 0.050] p<1e-03 | 0.003 [0.001, 0.005] p<1e-03 | -0.005 [-0.015, 0.004] p=0.35 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.215 | 0.148 | 0.121 | 0.103 |
| 20,000 | 0.398 | 0.124 | 0.273 | 0.105 |
| 40,000 | 0.612 | 0.096 | 0.457 | 0.098 |

### Jurkat_leChROseq (Jurkat leChRO-seq, new_cell)

n_pos = 2,600, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.031**; CpG-island rows 0.014 (n_pos 1,212), non-CpG-island rows 0.017 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.063 | 0.227 [0.207, 0.251] | 0.219 [0.203, 0.235] | 0.089 [0.081, 0.099] | 0.725 [0.714, 0.735] |
| pydreg, −log10 p over all candidates | 0.043 | 0.215 [0.196, 0.236] | 0.208 [0.193, 0.224] | 0.082 [0.075, 0.088] | 0.724 [0.713, 0.734] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.026 | 0.165 [0.147, 0.183] | 0.165 [0.150, 0.181] | 0.071 [0.066, 0.078] | 0.737 [0.727, 0.746] |
| Δ pydreg(SVR) − PINTS | 0.037 | 0.061 [0.045, 0.083] p<1e-03 | 0.054 [0.041, 0.069] p<1e-03 | 0.018 [0.012, 0.025] p<1e-03 | -0.011 [-0.022, -0.002] p=0.024 |
| Δ pydreg(−log p) − PINTS | 0.017 | 0.049 [0.032, 0.070] p<1e-03 | 0.043 [0.029, 0.057] p<1e-03 | 0.010 [0.005, 0.014] p<1e-03 | -0.012 [-0.023, -0.002] p=0.012 |
| GC-matched: pydreg (SVR) |  | 0.079 [0.066, 0.091] | 0.080 [0.070, 0.089] | 0.044 [0.040, 0.049] | 0.632 [0.621, 0.643] |
| GC-matched: PINTS |  | 0.051 [0.041, 0.061] | 0.056 [0.048, 0.064] | 0.041 [0.037, 0.045] | 0.635 [0.625, 0.645] |
| GC-matched: Δ pydreg − PINTS |  | 0.028 [0.014, 0.043] p<1e-03 | 0.024 [0.014, 0.035] p<1e-03 | 0.003 [0.001, 0.006] p=0.012 | -0.003 [-0.013, 0.006] p=0.496 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.045 [0.023, 0.065] p<1e-03 | 0.036 [0.020, 0.053] p<1e-03 | 0.004 [0.002, 0.008] p<1e-03 | 0.023 [0.009, 0.036] p=0.002 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.012 [-0.004, 0.027] p=0.16 | 0.011 [-0.000, 0.022] p=0.052 | -0.001 [-0.002, 0.001] p=0.192 | -0.027 [-0.040, -0.014] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.170 | 0.130 | 0.078 | 0.090 |
| 20,000 | 0.333 | 0.115 | 0.114 | 0.092 |
| 40,000 | 0.546 | 0.094 | 0.208 | 0.094 |

### Jurkat_CAGE (Jurkat CAGE, ood)

n_pos = 2,600, n_neg = 82,472, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.031**; CpG-island rows 0.014 (n_pos 1,212), non-CpG-island rows 0.017 (n_pos 1,388).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.007 | 0.081 [0.070, 0.094] | 0.067 [0.059, 0.076] | 0.031 [0.028, 0.033] | 0.515 [0.510, 0.521] |
| pydreg, −log10 p over all candidates | 0.010 | 0.081 [0.070, 0.094] | 0.068 [0.060, 0.077] | 0.031 [0.029, 0.034] | 0.515 [0.510, 0.521] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.002 | 0.095 [0.083, 0.106] | 0.082 [0.074, 0.091] | 0.032 [0.030, 0.035] | 0.522 [0.517, 0.528] |
| Δ pydreg(SVR) − PINTS | 0.005 | -0.013 [-0.024, -0.002] p=0.02 | -0.015 [-0.023, -0.008] p<1e-03 | -0.002 [-0.002, -0.001] p<1e-03 | -0.007 [-0.012, -0.002] p=0.012 |
| Δ pydreg(−log p) − PINTS | 0.008 | -0.013 [-0.024, -0.002] p=0.02 | -0.014 [-0.022, -0.007] p<1e-03 | -0.001 [-0.002, -0.000] p<1e-03 | -0.007 [-0.012, -0.002] p=0.016 |
| GC-matched: pydreg (SVR) |  | 0.007 [0.004, 0.011] | 0.008 [0.006, 0.011] | 0.029 [0.027, 0.031] | 0.445 [0.438, 0.452] |
| GC-matched: PINTS |  | 0.004 [0.002, 0.009] | 0.019 [0.014, 0.025] | 0.030 [0.027, 0.032] | 0.490 [0.483, 0.496] |
| GC-matched: Δ pydreg − PINTS |  | 0.003 [-0.002, 0.007] p=0.246 | -0.011 [-0.016, -0.007] p<1e-03 | -0.001 [-0.001, -0.001] p<1e-03 | -0.045 [-0.051, -0.038] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.007 [0.002, 0.012] p=0.004 | 0.001 [-0.003, 0.005] p=0.718 | -0.000 [-0.001, -0.000] p<1e-03 | -0.070 [-0.079, -0.061] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.009 [-0.020, -0.002] p=0.014 | -0.020 [-0.026, -0.014] p<1e-03 | -0.000 [-0.000, -0.000] p<1e-03 | -0.022 [-0.029, -0.015] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.026 | 0.030 | 0.018 | 0.033 |
| 20,000 | 0.059 | 0.042 | 0.077 | 0.069 |
| 40,000 | 0.078 | 0.050 | 0.079 | 0.069 |


## dREG auxiliary labels: chromhmm_active vs dhs_unmarked

* **Positives:** ENCODE Broad ChromHMM segments in states 1 (active promoter) and 4–5 (strong enhancer), as in Danko 2015 Fig 1b (K562 43,377; GM12878 41,415; median segment about 1.4 kb).
* **Negatives:** **Diagnostic only, not a headline contrast.** Duke ∩ UW DHSs with no H3K27ac and outside ChromHMM states 1–7 (K562 23,291, GM12878 15,556; CTCF/insulator-rich). Lacking H3K27ac does not mean untranscribed, so a call here may be a true weak element. Results on this set cannot distinguish specificity from sensitivity.
* **Notes:** Not in either dREG paper as a ROC: the papers report recall only, so the negatives are ours. No nascent or TSS data is used, so these labels are not circular for any library, GRO-cap included. H3K27ac peaks and ChromHMM segments are broad, which favours wider calls under the overlap rule. On hg38 cells (HCT116, A549) the DHS is one ENCODE DNase peak set, not Duke ∩ UW, and there is no ChromHMM: quiet excludes only DNase and H3K27ac. Jurkat (hg19) uses UW DNase rep1 ∩ rep2 with H3K27ac from either ChIP-Atlas replicate (GSM4551183/4), and has no ChromHMM.

### K562_PROseq_G1 (K562 PRO-seq, train)

n_pos = 43,377, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.651**; CpG-island rows 0.226 (n_pos 6,818), non-CpG-island rows 0.611 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.489 | 0.666 [0.657, 0.675] | 0.633 [0.625, 0.641] | 0.895 [0.892, 0.899] | 0.842 [0.838, 0.845] |
| pydreg, −log10 p over all candidates | 0.464 | 0.663 [0.655, 0.672] | 0.619 [0.610, 0.628] | 0.892 [0.888, 0.896] | 0.840 [0.837, 0.844] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.265 | 0.521 [0.500, 0.544] | 0.489 [0.472, 0.507] | 0.918 [0.914, 0.922] | 0.876 [0.871, 0.880] |
| Δ pydreg(SVR) − PINTS | 0.224 | 0.145 [0.123, 0.166] p<1e-03 | 0.144 [0.126, 0.160] p<1e-03 | -0.023 [-0.026, -0.019] p<1e-03 | -0.034 [-0.039, -0.030] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.199 | 0.142 [0.120, 0.164] p<1e-03 | 0.130 [0.112, 0.147] p<1e-03 | -0.026 [-0.029, -0.023] p<1e-03 | -0.035 [-0.040, -0.031] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.685 [0.675, 0.695] | 0.648 [0.641, 0.656] | 0.897 [0.893, 0.900] | 0.845 [0.841, 0.848] |
| GC-matched: PINTS |  | 0.543 [0.523, 0.562] | 0.507 [0.491, 0.524] | 0.921 [0.917, 0.925] | 0.881 [0.878, 0.885] |
| GC-matched: Δ pydreg − PINTS |  | 0.141 [0.124, 0.161] p<1e-03 | 0.141 [0.125, 0.157] p<1e-03 | -0.024 [-0.027, -0.021] p<1e-03 | -0.037 [-0.041, -0.033] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.091 [0.075, 0.106] p<1e-03 | 0.117 [0.098, 0.136] p<1e-03 | 0.045 [0.036, 0.054] p<1e-03 | 0.014 [0.011, 0.017] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.155 [0.137, 0.175] p<1e-03 | 0.147 [0.130, 0.164] p<1e-03 | -0.030 [-0.034, -0.027] p<1e-03 | -0.047 [-0.051, -0.044] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.221 | 0.999 | 0.150 | 0.989 |
| 20,000 | 0.362 | 0.996 | 0.229 | 0.983 |
| 40,000 | 0.491 | 0.989 | 0.318 | 0.975 |

### K562_GROseq_G2 (K562 GRO-seq, train)

n_pos = 43,377, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.651**; CpG-island rows 0.226 (n_pos 6,818), non-CpG-island rows 0.611 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.439 | 0.498 [0.489, 0.507] | 0.487 [0.478, 0.495] | 0.817 [0.811, 0.822] | 0.735 [0.730, 0.739] |
| pydreg, −log10 p over all candidates | 0.438 | 0.498 [0.489, 0.507] | 0.485 [0.477, 0.494] | 0.816 [0.810, 0.822] | 0.734 [0.730, 0.739] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.299 | 0.548 [0.530, 0.566] | 0.485 [0.475, 0.495] | 0.836 [0.831, 0.842] | 0.763 [0.759, 0.768] |
| Δ pydreg(SVR) − PINTS | 0.140 | -0.050 [-0.069, -0.030] p<1e-03 | 0.002 [-0.008, 0.013] p=0.728 | -0.020 [-0.022, -0.018] p<1e-03 | -0.028 [-0.031, -0.026] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.139 | -0.050 [-0.069, -0.030] p<1e-03 | 0.001 [-0.010, 0.011] p=0.946 | -0.020 [-0.022, -0.019] p<1e-03 | -0.029 [-0.031, -0.026] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.499 [0.490, 0.508] | 0.490 [0.481, 0.499] | 0.817 [0.811, 0.823] | 0.735 [0.731, 0.740] |
| GC-matched: PINTS |  | 0.567 [0.552, 0.575] | 0.494 [0.484, 0.504] | 0.838 [0.832, 0.843] | 0.765 [0.760, 0.770] |
| GC-matched: Δ pydreg − PINTS |  | -0.068 [-0.072, -0.051] p<1e-03 | -0.004 [-0.013, 0.005] p=0.366 | -0.021 [-0.022, -0.019] p<1e-03 | -0.029 [-0.032, -0.027] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.046 [0.032, 0.065] p<1e-03 | 0.099 [0.083, 0.115] p<1e-03 | 0.053 [0.044, 0.062] p<1e-03 | 0.016 [0.012, 0.020] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.084 [-0.089, -0.079] p<1e-03 | -0.023 [-0.031, -0.014] p<1e-03 | -0.028 [-0.031, -0.026] p<1e-03 | -0.038 [-0.041, -0.035] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.216 | 0.999 | 0.125 | 0.994 |
| 20,000 | 0.325 | 0.997 | 0.209 | 0.989 |
| 40,000 | 0.422 | 0.990 | 0.295 | 0.983 |

### K562_PROseq_G7 (K562 PRO-seq, holdout)

n_pos = 43,377, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.651**; CpG-island rows 0.226 (n_pos 6,818), non-CpG-island rows 0.611 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.472 | 0.529 [0.517, 0.539] | 0.516 [0.505, 0.526] | 0.828 [0.822, 0.833] | 0.751 [0.745, 0.756] |
| pydreg, −log10 p over all candidates | 0.470 | 0.529 [0.517, 0.539] | 0.515 [0.504, 0.525] | 0.827 [0.822, 0.832] | 0.751 [0.745, 0.756] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.321 | 0.555 [0.538, 0.573] | 0.516 [0.504, 0.528] | 0.858 [0.854, 0.863] | 0.793 [0.788, 0.798] |
| Δ pydreg(SVR) − PINTS | 0.151 | -0.026 [-0.049, -0.007] p=0.012 | -0.000 [-0.014, 0.012] p=0.928 | -0.031 [-0.033, -0.029] p<1e-03 | -0.043 [-0.046, -0.040] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.149 | -0.026 [-0.049, -0.007] p=0.012 | -0.001 [-0.015, 0.011] p=0.82 | -0.031 [-0.033, -0.029] p<1e-03 | -0.043 [-0.046, -0.040] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.530 [0.518, 0.541] | 0.520 [0.509, 0.530] | 0.828 [0.823, 0.833] | 0.752 [0.746, 0.757] |
| GC-matched: PINTS |  | 0.590 [0.572, 0.612] | 0.535 [0.524, 0.546] | 0.860 [0.855, 0.865] | 0.797 [0.791, 0.802] |
| GC-matched: Δ pydreg − PINTS |  | -0.060 [-0.085, -0.039] p<1e-03 | -0.015 [-0.027, -0.004] p=0.006 | -0.032 [-0.034, -0.030] p<1e-03 | -0.045 [-0.048, -0.043] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.050 [0.035, 0.063] p<1e-03 | 0.067 [0.054, 0.078] p<1e-03 | 0.026 [0.020, 0.031] p<1e-03 | 0.008 [0.006, 0.011] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.108 [-0.127, -0.082] p<1e-03 | -0.036 [-0.048, -0.024] p<1e-03 | -0.042 [-0.044, -0.039] p<1e-03 | -0.056 [-0.059, -0.053] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.216 | 0.999 | 0.165 | 0.997 |
| 20,000 | 0.325 | 0.997 | 0.221 | 0.993 |
| 40,000 | 0.435 | 0.992 | 0.326 | 0.983 |

### GM12878_GROseq (GM12878 GRO-seq, holdout)

n_pos = 41,415, n_neg = 15,556, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.727**; CpG-island rows 0.311 (n_pos 7,034), non-CpG-island rows 0.688 (n_pos 34,381).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.502 | 0.690 [0.678, 0.701] | 0.651 [0.641, 0.661] | 0.924 [0.922, 0.927] | 0.853 [0.849, 0.857] |
| pydreg, −log10 p over all candidates | 0.465 | 0.689 [0.677, 0.699] | 0.637 [0.626, 0.648] | 0.922 [0.919, 0.924] | 0.851 [0.848, 0.855] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.367 | 0.650 [0.634, 0.672] | 0.601 [0.587, 0.616] | 0.945 [0.943, 0.947] | 0.888 [0.884, 0.891] |
| Δ pydreg(SVR) − PINTS | 0.134 | 0.040 [0.020, 0.055] p<1e-03 | 0.050 [0.036, 0.063] p<1e-03 | -0.021 [-0.022, -0.019] p<1e-03 | -0.035 [-0.038, -0.032] p<1e-03 |
| Δ pydreg(−log p) − PINTS | 0.098 | 0.039 [0.018, 0.053] p<1e-03 | 0.036 [0.021, 0.050] p<1e-03 | -0.023 [-0.025, -0.021] p<1e-03 | -0.036 [-0.039, -0.034] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.700 [0.690, 0.710] | 0.660 [0.650, 0.669] | 0.925 [0.922, 0.927] | 0.854 [0.851, 0.858] |
| GC-matched: PINTS |  | 0.673 [0.654, 0.687] | 0.619 [0.604, 0.633] | 0.946 [0.944, 0.948] | 0.891 [0.887, 0.895] |
| GC-matched: Δ pydreg − PINTS |  | 0.028 [0.014, 0.046] p<1e-03 | 0.041 [0.028, 0.055] p<1e-03 | -0.021 [-0.023, -0.020] p<1e-03 | -0.037 [-0.040, -0.034] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.050 [0.040, 0.060] p<1e-03 | 0.061 [0.044, 0.077] p<1e-03 | 0.016 [0.009, 0.024] p<1e-03 | 0.007 [0.004, 0.010] p=0.002 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | 0.023 [0.007, 0.039] p=0.01 | 0.030 [0.015, 0.046] p<1e-03 | -0.030 [-0.032, -0.028] p<1e-03 | -0.047 [-0.050, -0.044] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.224 | 0.999 | 0.171 | 0.997 |
| 20,000 | 0.368 | 0.996 | 0.262 | 0.996 |
| 40,000 | 0.508 | 0.992 | 0.360 | 0.990 |

### K562_GROcap (K562 GRO-cap, ood)

n_pos = 43,377, n_neg = 23,291, promoter seeds = 1. AUPRC baseline (positive fraction) = **0.651**; CpG-island rows 0.226 (n_pos 6,818), non-CpG-island rows 0.611 (n_pos 36,559).

| run | TPR @ 1 % FPR | TPR @ 5 % FPR | pAUROC (FPR ≤ 0.1) | AUPRC | AUROC (suppl.) |
|---|---|---|---|---|---|
| pydreg, SVR score over all candidates (raw.peak) | 0.494 | 0.661 [0.648, 0.672] | 0.617 [0.608, 0.625] | 0.879 [0.875, 0.883] | 0.821 [0.817, 0.825] |
| pydreg, −log10 p over all candidates | 0.497 | 0.659 [0.647, 0.671] | 0.606 [0.597, 0.615] | 0.875 [0.871, 0.879] | 0.820 [0.815, 0.823] |
| PINTS, −log10 p over all rebuilt divergent pairs | 0.513 | 0.674 [0.664, 0.685] | 0.645 [0.636, 0.655] | 0.898 [0.895, 0.902] | 0.847 [0.843, 0.850] |
| Δ pydreg(SVR) − PINTS | -0.019 | -0.013 [-0.023, -0.006] p<1e-03 | -0.029 [-0.034, -0.023] p<1e-03 | -0.019 [-0.020, -0.018] p<1e-03 | -0.026 [-0.028, -0.024] p<1e-03 |
| Δ pydreg(−log p) − PINTS | -0.017 | -0.016 [-0.023, -0.007] p<1e-03 | -0.039 [-0.045, -0.034] p<1e-03 | -0.023 [-0.025, -0.022] p<1e-03 | -0.027 [-0.029, -0.025] p<1e-03 |
| GC-matched: pydreg (SVR) |  | 0.672 [0.664, 0.679] | 0.628 [0.620, 0.635] | 0.880 [0.876, 0.884] | 0.823 [0.819, 0.827] |
| GC-matched: PINTS |  | 0.693 [0.683, 0.702] | 0.659 [0.651, 0.668] | 0.900 [0.896, 0.903] | 0.849 [0.845, 0.853] |
| GC-matched: Δ pydreg − PINTS |  | -0.021 [-0.030, -0.012] p<1e-03 | -0.031 [-0.036, -0.026] p<1e-03 | -0.019 [-0.021, -0.018] p<1e-03 | -0.026 [-0.028, -0.025] p<1e-03 |
| CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.004 [-0.010, 0.001] p=0.11 | -0.009 [-0.016, -0.002] p=0.008 | -0.010 [-0.014, -0.005] p<1e-03 | -0.006 [-0.008, -0.003] p<1e-03 |
| non-CpG-island positives, GC-matched: Δ pydreg − PINTS |  | -0.035 [-0.044, -0.027] p<1e-03 | -0.037 [-0.043, -0.032] p<1e-03 | -0.024 [-0.026, -0.023] p<1e-03 | -0.030 [-0.032, -0.028] p<1e-03 |

Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):

| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |
|---|---|---|---|---|
| 10,000 | 0.217 | 0.999 | 0.191 | 0.999 |
| 20,000 | 0.369 | 0.997 | 0.308 | 0.997 |
| 40,000 | 0.510 | 0.987 | 0.445 | 0.994 |


## Label-build checks

* `check/nonenh_frac_overlapping_hg38_promoters` = 0.032
* `check/k562_positive_vs_K562.positive_jaccard` = 1.0
* `check/k562_negative_rebuilt_hmm_jaccard` = 0.9618
* `check/k562_negative_rebuilt_grocap_signal_jaccard` = 0.878
* `check/K562_dhs_both_frac_h3k27ac` = 0.4516
* `check/GM12878_dhs_both_frac_h3k27ac` = 0.6225
* `check/HCT116_dnase_frac_h3k27ac` = 0.4881
* `check/A549_dnase_frac_h3k27ac` = 0.5886
* `check/HCT116_pints_tss_frac_dnase` = 0.5361
* `check/A549_pints_tss_frac_dnase` = 0.7663
* `check/K562_pints_tss_frac_dnase` = 0.8715
* `check/GM12878_hidra_groups_inactive` = 2226471
* `check/GM12878_hidra_groups_active_lfc0.5` = 794021
* `check/GM12878_enh_frac_in_hidra_active` = 0.9884
* `check/GM12878_hidra_inactive_regions_before_dhs` = 1238308
* `check/Jurkat_dhs_both_frac_h3k27ac` = 0.1601
* `check/Jurkat_dhs_both_frac_h3k27ac_reproduced` = 0.1319
* `check/K562_dhs_frac_tss_silent` = 0.2276
* `check/K562_dhs_tss_silent_source` = K562_GROcap
* `check/GM12878_dhs_frac_tss_silent` = 0.1797
* `check/GM12878_dhs_tss_silent_source` = GM12878_GROcap
* `check/HCT116_dhs_frac_tss_silent` = 0.1748
* `check/HCT116_dhs_tss_silent_source` = HCT116_PROcap
* `check/A549_dhs_frac_tss_silent` = 0.4436
* `check/A549_dhs_tss_silent_source` = A549_PROcap
* `check/Jurkat_dhs_frac_tss_silent` = 0.6339
* `check/Jurkat_dhs_tss_silent_source` = Jurkat_CAGE
* `check/Jurkat_pints_tss_frac_dnase` = 0.8826
