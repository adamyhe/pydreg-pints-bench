"""Verbose phase-1 report: label definitions, sizes, and every metric with bootstrap CIs.

Reads results/bootstrap/*/ci.tsv (falls back to results/summary point estimates if absent) and
writes results/summary/report.md (also printed).
"""
import glob
import json

import numpy as np
import pandas as pd

from ppbench.calls import PROJ, RESULTS

LIBS = pd.read_csv(PROJ / "config" / "libraries.tsv", sep="\t").set_index("lib_id")
MANIFEST = json.loads((PROJ / "data" / "labels" / "built" / "manifest.json").read_text())

LABELS = {
    ("pints_paper", "assay_filtered"): (
        "PINTS paper (Yao 2022) labels, with the paper's read filter",
        "**K562:** CRISPR-validated enhancers (Supp. Table 2, CRISPRi/CRISPRa rows; 801 after hg38→hg19 liftOver and "
        "blacklist) **plus** an equal number of promoters drawn from a pool (5 draws = 5 seeds). Pool: GENCODE v24 "
        "protein-coding TSS ± 500 bp with K562 CAGE reads on both strands within TSS ± 1 kb, FPKM ≥ 1 in both ENCODE "
        "total and polyA RNA-seq, excluding genes overlapping another gene (7,504 in hg19). **GM12878:** the 3,544 "
        "HiDRA STARR-seq enhancers of Supp. Table 2 (3,497 in hg19), plus promoters from the same recipe with "
        "GM12878 CAGE and ENCODE RNA-seq (pool 7,149).",
        "**K562:** Supp. Table 3 non-enhancers (CRISPR/MPRA/STARR-tested, inactive): ≥ 250 bp, merged, minus "
        "TSS ± 500 of all transcripts (5,952 in hg19). **GM12878** (ours; Supp. Table 3 is K562 only): HiDRA "
        "(GSE104001) fragment groups with ≥ 30 plasmid reads and library-normalised log2(RNA/plasmid) ≤ 0, "
        "merged, minus any group with log2FC ≥ 0.5 and the positives ± 250 bp. Kept if ≥ 250 bp, touching a "
        "Duke/UW DHS and clear of TSS ± 500 (75,311 in hg19). These are accessible-but-inactive hard negatives.",
        "Positives and negatives are both restricted to elements with ≥ 1 read on **each** strand in the library "
        "under test (as in PINTS_analysis). This leaves few negatives on run-on data, and it favours whichever tool "
        "calls what the assay sees."),
    ("pints_paper", "unfiltered"): (
        "PINTS paper labels, no read filter",
        "As above.", "As above.",
        "No read filter. This is closer to asking whether each tool recovers the orthogonally validated elements."),
    ("dreg_paper", "negative"): (
        "dREG paper (Wang 2019) labels, original negatives",
        "Danko-lab hub `k562.positive.bed.gz`: Duke ∩ UW DHS ∩ GRO-cap HMM transcribed sites (27,798 merged).",
        "Hub `k562.negative.bed.gz`: genome outside Duke DHS, UW DHS and GRO-cap HMM. Tiled into 500 bp tiles; only "
        "tiles holding ≥ 1 dREG informative position in this library are kept (≤ 200k, seed 0).",
        "These are the dREG authors' own labels. The positives require GRO-cap, so they are circular on K562 "
        "GRO-cap."),
    ("dreg_paper", "negative_rebuilt"): (
        "dREG paper labels, negatives rebuilt with one recipe for K562 and GM12878",
        "K562: hub positives as above. GM12878: hub `GM12878.positive.bed.gz` (36,627).",
        "Genome minus (Duke ∪ UW DHS ∪ nonzero GRO-cap signal ± 250 bp ∪ positives), then tiled as above. On K562 "
        "this recipe matches the original negatives (Jaccard 0.88 with the signal proxy used here, 0.96 with the "
        "HMM).",
        "The same recipe in both cells, so the K562 → GM12878 cross-cell comparison is like for like."),
}
AUX_POS = {
    "dhs_h3k27ac": "Duke ∩ UW DHS overlapping an ENCODE Broad H3K27ac peak, the \"DHSs marked by H3K27ac\" of "
                   "Wang 2019 Fig 1B (K562 34,763; GM12878 55,925).",
    "dhs_h3k27ac_distal": "As dhs_h3k27ac, but > 1 kb from any GENCODE v24 TSS, i.e. putative enhancers "
                          "(K562 20,553; GM12878 28,215).",
    "chromhmm_active": "ENCODE Broad ChromHMM segments in states 1 (active promoter) and 4–5 (strong enhancer), "
                       "as in Danko 2015 Fig 1b (K562 43,377; GM12878 41,415; median segment about 1.4 kb).",
}
AUX_NEG = {
    "quiet": "Genome minus (any Duke/UW DHS ∪ H3K27ac ∪ ChromHMM states 1–8), in 500 bp tiles holding ≥ 1 dREG "
             "informative position (≤ 200k), like the Danko 2015 tuning negatives.",
    "dhs_unmarked": "**Diagnostic only, not a headline contrast.** Duke ∩ UW DHSs with no H3K27ac and outside "
                    "ChromHMM states 1–7 (K562 23,291, GM12878 15,556; CTCF/insulator-rich). Lacking H3K27ac does not "
                    "mean untranscribed, so a call here may be a true weak element. Results on this set cannot "
                    "distinguish specificity from sensitivity.",
}
for _p, _pd in AUX_POS.items():
    for _n, _nd in AUX_NEG.items():
        LABELS[("dreg_aux", f"{_p}|{_n}")] = (
            f"dREG auxiliary labels: {_p} vs {_n}", _pd, _nd,
            "Not in either dREG paper as a ROC: the papers report recall only, so the negatives are ours. No "
            "nascent or TSS data is used, so these labels are not circular for any library, GRO-cap included. "
            "H3K27ac peaks and ChromHMM segments are broad, which favours wider calls under the overlap rule.")

TSS_POS = {
    "pints_tss": "PINTS default bidirectional calls on the matched-cell TSS library: ENCODE PRO-cap (hg38; K562 "
                 "lifted to hg19 for the K562 run-on libraries), or FANTOM5 CAGE for Jurkat (hg19).",
    "pints_tss_dhs": "pints_tss overlapping the ENCODE DNase peaks of that cell: the dREG-paper recipe (DHS ∩ TSS-assay "
                     "transcription), with PINTS-on-PRO-cap in place of the GRO-cap HMM.",
    "pints_tss_distal": "pints_tss > 1 kb from any GENCODE v24 TSS.",
}
for _p, _pd in TSS_POS.items():
    LABELS[("tss_pints", f"{_p}|tss_silent")] = (
        f"PINTS-on-TSS labels: {_p} vs tss_silent", _pd,
        "Genome minus (ENCODE DNase ∪ PRO-cap signal ± 250 bp ∪ positives), in informative 500 bp tiles (≤ 200k). "
        "For K562 hg19 the negatives are the dREG hub negatives minus the lifted positives ± 250 bp.",
        "Applied only to libraries other than the PRO-cap library the labels came from. The positives are PINTS "
        "calls, so they favour PINTS's element boundaries and its sensitivity profile.")
_HG38_AUX = (" On hg38 cells (HCT116, A549) the DHS is one ENCODE DNase peak set, not Duke ∩ UW, and there is no "
             "ChromHMM: quiet excludes only DNase and H3K27ac. Jurkat (hg19) uses UW DNase rep1 ∩ rep2 with H3K27ac "
             "from either ChIP-Atlas replicate (GSM4551183/4), and has no ChromHMM.")
for _k in [k for k in LABELS if k[0] == "dreg_aux"]:
    t, a, b, c = LABELS[_k]
    LABELS[_k] = (t, a, b, c + _HG38_AUX)

_SIL = ("DHSs with no read in an independent TSS assay within ± 250 bp (K562: GRO-cap + HMM; GM12878: GRO-cap; "
        "HCT116/A549: PRO-cap; Jurkat: CAGE). These are accessible-but-untranscribed hard negatives (K562 17,523; "
        "GM12878 16,074; HCT116 11,691; A549 62,494 and Jurkat 62,150, the last two inflated by shallow TSS assays).")
_SILN = ("Scored on 3′ run-on libraries only. Positives overlapping a negative are dropped. This replaces "
         "dhs_unmarked as the specificity test.")
LABELS[("dreg_paper", "dhs_tss_silent")] = ("dREG paper positives vs TSS-silent DHSs",
    "dREG hub positives (DHS ∩ GRO-cap HMM; K562, GM12878).", _SIL, _SILN)
for _p, _pd in AUX_POS.items():
    LABELS[("dreg_aux", f"{_p}|dhs_tss_silent")] = (f"dREG auxiliary labels: {_p} vs TSS-silent DHSs", _pd, _SIL, _SILN)
for _p, _pd in TSS_POS.items():
    LABELS[("tss_pints", f"{_p}|dhs_tss_silent")] = (f"PINTS-on-TSS labels: {_p} vs TSS-silent DHSs", _pd, _SIL,
        _SILN + " The negatives come from the same TSS library as the positives, so the two sets are disjoint by "
        "construction.")

RUNS = {
    "pydreg_score": "pydreg, SVR score over all candidates (raw.peak)",
    "pydreg_logp": "pydreg, −log10 p over all candidates",
    "pints_logp": "PINTS, −log10 p over all rebuilt divergent pairs",
    "pydreg_score - pints": "Δ pydreg(SVR) − PINTS",
    "pydreg_logp - pints": "Δ pydreg(−log p) − PINTS",
    "pydreg_score|gcw": "GC-matched: pydreg (SVR)",
    "pints_logp|gcw": "GC-matched: PINTS",
    "pydreg_score - pints|gcw": "GC-matched: Δ pydreg − PINTS",
    "pydreg_score - pints|cgi_gcw": "CpG-island positives, GC-matched: Δ pydreg − PINTS",
    "pydreg_score - pints|noncgi_gcw": "non-CpG-island positives, GC-matched: Δ pydreg − PINTS",
    "pydreg_score_paper_default": "PINTS rule: pydreg default calls only (as in the paper)",
    "pydreg_score_paper_full": "PINTS rule: pydreg all candidates",
    "pints_logp_paper_full": "PINTS rule: PINTS",
    "paper: pydreg_default - pints": "Δ PINTS rule, pydreg default − PINTS",
    "paper: pydreg_full - pints": "Δ PINTS rule, pydreg full − PINTS",
}
METRIC_COLS = ["tpr_at_fpr0.01", "tpr_at_fpr0.05", "pauroc_fpr0.1", "auprc", "auroc", "pints_paper_auc"]
COL_NAMES = {"tpr_at_fpr0.01": "TPR @ 1 % FPR", "tpr_at_fpr0.05": "TPR @ 5 % FPR", "pauroc_fpr0.1": "pAUROC (FPR ≤ 0.1)",
             "auprc": "AUPRC", "auroc": "AUROC (suppl.)", "pints_paper_auc": "PINTS-rule AUC (replication)"}
BUDGETS_SHOWN = (10_000, 20_000, 40_000)


def eval_point(metric: str) -> pd.DataFrame:
    """As-evaluated point estimates of one metric from results/eval (seed mean), in ci.tsv layout."""
    m = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(str(RESULTS / "eval" / "*" / "metrics.tsv"))])
    m = m[(m.metric == metric) & (m.flavour == "full")]
    m["run"] = m.tool + "_" + m["rank"]
    g = m.groupby(["lib", "family", "variant", "run"], dropna=False).value.mean().reset_index()
    g = g.rename(columns={"value": "estimate"}).assign(metric=metric)
    import ppbench.bootstrap as Bt
    extra = []
    for key, d in g.groupby(["lib", "family", "variant"], dropna=False):
        v = {(r.run, metric): r.estimate for r in d.itertuples()}
        for (run, _), x in Bt._add_diffs(v).items():
            if run in Bt.DIFFS:
                extra.append({"lib": key[0], "family": key[1], "variant": key[2], "run": run, "metric": metric,
                              "estimate": x})
    return pd.concat([g, pd.DataFrame(extra)], ignore_index=True)


def budget_table(lib, family, variant) -> list[str]:
    f = RESULTS / "eval" / lib / "budget_curves.tsv.gz"
    if not f.exists():
        return []
    b = pd.read_csv(f, sep="\t")
    b = b[(b.family == family) & (b.variant == variant) & (b.seed <= 0) & b.budget.isin(BUDGETS_SHOWN)
          & (((b.tool == "pydreg") & (b["rank"] == "score")) | (b.tool == "pints"))]
    if b.empty:
        return []
    out = ["Matched call budgets (top-N merged calls; recall / labelled precision; seed 0, no CI):", "",
           "| budget | pydreg recall | pydreg precision | PINTS recall | PINTS precision |", "|---|---|---|---|---|"]
    for n in BUDGETS_SHOWN:
        r = {t: b[(b.budget == n) & (b.tool == t)] for t in ("pydreg", "pints")}
        if any(x.empty for x in r.values()):
            continue
        out.append(f"| {n:,} | " + " | ".join(f"{r[t].recall.iloc[0]:.3f} | {r[t].precision.iloc[0]:.3f}"
                                              for t in ("pydreg", "pints")) + " |")
    return out + [""]
CIRCULAR = {("dreg_paper", "K562_GROcap")}


def fmt(r, scheme):
    if r is None:
        return ""
    est = f"{r['estimate']:.3f}"
    if scheme is None or pd.isna(r.get("lo", np.nan)):
        return est
    s = f"{est} [{r['lo']:.3f}, {r['hi']:.3f}]"
    if pd.notna(r.get("p_boot")):
        s += f" p={max(r['p_boot'], 1 / r['n_boot']):.3g}" if r["p_boot"] > 0 else f" p<{1 / r['n_boot']:.0e}"
    return s


def point_estimates() -> pd.DataFrame:
    """Same layout as ci.tsv, built from results/summary/metrics_by_seed_mean.tsv (no CIs)."""
    m = pd.read_csv(RESULTS / "summary" / "metrics_by_seed_mean.tsv", sep="\t")
    m = m[m.family.isin(["pints_paper", "dreg_paper", "dreg_aux", "tss_pints"])]
    rows = []
    for _, r in m.iterrows():
        if r.flavour == "full" and r.metric in METRIC_COLS:
            run = f"{r.tool}_{r['rank']}"
        elif r.metric == "pints_paper_auc":
            run = f"{r.tool}_{r['rank']}_{r.flavour}"
        else:
            continue
        rows.append({"lib": r.lib, "family": r.family, "variant": r.variant, "run": run,
                     "metric": r.metric, "estimate": r["mean"]})
    df = pd.DataFrame(rows)
    n = m[(m.metric.isin(["n_pos", "n_neg"])) & (m.tool == "pints")]
    n = n.pivot_table(index=["lib", "family", "variant"], columns="metric", values="mean").reset_index()
    df = df.merge(n, on=["lib", "family", "variant"], how="left")
    import ppbench.bootstrap as Bt
    extra = []
    for key, g in df.groupby(["lib", "family", "variant"]):
        v = {(r.run, r.metric): r.estimate for r in g.itertuples()}
        for (run, metric), x in Bt._add_diffs(v).items():
            if run in Bt.DIFFS:
                extra.append({**g.iloc[0].to_dict(), "run": run, "metric": metric, "estimate": x})
    df = pd.concat([df, pd.DataFrame(extra)], ignore_index=True)
    df["n_seeds"] = np.where(df.family == "pints_paper", 5, 1)
    df[["n_pos", "n_neg"]] = df[["n_pos", "n_neg"]].round().astype(int)
    return df


def main():
    parts = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(RESULTS / "bootstrap" / "*" / "ci.tsv")))]
    ci = pd.concat(parts, ignore_index=True) if parts else None
    lines = ["# Phase 1 detailed report", "",
             "pydreg 1.0.0 vs pyPINTS 1.3.0. Each cell gives the estimate [95 % CI], and Δ rows add a two-sided "
             "paired-bootstrap p.",
             "",
             "* CIs are paired bootstrap intervals (1,000 replicates; `src/ppbench/bootstrap.py`). **block5Mb** "
             "resamples 5 Mb genomic blocks and is the headline interval. **element** is stratified by label and "
             "is in `results/bootstrap/*/ci.tsv`.",
             "* For PINTS-paper labels, each replicate also redraws the promoter sample (5 seeds), and the estimate "
             "is the mean over seeds.",
             "* **Headline metrics:** TPR at 1 % and 5 % FPR, and pAUROC (FPR ≤ 0.1, normalised). These sit in "
             "the low-FPR range where genome-wide calling operates, and do not depend on the class ratio.",
             "* **AUPRC** is shown with its baseline, the positive fraction, which varies widely between label sets "
             "because negative counts are capped or filtered by construction. Compare the tools within a set; do "
             "not compare AUPRC across sets.",
             "* **AUROC** is supplementary only. It weights every FPR equally, and genome-wide only low FPRs matter.",
             "* TPR @ 1 % FPR cells without an interval are as-evaluated point estimates: bootstraps run before "
             "that metric was added.",
             "* All metrics are element-level and symmetric: a labelled element takes the max score of any call "
             "touching it.",
             "* GC-matched rows reweight negatives, in every replicate, so their GC histogram (20 bins) matches the "
             "positives'. CpG-island-like = GC > 0.5 and CpG o/e > 0.6; the non-CpG-island rows isolate "
             "enhancer-like positives.",
             "* `pints_paper_auc` uses the PINTS_analysis rule: TP needs ≥ 50 % of a call inside a positive, and FP "
             "is any touch.",
             "", "## Libraries", "",
             "| lib | cell | assay | signal | role for dREG | notes |", "|---|---|---|---|---|---|"]
    for lib, r in LIBS.iterrows():
        lines.append(f"| {lib} | {r.cell} | {r.assay} | {r.signal} | {r.dreg_role} | {r.notes} |")
    lines += ["", "dREG roles: train = in dREG's training data (optimistic for pydreg), holdout = held out by the dREG "
              "paper, new_cell = cell line dREG never saw, ood = TSS assay outside dREG's design."]

    if ci is None:
        ci, scheme = point_estimates(), None
        lines[2] = "pydreg 1.0.0 vs pyPINTS 1.3.0. **Point estimates only; bootstrap CIs not run yet.**"
    else:
        scheme = "block5Mb"
        ci = ci[ci.scheme == scheme]
        if not (ci.metric == "tpr_at_fpr0.01").any():
            # bootstraps predating TPR@1% FPR: show the as-evaluated point estimate (no CI) for plain runs
            n = ci.drop_duplicates(["lib", "family", "variant"])[["lib", "family", "variant", "n_pos", "n_neg",
                                                                 "n_seeds"] + (["n_pos_cgi"] if "n_pos_cgi" in ci else [])]
            ci = pd.concat([ci, eval_point("tpr_at_fpr0.01").merge(n, on=["lib", "family", "variant"])],
                           ignore_index=True)
    order = sorted(LABELS, key=lambda k: k[1].endswith("|dhs_unmarked"))  # diagnostic sets last
    for (family, variant) in order:
        title, pos, neg, note = LABELS[(family, variant)]
        sub = ci[(ci.family == family) & (ci.variant == variant)]
        if sub.empty:
            continue
        lines += ["", f"## {title}", "", f"* **Positives:** {pos}", f"* **Negatives:** {neg}", f"* **Notes:** {note}", ""]
        for lib in [l for l in LIBS.index if l in set(sub.lib)]:
            s = sub[sub.lib == lib]
            n = s.iloc[0]
            tag = " — **circular, do not interpret**" if (family, lib) in CIRCULAR else ""
            npos, nneg = int(n.n_pos), int(n.n_neg)
            ncgi = int(n.n_pos_cgi) if "n_pos_cgi" in n and pd.notna(n.n_pos_cgi) else None
            base = npos / (npos + nneg)
            bl = f"AUPRC baseline (positive fraction) = **{base:.3f}**"
            if ncgi is not None:
                bl += (f"; CpG-island rows {ncgi / (ncgi + nneg):.3f} (n_pos {ncgi:,}), non-CpG-island rows "
                       f"{(npos - ncgi) / (npos - ncgi + nneg):.3f} (n_pos {npos - ncgi:,})")
            lines += [f"### {lib} ({LIBS.loc[lib, 'cell']} {LIBS.loc[lib, 'assay']}, {LIBS.loc[lib, 'dreg_role']}){tag}",
                      "", f"n_pos = {npos:,}, n_neg = {nneg:,}, promoter seeds = {n.n_seeds}. {bl}.", ""]
            cols = [c for c in METRIC_COLS if c in set(s.metric)]
            lines += ["| run | " + " | ".join(COL_NAMES[c] for c in cols) + " |", "|---|" + "---|" * len(cols)]
            for run, label in RUNS.items():
                rr = s[s.run == run]
                if rr.empty:
                    continue
                by = {m: rr[rr.metric == m].iloc[0].to_dict() for m in cols if (rr.metric == m).any()}
                lines.append(f"| {label} | " + " | ".join(fmt(by.get(m), scheme) for m in cols) + " |")
            lines.append("")
            lines += budget_table(lib, family, variant)
    lines += ["", "## Label-build checks", ""] + [f"* `{k}` = {v}" for k, v in MANIFEST.items() if k.startswith("check/")]
    text = "\n".join(lines) + "\n"
    out = RESULTS / "summary" / "report.md"
    out.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
