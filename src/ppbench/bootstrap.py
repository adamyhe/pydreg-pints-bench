"""Paired bootstrap confidence intervals for the phase-1 metrics.

For each library and label set, every labelled element is scored once per tool run (the same
element scores as ``evaluate``). Replicates then reweight elements instead of re-running overlaps:

* ``element``  - stratified: positives and negatives resampled separately, with replacement.
* ``block5Mb`` - genomic block bootstrap: 5 Mb blocks (positives + negatives together) resampled
  with replacement. Adjacent 500 bp negative tiles and clustered enhancers are not independent,
  so this is the more conservative interval.

Resampling is *paired*: within a replicate all tools share the same weights, so the interval on
(pydreg - PINTS) accounts for the shared label noise. For PINTS-paper labels, each replicate also
draws one of the N_SEEDS promoter samples (replicate b uses seed b mod N_SEEDS), so promoter
sampling variance is part of the interval.

The PINTS-paper scoring rule is bootstrapped through per-element indicator matrices
(TP-eligible / touched, per threshold), which reproduce ``metrics.pints_paper_roc`` exactly at unit
weights.
"""
from __future__ import annotations

import json
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from . import calls as C
from . import metrics as M
from .evaluate import aux_label_sets, dreg_label_sets, tss_label_sets, libraries, pints_label_sets, tool_runs
from .gc import seq_stats
from .genome import BED3

BOOT = C.RESULTS / "bootstrap"
BLOCK = 5_000_000
METRICS = ("auroc", "auprc", "tpr_at_fpr0.01", "tpr_at_fpr0.05", "pauroc_fpr0.1")
GC_BINS = 20
COND_RUNS = ("pydreg_score", "pints_logp")  # runs also scored GC-matched and split by CpG-island status
# condition suffixes: gcw = negatives reweighted to the positives' GC histogram (per replicate);
# cgi_gcw / noncgi_gcw = only CpG-island-like (GC > 0.5, CpG o/e > 0.6) / other positives, negatives
# GC-matched to that subset.
CONDS = ("gcw", "cgi_gcw", "noncgi_gcw")


# ------------------------------------------------------------------ weighted ROC metrics
class Ranked:
    """Element scores pre-ranked once, so each replicate is a bincount + cumsum."""

    def __init__(self, pos_scores, neg_scores):
        s = M._finite(np.r_[pos_scores, neg_scores])
        self.uniq, inv = np.unique(s, return_inverse=True)
        self.ip, self.ineg = inv[: len(pos_scores)], inv[len(pos_scores):]

    def metrics(self, wp, wn) -> dict:
        K = len(self.uniq)
        P = np.bincount(self.ip, wp, minlength=K)[::-1]
        N = np.bincount(self.ineg, wn, minlength=K)[::-1]
        tp, fp = np.cumsum(P), np.cumsum(N)
        if tp[-1] <= 0 or fp[-1] <= 0:
            return {k: np.nan for k in METRICS}
        tpr, fpr = np.r_[0.0, tp / tp[-1]], np.r_[0.0, fp / fp[-1]]
        prec = tp / np.maximum(tp + fp, 1e-300)
        out = {"auroc": float(np.trapezoid(tpr, fpr)),
               "auprc": float(np.sum(P / tp[-1] * prec))}  # sklearn average_precision (step)
        out["tpr_at_fpr0.01"] = float(np.interp(0.01, fpr, tpr))
        out["tpr_at_fpr0.05"] = float(np.interp(0.05, fpr, tpr))
        m = fpr <= 0.1
        out["pauroc_fpr0.1"] = float(np.trapezoid(np.r_[tpr[m], np.interp(0.1, fpr, tpr)],
                                                  np.r_[fpr[m], 0.1]) / 0.1)
        return out


# ------------------------------------------------------------------ PINTS-paper rule
def pp_indicators(pos, neg, calls, thresholds):
    """Per element and threshold: positive TP-eligible (>= 50 % of a clean call inside it),
    positive touched, negative touched. Mirrors metrics._pp_point."""
    from .genome import overlaps_any
    import bioframe as bf

    scope = pd.concat([pos[BED3], neg[BED3]], ignore_index=True)
    calls = calls[overlaps_any(calls, scope)]
    posr = pos[BED3].reset_index(drop=True).assign(_i=lambda d: np.arange(len(d)))
    T = len(thresholds)
    tpi, tpos, tneg = (np.zeros((len(pos), T), bool), np.zeros((len(pos), T), bool),
                       np.zeros((len(neg), T), bool))
    for j, t in enumerate(thresholds):
        c = calls[calls.score >= t]
        if len(c) == 0:
            continue
        clean = c[~overlaps_any(c, neg)]
        if len(clean):
            ov = bf.overlap(posr, clean[BED3], how="inner", return_overlap=True, suffixes=("", "_c"))
            frac = (ov.overlap_end - ov.overlap_start) / (ov.end_c - ov.start_c)
            tpi[ov.loc[frac >= 0.5, "_i"].unique(), j] = True
        tpos[:, j] = overlaps_any(pos, c)
        tneg[:, j] = overlaps_any(neg, c)
    return tpi, tpos, tneg


def pp_auc(ind, wp, wn) -> float:
    tpi, tpos, tneg = ind
    TP = wp @ tpi
    FN = wp @ ~tpos
    FP = wn @ tneg
    TN = wn.sum() - FP
    ok = TP > 0
    tpr = np.where(ok, TP / np.maximum(TP + FN, 1e-300), 0.0)
    fpr = np.where(ok, FP / np.maximum(FP + TN, 1e-300), 0.0)
    fpr, tpr = np.r_[0.0, 1.0, fpr], np.r_[0.0, 1.0, tpr]
    o = np.lexsort((tpr, fpr))
    return float(np.trapezoid(tpr[o], fpr[o]))


# ------------------------------------------------------------------ label-set assembly
def _blocks(df):
    return (df.chrom.astype(str) + ":" + (df.start // BLOCK).astype(str)).to_numpy()


def build_sets(lib: str):
    """{(family, variant): [per-seed dict(pos, neg, ranked{run}, pp{run})]}."""
    row = libraries().set_index("lib_id").loc[lib].copy()
    row["lib_id"] = lib
    runs = tool_runs(lib)
    full = {f"{t}_{r}": f for t, r, f, _ in runs}
    dflt = {f"{t}_{r}": d for t, r, _, d in runs}
    out, tables = {}, []

    def one(family, variant, seed, pos, neg, replica):
        pos, neg = pos[BED3].reset_index(drop=True), neg[BED3].reset_index(drop=True)
        d = {"pos": pos, "neg": neg, "ranked": {}, "pp": {},
             "blk_pos": _blocks(pos), "blk_neg": _blocks(neg)}
        tab = pd.concat([pos.assign(label=1), neg.assign(label=0)], ignore_index=True)
        st = seq_stats(tab, row.genome)
        gc = st.gc.to_numpy(float)
        gc = np.where(np.isfinite(gc), gc, np.nanmedian(gc))
        edges = np.unique(np.quantile(gc, np.linspace(0, 1, GC_BINS + 1)))
        b = np.clip(np.searchsorted(edges, gc, side="right") - 1, 0, max(len(edges) - 2, 0))
        d["gc_bin_pos"], d["gc_bin_neg"], d["gc_k"] = b[: len(pos)], b[len(pos):], max(len(edges) - 1, 1)
        d["cgi_pos"] = ((st.gc > 0.5) & (st.cpg_oe > 0.6)).to_numpy()[: len(pos)]
        tab["gc"], tab["cpg_oe"] = st.gc.to_numpy(), st.cpg_oe.to_numpy()
        for run, calls in full.items():
            ps, ns = M.element_scores(pos, calls), M.element_scores(neg, calls)
            d["ranked"][run] = Ranked(ps, ns)
            tab[run] = np.r_[ps, ns]
        if replica:
            if "pydreg_score" in full:
                for fl, df in (("paper_default", dflt["pydreg_score"]), ("paper_full", full["pydreg_score"])):
                    d["pp"][f"pydreg_score_{fl}"] = pp_indicators(pos, neg, df, M.linspace_thresholds(df))
            if "pints_logp" in full:
                f = full["pints_logp"]
                d["pp"]["pints_logp_paper_full"] = pp_indicators(pos, neg, f, M.pints_pvalue_thresholds(f))
        out.setdefault((family, variant), []).append(d)
        tables.append(tab.assign(family=family, variant=variant, seed=seed))

    for filtered in (True, False):
        variant = "assay_filtered" if filtered else "unfiltered"
        for seed, pos, neg in pints_label_sets(row, filtered) or []:
            one("pints_paper", variant, seed, pos, neg, replica=True)
    for neg_name, pos, neg in dreg_label_sets(row) or []:
        one("dreg_paper", neg_name, 0, pos, neg, replica=False)
    for variant, pos, neg in aux_label_sets(row):
        one("dreg_aux", variant, 0, pos, neg, replica=False)
    for variant, pos, neg in tss_label_sets(row) or []:
        one("tss_pints", variant, 0, pos, neg, replica=False)
    return out, (pd.concat(tables, ignore_index=True) if tables else None)


# ------------------------------------------------------------------ bootstrap
def _weights(d, scheme, rng):
    npos, nneg = len(d["pos"]), len(d["neg"])
    if scheme == "element":
        return (np.bincount(rng.integers(0, npos, npos), minlength=npos).astype(float),
                np.bincount(rng.integers(0, nneg, nneg), minlength=nneg).astype(float))
    ub, inv = np.unique(np.r_[d["blk_pos"], d["blk_neg"]], return_inverse=True)
    cnt = np.bincount(rng.integers(0, len(ub), len(ub)), minlength=len(ub)).astype(float)
    w = cnt[inv]
    return w[:npos], w[npos:]


def gc_match(d, wp, wn):
    """Negative weights reweighted so the negatives' GC histogram matches the (weighted) positives'."""
    K = d["gc_k"]
    fp = np.bincount(d["gc_bin_pos"], wp, minlength=K)
    fn = np.bincount(d["gc_bin_neg"], wn, minlength=K)
    if fp.sum() <= 0 or fn.sum() <= 0:
        return None
    r = np.where(fn > 0, (fp / fp.sum()) / np.maximum(fn / fn.sum(), 1e-300), 0.0)
    w = wn * r[d["gc_bin_neg"]]
    return w * (wn.sum() / w.sum()) if w.sum() > 0 else None


def _evaluate(d, wp, wn) -> dict:
    """{(run, metric): value} for one replicate (or the point estimate at unit weights)."""
    v = {}
    for run, rk in d["ranked"].items():
        for k, x in rk.metrics(wp, wn).items():
            v[(run, k)] = x
    if "gc_bin_pos" in d:
        cgi = d["cgi_pos"]
        subsets = {"gcw": wp, "cgi_gcw": wp * cgi, "noncgi_gcw": wp * ~cgi}
        for cond, wps in subsets.items():
            wns = gc_match(d, wps, wn)
            if wns is None:
                continue
            for run in COND_RUNS:
                if run in d["ranked"]:
                    for k, x in d["ranked"][run].metrics(wps, wns).items():
                        v[(f"{run}|{cond}", k)] = x
    for run, ind in d["pp"].items():
        v[(run, "pints_paper_auc")] = pp_auc(ind, wp, wn)
    return v


DIFFS = {  # label: (minuend, subtrahend) runs
    "pydreg_score - pints": ("pydreg_score", "pints_logp"),
    "pydreg_logp - pints": ("pydreg_logp", "pints_logp"),
    "paper: pydreg_default - pints": ("pydreg_score_paper_default", "pints_logp_paper_full"),
    "paper: pydreg_full - pints": ("pydreg_score_paper_full", "pints_logp_paper_full"),
}


def _add_diffs(v: dict) -> dict:
    out = dict(v)
    for label, (a, b) in DIFFS.items():
        for (run, k), x in v.items():
            base, _, cond = run.partition("|")
            if base == a and (b + ("|" + cond if cond else ""), k) in v:
                out[(label + ("|" + cond if cond else ""), k)] = x - v[(b + ("|" + cond if cond else ""), k)]
    return out


def _is_diff(run: str) -> bool:
    return run.partition("|")[0] in DIFFS


def bootstrap_library(lib: str, n_boot: int = 1000, seed: int = 0, outdir: Path | None = None) -> Path:
    outdir = Path(outdir) if outdir else BOOT / lib
    outdir.mkdir(parents=True, exist_ok=True)
    sets, table = build_sets(lib)
    if table is not None:
        table.to_csv(outdir / "element_scores.tsv.gz", sep="\t", index=False)
    rows = []
    for (family, variant), per_seed in sets.items():
        pts = [_add_diffs(_evaluate(d, np.ones(len(d["pos"])), np.ones(len(d["neg"])))) for d in per_seed]
        point = {k: float(np.nanmean([p[k] for p in pts])) for k in pts[0]}
        for scheme in ("element", "block5Mb"):
            rng = np.random.default_rng([seed, zlib.crc32(f"{family}|{variant}|{scheme}".encode())])
            reps = []
            for b in range(n_boot):
                d = per_seed[b % len(per_seed)]
                reps.append(_add_diffs(_evaluate(d, *_weights(d, scheme, rng))))
            for key, est in point.items():
                x = np.array([r.get(key, np.nan) for r in reps], float)
                x = x[np.isfinite(x)]
                is_diff = _is_diff(key[0])
                rows.append({
                    "lib": lib, "family": family, "variant": variant, "scheme": scheme,
                    "run": key[0], "metric": key[1], "estimate": est,
                    "lo": np.quantile(x, 0.025) if len(x) else np.nan,
                    "hi": np.quantile(x, 0.975) if len(x) else np.nan,
                    "se": x.std(ddof=1) if len(x) > 1 else np.nan,
                    # two-sided bootstrap p for diff == 0 (resolution 1/n_boot)
                    "p_boot": (min(1.0, 2 * min((x <= 0).mean(), (x >= 0).mean())) if is_diff and len(x)
                               else np.nan),
                    "n_boot": len(x), "n_pos": len(per_seed[0]["pos"]), "n_neg": len(per_seed[0]["neg"]),
                    "n_seeds": len(per_seed), "n_pos_cgi": int(per_seed[0]["cgi_pos"].sum()),
                })
    pd.DataFrame(rows).to_csv(outdir / "ci.tsv", sep="\t", index=False)
    (outdir / "done.json").write_text(json.dumps({"lib": lib, "n_boot": n_boot, "n_rows": len(rows)}))
    return outdir
