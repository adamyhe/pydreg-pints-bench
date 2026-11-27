"""Evaluation metrics.

Primary (symmetric) protocol
  Each labelled element gets the max score of any call overlapping it by >= 1 bp; elements no call
  touches get -inf. ROC/PR are computed over labelled elements from those scores, using each
  tool's *full* candidate set, so both tools are swept over their whole range with the same
  overlap rule for positives and negatives.

PINTS-paper replica (``pints_paper_roc``)
  Re-implements PINTS_analysis/utils_roc.py ``evaluation_core_element`` exactly:
  TP = positives overlapped by calls covering >= 50 % of the call's length (bedtools -F 0.5),
       counting only calls that touch no negative;
  FP = negatives touched (>= 1 bp) by any call;  FN/TN analogously;  calls restricted to those
  touching a labelled element; ROC closed with (0,0) and (1,1) and integrated by trapezoid.
"""
from __future__ import annotations

import bioframe as bf
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve

from .genome import BED3, overlaps_any

NEG_INF = -np.inf


def element_scores(elements: pd.DataFrame, calls: pd.DataFrame) -> np.ndarray:
    """Max overlapping call score per element (-inf if untouched)."""
    el = elements[BED3].reset_index(drop=True).copy()
    el["_i"] = np.arange(len(el))
    out = np.full(len(el), NEG_INF)
    if len(calls) == 0 or len(el) == 0:
        return out
    hit = bf.overlap(el, calls[BED3 + ["score"]], how="inner", suffixes=("", "_c"))
    if len(hit):
        best = hit.groupby("_i")["score_c"].max()
        out[best.index.to_numpy()] = best.to_numpy()
    return out


def _finite(scores: np.ndarray) -> np.ndarray:
    s = scores.astype(float).copy()
    fin = np.isfinite(s)
    lo = s[fin].min() - 1.0 if fin.any() else 0.0
    hi = s[fin].max() + 1.0 if fin.any() else 1.0
    s[s == NEG_INF] = lo
    s[s == np.inf] = hi
    return s


def roc_metrics(pos_scores: np.ndarray, neg_scores: np.ndarray) -> dict:
    y = np.r_[np.ones(len(pos_scores)), np.zeros(len(neg_scores))]
    s = _finite(np.r_[pos_scores, neg_scores])
    fpr, tpr, _ = roc_curve(y, s)
    out = {"n_pos": len(pos_scores), "n_neg": len(neg_scores),
           "auroc": roc_auc_score(y, s), "auprc": average_precision_score(y, s),
           "pos_called": float(np.isfinite(pos_scores).mean() if len(pos_scores) else np.nan),
           "neg_called": float(np.isfinite(neg_scores).mean() if len(neg_scores) else np.nan)}
    for f in (0.01, 0.05, 0.1):
        out[f"tpr_at_fpr{f:g}"] = float(np.interp(f, fpr, tpr))
    m = fpr <= 0.1
    out["pauroc_fpr0.1"] = float(np.trapezoid(np.r_[tpr[m], np.interp(0.1, fpr, tpr)],
                                              np.r_[fpr[m], 0.1]) / 0.1)
    return out


def roc_curve_table(pos_scores, neg_scores, n=201) -> pd.DataFrame:
    y = np.r_[np.ones(len(pos_scores)), np.zeros(len(neg_scores))]
    s = _finite(np.r_[pos_scores, neg_scores])
    fpr, tpr, _ = roc_curve(y, s)
    grid = np.linspace(0, 1, n)
    return pd.DataFrame({"fpr": grid, "tpr": np.interp(grid, fpr, tpr)})


def operating_point(pos: pd.DataFrame, neg: pd.DataFrame, calls: pd.DataFrame) -> dict:
    tp = overlaps_any(pos, calls).sum()
    fp = overlaps_any(neg, calls).sum()
    return {"op_tpr": tp / max(len(pos), 1), "op_fpr": fp / max(len(neg), 1),
            "op_precision": tp / max(tp + fp, 1), "op_n_calls": len(calls)}


def _n_merged(calls: pd.DataFrame) -> int:
    return len(bf.merge(calls[BED3])) if len(calls) else 0


def matched_budget(pos: pd.DataFrame, neg: pd.DataFrame, calls: pd.DataFrame, budgets) -> pd.DataFrame:
    """Threshold each tool so it reports ~N *merged* calls (PINTS pairs overlap; pydreg peaks don't),
    then recall over positives and labelled precision (elements)."""
    scores = np.sort(calls.score.to_numpy())[::-1]
    rows = []
    for n in budgets:
        lo, hi = 0, len(scores) - 1  # binary search on rank k: smallest k with >= n merged calls
        if len(scores) == 0 or _n_merged(calls) < n:
            k = len(scores) - 1
        else:
            while lo < hi:
                mid = (lo + hi) // 2
                if _n_merged(calls[calls.score >= scores[mid]]) >= n:
                    hi = mid
                else:
                    lo = mid + 1
            k = lo
        top = calls[calls.score >= scores[k]] if len(scores) else calls
        tp = overlaps_any(pos, top).sum()
        fp = overlaps_any(neg, top).sum()
        rows.append({"budget": int(n), "n_used": _n_merged(top), "threshold": float(scores[k]) if len(scores) else np.nan,
                     "recall": tp / max(len(pos), 1), "fpr": fp / max(len(neg), 1),
                     "precision": tp / max(tp + fp, 1)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ PINTS-paper replica
def _pp_point(pos, neg, calls):
    if len(calls) == 0:
        return 0.0, 1.0
    clean_calls = calls[~overlaps_any(calls, neg)]
    if len(clean_calls):
        ov = bf.overlap(pos[BED3].reset_index(drop=True).assign(_i=lambda d: np.arange(len(d))),
                        clean_calls[BED3], how="inner", return_overlap=True, suffixes=("", "_c"))
        frac = (ov.overlap_end - ov.overlap_start) / (ov.end_c - ov.start_c)
        tp = ov.loc[frac >= 0.5, "_i"].nunique()
    else:
        tp = 0
    fp = overlaps_any(neg, calls).sum()
    fn = (~overlaps_any(pos, calls)).sum()
    tn = len(neg) - fp
    if tp == 0:
        return 0.0, 1.0
    return tp / (tp + fn), tn / (fp + tn)


def pints_paper_roc(pos, neg, calls, thresholds) -> tuple[float, pd.DataFrame]:
    scope = pd.concat([pos[BED3], neg[BED3]], ignore_index=True)
    calls = calls[overlaps_any(calls, scope)]
    pts = [(1.0 - spec, rec, t) for t in thresholds
           for rec, spec in [_pp_point(pos, neg, calls[calls.score >= t])]]
    df = pd.DataFrame(pts, columns=["fpr", "tpr", "threshold"])
    xy = pd.concat([pd.DataFrame({"fpr": [0.0, 1.0], "tpr": [0.0, 1.0]}), df[["fpr", "tpr"]]])
    xy = xy.sort_values(["fpr", "tpr"])
    return float(np.trapezoid(xy.tpr, xy.fpr)), df


def linspace_thresholds(calls: pd.DataFrame, n=101):
    s = calls.score[np.isfinite(calls.score)]
    return np.linspace(s.min(), s.max(), n) if len(s) else np.array([0.0])


def pints_pvalue_thresholds(calls: pd.DataFrame, n=101):
    """`ROCProfiler._unbalanced_range_for_pqs`, mapped onto -log10(p) scores."""
    pmin = 10 ** (-calls.score[np.isfinite(calls.score)].max()) if len(calls) else 1e-10
    n_small = int(n * 0.25)
    n_exp = min(int(n * 0.25), 51)
    small = np.linspace(pmin, 0.1, n_small).tolist()
    exp = [10.0 ** i for i in range(-n_exp, -1)]
    n_gen = n - len(set(small + exp)) + 1
    general = np.linspace(0.1, 1, n_gen).tolist()
    p = np.array(sorted(set(small + exp + general)))
    return -np.log10(np.clip(p, 1e-300, 1))
