import numpy as np
import pandas as pd
import pytest

from ppbench import calls as C
from ppbench import metrics as M
from ppbench.evaluate import _tile
from ppbench.genome import DATA, read_bed


def bed(rows, extra=()):
    return pd.DataFrame(rows, columns=["chrom", "start", "end", *extra])


def test_pairing_follows_pints_windows():
    cand = pd.DataFrame([
        # + anchor at 1000-1050; '-' partner must overlap [700, 1050)
        ("chr1", 1000, 1050, "+", 1e-8, 1e-6, 1020),
        ("chr1", 800, 850, "-", 0.2, 0.5, 820),      # upstream within 300 -> pairs
        # isolated '+' peak, no opposite partner within 300 bp
        ("chr1", 5000, 5050, "+", 1e-12, 1e-10, 5020),
        ("chr1", 5400, 5450, "-", 1e-3, 1e-2, 5420),  # downstream of + anchor: not in + window
    ], columns=["chrom", "start", "end", "strand", "pval", "qval", "summit"])
    el = C.pair_pints_candidates(cand)
    # pair 1: + anchor (p=1e-8) and - anchor (p=0.2) merge into one element [800,1050), score 8
    first = el[el.start == 800].iloc[0]
    assert first.end == 1050 and first.score == pytest.approx(8)
    # the '-' peak at 5400 looks downstream (+300) -> [5400, 5750): no '+' there; '+' at 5000 looks
    # upstream [4700, 5050): no '-' there. So no element near 5000.
    assert not ((el.start >= 4600) & (el.start <= 5500)).any()


def test_element_scores_and_roc():
    pos = bed([("chr1", 0, 100), ("chr1", 1000, 1100), ("chr1", 5000, 5100)])
    neg = bed([("chr1", 2000, 2100), ("chr1", 3000, 3100)])
    calls = bed([("chr1", 50, 60, 5.0), ("chr1", 1090, 1200, 3.0), ("chr1", 2050, 2060, 1.0)], extra=["score"])
    ps, ns = M.element_scores(pos, calls), M.element_scores(neg, calls)
    assert list(ps[:2]) == [5.0, 3.0] and ps[2] == -np.inf
    assert ns[0] == 1.0 and ns[1] == -np.inf
    r = M.roc_metrics(ps, ns)
    assert 0.5 < r["auroc"] <= 1.0 and r["pos_called"] == pytest.approx(2 / 3)


def test_pints_paper_rule_penalises_wide_calls():
    pos = bed([("chr1", 1000, 1300)])
    neg = bed([("chr1", 5000, 5300)])
    narrow = bed([("chr1", 1050, 1250, 1.0)], extra=["score"])
    wide = bed([("chr1", 900, 1700, 1.0)], extra=["score"])       # only 300/800 inside positive
    assert M._pp_point(pos, neg, narrow)[0] == 1.0
    assert M._pp_point(pos, neg, wide)[0] == 0.0                   # fails -F 0.5 although it overlaps
    touching = bed([("chr1", 1100, 1200, 1.0), ("chr1", 5299, 5400, 1.0)], extra=["score"])
    rec, spec = M._pp_point(pos, neg, touching)
    assert rec == 1.0 and spec == 0.0                              # 1-bp graze of a negative = FP


def test_tiling():
    t = _tile(bed([("chr1", 0, 1250), ("chr2", 100, 400), ("chr3", 10, 1010)]))
    assert t.values.tolist() == [["chr1", 0, 500], ["chr1", 500, 1000], ["chr3", 10, 510], ["chr3", 510, 1010]]


def test_pvalue_grid_monotone():
    calls = bed([("chr1", 0, 10, 12.0), ("chr1", 20, 30, 0.5)], extra=["score"])
    g = M.pints_pvalue_thresholds(calls)
    # like ROCProfiler._unbalanced_range_for_pqs: linspace from the smallest p, a fixed 1e-25..1e-2
    # ladder, then 0.1..1 -- so the grid reaches p=1 (score 0) and covers the smallest observed p
    assert abs(g.min()) < 1e-12 and np.isclose(g, 12.0).any() and g.max() == pytest.approx(25.0)


@pytest.mark.skipif(not (DATA / "labels/pints_supp/41587_2022_1211_MOESM3_ESM.xlsx").exists(), reason="no supp")
def test_supp_tables_parse():
    from ppbench.labels import _supp_table
    e = _supp_table("SupplementaryTable2")
    ne = _supp_table("SupplementaryTable3")
    k = e[e["Cell line"] == "K562"]
    assert k["Supported by"].str.contains("CRISPR").sum() == 803
    assert len(ne) == 7097 and set(ne.columns[:3]) == {"chrom", "start", "end"}


@pytest.mark.skipif(not (DATA / "labels/dreg_ftp/k562.negative.bed.gz").exists(), reason="no hub labels")
def test_hub_labels_parse_with_r_headers():
    neg = read_bed(DATA / "labels/dreg_ftp/k562.negative.bed.gz")
    pos = read_bed(DATA / "labels/dreg_ftp/k562.positive.bed.gz")
    assert neg.chrom.iloc[0].startswith("chr") and len(neg) == 224864
    assert len(pos) == 158111


def test_bootstrap_unit_weights_match_point_metrics():
    from ppbench import bootstrap as Bt
    rng = np.random.default_rng(1)
    ps = np.r_[rng.normal(1, 1, 300), np.full(20, -np.inf)]
    ns = np.r_[rng.normal(0, 1, 500), np.full(200, -np.inf)]
    ps[:30] = np.round(ps[:30])  # ties
    ref = M.roc_metrics(ps, ns)
    got = Bt.Ranked(ps, ns).metrics(np.ones(len(ps)), np.ones(len(ns)))
    for k in Bt.METRICS:
        assert got[k] == pytest.approx(ref[k], abs=1e-9), k
    # integer weights == explicit duplication
    wp = rng.integers(0, 3, len(ps)).astype(float)
    wn = rng.integers(0, 3, len(ns)).astype(float)
    ref2 = M.roc_metrics(np.repeat(ps, wp.astype(int)), np.repeat(ns, wn.astype(int)))
    got2 = Bt.Ranked(ps, ns).metrics(wp, wn)
    for k in Bt.METRICS:
        assert got2[k] == pytest.approx(ref2[k], abs=1e-9), k


def test_bootstrap_pints_rule_matches_replica():
    from ppbench import bootstrap as Bt
    rng = np.random.default_rng(2)
    starts = np.sort(rng.choice(np.arange(0, 2_000_000, 1000), 400, replace=False))
    el = bed([("chr1", s, s + 400) for s in starts])
    pos, neg = el.iloc[::2].reset_index(drop=True), el.iloc[1::2].reset_index(drop=True)
    cs = rng.choice(np.arange(0, 2_000_000, 50), 600, replace=False)
    w = rng.integers(50, 900, 600)
    calls = bed([("chr1", int(s), int(s + x), float(v)) for s, x, v in zip(cs, w, rng.random(600))],
                extra=["score"])
    th = M.linspace_thresholds(calls, 21)
    auc, _ = M.pints_paper_roc(pos, neg, calls, th)
    ind = Bt.pp_indicators(pos, neg, calls, th)
    assert Bt.pp_auc(ind, np.ones(len(pos)), np.ones(len(neg))) == pytest.approx(auc, abs=1e-12)


def test_gc_match_equalises_histograms():
    from ppbench import bootstrap as Bt
    rng = np.random.default_rng(3)
    d = {"gc_bin_pos": rng.integers(10, 20, 500), "gc_bin_neg": rng.integers(0, 20, 3000), "gc_k": 20}
    wp, wn = np.ones(500), rng.integers(0, 3, 3000).astype(float)
    w = Bt.gc_match(d, wp, wn)
    hp = np.bincount(d["gc_bin_pos"], wp, minlength=20) / wp.sum()
    hn = np.bincount(d["gc_bin_neg"], w, minlength=20) / w.sum()
    assert np.allclose(hp, hn) and w.sum() == pytest.approx(wn.sum())
    # diffs propagate condition suffixes
    v = {("pydreg_score|gcw", "auroc"): 0.8, ("pints_logp|gcw", "auroc"): 0.7}
    assert Bt._add_diffs(v)[("pydreg_score - pints|gcw", "auroc")] == pytest.approx(0.1)
