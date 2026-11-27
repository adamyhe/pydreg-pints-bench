"""Per-library evaluation: assemble assay-matched labels, load calls, compute metrics."""
from __future__ import annotations

import json
from pathlib import Path

import bioframe as bf
import numpy as np
import pandas as pd
from figwig import BigWigReader

from . import calls as C
from . import metrics as M
from .genome import BED3, DATA, overlaps_any, read_bed
from .labels import BUILT

EVAL = C.RESULTS / "eval"
N_SEEDS = 5
NEG_TILE = 500
MAX_NEG_TILES = 200_000
BUDGETS = [5_000, 10_000, 20_000, 30_000, 40_000, 60_000, 80_000, 100_000]


def libraries() -> pd.DataFrame:
    return pd.read_csv(C.PROJ / "config" / "libraries.tsv", sep="\t")


def _lbl(family, cell, name, genome):
    p = BUILT / family / cell / f"{name}.{genome}.bed.gz"
    return read_bed(p)[BED3] if p.exists() else None


# ------------------------------------------------------------------ read support (PINTS filter)
def strand_sums(el: pd.DataFrame, lib: str) -> tuple[np.ndarray, np.ndarray]:
    """Total |signal| on each strand over each element, from the library's bigWigs."""
    el = el.reset_index(drop=True)
    L = (el.end - el.start).to_numpy()
    out = [np.zeros(len(el)), np.zeros(len(el))]
    for k, strand in enumerate(("pl", "mn")):
        rd = BigWigReader(DATA / "libs" / f"{lib}.{strand}.bw")
        width_bin = np.ceil(L / 512).astype(int) * 512
        for w in np.unique(width_bin):
            idx = np.flatnonzero(width_bin == w)
            vals = rd.read(el.chrom.to_numpy()[idx], el.start.to_numpy()[idx], int(w), oob=0.0)
            mask = np.arange(w)[None, :] < L[idx][:, None]
            out[k][idx] = np.abs(np.nan_to_num(vals) * mask).sum(axis=1)
    return out[0], out[1]


def both_strands(el: pd.DataFrame, lib: str) -> pd.DataFrame:
    if len(el) == 0:
        return el
    pl, mn = strand_sums(el, lib)
    return el.reset_index(drop=True)[(pl > 0) & (mn > 0)]


# ------------------------------------------------------------------ label sets per library
def pints_label_sets(row, filtered: bool):
    """PINTS-paper positives/negatives for this library; yields (seed, pos, neg)."""
    g, cell, lib = row.genome, row.cell, row.lib_id
    enh = _lbl("pints", cell, "enhancers_crispr", g)  # K562: CRISPR-validated
    if enh is None:
        enh = _lbl("pints", cell, "enhancers_hidra", g)  # GM12878: HiDRA STARR-seq
    if enh is None:
        return
    ne = _lbl("pints", cell, "nonenhancers", g)
    pool = _lbl("pints", cell, "promoter_pool", g)
    if filtered:
        enh, ne, pool = (both_strands(x, lib) for x in (enh, ne, pool))
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed)
        take = rng.choice(len(pool), size=min(len(enh), len(pool)), replace=False)
        pos = bf.sort_bedframe(pd.concat([enh, pool.iloc[take]], ignore_index=True))
        yield seed, pos, ne


def dreg_label_sets(row):
    g, cell, lib = row.genome, row.cell, row.lib_id
    pos = _lbl("dreg", cell, "positive", g)
    if pos is None:
        return
    infp = None
    for neg_name in ("negative", "negative_rebuilt"):
        regions = _lbl("dreg", cell, neg_name, g)
        if regions is None:
            continue
        if infp is None:
            infp = C.load_infp(lib)
        yield neg_name, pos, informative_tiles(regions, infp)
    sil = silent_negatives(row)
    if sil is not None:
        yield "dhs_tss_silent", _drop_overlap(pos, sil), sil


def informative_tiles(regions, infp):
    """500 bp tiles of negative regions holding >= 1 dREG informative position (the dREG paper's
    evaluation universe), capped at MAX_NEG_TILES (seed 0)."""
    tiles = pd.concat([_tile(r) for r in _chunks(bf.sort_bedframe(regions))], ignore_index=True)
    tiles = tiles[overlaps_any(tiles, infp)]
    if len(tiles) > MAX_NEG_TILES:
        tiles = tiles.sample(MAX_NEG_TILES, random_state=0)
    return bf.sort_bedframe(tiles)


AUX_POS = ("dhs_h3k27ac", "dhs_h3k27ac_distal", "chromhmm_active")
AUX_NEG = ("quiet", "dhs_unmarked")


def aux_label_sets(row):
    """dREG-paper auxiliary labels (labels.build_dreg_aux); yields ('pos|neg', pos, neg)."""
    g, cell, lib = row.genome, row.cell, row.lib_id
    negs = {}
    for n in AUX_NEG:
        x = _lbl("dreg_aux", cell, f"neg_{n}", g)
        if x is not None:
            negs[n] = informative_tiles(x, C.load_infp(lib)) if n == "quiet" else x
    sil = silent_negatives(row)
    for pn in AUX_POS:
        pos = _lbl("dreg_aux", cell, f"pos_{pn}", g)
        if pos is None:
            continue
        for n, neg in negs.items():
            yield f"{pn}|{n}", pos, neg
        if sil is not None:
            yield f"{pn}|dhs_tss_silent", _drop_overlap(pos, sil), sil


def silent_negatives(row):
    """labels.build_dhs_tss_silent for this cell, or None on TSS-assay (5') libraries, where it is circular."""
    if row.signal == "5p":
        return None
    return _lbl("silent", row.cell, "dhs_tss_silent", row.genome)


def _drop_overlap(pos, neg):
    return pos[~overlaps_any(pos, neg)].reset_index(drop=True)


def _chunks(df, n=200_000):
    for i in range(0, len(df), n):
        yield df.iloc[i:i + n]


def _tile(regions: pd.DataFrame) -> pd.DataFrame:
    r = regions[(regions.end - regions.start) >= NEG_TILE]
    n = ((r.end - r.start) // NEG_TILE).to_numpy()
    starts = np.repeat(r.start.to_numpy(), n) + NEG_TILE * (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    return pd.DataFrame({"chrom": np.repeat(r.chrom.to_numpy(), n), "start": starts, "end": starts + NEG_TILE})


TSS_POS = ("pints_tss", "pints_tss_dhs", "pints_tss_distal")


def tss_label_sets(row):
    """PINTS-on-matched-PRO-cap labels (labels.build_tss_pints); skipped on the source library."""
    from .labels import TSS_LIB
    g, cell, lib = row.genome, row.cell, row.lib_id
    if TSS_LIB.get(cell) == lib:
        return
    neg = _lbl("tss_pints", cell, "neg_tss_silent", g)
    if neg is None:
        return
    neg = informative_tiles(neg, C.load_infp(lib))
    sil = silent_negatives(row)
    for pn in TSS_POS:
        pos = _lbl("tss_pints", cell, f"pos_{pn}", g)
        if pos is not None:
            yield f"{pn}|tss_silent", pos, neg
            if sil is not None:
                yield f"{pn}|dhs_tss_silent", _drop_overlap(pos, sil), sil


# ------------------------------------------------------------------ main
def tool_runs(lib):
    """(tool, rank, full_calls, default_calls, sweep_thresholds_for_replica)."""
    out = []
    for rank in ("score", "logp"):
        try:
            full = C.load_pydreg(lib, "full", rank)
            dflt = C.load_pydreg(lib, "default", rank)
            out.append(("pydreg", rank, full, dflt))
        except FileNotFoundError as e:
            print("missing pydreg output:", e)
            break
    try:
        out.append(("pints", "logp", C.load_pints(lib, "full"), C.load_pints(lib, "default")))
    except FileNotFoundError as e:
        print("missing PINTS output:", e)
    return out


def _record(rows, base, d):
    for k, v in d.items():
        rows.append({**base, "metric": k, "value": v})


def evaluate_library(lib: str) -> Path:
    row = libraries().set_index("lib_id").loc[lib]
    row = row.copy()
    row["lib_id"] = lib
    outdir = EVAL / lib
    outdir.mkdir(parents=True, exist_ok=True)
    runs = tool_runs(lib)
    metrics, curves, budgets, replica = [], [], [], []

    # call-set descriptors
    for tool, rank, full, dflt in runs:
        for fl, df in (("full", full), ("default", dflt)):
            w = (df.end - df.start)
            _record(metrics, {"lib": lib, "tool": tool, "rank": rank, "flavour": fl, "family": "calls",
                              "variant": "", "seed": -1},
                    {"n_calls": len(df), "width_median": float(w.median()) if len(df) else np.nan,
                     "width_q25": float(w.quantile(.25)) if len(df) else np.nan,
                     "width_q75": float(w.quantile(.75)) if len(df) else np.nan})
    # PINTS reconstruction check: anchors at the default cutoff vs PINTS's own bidirectional calls
    pints_runs = [r for r in runs if r[0] == "pints"]
    if pints_runs:
        full, dflt = pints_runs[0][2], pints_runs[0][3]
        a, b = bf.merge(full[BED3]), bf.merge(dflt[BED3])
        _record(metrics, {"lib": lib, "tool": "pints", "rank": "logp", "flavour": "check", "family": "calls",
                          "variant": "", "seed": -1},
                {"default_covered_by_full": float(overlaps_any(b, a).mean()) if len(b) else np.nan})

    def score_all(family, variant, seed, pos, neg, replica_ok):
        for tool, rank, full, dflt in runs:
            base = {"lib": lib, "tool": tool, "rank": rank, "flavour": "full", "family": family,
                    "variant": variant, "seed": seed}
            ps, ns = M.element_scores(pos, full), M.element_scores(neg, full)
            _record(metrics, base, M.roc_metrics(ps, ns))
            _record(metrics, {**base, "flavour": "default"}, M.operating_point(pos, neg, dflt))
            if seed <= 0:
                c = M.roc_curve_table(ps, ns)
                curves.append(c.assign(**base))
                budgets.append(M.matched_budget(pos, neg, full, BUDGETS).assign(**base))
            if replica_ok:
                # (a) PINTS-paper protocol as run: dREG = default calls + linspace(SVR); PINTS = full + p grid
                # (b) same rule, but dREG also given its full candidate set
                if tool == "pydreg" and rank == "score":
                    for fl, df in (("paper_default", dflt), ("paper_full", full)):
                        auc, pts = M.pints_paper_roc(pos, neg, df, M.linspace_thresholds(df))
                        _record(metrics, {**base, "flavour": fl}, {"pints_paper_auc": auc})
                        if seed <= 0:
                            replica.append(pts.assign(**{**base, "flavour": fl}))
                elif tool == "pints":
                    auc, pts = M.pints_paper_roc(pos, neg, full, M.pints_pvalue_thresholds(full))
                    _record(metrics, {**base, "flavour": "paper_full"}, {"pints_paper_auc": auc})
                    if seed <= 0:
                        replica.append(pts.assign(**{**base, "flavour": "paper_full"}))

    for filtered in (True, False):
        variant = "assay_filtered" if filtered else "unfiltered"
        for seed, pos, neg in pints_label_sets(row, filtered) or []:
            score_all("pints_paper", variant, seed, pos, neg, replica_ok=True)
    for neg_name, pos, neg in dreg_label_sets(row) or []:
        score_all("dreg_paper", neg_name, 0, pos, neg, replica_ok=False)
    for variant, pos, neg in aux_label_sets(row):
        score_all("dreg_aux", variant, 0, pos, neg, replica_ok=False)
    for variant, pos, neg in tss_label_sets(row) or []:
        score_all("tss_pints", variant, 0, pos, neg, replica_ok=False)

    pd.DataFrame(metrics).to_csv(outdir / "metrics.tsv", sep="\t", index=False)
    for name, parts in (("roc_curves", curves), ("budget_curves", budgets), ("pints_paper_roc_points", replica)):
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(outdir / f"{name}.tsv.gz", sep="\t", index=False)
    (outdir / "done.json").write_text(json.dumps({"lib": lib, "n_metric_rows": len(metrics)}))
    return outdir
