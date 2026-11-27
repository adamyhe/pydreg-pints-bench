"""Load peak calls from pydreg and PINTS into a common table.

Every call table has columns chrom, start, end, score (higher = more confident) and is restricted
to canonical chromosomes. Two flavours per tool:

* ``default`` - the tool's own final calls at its default significance cutoff.
* ``full``    - every candidate the tool scored, so a threshold sweep can reach high recall:
    - pydreg: ``*.dREG.raw.peak.bed.gz`` (all candidate peaks before FDR filtering);
    - PINTS:  divergent pairs rebuilt from ``*_peaks.gtf`` (written with ``--debug``), which holds
      every per-strand candidate peak with its p-value. This mirrors how the PINTS paper built
      its own ROC (``ProfilerPINTS``: threshold per-strand p-values, then re-pair opposite peaks
      within 300 bp).
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import bioframe as bf
import numpy as np
import pandas as pd

from .genome import BED3, CHROMS, PROJ, _open, read_bed

RESULTS = PROJ / "results"
CLOSE = 300  # PINTS --close-threshold default and the PINTS-paper pairing distance


def _neglog10(p):
    return -np.log10(np.clip(np.asarray(p, float), 1e-300, 1.0))


def _finish(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df.chrom.isin(CHROMS) & (df.end > df.start)]
    return bf.sort_bedframe(df).reset_index(drop=True)


# ---------------------------------------------------------------- pydreg
def pydreg_dir(lib: str) -> Path:
    return RESULTS / "calls" / lib / "pydreg"


def load_pydreg(lib: str, flavour: str = "full", rank: str = "score") -> pd.DataFrame:
    """rank: 'score' = max SVR score (what the PINTS paper swept); 'logp' = -log10 raw peak p-value."""
    d = pydreg_dir(lib)
    if flavour == "full":
        df = read_bed(d / f"{lib}.dREG.raw.peak.bed.gz",
                      names=BED3 + ["svr", "prob", "smooth_mode", "original_mode", "centroid"])
    else:
        df = read_bed(d / f"{lib}.dREG.peak.full.bed.gz", names=BED3 + ["svr", "prob", "center"])
    df["svr"] = df.svr.astype(float)
    df["prob"] = df.prob.astype(float)
    df["score"] = df.svr if rank == "score" else _neglog10(df.prob)
    return _finish(df[BED3 + ["score", "svr", "prob"]])


def load_infp(lib: str) -> pd.DataFrame:
    """dREG informative positions for a library (the dREG paper's evaluation universe).
    pydreg 1.0 writes them as a bigWig (``*.dREG.infp.bw``); older runs may have a bed.gz."""
    bed = pydreg_dir(lib) / f"{lib}.dREG.infp.bed.gz"
    if bed.exists():
        return _finish(read_bed(bed)[BED3])
    import os, subprocess, tempfile
    bw = pydreg_dir(lib) / f"{lib}.dREG.infp.bw"
    with tempfile.NamedTemporaryFile(suffix=".bg", dir=os.environ.get("TMPDIR")) as tf:
        subprocess.run(["/programs/kentUtils/bin/bigWigToBedGraph", str(bw), tf.name], check=True)
        df = pd.read_csv(tf.name, sep="\t", header=None, names=BED3 + ["v"], usecols=[0, 1, 2])
    return _finish(df)


# ---------------------------------------------------------------- PINTS
def pints_dir(lib: str) -> Path:
    return RESULTS / "calls" / lib / "pints"


def _one(pattern: str) -> str:
    hits = sorted(glob.glob(pattern))
    hits = [h for h in hits if "_sig_" not in h]
    if not hits:
        raise FileNotFoundError(pattern)
    return hits[0]


def load_pints_default(lib: str) -> pd.DataFrame:
    """PINTS final bidirectional elements. They carry no numeric score, so score = +inf
    (usable only as an operating point, not for a sweep)."""
    f = _one(str(pints_dir(lib) / f"{lib}*_bidirectional_peaks.bed*"))
    df = read_bed(f, names=BED3 + ["confidence", "tss_fwd", "tss_rev"])
    df["score"] = np.inf
    return _finish(df)


_ATTR = {k: re.compile(rf"\b{k} ([^;]+)") for k in ("pval", "qval", "summit")}


def load_pints_candidates(lib: str) -> pd.DataFrame:
    """All per-strand PINTS candidate peaks (BED coords) with pval/qval/summit."""
    f = _one(str(pints_dir(lib) / f"{lib}*_peaks.gtf*"))
    rows = []
    with _open(f) as fh:
        for line in fh:
            it = line.rstrip("\n").split("\t")
            if len(it) < 9:
                continue
            attr = it[8]
            rows.append((it[0], int(it[3]) - 1, int(it[4]), it[6],
                         *(float(_ATTR[k].search(attr).group(1)) for k in ("pval", "qval", "summit"))))
    df = pd.DataFrame(rows, columns=BED3 + ["strand", "pval", "qval", "summit"])
    return _finish(df)


def pair_pints_candidates(cand: pd.DataFrame, close: int = CLOSE) -> pd.DataFrame:
    """Rebuild divergent/bidirectional elements from per-strand candidates, following
    pints.calling_engine.merge_opposite_peaks: an anchor peak on '+' pairs with an opposite-strand
    candidate overlapping [start-close, end); a '-' anchor with one overlapping [start, end+close).
    The partner is the lowest-p overlapping candidate. An element exists at p-threshold t iff its
    anchor has p <= t, so element score = -log10(anchor p).

    Pairs are NOT merged across scores: merging weak pairs into strong ones would widen every
    element at every threshold (on K562 GRO-cap, median 429 bp vs PINTS's own 179 bp). Overlap
    queries are unaffected by merging; budget counts merge at the chosen threshold
    (metrics.matched_budget). This matches PINTS_analysis ProfilerPINTS, which concatenated the
    pl/mn pairs per threshold without merging. Exact duplicate pairs keep the max score."""
    out = []
    for strand, opp in (("+", "-"), ("-", "+")):
        a = cand[cand.strand == strand].reset_index(drop=True).copy()
        b = cand[cand.strand == opp].reset_index(drop=True)
        if len(a) == 0 or len(b) == 0:
            continue
        a["a_id"] = np.arange(len(a))
        w = a[BED3 + ["a_id"]].copy()
        if strand == "+":
            w["start"] = (w.start - close).clip(lower=0)
        else:
            w["end"] = w.end + close
        hits = bf.overlap(w, b[BED3 + ["pval"]], how="inner", suffixes=("", "_b"))
        if len(hits) == 0:
            continue
        hits = hits.sort_values("pval_b").drop_duplicates("a_id")  # lowest-p partner
        anchors = a.set_index("a_id")[["start", "end", "pval"]].add_suffix("_a")
        hits = hits.join(anchors, on="a_id")
        el = pd.DataFrame({
            "chrom": hits.chrom.values,
            "start": np.minimum(hits.start_a.values, hits.start_b.values.astype(np.int64)),
            "end": np.maximum(hits.end_a.values, hits.end_b.values.astype(np.int64)),
            "score": _neglog10(hits.pval_a.values),
        })
        out.append(el)
    if not out:
        return pd.DataFrame(columns=BED3 + ["score"])
    el = pd.concat(out, ignore_index=True)
    el = el.groupby(BED3, as_index=False)["score"].max()
    return _finish(el)


def load_pints(lib: str, flavour: str = "full") -> pd.DataFrame:
    if flavour == "default":
        return load_pints_default(lib)
    return pair_pints_candidates(load_pints_candidates(lib))


def load_calls(tool: str, lib: str, flavour: str, rank: str = "score") -> pd.DataFrame:
    if tool == "pydreg":
        return load_pydreg(lib, flavour, rank)
    if tool == "pints":
        return load_pints(lib, flavour)
    raise ValueError(tool)
