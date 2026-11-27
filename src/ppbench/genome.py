"""Paths, chromosome sets, BED I/O, blacklist and liftOver helpers.

Interval arithmetic uses bioframe (pandas); no pybedtools / pyBigWig.
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import bioframe as bf
import numpy as np
import pandas as pd

PROJ = Path(os.environ.get("PPBENCH_PROJ", Path(__file__).resolve().parents[2]))
DATA = PROJ / "data"
REF = DATA / "ref"


def _site_env() -> dict:
    """config/site.env (KEY=VALUE, shell-style) overlaid by the process environment."""
    env, f = {}, PROJ / "config" / "site.env"
    if f.exists():
        for line in f.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                k, v = line.split("=", 1)
                env[k.strip()] = os.path.expandvars(v.strip().strip("'\""))
    env.update({k: v for k, v in os.environ.items() if v})
    return env


SITE = _site_env()


def kent(tool: str) -> str:
    """Path to a kentUtils binary: $KENT_BIN, else PATH, else /programs/kentUtils/bin, else tools/kent."""
    if SITE.get("KENT_BIN"):
        return str(Path(SITE["KENT_BIN"]) / tool)
    site = Path("/programs/kentUtils/bin") / tool
    return shutil.which(tool) or (str(site) if site.exists() else str(PROJ / "tools" / "kent" / tool))


LIFTOVER = SITE.get("LIFTOVER_BIN") or kent("liftOver")
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
BED3 = ["chrom", "start", "end"]


def chrom_sizes(genome: str) -> pd.Series:
    s = pd.read_csv(REF / f"{genome}.chrom.sizes", sep="\t", header=None, names=["chrom", "size"])
    s = s[s.chrom.isin(CHROMS)]
    return s.set_index("chrom")["size"]


def genome_view(genome: str) -> pd.DataFrame:
    s = chrom_sizes(genome)
    return pd.DataFrame({"chrom": s.index, "start": 0, "end": s.values})


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _sniff_sep(path) -> str:
    """Tab if the first data line has a tab (fields may then contain spaces, e.g. the ENCODE
    blacklist's 'Low Mappability'); otherwise whitespace (R write.table output on the dREG hub)."""
    with _open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                return "\t" if "\t" in line else r"\s+"
    return "\t"


def read_bed(path, names=None, sep=None) -> pd.DataFrame:
    """Read a BED-like file; tolerates header lines ("V1 V2 V3", "track", "#")."""
    sep = sep or _sniff_sep(path)
    with _open(path) as fh:
        df = pd.read_csv(fh, sep=sep, header=None, comment="#", dtype=str,
                         engine="c" if sep == "\t" else "python")
    df = df[pd.to_numeric(df[1], errors="coerce").notna()].copy()
    df[1] = df[1].astype(np.int64)
    df[2] = df[2].astype(np.int64)
    cols = list(BED3) + [f"c{i}" for i in range(3, df.shape[1])]
    if names:
        cols[: len(names)] = names
    df.columns = cols
    return df.reset_index(drop=True)


def write_bed(df: pd.DataFrame, path, cols=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df[cols] if cols else df
    out.to_csv(path, sep="\t", header=False, index=False, compression="gzip" if str(path).endswith(".gz") else None)


def clean(df: pd.DataFrame, genome: str, blacklist: str = "drop") -> pd.DataFrame:
    """Keep canonical chromosomes, drop empty / out-of-range intervals; sort.

    blacklist: 'drop' removes elements touching the ENCODE blacklist (for element sets);
    'subtract' cuts the blacklist out of large regions (for negative-region sets); 'none'.
    """
    sizes = chrom_sizes(genome)
    df = df[df.chrom.isin(CHROMS)].copy()
    df = df[(df.end > df.start) & (df.start >= 0)]
    df = df[df.end <= df.chrom.map(sizes)]
    if blacklist != "none":
        bl = read_bed(REF / f"{genome}-blacklist.v2.bed.gz")[BED3]
        if blacklist == "drop":
            df = df[~overlaps_any(df, bl)]
        else:
            df = bf.subtract(df[BED3], bl)[BED3]
    return bf.sort_bedframe(df).reset_index(drop=True)


def overlaps_any(df: pd.DataFrame, other: pd.DataFrame) -> np.ndarray:
    """Boolean mask: does each interval of `df` overlap (>=1 bp) any interval of `other`?"""
    if len(df) == 0 or len(other) == 0:
        return np.zeros(len(df), bool)
    cnt = bf.count_overlaps(df[BED3].reset_index(drop=True), other[BED3])
    return cnt["count"].to_numpy() > 0


def liftover(df: pd.DataFrame, src: str, dst: str, min_match: float = 0.95) -> pd.DataFrame:
    """Lift intervals with UCSC liftOver; extra columns are carried through. Unmapped rows are dropped."""
    if src == dst:
        return df.copy()
    chain = REF / f"{src}To{dst[0].upper()}{dst[1:]}.over.chain.gz"
    tmp = df.reset_index(drop=True).copy()
    tmp["_id"] = np.arange(len(tmp))
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
        inp, out, un = (Path(td) / n for n in ("in.bed", "out.bed", "un.bed"))
        tmp[BED3 + ["_id"]].to_csv(inp, sep="\t", header=False, index=False)
        subprocess.run([LIFTOVER, f"-minMatch={min_match}", str(inp), str(chain), str(out), str(un)],
                       check=True, capture_output=True)
        lifted = pd.read_csv(out, sep="\t", header=None, names=BED3 + ["_id"])
    rest = tmp.drop(columns=BED3).set_index("_id")
    lifted = lifted.join(rest, on="_id").drop(columns="_id")
    return lifted


def merge(df: pd.DataFrame, min_dist: int = 0) -> pd.DataFrame:
    if len(df) == 0:
        return df[BED3].copy()
    return bf.merge(df[BED3], min_dist=min_dist)[BED3]
