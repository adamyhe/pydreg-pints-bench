"""GC / CpG content of intervals from UCSC .2bit genomes (kentUtils twoBitToFa; no pysam/pybedtools)."""
from __future__ import annotations

import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

from .genome import BED3, REF

TWOBIT2FA = "/programs/kentUtils/bin/twoBitToFa"


def seq_stats(df: pd.DataFrame, genome: str) -> pd.DataFrame:
    """Per interval: gc (fraction G+C of non-N bases), cpg_oe (observed/expected CpG), n_frac."""
    iv = df[BED3].drop_duplicates().reset_index(drop=True)
    iv["name"] = np.arange(len(iv)).astype(str)
    gc = np.full(len(iv), np.nan)
    oe = np.full(len(iv), np.nan)
    nf = np.full(len(iv), np.nan)
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
        bed, fa = f"{td}/in.bed", f"{td}/out.fa"
        iv[BED3 + ["name"]].to_csv(bed, sep="\t", header=False, index=False)
        subprocess.run([TWOBIT2FA, "-noMask", f"-bed={bed}", str(REF / f"{genome}.2bit"), fa], check=True)
        i, buf = None, []

        def flush():
            if i is None:
                return
            s = "".join(buf).upper()
            L = len(s)
            n = s.count("N")
            c, g = s.count("C"), s.count("G")
            eff = L - n
            if eff > 0:
                gc[i] = (c + g) / eff
                exp = c * g / eff
                oe[i] = s.count("CG") / exp if exp > 0 else 0.0
            nf[i] = n / L if L else np.nan

        with open(fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    i, buf = int(line[1:].split(":")[0].split()[0]), []
                else:
                    buf.append(line.strip())
            flush()
    iv["gc"], iv["cpg_oe"], iv["n_frac"] = gc, oe, nf
    return df[BED3].merge(iv.drop(columns="name"), on=BED3, how="left")
