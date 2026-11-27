"""GC diagnostic over every evaluated label set (results/bootstrap/*/element_scores.tsv.gz, seed 0).

Per (library, family, variant):
  gc_pos / gc_neg      median GC of positives / negatives
  cgi_pos / cgi_neg    fraction CpG-island-like (GC > 0.5 and CpG o/e > 0.6)
  auroc_gc_only        how well GC alone separates positives from negatives (0.5 = matched)
  <tool>_auroc / _gcw  tool AUROC as evaluated / with negatives reweighted to the positives' GC histogram
                       (20 pooled-quantile GC bins; importance weights, no resampling)
Writes results/summary/gc_diagnostic.tsv.
"""
import glob
import sys

import numpy as np
import pandas as pd

from ppbench.bootstrap import Ranked
from ppbench.calls import RESULTS
from ppbench.evaluate import libraries
from ppbench.gc import seq_stats

RUNS = ["pydreg_score", "pints_logp"]
NBINS = 20


def gc_weights(gc_pos, gc_neg, nbins=NBINS):
    edges = np.unique(np.quantile(np.r_[gc_pos, gc_neg], np.linspace(0, 1, nbins + 1)))
    bp = np.clip(np.searchsorted(edges, gc_pos, side="right") - 1, 0, len(edges) - 2)
    bn = np.clip(np.searchsorted(edges, gc_neg, side="right") - 1, 0, len(edges) - 2)
    fp = np.bincount(bp, minlength=len(edges) - 1) / len(bp)
    fn = np.bincount(bn, minlength=len(edges) - 1) / len(bn)
    w = np.where(fn[bn] > 0, fp[bn] / np.maximum(fn[bn], 1e-12), 0.0)
    return w * len(w) / w.sum(), float((fp[fn == 0]).sum())  # weights; positive mass in bins with no negatives


def main():
    libs = libraries().set_index("lib_id")
    tabs = []
    for f in sorted(glob.glob(str(RESULTS / "bootstrap" / "*" / "element_scores.tsv.gz"))):
        lib = f.split("/")[-2]
        t = pd.read_csv(f, sep="\t")
        t = t[t.seed <= 0]
        t["lib"], t["genome"] = lib, libs.loc[lib, "genome"]
        tabs.append(t)
    t = pd.concat(tabs, ignore_index=True)
    parts = []
    for g, sub in t.groupby("genome"):
        print("GC for", g, sub[["chrom", "start", "end"]].drop_duplicates().shape[0], "intervals", file=sys.stderr)
        parts.append(sub.merge(seq_stats(sub, g).drop_duplicates(["chrom", "start", "end"]),
                               on=["chrom", "start", "end"], how="left"))
    t = pd.concat(parts, ignore_index=True)
    t["cgi"] = (t.gc > 0.5) & (t.cpg_oe > 0.6)
    rows = []
    for (lib, fam, var), d in t.groupby(["lib", "family", "variant"]):
        d = d[d.gc.notna()]
        p, n = d[d.label == 1], d[d.label == 0]
        w, lost = gc_weights(p.gc.to_numpy(), n.gc.to_numpy())
        r = {"lib": lib, "family": fam, "variant": var, "n_pos": len(p), "n_neg": len(n),
             "gc_pos": p.gc.median(), "gc_neg": n.gc.median(), "cgi_pos": p.cgi.mean(), "cgi_neg": n.cgi.mean(),
             "auroc_gc_only": Ranked(p.gc.to_numpy(), n.gc.to_numpy()).metrics(np.ones(len(p)), np.ones(len(n)))["auroc"],
             "pos_mass_unmatched": lost}
        for run in RUNS:
            rk = Ranked(p[run].to_numpy(float), n[run].to_numpy(float))
            r[f"{run}_auroc"] = rk.metrics(np.ones(len(p)), np.ones(len(n)))["auroc"]
            r[f"{run}_auroc_gcw"] = rk.metrics(np.ones(len(p)), w)["auroc"]
            r[f"{run}_auprc"] = rk.metrics(np.ones(len(p)), np.ones(len(n)))["auprc"]
            r[f"{run}_auprc_gcw"] = rk.metrics(np.ones(len(p)), w)["auprc"]
        r["diff_auroc"] = r["pydreg_score_auroc"] - r["pints_logp_auroc"]
        r["diff_auroc_gcw"] = r["pydreg_score_auroc_gcw"] - r["pints_logp_auroc_gcw"]
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "summary" / "gc_diagnostic.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 300)
    cols = ["lib", "family", "variant", "gc_pos", "gc_neg", "cgi_pos", "cgi_neg", "auroc_gc_only",
            "pydreg_score_auroc", "pydreg_score_auroc_gcw", "pints_logp_auroc", "pints_logp_auroc_gcw",
            "diff_auroc", "diff_auroc_gcw"]
    print(out[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
