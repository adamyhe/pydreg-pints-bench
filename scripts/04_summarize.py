"""Aggregate results/eval/*/metrics.tsv into results/summary/*.tsv and a markdown digest."""
import glob

import pandas as pd

from ppbench.calls import RESULTS

KEY = ["lib", "family", "variant", "tool", "rank", "flavour", "metric"]


def main():
    parts = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(RESULTS / "eval" / "*" / "metrics.tsv")))]
    if not parts:
        raise SystemExit("no metrics found")
    m = pd.concat(parts, ignore_index=True)
    m["value"] = pd.to_numeric(m["value"], errors="coerce")
    agg = m.groupby(KEY, dropna=False)["value"].agg(["mean", "std", "count"]).reset_index()
    out = RESULTS / "summary"
    out.mkdir(parents=True, exist_ok=True)
    agg.to_csv(out / "metrics_by_seed_mean.tsv", sep="\t", index=False)

    headline = agg[agg.metric.isin(["auroc", "auprc", "pauroc_fpr0.1", "tpr_at_fpr0.05", "pints_paper_auc",
                                    "op_tpr", "op_fpr", "op_precision", "op_n_calls"])]
    wide = headline.pivot_table(index=["family", "variant", "lib", "metric"],
                                columns=["tool", "rank", "flavour"], values="mean")
    wide.to_csv(out / "headline_wide.tsv", sep="\t")
    with open(out / "headline.md", "w") as fh:
        fh.write(wide.round(3).to_string())
    print(wide.round(3).to_string())


if __name__ == "__main__":
    main()
