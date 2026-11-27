"""Paired bootstrap CIs for one or more libraries: python scripts/05_bootstrap.py [--n-boot N] LIB [LIB...]"""
import argparse

from ppbench.bootstrap import bootstrap_library

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("libs", nargs="+")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--outdir", default=None, help="override output dir (single library; for tests)")
    a = ap.parse_args()
    for lib in a.libs:
        print("wrote", bootstrap_library(lib, a.n_boot, outdir=a.outdir))
