"""Evaluate pydreg and PINTS calls for one library: python scripts/03_evaluate.py LIB_ID"""
import sys

from ppbench.evaluate import evaluate_library

if __name__ == "__main__":
    for lib in sys.argv[1:]:
        print("wrote", evaluate_library(lib))
