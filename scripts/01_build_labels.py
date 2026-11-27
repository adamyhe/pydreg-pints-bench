"""Build the phase-1 PINTS-paper and dREG-paper label sets (data/labels/built/)."""
import json

from ppbench.labels import build_all

if __name__ == "__main__":
    print(json.dumps(build_all(), indent=2))
