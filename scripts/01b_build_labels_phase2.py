"""Phase-2 labels (hg38 dreg_aux for HCT116/A549; PINTS-on-PRO-cap labels). Run after PINTS."""
import json

from ppbench.labels import build_phase2

if __name__ == "__main__":
    m = build_phase2()
    print(json.dumps({k: v for k, v in m.items() if "hg38" in k or "tss_pints" in k or "check/" in k}, indent=2))
