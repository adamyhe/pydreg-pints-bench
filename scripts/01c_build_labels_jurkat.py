"""Phase-3 (Jurkat) labels. The tss_pints part needs PINTS calls on Jurkat_CAGE."""
import json

from ppbench.labels import BUILT, build_jurkat

if __name__ == "__main__":
    mf = BUILT / "manifest.json"
    m = json.loads(mf.read_text())
    new = {}
    build_jurkat(new)
    m.update(new)
    mf.write_text(json.dumps(m, indent=2))
    print(json.dumps(new, indent=1))
