#!/bin/bash
# Make sure the kentUtils programs the pipeline uses are available. If they are not on PATH (or at
# KENT_BIN / /programs/kentUtils/bin), download UCSC's static Linux x86_64 builds into tools/kent/.
# UCSC binaries are free for academic and non-profit use; commercial use needs a licence (genome.ucsc.edu).
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$PROJ/scripts/site.sh"
TOOLS="bigWigToBedGraph bedGraphToBigWig bigWigMerge bigWigInfo twoBitToFa liftOver"
missing=""
for t in $TOOLS; do [[ -x "$KENT_BIN/$t" ]] || missing+=" $t"; done
if [[ -z "$missing" ]]; then echo "kentUtils: all present in $KENT_BIN"; exit 0; fi
if [[ "$KENT_BIN" != "$PROJ/tools/kent" ]]; then
    echo "kentUtils missing from KENT_BIN=$KENT_BIN:$missing" >&2
    echo "Install them there, or unset KENT_BIN to download into tools/kent." >&2
    exit 1
fi
[[ "$(uname -sm)" == "Linux x86_64" ]] || { echo "auto-download supports Linux x86_64 only; install kentUtils and set KENT_BIN" >&2; exit 1; }
mkdir -p "$KENT_BIN"
for t in $missing; do
    echo "fetch $t"
    curl -fsSL --retry 5 -o "$KENT_BIN/$t.part" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$t"
    chmod +x "$KENT_BIN/$t.part"; mv "$KENT_BIN/$t.part" "$KENT_BIN/$t"
done
echo "kentUtils installed in $KENT_BIN"
