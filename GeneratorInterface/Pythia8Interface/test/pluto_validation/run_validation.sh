#!/bin/bash
# PlutoDecayer vs. theory papers: generate MC, overlay on the reference curves.
# Run from anywhere inside a CMSSW work area (cmsenv done):
#   ./run_validation.sh [output_dir]        # default: ./output
# Produces output/paper_compare_1511.png and output/paper_compare_0705.png.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$HERE/output}"
mkdir -p "$OUT"
: "${CMSSW_BASE:?run cmsenv first (need CMSSW_BASE and root-config)}"

echo "[1/3] compiling compare_papers.cc (uses the headers in GeneratorInterface/Pythia8Interface/interface)"
g++ -std=c++17 -O2 -I"$CMSSW_BASE/src" $(root-config --cflags) \
    "$HERE/compare_papers.cc" $(root-config --libs) -o "$OUT/compare_papers"

echo "[2/3] generating MC histograms (~30 s)"
"$OUT/compare_papers" "$OUT/pluto_mc.csv"

echo "[3/3] merging with reference curves and plotting"
python3 "$HERE/make_comparison_data.py" "$OUT/pluto_mc.csv" "$OUT/comparison_data.json"
python3 "$HERE/plot_comparison_png.py" "$OUT/comparison_data.json" "$OUT"
