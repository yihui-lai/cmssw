# PlutoDecayer shape validation against theory papers

Checks that the decay-kinematics models in `PlutoDecayer` reproduce the
published theoretical spectra of eta/eta' -> l+l-gamma, l+l-l+l-, l+l-pi+pi-.
No CMSSW build, `cmsRun` or Pythia is needed: `compare_papers.cc` calls the
header-only samplers in `GeneratorInterface/Pythia8Interface/interface/`
directly (only ROOT's `TLorentzVector` is used), so it always tests the
physics that is currently in the tree.

## Run

```
cmsenv                       # any CMSSW area containing this package
cd GeneratorInterface/Pythia8Interface/test/pluto_validation
./run_validation.sh          # ~30 s; results in ./output/
```

Outputs (`output/`, not tracked by git):

| file | content |
|---|---|
| `paper_compare_1511.png` | 12 panels vs. arXiv:1511.04916: eta/eta' -> `2mugamma`/`2mu2e`/`4e`/`4mu` |
| `paper_compare_0705.png` | 4 panels vs. arXiv:0705.0954 Fig. 7: eta/eta' -> `2e2pi`/`2mu2pi` |
| `comparison_data.json`, `pluto_mc.csv` | reference and MC curves behind the plots |

Line = reference curve digitized from the paper; dots = `PlutoDecayer` MC.
`PlutoDecayer` predicts **shapes only** (forced exclusive decays, no branching
fractions), so each MC histogram is multiplied by one constant to match the
reference area; the constants are printed by `make_comparison_data.py`.
Compare shapes, not heights.

## Files

| file | role |
|---|---|
| `run_validation.sh` | driver: compile, generate, merge, plot |
| `compare_papers.cc` | generates the 16 MC histograms (40k events each, fixed seed) |
| `make_comparison_data.py` | reads `csv/` + MC, area-matches, writes the JSON; holds the series -> csv file/column mapping |
| `plot_comparison_png.py` | draws the PNGs (PIL only, no matplotlib) |
| `csv/` | reference curves digitized from the papers' figures (one file per curve) |

## What to expect (current physics)

- eta channels and eta' `2mugamma`/`2mu2e`: shapes agree with the papers,
  including the omega peak (~0.78 GeV) and the rho/omega bump in eta'
  `2mu2e` (implemented in `interface/PlutoEtaPrimeTFF.h`).
- eta' `4e`/`4mu` use a simpler rho-pole model; the small omega bump in
  the paper's eta' `4e` curve is not reproduced.
- Just above the omega peak (0.80-0.85 GeV) the eta' `2mugamma` MC is 25-40%
  below the reference (known residual, see `doc/PlutoDecayer.md`).
- The steep, nearly vertical strokes at the left edge of some reference
  curves are the real threshold rise digitized at finite resolution.
- Peak-position ratios are a poor statistic for the electron channels, which
  peak at threshold; judge them by the overlay.

## Adding a channel

Add the CSV to `csv/`, add its line to `CHANNELS` in `make_comparison_data.py`,
add a histogram block with the same series name in `compare_papers.cc`, and
add a panel entry to `P1511` or `P0705` in `plot_comparison_png.py`.
