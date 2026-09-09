# Pluto external-decayer integration

This branch includes the adapter/CMake project and model source under
`test/externalPlugins/pluto`, its upstream compatibility patch and standalone
tests, and a five-event CMSSW config. It is independent of the jet-flavour and
more-showers branches. No new CMSSW hadronizer API is needed: Pythia loads the
external DecayHandler plugin through `Init:plugins`.

The library is built externally, not by the normal SCRAM package build. Use
`DickyChant/cmsdist:integration/pluto` with the source bundle, or build the
included CMake project against Pluto 6.2.2 commit
`a242b7fc1246d4ee639dd64e9afae7f8eae1fd30`, ROOT and Pythia 8.317. Apply
`externalPlugins/patches/pluto6-cmssw-rng-lifetime.patch` to the upstream Pluto
source first. The cmsdist source URL remains a local candidate, not a centrally
published source release. Adapter/model source is included here for review.

With `libCMSPluto.so`, Pluto's shared library/dictionary and dependencies
available through the external SCRAM tool or runtime library path:

```bash
cmsRun GeneratorInterface/Pythia8Interface/test/plutoDecayer_cfg.py parent=eta mode=2mu2e
cmsRun GeneratorInterface/Pythia8Interface/test/plutoDecayer_cfg.py parent=etaprime mode=4mu
```

Both parents support 4e, 2mu2e, 4mu, 2e2pi and 2mu2pi. These are forced exclusive
finite-mass LO constant-form-factor models, not branching-rate predictions.
Identical four-lepton channels include exchange interference; mixed leptons do
not have identical-lepton exchange. Tracks mean charged pions, with magnetic-only
dynamics. Nonconstant form factors/radiative corrections are not supplied.
The native electron+pion model is quarantined; do not silently enable it.
See the included model notes and tests for formulas, assumptions and checks.
