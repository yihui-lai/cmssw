"""Read five-event output: python3 checkPluginSmoke.py FILE [dire]."""
import math
import sys
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.FWLiteEnabler.enable()
ROOT.gSystem.Load("libDataFormatsFWLite")
ROOT.gSystem.Load("libSimDataFormatsGeneratorProducts")
ROOT.gSystem.Load("libDataFormatsHepMCCandidate")
from DataFormats.FWLite import Events, Handle

info = Handle("io_v1::GenEventInfoProduct")
particles = Handle("std::vector<reco::GenParticle>")
weights = []
for event in Events(sys.argv[1]):
    event.getByLabel("generator", info)
    event.getByLabel("genParticles", particles)
    assert particles.product().size() > 0
    weights.append(info.product().weight())
assert len(weights) == 5 and all(math.isfinite(w) for w in weights), weights
if len(sys.argv) > 2 and sys.argv[2] == "dire":
    assert any(abs(w - 1.) > 1.e-8 for w in weights), weights
print("Five events with persisted genParticles and finite generator weights PASS:", weights)
