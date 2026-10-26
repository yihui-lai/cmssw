# Example and test of the exclusive eta/eta' decayer (interface/EtaDecayer.h): minimum-bias events in which every
# eta decays to mu+ mu- e+ e-, run as in production through ExternalGeneratorFilter. The decayer is enabled by the
# EtaDecayer:* Pythia8 settings, see EtaDecayer.h. The PythiaDauVFilter keeps events with an eta -> mu- mu+ e- e+,
# so that test_EtaDecayer.sh can check that the decayer is active.
# For P -> mu+ mu- A or pi+ pi- A, also define the boson A as a Pythia particle, e.g. a 0.1 GeV dark photon
# decaying to e+ e-:
#   '9000001:all = A void 3 0 0 0.1 0. 0. 0. 0.', '9000001:addChannel = 1 1. 0 11 -11',
#   'EtaDecayer:mode = 2muA', 'EtaDecayer:A_id = 9000001', 'EtaDecayer:A_type = darkphoton'
# Run with onePerEvent=1 to force one eta per event (EtaDecayer:onePerEvent = on) instead of every eta.
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from GeneratorInterface.Core.ExternalGeneratorFilter import ExternalGeneratorFilter

options = VarParsing()
options.register("onePerEvent", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "force one eta per event")
options.parseArguments()

process = cms.Process("TEST")
process.source = cms.Source("EmptySource")
process.maxEvents.input = 10

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.generator = ExternalGeneratorFilter(cms.EDFilter("Pythia8GeneratorFilter",
    comEnergy = cms.double(13600.),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring(
            'SoftQCD:nonDiffractive = on',
            'EtaDecayer:filter = on',
            'EtaDecayer:parent = 221',
            'EtaDecayer:mode = 2mu2e',
            'EtaDecayer:onePerEvent = ' + ('on' if options.onePerEvent else 'off'),
        ),
        parameterSets = cms.vstring('processParameters')
    )
))

process.decayfilter = cms.EDFilter("PythiaDauVFilter",
    ParticleID = cms.untracked.int32(221),
    NumberDaughters = cms.untracked.int32(4),
    DaughterIDs = cms.untracked.vint32(13, -13, 11, -11),
    MinPt = cms.untracked.vdouble(0., 0., 0., 0.),
    MinEta = cms.untracked.vdouble(-99., -99., -99., -99.),
    MaxEta = cms.untracked.vdouble(99., 99., 99., 99.)
)

process.options.wantSummary = True
process.p = cms.Path(process.generator * process.decayfilter)
