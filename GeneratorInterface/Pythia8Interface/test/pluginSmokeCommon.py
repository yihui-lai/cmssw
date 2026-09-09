"""Five-event, single-stream Pythia plugin fixture; no jet-flavour dependency."""
import FWCore.ParameterSet.Config as cms


def make_process(commands, events, output, electron_positron=False):
    process = cms.Process("PLUGINTEST")
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.options = cms.untracked.PSet(
        numberOfThreads=cms.untracked.uint32(1),
        numberOfStreams=cms.untracked.uint32(1))
    process.source = cms.Source("EmptySource")
    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(events))
    process.RandomNumberGeneratorService = cms.Service(
        "RandomNumberGeneratorService",
        generator=cms.PSet(initialSeed=cms.untracked.uint32(12345),
                           engineName=cms.untracked.string("HepJamesRandom")))
    process.generator = cms.EDFilter(
        "Pythia8GeneratorFilter", comEnergy=cms.double(91.2 if electron_positron else 13600.),
        maxEventsToPrint=cms.untracked.int32(0),
        pythiaPylistVerbosity=cms.untracked.int32(0),
        pythiaHepMCVerbosity=cms.untracked.bool(False),
        PythiaParameters=cms.PSet(parameterSets=cms.vstring("processParameters"),
                                 processParameters=cms.vstring(*commands)))
    if electron_positron:
        process.generator.ElectronPositronInitialState = cms.bool(True)
    process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
    process.genParticles.src = "generator:unsmeared"
    process.path = cms.Path(process.generator * process.genParticles)
    process.out = cms.OutputModule(
        "PoolOutputModule", fileName=cms.untracked.string(output),
        outputCommands=cms.untracked.vstring("drop *", "keep *_generator_*_*", "keep *_genParticles_*_*"))
    process.end = cms.EndPath(process.out)
    return process
