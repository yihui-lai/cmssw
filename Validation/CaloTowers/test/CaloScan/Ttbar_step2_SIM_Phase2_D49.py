import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('TEST',Phase2C9)

### RANDOM setting (change last digit(s) to make runs different !)
process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
# ---- D49 workflow
####process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load("SimG4Core.Application.g4SimHits_cfi")
process.load("SimCalorimetry.Configuration.ecalDigiSequence_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load("SimGeneral.PileupInformation.AddPileupSummary_cfi")
process.load('Configuration.StandardSequences.Digi_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
#process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# ---- D49 workflow
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.g4SimHits.UseMagneticField = False
process.load("DQMServices.Core.DQMStore_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")

# import of standard configurations
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
#process.load("SimGeneral.MixingModule.mixNoPU_cfi")
#process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.Digi_cff')
#process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
#process.load('Configuration.StandardSequences.SimL1Emulator_cff')
#process.load('Configuration.StandardSequences.DigiToRaw_cff')
#process.load('HLTrigger.Configuration.HLT_Fake2_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load('Configuration/StandardSequences/RawToDigi_cff')

#process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("DQMServices.Core.DQMStore_cfi")
#process.load("DQMServices.Components.MEtoEDMConverter_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100) 
)
# Input source
process.source = cms.Source("PoolSource",
    firstEvent = cms.untracked.uint32(1), 
    fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/work/y/yilai/phase2_zs/CMSSW_11_3_0_pre5/src/Validation/CaloTowers/test/CaloScan/Gen_ttbar/pi50_INPUTXXX.root',
    )
) 

process.FEVT = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring('drop *', 'keep *_MEtoEDMConverter_*_*'),
     splitLevel = cms.untracked.int32(0),
     fileName = cms.untracked.string("output.root")
)

process.VtxSmeared.SigmaX = 0.00001
process.VtxSmeared.SigmaY = 0.00001
process.VtxSmeared.SigmaZ = 0.00001

process.load("Validation.HcalHits.HcalSimHitsValidation_cfi")
process.HcalSimHitsAnalyser.outputFile = cms.untracked.string('HcalSimHitsValidation.root')

process.load("Validation.HcalDigis.HcalDigisParam_cfi")
process.hcaldigisAnalyzer.outputFile = cms.untracked.string('HcalDigisValidationRelVal.root')

process.load("Validation.HcalRecHits.HcalRecHitParam_cfi")

process.load("Validation.CaloTowers.CaloTowersParam_cfi")
process.calotowersAnalyzer.outputFile = cms.untracked.string('CaloTowersValidationRelVal.root')

process.load("Validation.RecoParticleFlow.pfClusterValidation_cfi")
process.load("Validation.RecoMET.METValidation_cfi")

#------------- CUSTOMIZATION - replace hbhereco with hbheprereco (everywhere)
delattr(process,"hbhereco")
process.hbhereco = process.hbheprereco.clone()
process.hcalLocalRecoTask = cms.Task(process.hbhereco,process.hfprereco,process.hfreco,process.horeco,process.zdcreco)

#---------- PATH
# -- NB: for vertex smearing the Label should be: "unsmeared" 
# for GEN produced since 760pre6, for older GEN - just "": 


process.VtxSmeared.src = cms.InputTag("generator", "")
process.generatorSmeared = cms.EDProducer("GeneratorSmearedProducer")
process.g4SimHits.Generator.HepMCProductLabel = cms.InputTag('VtxSmeared')

#watch which is running
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#--- ZS (re-)set 
#process.simHcalDigis.HBlevel = 8
#process.simHcalDigis.HElevel = FILE_NAME
#process.simHcalDigis.HOlevel = 21
#process.simHcalDigis.HFlevel = -99
#process.simHcalDigis.useConfigZSvalues = 1
#--- ZS reset

#process.es_ascii = cms.ESSource('HcalTextCalibrations',
#   input = cms.VPSet(
#     cms.PSet(
#       object = cms.string('ZSThresholds'),
#       file   = cms.FileInPath('FILE_NAME')
#       file   = cms.FileInPath('HcalZSThresholds_2024_v1.0_mc_mean_0.5sigma.txt')
#       file   = cms.FileInPath('HcalZSThresholds_2024_v1.0_mc_mean_1.0sigma.txt')
#       file   = cms.FileInPath('HcalZSThresholds_2024_v1.0_mc_mean_1.5sigma.txt')
#       file   = cms.FileInPath('HcalZSThresholds_2024_v1.0_mc_mean_2.0sigma.txt')
#       file   = cms.FileInPath('HcalZSThresholds_2024_v1.0_mc_mean_2.5sigma.txt')
#       file   = cms.FileInPath('HcalZSThresholds_2024_v1.0_mc_mean_3.0sigma.txt')
#       file   = cms.FileInPath('HcalZSThresholds_2024_v1.0_mc_mean_30sigma.txt')
#     ),
#   )
#)
#process.es_prefer = cms.ESPrefer('HcalTextCalibrations','es_ascii' )


#process.hbhereco.algorithm.useM3 = cms.bool(False)
#process.hbhereco.algorithm.useMahi = cms.bool(False)
#process.calotowermaker.HBThreshold = cms.double(-1000.)
#process.calotowermaker.HBThreshold1 = cms.double(-1000.)
#process.calotowermaker.HBThreshold2 = cms.double(-1000.)

# For PU
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
process.mix.input.nbPileupEvents.averageNumber = cms.double(200.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-3)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/2d25844c-b7c5-4ada-8de6-71d396ff25f1.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/44cbf8cf-a514-4795-9a3f-25dd9e877ff5.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/4fa4dfb0-4631-4ff3-9907-70ca007b0dce.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/50852062-f221-4eb1-8851-458c5aca0060.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/52aa2962-c03a-459a-9983-ed600f591de7.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/68682f83-a6bd-45f1-9685-c2b6c5fb8e94.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/9edb9a7d-c16a-4edc-a49c-850e6f7a5baa.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/a0ec3c74-57b5-46fe-b0c4-5d8f66dc1443.root','/store/relval/CMSSW_11_2_0_pre8/RelValMinBias_14TeV/GEN-SIM/112X_mcRun4_realistic_v3_2026D49noPU-v1/00000/ce8f16fe-caa2-4c8c-a931-4f8129d071fc.root',])
process.mix.digitizers = cms.PSet(process.theDigitizersValid)


process.p = cms.Path(
  process.VtxSmeared *
 process.generatorSmeared *
  process.g4SimHits *
 process.mix *
 process.ecalDigiSequence * 
 process.hcalDigiSequence *
 process.addPileupInfo *
 process.bunchSpacingProducer *
 process.ecalPacker *
 process.esDigiToRaw *
 process.hcalRawData *
 process.rawDataCollector *
 process.ecalDigis *
 process.ecalPreshowerDigis *
 process.hcalDigis *
 process.castorDigis *
 process.calolocalreco *
 process.caloTowersRec *
 process.hcalnoise *
 process.particleFlowCluster *
#
 process.HcalSimHitsAnalyser *
 process.hcaldigisAnalyzer *
 process.calotowersAnalyzer *
 process.hcalRecoAnalyzer *
 process.pfClusterValidation * 
 process.caloMet *
 process.metAnalyzer *
 process.pfMetAnalyzer *
 process.genMetTrueAnalyzer *
#
 process.MEtoEDMConverter
)

process.outpath = cms.EndPath(process.FEVT)

