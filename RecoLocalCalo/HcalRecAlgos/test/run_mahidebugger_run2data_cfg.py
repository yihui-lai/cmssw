# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3 --conditions auto:run2_data_promptlike -s RAW2DIGI,L1Reco,RECO --runUnscheduled --process reRECO --data --era Run2_2017 --eventcontent RECO --hltProcess reHLT --scenario pp --datatier RECO --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2017 -n 100 --filein file:step2.root --fileout file:step3.root --no_exec
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
#from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Run3_cff import Run3

#process = cms.Process('reRECO',eras.Run2_2018)


options = VarParsing.VarParsing()
options.register('inputFile',
                 "file:/eos/cms/store/user/jaehyeok/store_data_Run2018A_JetHT_RAW_v1_000_316_944_00000_2EF71399-6B64-E811-9B71-FA163ED59971.root", # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "input file")
options.register('maxEvents',
                 10, # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "maxEvents")
options.register('outputFile',
                  "mahidebugger.root",       
#                 "file:/eos/cms/store/user/jaehyeok/store_data_Run2018A_JetHT_RAW_v1_000_316_944_00000_2EF71399-6B64-E811-9B71-FA163ED59971.root", # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "output file")
options.parseArguments()


process = cms.Process('reRECO',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )


process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring('file:{0}'.format(options.inputFile)),
#'file:/eos/user/y/yilai/OOTPU/phase_scan_off_TimeSlew/Job88_timePhase6_syncPhase1_OOTPU_step2.root'),
#'file:/afs/cern.ch/user/y/yilai/OOTPU/CMSSW_11_1_0_pre8/src/11634.0_TTbar_14TeV+TTbar_14TeV_TuneCP5_2021_GenSimFull+DigiFull_2021+RecoFull_2021+HARVESTFull_2021+ALCAFull_2021/step3_inMINIAODSIM.root'),
#'file:/eos/cms/store/user/jaehyeok/store_data_Run2018A_JetHT_RAW_v1_000_316_944_00000_2EF71399-6B64-E811-9B71-FA163ED59971.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3.root'),
    outputCommands = process.RECOEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('file:{0}'.format(options.outputFile)) #"mahidebugger.root")
    #fileName = cms.string("mahidebugger_8p_pickevents.root")
    )

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_dataRun2_v12', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

process.hbheprereco.saveInfos    = cms.bool(True)

process.load("RecoLocalCalo.HcalRecAlgos.mahiDebugger_cfi")


# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.hcalDigis)
process.reconstruction_step = cms.Path(process.hbheprereco)
process.flat_step = cms.Path(process.mahiDebugger)
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.L1Reco_step = cms.Path(process.L1Reco)
#process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.flat_step)
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOoutput_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2017

#call to customisation function customisePostEra_Run2_2017 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2017(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion



#use shape 208
process.es_prefer = cms.ESPrefer('HcalTextCalibrations','es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
     input = cms.VPSet(
        cms.PSet(
           object = cms.string('RecoParams'),
           file   = cms.FileInPath('HcalRecoParams_2021_shape208.txt')
        )
     )
)




# print process time
#process.FastTimerService = cms.Service("FastTimerService",
#    dqmLumiSectionsRange = cms.untracked.uint32(2500),
#    dqmMemoryRange = cms.untracked.double(1000000.0),
#    dqmMemoryResolution = cms.untracked.double(5000.0),
#    dqmModuleMemoryRange = cms.untracked.double(100000.0),
#    dqmModuleMemoryResolution = cms.untracked.double(500.0),
#    dqmModuleTimeRange = cms.untracked.double(40.0),
#    dqmModuleTimeResolution = cms.untracked.double(0.2),
#    dqmPath = cms.untracked.string('HLT/TimerService'),
#    dqmPathMemoryRange = cms.untracked.double(1000000.0),
#    dqmPathMemoryResolution = cms.untracked.double(5000.0),
#    dqmPathTimeRange = cms.untracked.double(100.0),
#    dqmPathTimeResolution = cms.untracked.double(0.5),
#    dqmTimeRange = cms.untracked.double(2000.0),
#    dqmTimeResolution = cms.untracked.double(5.0),
#    enableDQM = cms.untracked.bool(True),
#    enableDQMTransitions = cms.untracked.bool(False),
#    enableDQMbyLumiSection = cms.untracked.bool(True),
#    enableDQMbyModule = cms.untracked.bool(False),
#    enableDQMbyPath = cms.untracked.bool(False),
#    enableDQMbyProcesses = cms.untracked.bool(True),
#    printEventSummary = cms.untracked.bool(False),
#    printJobSummary = cms.untracked.bool(True),
#    printRunSummary = cms.untracked.bool(True)
#    )
#process.options = cms.untracked.PSet(
#    numberOfStreams = cms.untracked.uint32(0),
#    numberOfThreads = cms.untracked.uint32(1),
#    sizeOfStackForThreadsInKB = cms.untracked.uint32(10240),
#    wantSummary = cms.untracked.bool(True)
#    )

# dump config
dump = file('dump.py', 'w')
dump.write( process.dumpPython() )
dump.close()
