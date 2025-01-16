import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## GT to be used
#options.register('globalTag','124X_dataRun3_Prompt_v4',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','112X_dataRun3_Prompt_v2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','124X_dataRun3_PromptAnalysis_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloo be used');
#options.register('globalTag','140X_dataRun3_Prompt_v2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','133X_mcRun3_2024_realistic_v10',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabbe used');
#options.register('globalTag','94X_mc2017_realistic_v11',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag toe used');
#options.register('globalTag','106X_dataRun2_v36',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
options.register('globalTag','94X_dataRun2_ReReco_EOY17_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tae used');

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

## outputFile Name
#options.register('outputFileName','ku_QCD_AOD_diag_140_gammares_v12.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');
options.register('outputFileName','ku_Met_17D_MiniAOD_diag_140_gammares_v12.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name in cmsRun');

options.register('doTwoTier',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to do twotier processing');
options.register('doDiag',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to store diagnostic info');

## parsing command line arguments
options.parseArguments()

## Define the CMSSW process
from Configuration.StandardSequences.Eras import eras
process = cms.Process(options.processName,eras.Run3)

## Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')#
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')#
#process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')#

#process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")#
#process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")#
#process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi")#
#process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")#
#process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")#
#process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
#process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('RecoLocalCalo.EcalRecProducers.ecalCPUUncalibRecHitProducer_cfi')#
#process.load("RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi")
#process.load('GammaResTool.GammaResTool.jwk_ku_ecalLocalRecoSequence_cff')
#process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff' )

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')#

process.load('Configuration.StandardSequences.EndOfProcess_cff')#

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")#
#process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## Define the input source
#eventList = open(options.rlelist,'r')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'file:jwk_reco_data_DIGI2RAW.root'),

        #'root://cms-xrd-global.cern.ch//store/mc/Run3Winter24MiniAOD/DYto2L-4Jets_MLL-50_1J_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/133X_mcRun3_2024_realistic_v10-v2/2830000/000a0b08-4970-4a08-bbd1-69c4ae918e66.root',
        #'file:967aebe0-e567-4139-9f91-d9e67f6b2ace.root'
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/43f9417c-f2bc-4e2c-a114-ac6d5c5b7052.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/444fe16d-0217-4915-92ca-101fdad77998.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/67beed64-8d00-4f97-93b2-72292c720443.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/6f38df6b-6ddb-400d-85ba-e8da46fb9500.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/967aebe0-e567-4139-9f91-d9e67f6b2ace.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/9aaff3a9-b4d8-4025-9dfa-a64751429a80.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/9ba1c87f-217b-46b0-b44b-7a6bc0a2ba17.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/b15ef1a9-d310-4a98-9da3-edbd85b93b01.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/b8a7f31d-8c9c-4254-99aa-38fd0d6d638f.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/d868dddd-3fff-456a-bc1f-09fba3339620.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/e0bd624f-b25d-48f0-890c-1340ea1eb666.root',
        #'root://cms-xrd-global.cern.ch//eos/cms/tier0/store/data/Run2024E/EGamma0/MINIAOD/PromptReco-v2/000/381/384/00000/f7b05c22-6939-429c-9a83-51cf69cdedc1.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2018B/EGamma/MINIAOD/15Feb2022_UL2018-v1/2560000/11B9426C-D91F-304B-87A1-054DD706ED4C.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2018B/EGamma/MINIAOD/UL2018_MiniAODv2-v1/50000/C2FAC52C-9570-454C-B93D-865F3502E49A.root',

        'root://cms-xrd-global.cern.ch//store/data/Run2017D/MET/MINIAOD/17Nov2017-v1/60000/C2906A84-83EC-E711-AA33-008CFAC93BD8.root',

	),
#    secondaryFileNames = cms.untracked.vstring(

#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/116043a5-e71e-4df8-8014-34fdb3391301.root'

#	),
    #eventsToProcess = cms.untracked.VEventRange(eventList)
)


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFileName))
				   
## Decide which label to use for MET Flags
#if   options.isMC : triggerFlagsProcess = "PAT"
#else              : triggerFlagsProcess = "RECO"
triggerFlagsProcess = "RECO"

## generate track collection at miniAOD
from PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi import unpackedTracksAndVertices
process.unpackedTracksAndVertices = unpackedTracksAndVertices.clone()

# Make the tree 
process.tree = cms.EDAnalyzer("GammaResTool",
   ## additional collections
   #triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
   #triggerObjects = cms.InputTag("slimmedPatTrigger"),
   ## met filters
   #inputFlags       = cms.string(options.inputFlags),
   #triggerFlags     = cms.InputTag("TriggerResults", "", triggerFlagsProcess),
   #ecalBadCalibFlag = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),			      
   ## tracks
   #tracks = cms.InputTag("unpackedTracksAndVertices"),
   ## vertices
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   ## rho
   #rho = cms.InputTag("fixedGridRhoFastjetAll"), #fixedGridRhoAll
   ## METs
   #mets = cms.InputTag("slimmedMETsModifiedMET"),
   #mets = cms.InputTag("slimmedMETs"),
   ## jets
   #jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
   #jets = cms.InputTag("slimmedJets"),
   ## electrons
   electrons = cms.InputTag("slimmedElectrons"),
   ## muons
   #muons = cms.InputTag("slimmedMuons"),
   ## photons
   gedPhotons = cms.InputTag("slimmedPhotons"),
   ootPhotons = cms.InputTag("slimmedOOTPhotons"),
   ## ecal recHits
   recHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
   recHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits"),

   ## do two tier reconstruction of second rechit collection
   doTwoTier = cms.bool(options.doTwoTier), 
   doDiag = cms.bool(options.doDiag),

   ## ecal Multi rechits

   kuCCStcRecHitsEB = cms.InputTag("none", "none"),
   kuCCStcRecHitsEE = cms.InputTag("none", "none"),

   kuRtStcRecHitsEB = cms.InputTag("none", "none"),
   kuRtStcRecHitsEE = cms.InputTag("none", "none"),

   ku_uncalibratedRecHitsEB = cms.InputTag("none","none"),
   ku_uncalibratedRecHitsEE = cms.InputTag("none","none"),

)


# Set up the path
process.tree_step = cms.EndPath(
	process.unpackedTracksAndVertices +
	process.tree
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.content_step = cms.Path(process.content)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule( process.tree_step )

### Extra bits from other configs
process.options = cms.untracked.PSet(
    #numberOfThreads = cms.untracked.uint32(4),
    #numberOfStreams = cms.untracked.uint32(4),
    #####SkipEvent = cms.untracked.vstring('ProductNotFound'),
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    #wantSummary = cms.untracked.bool(True)
)

#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#if options.runUnscheduled : 
#	process = convertToUnscheduled(process)

#from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
#if options.deleteEarly :
#	process = customiseEarlyDelete(process)
