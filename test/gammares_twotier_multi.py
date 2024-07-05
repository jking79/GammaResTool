import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## GT to be used
#options.register('globalTag','124X_dataRun3_Prompt_v4',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','112X_dataRun3_Prompt_v2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','124X_dataRun3_PromptAnalysis_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','140X_dataRun3_Prompt_Candidate_2024_05_31_21_23_47',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','140X_dataRun3_Prompt_v2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
options.register('globalTag','140X_dataRun3_Candidate_2024_06_26_20_41_23',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

## outputFile Name
options.register('outputFileName','ku_twotier_140_gammares_punt.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');

options.register('doTwoTier',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to do twotier processing');
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
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")#
process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")#
#process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
#process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('RecoLocalCalo.EcalRecProducers.ecalCPUUncalibRecHitProducer_cfi')#
#process.load("RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi")
process.load('GammaResTool.GammaResTool.jwk_ku_ecalMultiFitUncalRecHit_cff')
process.load('GammaResTool.GammaResTool.jwk_ku_cc_ecalRecHit_cff')
process.load('GammaResTool.GammaResTool.jwk_ku_ecalRecHit_cff')
process.load('GammaResTool.GammaResTool.jwk_ku_ecalLocalRecoSequence_cff')
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

        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/MINIAOD/PromptReco-v1/000/369/844/00000/1099c5f3-61ca-4f86-bf1b-8215df8b2e7d.root'

	),
    secondaryFileNames = cms.untracked.vstring(

        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/26d660ce-90ee-4dd5-8daa-1ca6afc40e13.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/1d5a7dc7-49b1-4800-bdb7-cab9cdf39687.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/1ca84eb6-c58d-47d7-abc8-6a7cf6b36db2.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/1d8f9009-0764-42a5-b76e-7151e158670b.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/328eeea0-73a7-4833-8c75-e9de909a5a25.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/43e60438-0172-4bcb-a4db-2aba048fea16.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/4a077359-a254-4540-83af-56100ec4e947.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/4aa408ef-5eec-46a2-b25f-ecb366b7bdc3.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/53102e06-1229-4241-928d-bd93ec57dc5e.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/9290b2f8-f608-4d8e-b499-810b70bb27ce.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/a3211580-8c85-4818-bb74-6697255a5458.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/b1d8463a-e508-470e-a3bb-dbad2eeabaa4.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/ec871711-0fca-4a2d-b2f6-24c3819d8a83.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/ee992f90-9537-4e59-bc99-a207cf5d3d47.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/f062b321-c0a3-4554-b9fb-25de47e3722a.root',
        'root://xrootd-cms.infn.it///store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/ff999bba-3944-4e11-8309-19ae3714f6b6.root',

	),
    #eventsToProcess = cms.untracked.VEventRange(eventList)
)
#Successfully opened file root://cmsxrootd.fnal.gov//store/data/Run2023D/EGamma0/RAW/v1/000/369/844/00000/4aa408ef-5eec-46a2-b25f-ecb366b7bdc3.root

## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
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
   triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
   triggerObjects = cms.InputTag("slimmedPatTrigger"),
   ## met filters
   #inputFlags       = cms.string(options.inputFlags),
   triggerFlags     = cms.InputTag("TriggerResults", "", triggerFlagsProcess),
   ecalBadCalibFlag = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),			      
   ## tracks
   tracks = cms.InputTag("unpackedTracksAndVertices"),
   ## vertices
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   ## rho
   rho = cms.InputTag("fixedGridRhoFastjetAll"), #fixedGridRhoAll
   ## METs
   #mets = cms.InputTag("slimmedMETsModifiedMET"),
   mets = cms.InputTag("slimmedMETs"),
   ## jets
   #jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
   jets = cms.InputTag("slimmedJets"),
   ## electrons
   electrons = cms.InputTag("slimmedElectrons"),
   ## muons
   muons = cms.InputTag("slimmedMuons"),
   ## photons
   gedPhotons = cms.InputTag("slimmedPhotons"),
   ootPhotons = cms.InputTag("slimmedOOTPhotons"),
   ## ecal recHits
   recHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
   recHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits"),

   ## do two tier reconstruction of second rechit collection
   doTwoTier = cms.bool(options.doTwoTier), 
   doDiag = cms.bool(options.doDiag),

   ## ecal kuRecHits
   #kuRecHitsEB = cms.InputTag("kuEcalRecHit", "kuRecHitsEB"),
   #kuRecHitsEE = cms.InputTag("kuEcalRecHit", "kuRecHitsEE"),

   #kuStcRecHitsEB = cms.InputTag("kuStcEcalRecHit", "kuStcRecHitsEB"),
   #kuStcRecHitsEE = cms.InputTag("kuStcEcalRecHit", "kuStcRecHitsEE"),
   #kuStcRecHitsEB = cms.InputTag("kuStcEcalLHCRecHit", "kuStcRecHitsEB"),
   #kuStcRecHitsEE = cms.InputTag("kuStcEcalLHCRecHit", "kuStcRecHitsEE"),

   ## ecal Multi rechits

   ##kuWtRecHitsEB = cms.InputTag("kuWtEcalRecHit", "kuWtRecHitsEB"),
   ##kuWtRecHitsEE = cms.InputTag("kuWtEcalRecHit", "kuWtRecHitsEE"),

   #kuWtStcRecHitsEB = cms.InputTag("kuWtStcEcalRecHit", "kuWtStcRecHitsEB"),
   #kuWtStcRecHitsEE = cms.InputTag("kuWtStcEcalRecHit", "kuWtStcRecHitsEE"),

   #kuCCStcRecHitsEB = cms.InputTag("none", "none"),
   #kuCCStcRecHitsEE = cms.InputTag("none", "none"),
   #kuCCStcRecHitsEB = cms.InputTag("kuCCStcEcalRecHit", "kuCCStcRecHitsEB"),
   #kuCCStcRecHitsEE = cms.InputTag("kuCCStcEcalRecHit", "kuCCStcRecHitsEE"),
   kuCCStcRecHitsEB = cms.InputTag("kuCCEcalRecHit", "kuCCRecHitsEB"),
   kuCCStcRecHitsEE = cms.InputTag("kuCCEcalRecHit", "kuCCRecHitsEE"),
   #kuCCStcRecHitsEB = cms.InputTag("kuCCNativeEcalRecHit", "kuCCRecHitsEB"),
   #kuCCStcRecHitsEE = cms.InputTag("kuCCNativeEcalRecHit", "kuCCRecHitsEE"),


   #kuRtStcRecHitsEB = cms.InputTag("kuStcEcalLHCRecHit", "kuStcRecHitsEB"),
   #kuRtStcRecHitsEE = cms.InputTag("kuStcEcalLHCRecHit", "kuStcRecHitsEE"),
   kuRtStcRecHitsEB = cms.InputTag("none", "none"),
   kuRtStcRecHitsEE = cms.InputTag("none", "none"),

   ## ecal uncalib recHits
   #ku_uncalibratedRecHitsEB = cms.InputTag("kuMfootEcalMultiFitUncalibRecHit","kuMfootEcalUncalibRecHitsEB"),
   #ku_uncalibratedRecHitsEE = cms.InputTag("kuMfootEcalMultiFitUncalibRecHit","kuMfootEcalUncalibRecHitsEE"),

   #uncalibratedRecHitsEB = cms.InputTag("ecalMultiFitUncalibRecHitBase","EcalUncalibRecHitsBaseEB"),
   #uncalibratedRecHitsEE = cms.InputTag("ecalMultiFitUncalibRecHitBase","EcalUncalibRecHitsBaseEE"),

   #uncalibratedRecHitsEB = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEB"),
   #uncalibratedRecHitsEE = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEE"),

   ku_uncalibratedRecHitsEB = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEB"),
   ku_uncalibratedRecHitsEE = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEE"),
   #ku_uncalibratedRecHitsEB = cms.InputTag("kuCCNativeEcalMultiFitUncalibRecHit","EcalUncalibRecHitsEB"),
   #ku_uncalibratedRecHitsEE = cms.InputTag("kuCCNativeEcalMultiFitUncalibRecHit","EcalUncalibRecHitsEE"),

   ## digis
   #EBdigiCollection = cms.InputTag("ecalDigis","ebDigis"),
   #EEdigiCollection = cms.InputTag("ecalDigis","eeDigis"),

)


# Set up the path
#process.treePath = cms.Path(
process.tree_step = cms.EndPath(
	process.unpackedTracksAndVertices +
	process.tree
)

process.ecalDigis = process.ecalEBunpacker.clone()
process.ecalDigis.InputLabel = cms.InputTag('rawDataCollector')
process.digiPath = cms.Path( process.ecalDigis )
process.bunchSpacing = cms.Path( process.bunchSpacingProducer )
#process.beamSpot = cms.Path( process.offlineBeamSpot )

process.jwk_calolocalreco = cms.Sequence(
		##process.ku_min_ecalLocalRecoSequence
               	##process.ku_multi_ecalLocalRecoSequence
                #process.kucc_only_ecalLocalRecoSequence
               	#process.ku_reduced_multi_ecalLocalRecoSequence
                #process.ku_spike_nomulti_ecalLocalRecoSequence # vary on reduced
                process.ku_cc_gt_ecalLocalRecoSequence
                #process.ku_cc_native_ecalLocalRecoSequence # for use only with gt that have cc time calibrations
               	##process.ku_ecalLocalRecoSequence
               	##process.ecalLocalRecoSequence
		##process.hcalLocalRecoSequence
)

process.jwk_localreco = cms.Sequence(
		#process.bunchSpacingProducer+
		#process.trackerlocalreco+
		#process.muonlocalreco+
		process.jwk_calolocalreco
		#process.castorreco
)

process.jwk_highlevelreco = cms.Sequence(
		#process.egammaHighLevelRecoPrePF*
                #process.particleFlowReco*
                #process.egammaHighLevelRecoPostPF*
                #process.muoncosmichighlevelreco*
                #process.muonshighlevelreco *
                #process.particleFlowLinks*
                #process.jetHighLevelReco*
                #process.metrecoPlusHCALNoise*
                #process.btagging*
                #process.recoPFMET*
                #process.PFTau*
                #process.reducedRecHits #*
                #process.cosmicDCTracksSeq
)

process.jwk_reconstruction = cms.Sequence(
		#process.localreco*
                process.jwk_localreco
		#process.globalreco*
		#process.jwk_highlevelreco*
		#process.logErrorHarvester
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.content_step = cms.Path(process.content)

#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.ecalraw2digi_step = cms.Path(process.jwk_digisunpacker)
#process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.jwk_reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
		process.bunchSpacing,
                #process.beamSpot,
		process.digiPath,
		#process.raw2digi_step,
		#process.ecalraw2digi_step,
      	        #process.L1Reco_step,
		process.reconstruction_step,
		#process.content_step,
		#process.endjob_step,
		process.tree_step
)

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
