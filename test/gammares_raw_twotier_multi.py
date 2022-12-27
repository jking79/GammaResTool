import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## GT to be used
#options.register('globalTag','124X_dataRun3_Prompt_v4',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','112X_dataRun3_Prompt_v2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
options.register('globalTag','124X_dataRun3_PromptAnalysis_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

## outputFile Name
options.register('outputFileName','ku_KUCC_only_126_gammares_v2.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');

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
process.load('GammaResTool.gammaResTool.jwk_ku_ecalLocalRecoSequence_cff')
#process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff' )

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')#

process.load('Configuration.StandardSequences.EndOfProcess_cff')#

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")#
#process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## Define the input source
#eventList = open(options.rlelist,'r')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'file:jwk_reco_data_DIGI2RAW.root'),
		#'/store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/355/809/00000/60ef0541-b8f5-479b-becd-4fdbd0e0599b.root'
    	'/store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/355/892/00000/08458ea6-e15b-46ff-a71d-6f61d4c6d288.root'
		#'/store/data/Run2022F/EGamma/MINIAOD/PromptReco-v1/000/362/154/00000/0392b31a-851b-4c61-98f5-cfb1598aef5f.root'
        #'/store/data/Run2018A/EGamma/MINIAOD/17Sep2018-v2/100000/10E13819-4C42-FA4F-B6A1-53A26A2F388A.root'
	),
    secondaryFileNames = cms.untracked.vstring(

#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/26201159-BA64-E811-B349-02163E019F92.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/28389948-BA64-E811-B1DE-FA163E0639A2.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/36DF98ED-BA64-E811-8926-02163E01A036.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/3ED17F49-BA64-E811-B3AA-FA163E649742.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/4E09C5E3-BA64-E811-B4A6-FA163EE8E7AA.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/8A6280FC-BA64-E811-A47B-02163E01A029.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/A8557BB9-BA64-E811-9F57-FA163E1D3FC8.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/F03848F8-BA64-E811-A717-FA163E4907DA.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00001/8434C310-BB64-E811-BEE4-FA163EF0B127.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/944/00000/42332B35-9064-E811-8D6C-FA163EBE4E78.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/985/00000/123ED0DD-B764-E811-B2CD-02163E01A0B0.root',
#        '/store/data/Run2018A/EGamma/RAW/v1/000/316/985/00000/24D8862F-B864-E811-8087-FA163EC5FAA0.root',


        '/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/f966971d-944d-4b72-8b31-f23914a64695.root',
        '/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/f9fd7d6f-d376-4dd8-a456-065dfa315842.root',
        '/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/faa18a05-0c6b-4319-b684-03a7b4953cc8.root',
        '/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/fb874daf-4226-4e33-8041-d0e657517b87.root',
        '/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/ffeab3bd-21aa-491f-8813-5e40e430dcd8.root'
		
#		'/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01df3a20-4381-4103-8b98-5fd77e266823.root',

#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0010965c-a872-454b-9b4f-58ac78517be6.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/003464aa-8285-40a1-bdfc-623b8c3c9c11.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/00f12706-41f5-4de7-b6ec-3d5642ce1079.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01367a50-7548-46f7-a1ab-21ee2b461c09.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/013c0377-784e-4fec-a671-2c0001950dda.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01754c90-e053-4e3a-b46d-9297dd14d1ac.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01df0966-fd10-4513-a227-ba0116eaf323.root',
#        #'/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01df3a20-4381-4103-8b98-5fd77e266823.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/02e2a278-bb19-4bfa-82c9-ec160adce841.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/035b02fb-dbe9-4889-a46d-0dc40acc49a9.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/03b7d404-2afd-459a-8c4c-88be935fb000.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/03f4c919-314b-49d4-a979-13596d55997a.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0416a8d1-cc81-41d3-aa6a-e3850a174fbf.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/04f1bea2-4572-4089-a1cf-f0b98ebb58fc.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/055fa8c0-2f83-4ecf-b3e9-bbbdf8a6b701.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/05679cb3-468f-4ed3-a795-9c94c9027b43.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/06013a2c-d014-4b40-b6bf-391617ce1897.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/07025fe3-45e6-49de-8bf9-2011f65cef75.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/076f7b4c-449f-4d7d-9f64-34a9cb073793.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/081d0f90-4314-4b2e-abf9-21b8eef00af4.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/08b1c8b6-b8eb-4ee5-99df-6a9bb3d35a8c.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/08fc9701-09a4-4423-bf54-7573634e2f8b.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/09614c60-7af4-49b3-b16f-8373bd8e165c.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0a3b1163-529d-4887-ba78-0c7677ba724a.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0a430a98-5f17-44e7-bfc3-b67e55b5fbc2.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0a685ecd-5220-43e4-b1ad-3d9d2ee3208a.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0b62d850-f979-4aa3-bb96-7bdf98e35529.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0bfd60c7-fa0e-454d-bfb9-4e7508d50eb6.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0bfdec2a-098d-461d-8570-f03681ebf7c1.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0c3188cd-249c-48e0-b5c7-cbf72c852e74.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0cd4d5ca-0244-474e-91e2-dd7137ffb2c7.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0d347e7b-bc3f-4e3e-8328-580874a53c32.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0d798a37-8423-427e-9049-07b6911fdbe5.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0d9b834e-a110-42ff-afba-0ef75ff2d58b.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0e9d6c19-31cd-4f7d-abb3-f94db32c6234.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0f0ceb8a-eeaf-4fba-a88a-11d3fcf69ac8.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0f25b562-9ff3-4e92-926a-62ae92fef387.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0f28ad51-b4a2-4a63-a72c-e4e28756c43a.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0f7c6b44-aaf6-461c-b690-85dea88ca036.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0facce5e-4b5a-45f3-b72c-5b96be3ede0b.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0fba0679-fac9-48b6-94cb-fc206713ad21.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0fc2731d-80d2-4ba1-8154-3dcad3b615d9.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0fc889c2-c932-4e08-97b7-2bba83f0484a.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/10a2f010-1f2b-4081-8c45-05538991b453.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/10f5c8e7-41dc-40fd-b208-d8bd514b2db2.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/1119be7c-35cc-4f5a-bb33-043b2b121f03.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/11368d45-2c64-46b2-aae5-cc35b1e18d42.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/116043a5-e71e-4df8-8014-34fdb3391301.root'

	),
    #eventsToProcess = cms.untracked.VEventRange(eventList)
)


## How many events to process
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

   kuCCStcRecHitsEB = cms.InputTag("kuCCStcEcalRecHit", "kuCCStcRecHitsEB"),
   kuCCStcRecHitsEE = cms.InputTag("kuCCStcEcalRecHit", "kuCCStcRecHitsEE"),
   #kuCCStcRecHitsEB = cms.InputTag("kuCCStcEcalLHCRecHit", "kuCCStcRecHitsEB"),
   #kuCCStcRecHitsEE = cms.InputTag("kuCCStcEcalLHCRecHit", "kuCCStcRecHitsEE"),

   ## ecal uncalib recHits
   #ku_uncalibratedRecHitsEB = cms.InputTag("kuMfootEcalMultiFitUncalibRecHit","kuMfootEcalUncalibRecHitsEB"),
   #ku_uncalibratedRecHitsEE = cms.InputTag("kuMfootEcalMultiFitUncalibRecHit","kuMfootEcalUncalibRecHitsEE"),

   #uncalibratedRecHitsEB = cms.InputTag("ecalMultiFitUncalibRecHitBase","EcalUncalibRecHitsBaseEB"),
   #uncalibratedRecHitsEE = cms.InputTag("ecalMultiFitUncalibRecHitBase","EcalUncalibRecHitsBaseEE"),

   #uncalibratedRecHitsEB = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEB"),
   #uncalibratedRecHitsEE = cms.InputTag("kuEcalMultiFitUncalibRecHit","kuEcalUncalibRecHitsEE"),

   ku_uncalibratedRecHitsEB = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEB"),
   ku_uncalibratedRecHitsEE = cms.InputTag("kuCCEcalMultiFitUncalibRecHit","kuCCEcalUncalibRecHitsEE"),

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
                process.kucc_only_ecalLocalRecoSequence
               	#process.ku_reduced_multi_ecalLocalRecoSequence
                #process.ku_spike_multi_ecalLocalRecoSequence # vary on reduced
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
		process.endjob_step,
		process.tree_step
)

### Extra bits from other configs
process.options = cms.untracked.PSet(
    #numberOfThreads = cms.untracked.uint32(4),
    #numberOfStreams = cms.untracked.uint32(4),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    #wantSummary = cms.untracked.bool(True)
)

#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#if options.runUnscheduled : 
#	process = convertToUnscheduled(process)

#from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
#if options.deleteEarly :
#	process = customiseEarlyDelete(process)
