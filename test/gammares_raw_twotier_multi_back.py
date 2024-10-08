import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## GT to be used
##options.register('globalTag','124X_dataRun3_Prompt_v4',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabltag to be used');
##options.register('globalTag','112X_dataRun3_Prompt_v2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabltag to be used');
#options.register('globalTag','130X_dataRun3_Candidate_2023_08_08_21_30_44',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','130X_dataRun3_Prompt_v2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
options.register('globalTag','132X_dataRun3_Express_v4',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

## outputFile Name
options.register('outputFileName','ku_tt_131diag_126_gammares_v13.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');
#options.register('outputFileName','ku_KUCC_tt_R2018A_126_gammares_v11.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');

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
process.load('GammaResTool.GammaResTool.jwk_ku_ecalLocalRecoSequence_cff')
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

        #'/store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/355/809/00000/60ef0541-b8f5-479b-becd-4fdbd0e0599b.root',

        '/store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/355/828/00000/4c479772-def0-43f1-b615-9a9fc965d36d.root',

        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/081/00000/1544045f-d9c3-47a4-b31c-e08bc488f0dc.root',
        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/081/00000/90a1f768-e33b-46df-adc7-c2e6318b36d5.root',
        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/081/00000/e456709e-b59f-41d3-8d92-6b88c69ffc18.root',
        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/082/00000/9948e5c1-f309-4fe1-bdab-9c61b83e4ac5.root',
        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/082/00000/ddee20fb-56fe-45a9-a7c8-aa38850e2592.root',
        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/086/00000/d805d8b5-173a-41d5-beae-1b15532e6c85.root',
        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/087/00000/cc2657bd-87a2-43c5-b904-725f17216d25.root',
        #'/store/data/Run2023F/ZeroBias/RAW/v1/000/373/088/00000/03e9d5a0-2d39-4090-aba5-e378e256c09a.root',

     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/04f24116-3e5e-4d66-a029-fd2c2fcfdaca.root',
     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/0c9575ba-8a64-4537-8dee-83af6ed02b59.root',
     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/1c78c82f-99e9-4fd8-8551-74591e745952.root',
     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/2dc15752-63dd-49fb-8f29-f278b13dd08a.root',
     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/3885916c-9890-40e8-bb6c-70863d9b5d6e.root',
     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/3cf46990-85b4-40b1-950b-338179cf004f.root',
     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/3d3fa3d5-8419-41b7-9218-b8bc177e7743.root',
     #   '/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/41a3b249-7bab-4cf9-8138-03ba4a490ee1.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/534a9b86-9915-4962-8f89-38abd1d427e2.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/5350f914-ce7c-4e48-8d19-e7ffbed52ebd.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/58a05630-45b8-40d4-b237-aa5750cda69a.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/5b75b186-c7aa-4de3-ba11-f83d5cd03557.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/5e6cca5b-b3c6-4cae-a00a-d770e8f8bf30.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/73224f95-f4c7-4b75-b9b7-1c861f2eac08.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/751087b8-d58e-4cd2-9d38-2dff67d09f56.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/754eb530-3522-443c-8371-1157d31406e1.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/7642825e-868c-498c-96bf-1fd7327f2986.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/7cb8eeb9-50a7-4f2c-8b33-63392802bf5d.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/8935c147-5a6f-45cb-aec6-0f1c7ac851e6.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/8c90f76b-af6d-4c0c-874f-efb78acc912a.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/8eff40dc-371f-4318-ac46-00e6b43e55c0.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/91962d78-7162-4378-9975-fbd57b0c048c.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/92811b2b-365e-4700-8c21-2f4b3fcd93e4.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/959235f3-7230-4822-892a-3a504f327fb7.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/98d1b86e-a754-4715-adcc-34b3e0ffbfcc.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/9ca091dc-e01f-4ed3-b2f8-794d11e192b1.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/aee44530-0734-491d-a97b-b8de9d25bd86.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/c5e8af1e-2d5e-47d8-8623-9fb825dd94da.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/cf5361c2-b1e4-45df-aa73-46f3bb896780.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/da300b9a-860a-4ccb-a76b-60cd627ae706.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/dbb37251-0b95-4d73-a262-928ae2fc0d39.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/de119dc0-c657-4520-b1ed-4768725908da.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/e4109b21-9846-408d-820f-5d1a56ba97a4.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/e4eb9fb4-0fc8-43e3-a029-42647cca246a.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/ef799316-8b8f-46bb-853e-379f0a0330de.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/f1632c55-5905-4251-83d1-4625503463f6.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/f2b69df0-29a9-4e1e-8ea3-f86a6b961680.root',
        #'/store/data/Run2023B/ZeroBias/RAW/v1/000/366/850/00000/f721444f-b367-4dfc-b0fc-8de7416ac50a.root',

        #'/store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/355/809/00000/60ef0541-b8f5-479b-becd-4fdbd0e0599b.root'
##i    	'/store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/355/892/00000/08458ea6-e15b-46ff-a71d-6f61d4c6d288.root'
	#'/store/data/Run2022F/EGamma/MINIAOD/PromptReco-v1/000/362/154/00000/0392b31a-851b-4c61-98f5-cfb1598aef5f.root'
        #'/store/data/Run2018A/EGamma/MINIAOD/17Sep2018-v2/100000/10E13819-4C42-FA4F-B6A1-53A26A2F388A.root'
	),
    secondaryFileNames = cms.untracked.vstring(

         #'/store/data/Run2022C/EGamma/RAW/v1/000/355/809/00000/214d5b97-90c6-44ae-85b7-cab7a364d56e.root',
         #'/store/data/Run2022C/EGamma/RAW/v1/000/355/828/00000/144d4d6e-ff92-431a-b488-ab82c3ca1109.root',
         #'/store/data/Run2022C/EGamma/RAW/v1/000/355/828/00000/8df27f43-9b13-4c81-b7ea-be23d9ab67d1.root',

        #'/store/data/Run2018A/EGamma/RAW/v1/000/316/758/00000/26201159-BA64-E811-B349-02163E019F92.root',
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


        #'/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/f966971d-944d-4b72-8b31-f23914a64695.root',
        #'/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/f9fd7d6f-d376-4dd8-a456-065dfa315842.root',
        #'/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/faa18a05-0c6b-4319-b684-03a7b4953cc8.root',
        #'/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/fb874daf-4226-4e33-8041-d0e657517b87.root',
##        '/store/data/Run2022C/EGamma/RAW/v1/000/355/892/00000/ffeab3bd-21aa-491f-8813-5e40e430dcd8.root'
		
#	'/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01df3a20-4381-4103-8b98-5fd77e266823.root',

#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/0010965c-a872-454b-9b4f-58ac78517be6.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/003464aa-8285-40a1-bdfc-623b8c3c9c11.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/00f12706-41f5-4de7-b6ec-3d5642ce1079.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01367a50-7548-46f7-a1ab-21ee2b461c09.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/013c0377-784e-4fec-a671-2c0001950dda.root',
#        '/store/data/Run2022F/EGamma/RAW/v1/000/362/154/00000/01754c90-e053-4e3a-b46d-9297dd14d1ac.root',

	),
    #eventsToProcess = cms.untracked.VEventRange(eventList)
)


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(47000))

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
 
   #recHitsEB = cms.InputTag("kuStcEcalRecHit", "kuStcRecHitsEB"),
   #recHitsEE = cms.InputTag("kuStcEcalRecHit", "kuStcRecHitsEE"),

   ## do two tier reconstruction of second rechit collection
   #doTwoTier = cms.bool(True), 
   doTwoTier = cms.bool(False),
   #doDiag = cms.bool(True),
   doDiag = cms.bool(False),

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

   #kuCCStcRecHitsEB = cms.InputTag("kuCCStcEcalRecHit", "kuCCStcRecHitsEB"),
   #kuCCStcRecHitsEE = cms.InputTag("kuCCStcEcalRecHit", "kuCCStcRecHitsEE"),
   kuCCStcRecHitsEB = cms.InputTag("kuCCEcalRecHit", "kuCCRecHitsEB"),
   kuCCStcRecHitsEE = cms.InputTag("kuCCEcalRecHit", "kuCCRecHitsEE"),
   #kuCCStcRecHitsEB = cms.InputTag("kuCCStcEcalLHCRecHit", "kuCCStcRecHitsEB"),
   #kuCCStcRecHitsEE = cms.InputTag("kuCCStcEcalLHCRecHit", "kuCCStcRecHitsEE"),

   kuRtStcRecHitsEB = cms.InputTag("kuStcEcalLHCRecHit", "kuStcRecHitsEB"),
   kuRtStcRecHitsEE = cms.InputTag("kuStcEcalLHCRecHit", "kuStcRecHitsEE"),

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
process.tree_step = cms.EndPath(
	process.unpackedTracksAndVertices +
	process.tree
)

process.ecalDigis = process.ecalEBunpacker.clone()
process.ecalDigis.InputLabel = cms.InputTag('rawDataCollector')
process.digiPath = cms.Path( process.ecalDigis )
process.bunchSpacing = cms.Path( process.bunchSpacingProducer )
#process.beamSpot = cms.Path( process.offlineBeamSpot )

process.jwk_localreco = cms.Sequence(
                #process.kucc_only_ecalLocalRecoSequence
                process.ku_cc_gt_ecalLocalRecoSequence 
                #process.ku_reduced_multi_ecalLocalRecoSequence
                #process.ku_spike_multi_ecalLocalRecoSequence # vary on reduced
                ##process.ku_reduced_flipped_ecalLocalRecoSequence
)

#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.ecalraw2digi_step = cms.Path(process.jwk_digisunpacker)
#process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.jwk_localreco)
#process.reconstruction_step = cms.Path( process.ecalLocalRecoSequence )
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
		process.bunchSpacing,
                #process.beamSpot,
		process.digiPath,
		#process.raw2digi_step,
		#process.ecalraw2digi_step,
      	        #process.L1Reco_step,
		process.reconstruction_step,
		process.endjob_step,
		#process.tree_step
)

### Extra bits from other configs
process.options = cms.untracked.PSet(
    #numberOfThreads = cms.untracked.uint32(4),
    #numberOfStreams = cms.untracked.uint32(4),
    ####SkipEvent = cms.untracked.vstring('ProductNotFound'),
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    #wantSummary = cms.untracked.bool(True)
)

