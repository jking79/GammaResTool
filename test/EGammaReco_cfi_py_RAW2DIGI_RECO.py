# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleMuPt1_cfi.py -s RAW2DIGI,RECO --conditions=MC_3XY_V10::All --eventcontent=AOD -n 10 --no_exec
import FWCore.ParameterSet.Config as cms



process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'file:SingleMuPt1_cfi.py_DIGI2RAW.root'),
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
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('EGammaReco_cfi_py_RAW2DIGI_RECO'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('EGammaReco_cfi_py_RAW2DIGI_RECO.root'),
    outputCommands = process.AODEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#140X_dataRun3_Candidate_2024_06_26_20_41_23
process.GlobalTag.globaltag = '140X_dataRun3_Candidate_2024_06_26_20_41_23'
#process.GlobalTag = GlobalTag(process.GlobalTag, 'MC_3XY_V10::All', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODoutput_step = cms.EndPath(process.AODoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.AODoutput_step)
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.AODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
