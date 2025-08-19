#! /usr/bin/env python

import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from httplib import HTTPException

def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = 'submit',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = 'myWorkingArea',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options

#-----------------------------------------------------------------------------------------------------------------------------

def subcrab( runs, events, reqmem ):

    options = getOptions()

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        # External files needed by CRAB
        inputDir     = ''

        #inputJSON     = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
        #inputJSON     = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
        #inputJSON     = 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'

        #inputJSON    = 'Cert_Collisions2022_355100_362760_Golden.json'
        #inputJSON    = 'Cert_Collisions2023_366442_370790_Golden.json'
        inputJSON    = 'Cert_Collisions2024_378981_386951_Golden.json'
        #inputJSON    = ''

        #inputJSON     = 'Collisions25_13p6TeV_Latest.json'

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        import CRABClient
        from CRABClient.UserUtilities import config
        config = config()

        config.General.workArea    = options.workArea
        #config.General.requestName = None

        config.JobType.pluginName  = 'Analysis'
        config.JobType.psetName    = 'gammares_notwotier_multi.py'
        #config.JobType.numCores    = 8
        #config.JobType.maxMemoryMB = 2250 #reqmem
        #config.JobType.maxJobRuntimeMin = 1500
        config.JobType.pyCfgParams = None
        #config.JobType.inputFiles  = [ inputDir+inputPaths , inputDir+inputFilters , inputDir+inputFlags ]

        config.section_("Site")
        config.Site.whitelist = ['T1_*','T2_US_*','T2_IT_*','T2_DE_*','T2_ES_*','T2_FR_*','T2_UK_*']
        config.Site.blacklist = ['T2_ES_IFCA']

        config.Data.inputDataset = None
        config.Data.partialDataset = True
        #config.Data.useParent      = True
	#config.Data.secondaryInputDataset = secInputPaths
        #config.Data.useParent      = False
        config.Data.lumiMask     = inputDir+inputJSON
        #config.Data.splitting    = 'LumiBased'
        config.Data.splitting    = 'EventAwareLumiBased'
        #config.Data.splitting    = 'Automatic'

        #config.Data.unitsPerJob  =  25000 # for auto job splitting
        #config.Data.runRange  = runs #'321122-321128'
        config.Data.unitsPerJob  =  400000 # unitsPerJob = 1000 for 321122-321128 and maxMemoryMB = 4000  on EventAwareLumiBased
        #config.Data.unitsPerJob  =  100000 # 2024F cc test rereco

        #config.Data.outputDatasetTag = 'reducedRAW_EGamma_ntuple'
	     
        config.JobType.allowUndistributedCMSSW = True
        config.Data.publication      = False
        config.Site.storageSite      = 'T3_US_FNALLPC'
        #config.Data.outLFNDirBase    = '/store/user/jaking/ecalTiming/'
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDataAndOpts = [

            # Dataset: Run2025

            #['/EGamma0/Run2025B-PromptReco-v1/MINIAOD'],
            #['/EGamma1/Run2025B-PromptReco-v1/MINIAOD'],
            #['/EGamma2/Run2025B-PromptReco-v1/MINIAOD'],
            #['/EGamma3/Run2025B-PromptReco-v1/MINIAOD'],

            #['/EGamma0/Run2025C-PromptReco-v1/MINIAOD'],
            #['/EGamma1/Run2025C-PromptReco-v1/MINIAOD'],
            #['/EGamma2/Run2025C-PromptReco-v1/MINIAOD'],
            #['/EGamma3/Run2025C-PromptReco-v1/MINIAOD'],

            #['/EGamma0/Run2025C-PromptReco-v2/MINIAOD'],
            #['/EGamma1/Run2025C-PromptReco-v2/MINIAOD'],
            #['/EGamma2/Run2025C-PromptReco-v2/MINIAOD'],
            #['/EGamma3/Run2025C-PromptReco-v2/MINIAOD'],

            #['/EGamma0/Run2025D-PromptReco-v1/MINIAOD'],
            #['/EGamma1/Run2025D-PromptReco-v1/MINIAOD'],
            #['/EGamma2/Run2025D-PromptReco-v1/MINIAOD'],
            #['/EGamma3/Run2025D-PromptReco-v1/MINIAOD'],

            # Dataset: Run2024

            ['/EGamma1/Run2024F-ECAL_CC_HCAL_DI-v3/MINIAOD'],

            ##['/EGamma0/Run2024A-PromptReco-v1/MINIAOD'], 	# 378927-378962 52.6M
            #['/EGamma0/Run2024B-PromptReco-v1/MINIAOD'],	# 378981-379350	0.54T
            #['/EGamma0/Run2024C-PromptReco-v1/MINIAOD'],	# 379413-379765	3.3T
            #['/EGamma0/Run2024D-PromptReco-v1/MINIAOD'],       # 380306-380933 8.7T 
            #['/EGamma0/Run2024E-PromptReco-v2/MINIAOD'], 
            ['/EGamma0/Run2024F-PromptReco-v1/MINIAOD'], 
            #['/EGamma0/Run2024G-PromptReco-v1/MINIAOD'], 
            #['/EGamma0/Run2024H-PromptReco-v1/MINIAOD'], 
            #['/EGamma0/Run2024I-PromptReco-v2/MINIAOD'], 
            #['/EGamma0/Run2024J-PromptReco-v1/MINIAOD'], 

            ##['/EGamma1/Run2024A-PromptReco-v1/MINIAOD'],	# 378919-378961	52.1M
	    #['/EGamma1/Run2024B-PromptReco-v1/MINIAOD'],	# 378981-379349	0.54T
            #['/EGamma1/Run2024C-PromptReco-v1/MINIAOD'],	# 379415-379774	3.4T
            #['/EGamma1/Run2024D-PromptReco-v1/MINIAOD'],       # 380306-380933 8.7T
            #['/EGamma1/Run2024E-PromptReco-v2/MINIAOD'], 
            ['/EGamma1/Run2024F-PromptReco-v1/MINIAOD'],
            #['/EGamma1/Run2024G-PromptReco-v1/MINIAOD'], 
            #['/EGamma1/Run2024H-PromptReco-v1/MINIAOD'], 
            #['/EGamma1/Run2024I-PromptReco-v2/MINIAOD'], 
            #['/EGamma1/Run2024J-PromptReco-v1/MINIAOD'], 

            # Dataset: Run2023


            # Dataset: Run2024


            # Dataset: /EGamma/Run2018-12Nov2019_UL2018-/MINIAOD

            #['/EGamma/Run2018A-15Feb2022_UL2018-v1/MINIAOD'],
            #['/EGamma/Run2018B-15Feb2022_UL2018-v1/MINIAOD'],
            #['/EGamma/Run2018C-15Feb2022_UL2018-v1/MINIAOD'],
            #['/EGamma/Run2018D-15Feb2022_UL2018-v1/MINIAOD'],

            # Dataset: /EGamma/Run2018-12Nov2019_UL2018-/AOD

            #['/EGamma/Run2018A-12Nov2019_UL2018-v2/MINIAOD'],
            #['/EGamma/Run2018B-12Nov2019_UL2018-v2/MINIAOD'],
            #['/EGamma/Run2018C-12Nov2019_UL2018-v2/MINIAOD'],
            #['/EGamma/Run2018D-12Nov2019_UL2018-v8/MINIAOD'],

            # Dataset: /DoubleEG/Run2016-21Feb2020_UL2016-/MINIAOD

            #['/DoubleEG/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'],
            #['/DoubleEG/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD'],
            #['/DoubleEG/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD'],
            #['/DoubleEG/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD'],
            #['/DoubleEG/Run2016F-21Feb2020_UL2016-v1/MINIAOD'],
            #['/DoubleEG/Run2016G-21Feb2020_UL2016-v1/MINIAOD'],
            #['/DoubleEG/Run2016H-21Feb2020_UL2016-v1/MINIAOD'],

            # Dataset: /DoubleEG/Run2017-09Aug2019_UL2017-/MINIAOD

            #['/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD'],
            #['/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD'],
            #['/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD'],
            #['/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD'],
            #['/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD'],

            #['/DoubleEG/Run2017B-09Aug2019_UL2017-v1/MINIAOD'],
            #['/DoubleEG/Run2017C-09Aug2019_UL2017-v1/MINIAOD'],
            #['/DoubleEG/Run2017D-09Aug2019_UL2017-v1/MINIAOD'],
            #['/DoubleEG/Run2017E-09Aug2019_UL2017-v1/MINIAOD'],
            #['/DoubleEG/Run2017F-09Aug2019_UL2017-v1/MINIAOD'],

	    ]
 
        for inDO in inputDataAndOpts:
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
	         #print( inDO[0] )
            primaryDataset = inDO[0].split('/')[1]
            #primaryDataset = (inDO[0].split('/')[1]).split('_Tune')[0]
            runEra         = inDO[0].split('/')[2]
            #runEra         = inDO[0].split('/')[2].split('-')[0]
            #runEra         = ((inDO[0].split('/')[2]).split('_')[0]+'_'+(inDO[0].split('/')[2]).split('_')[1]).split('-PU')[0] # MC
            #runEra         = (inDO[0].split('/')[2]).split('-133X')[0]
            dataset	   = inDO[0].split('/')[3]

            config.Data.allowNonValidInputDataset = True # speacial rereco pd

            #trial          = 'gammares_ratio_126_v2'
            #trial          = 'gammares_cc_140_v2'
            #trial          = 'gammares_ttcc_140_v5' # 24C and earlier only
            #trial          = 'gammares_mc'
            #trial          = 'gammares_llpana'
            #trial          = 'gammares_r24f_cctest'
            #trial          = 'gammares_cali'
            #trial          = 'gammares_24mc'
            #trial          = 'gammares_r24f_prompt'
            #trial          = 'gammares_r25_val'
            #trial          = 'gammares_ECAL_CC_HCAL_DI-v3'
            trial          = 'gammares_DPG_24'
            #trial          = 'gammares_ul18'

            #config.Data.outLFNDirBase    = "/store/user/jaking/ecalTiming/"+trial+"/"
            config.Data.outLFNDirBase    = "/store/group/lpcsusylep/jaking/ecalTiming/"+trial+"/"
            #config.General.requestName   = trial+"_"+primaryDataset+"_"+runEra+"_"+runs+"_"+dataset
            #config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_"+runs+"_dispho"
            #config.Data.outputDatasetTag = runEra+"_"+runs+"_"+dataset
            config.General.requestName   = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_request"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra



            ofn = 'outputFileName=output.root'
            dTTf = 'doTwoTier=False'
            dTTt = 'doTwoTier=True'
            dDt = 'doDiag=True'
            dDf = 'doDiag=False'
#>>>>>>>>>>>>>  Run2 UL 16/17/18
            #gtag = 'globalTag=106X_dataRun2_v36'
            #gtag = 'globalTag=106X_dataRun2_v20'
            #config.JobType.pyCfgParams   = [gtag,ofn,dTTf,dDt]
            #config.JobType.pyCfgParams   = [gtag, ofn,dTTf,dDt]
#>>>>>>>>>>>>>  2017 EOY 94X_dataRun2_ReReco_EOY17_v1
            #gtag = 'globalTag=94X_dataRun2_ReReco_EOY17_v1'
            #config.JobType.pyCfgParams   = [gtag, ofn,dTTf,dDt]
#>>>>>>>>>>>>>  2018 EOY
            #gtag = 'globalTag= 102X_dataRun2_Prompt_v11'
            #config.JobType.pyCfgParams   = [gtag, ofn,dTTf,dDt]
            ### MC 2017
            #gtag = 'globalTag=94X_mc2017_realistic_v12'
            #config.JobType.pyCfgParams = [gtag, ofn,dTTf,dDt]

#>>>>>>>>>>>> 2024 tested
            #gtag = 'globalTag=140X_dataRun3_Prompt_v3'
            #config.JobType.pyCfgParams   = [gtag, ofn,dTTf,dDt]
            #gtag = 'globalTag=141X_dataRun3_HLT_frozen_v2'
            #config.JobType.pyCfgParams   = [gtag, ofn,dTTf,dDt]

#>>>>>>>>>>>> 2025
            gtag = 'globalTag=150X_dataRun3_Prompt_v1'
            config.JobType.pyCfgParams   = [gtag, ofn,dTTf,dDt]

            ## MC Run3Winter24
            #config.JobType.pyCfgParams   = ['globalTag=133X_mcRun3_2024_realistic_v8', ofn,dTTf,dDt]
            #config.JobType.pyCfgParams   = ['globalTag=140X_dataRun3_Prompt_v2', ofn,dTTf,dDt]  # 2024 tested
            #gtag = 'globalTag=140X_dataRun3_Prompt_Candidate_2024_05_31_21_23_47'
            #config.JobType.pyCfgParams   = [gtag, ofn,dTTf,dDt]  # 5 GeV cali - EE cali bad
            ##config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v4',dTTf,dDf]  
            #config.JobType.pyCfgParams   = ['globalTag=112X_dataRun3_Prompt_v2', ofn,dTTf,dDf]]  # 2018A tested
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_PromptAnalysis_v1', ofn,dTTf,dDf] # ABCD
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v8', ofn,dTTf,dDf] # EF
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v10', ofn,dTTf,dDf] # G
            #config.JobType.pyCfgParams   = ['globalTag=130X_mcRun3_2023_realistic_v14', ofn,dTTf,dDt] # R3 DY summer23
            #config.JobType.pyCfgParams   = ['globalTag=133X_mcRun3_2024_realistic_v10', ofn,dTTf,dDt] # R3 DY winter24
            #config.JobType.pyCfgParams   = ['globalTag=133X_mcRun3_2024_realistic_v8', ofn,dTTf,dDt] #
            #config.JobType.pyCfgParams   = ['globalTag=102X_dataRun2_Prompt_v11', ofn,dTTf,dDt] #


            config.Data.inputDataset     = inDO[0]
            # Submit.
            try:
                print( "Submitting for input dataset %s for runs %s" % (inDO[0], runs) )
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
                os.system("rm -rf %s/crab_%s/inputs" % (config.General.workArea, config.General.requestName))
            #except HTTPException as hte:
            #    print( "Submission for input dataset %s failed: %s" % (inDO[0], hte.headers) )
            except ClientException as cle:
                print( "Submission for input dataset %s failed: %s" % (inDO[0], cle) )

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print( "-"*len(msg) )
            print( msg )
            print( "-"*len(msg) )
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers) )
            except ClientException as cle:
                print( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle) )

#--------------------------------------------------------------------------------------------------------------------------

def submit_run():

	#subcrab( "378981-379391","",2500)#	Run3 24B		
        #subcrab( "379415-380238","",2500)#	Run3 24C
        #subcrab( "380306-380947","",2500)#      Run3 24D

	#subcrab( "352400-358400","",2500)# Run3 22ABCD
        #subcrab( "359000-362200","",2500)# Run3 22EF
        #subcrab( "362300-362800","",2500)# Run3 22G

        #subcrab( "316241-316245","",2500) # very small test batch 18A
        #subcrab( "316000-316499","",2500) # 18A	
        #subcrab( "357101-357268","",2500)#22C?
	#subcrab( "360395-360415","",2500)#22F?

        #subcrab( "378981-386951","",2500)#R3_2024
        subcrab( "000000-999999","",2500)#MC - any

########################################################

submit_run()
