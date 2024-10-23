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
        #inputJSON    = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
        #inputJSON    = 'Cert_Collisions2024_eraC_Golden.json'
        inputJSON    = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
        #inputJSON    = ''
        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        import CRABClient
        from CRABClient.UserUtilities import config
        config = config()

        config.General.workArea    = options.workArea
        #config.General.requestName = None

        config.JobType.pluginName  = 'Analysis'
        config.JobType.psetName    = 'gammares_notwotier_multi_aod.py'
        #config.JobType.numCores    = 8
        #config.JobType.maxMemoryMB = 2250 #reqmem
        #config.JobType.maxJobRuntimeMin = 1500
        config.JobType.pyCfgParams = None
        #config.JobType.inputFiles  = [ inputDir+inputPaths , inputDir+inputFilters , inputDir+inputFlags ]

        config.Data.inputDataset = None
        #config.Data.useParent      = True
	#config.Data.secondaryInputDataset = secInputPaths
        config.Data.partialDataset = True
        #config.Data.useParent      = False
        #config.Data.lumiMask     = inputDir+inputJSON
        #config.Data.splitting    = 'LumiBased'
        config.Data.splitting    = 'EventAwareLumiBased'
        #config.Data.splitting    = 'Automatic'
        config.Data.unitsPerJob  =  25000 # MET 2017 AOD
        #config.Data.unitsPerJob  =  1500 # MC GMSB
        #config.Data.unitsPerJob  =  15000 # MC QCD

        #config.Data.unitsPerJob  =  300 # for auto job splitting
        #config.Data.runRange  = runs #'321122-321128'
        #config.Data.unitsPerJob  =  250000 # unitsPerJob = 1000 for 321122-321128 and maxMemoryMB = 4000  on EventAwareLumiBased


        #config.Data.outputDatasetTag = 'reducedRAW_EGamma_ntuple'
	     
        config.JobType.allowUndistributedCMSSW = True
        config.Data.publication      = False
        config.Site.storageSite      = 'T3_US_FNALLPC'
        #config.Data.outLFNDirBase    = '/store/user/jaking/ecalTiming/'
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDataAndOpts = [

            ##['/EGamma0/Run2024A-PromptReco-v1/MINIAOD'], 	# 378927-378962 52.6M
            #['/EGamma0/Run2024B-PromptReco-v1/MINIAOD'],	# 378981-379350	0.54T
            #['/EGamma0/Run2024C-PromptReco-v1/MINIAOD'],	# 379413-379765	3.3T
            #['/EGamma0/Run2024D-PromptReco-v1/MINIAOD'],       # 380306-380933 8.7T 

            ##['/EGamma1/Run2024A-PromptReco-v1/MINIAOD'],	# 378919-378961	52.1M
	    #['/EGamma1/Run2024B-PromptReco-v1/MINIAOD'],	# 378981-379349	0.54T
            #['/EGamma1/Run2024C-PromptReco-v1/MINIAOD'],	# 379415-379774	3.4T
            #['/EGamma1/Run2024D-PromptReco-v1/MINIAOD'],       # 380306-380933 8.7T

            #['/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
            #['/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'],
            ##['/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM'],
            #['/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
            #['/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
            #['/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'],
            #['/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
            #['/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_new_pmx_94X_mc2017_realistic_v11-v1/AODSIM'],
            #['/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'],
            #['/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'],

            #['/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],

            #['/EGamma/Run2018A-17Sep2018-v2/MINIAOD'], # for ABC change to GT 102X_dataRun2_v13

            #['/EGamma/Run2022A-PromptReco-v1/MINIAOD'], # 2022 ABCD Prompt 124X_dataRun3_PromptAnalysis_v1 352499 - 355062    
            #['/EGamma/Run2022B-PromptReco-v1/MINIAOD'], #						    355094 - 355769											 
            #['/EGamma/Run2022C-PromptReco-v1/MINIAOD'], #						    355809 - 357482
            #['/EGamma/Run2022D-PromptReco-v1/MINIAOD'], #						    357538 - 357733
            #['/EGamma/Run2022D-PromptReco-v2/MINIAOD'], #						    357734 - 357902
            #['/EGamma/Run2022D-PromptReco-v3/MINIAOD'], #						    358381
            #['/EGamma/Run2022E-PromptReco-v1/MINIAOD'], # 2022 EF Prompt 124X_dataRun3_Prompt_v8	    359090 - 360327
            #['/EGamma/Run2022F-PromptReco-v1/MINIAOD'], #						    360389 - 362167
            #['/EGamma/Run2022G-PromptReco-v1/MINIAOD']  # ? 2022 G Prompt 124X_dataRun3_Prompt_v10   	    362399 - 362760

#
#            ['/GMSB_L-100TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            #['/GMSB_L-100TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            #['/GMSB_L-150TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-150TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            #['/GMSB_L-200TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-200TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            #['/GMSB_L-250TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
#            ['/GMSB_L-250TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],

            #['/MET/Run2017B-17Nov2017-v1/AOD',''],#12.5
            #['/MET/Run2017C-17Nov2017-v1/AOD',''],#27.9
            #['/MET/Run2017D-17Nov2017-v1/AOD',''],#5.3
            #['/MET/Run2017E-17Nov2017-v1/AOD',''],#28.5
            #['/MET/Run2017F-17Nov2017-v1/AOD',''],#75.3

            # PD: /DoubleEG/Run2017*-17Nov2017-v1/AOD
            #['/DoubleEG/Run2017B-17Nov2017-v1/AOD',''],
            #['/DoubleEG/Run2017D-17Nov2017-v1/AOD',''],
            #['/DoubleEG/Run2017E-17Nov2017-v1/AOD',''],
            #['/DoubleEG/Run2017F-17Nov2017-v1/AOD',''],

            #[/DoubleEG/Run2017B-09Aug2019_UL2017-v1/AOD,''],
            #[/DoubleEG/Run2017B-17Nov2017-v1/AOD,''],
            #[/DoubleEG/Run2017B-22Jun2017-v1/AOD,''],
            #[/DoubleEG/Run2017B-23Jun2017-v1/AOD,''],
            #[/DoubleEG/Run2017C-09Aug2019_UL2017-v1/AOD,''],
            #[/DoubleEG/Run2017C-11Dec2018_AB-v1/AOD,''],
            #[/DoubleEG/Run2017C-11Dec2018_AC-v1/AOD,''],
            #[/DoubleEG/Run2017C-12Sep2017-v1/AOD,''],
            #[/DoubleEG/Run2017D-09Aug2019_UL2017-v1/AOD,''],
            #[/DoubleEG/Run2017D-17Nov2017-v1/AOD,''],
            #['/DoubleEG/Run2017E-09Aug2019_UL2017-v1/AOD',''],
            #[/DoubleEG/Run2017E-17Nov2017-v1/AOD,''],
            #[/DoubleEG/Run2017F-09Aug2019_UL2017-v1/AOD,''],
            #[/DoubleEG/Run2017F-17Nov2017-v1/AOD,''],

            #['/DoubleEG/Run2017F-09Aug2019_UL2017-v1/AOD',''],
            #['/DoubleEG/Run2017E-09Aug2019_UL2017-v1/AOD',''],

            ['/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',''],

	    ]
 
        for inDO in inputDataAndOpts:
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
	         #print( inDO[0] )
            primaryDataset = inDO[0].split('/')[1]
            #runEra         = inDO[0].split('/')[2]
            runEra         = ((inDO[0].split('/')[2]).split('_')[0]+'_'+(inDO[0].split('/')[2]).split('_')[1]).split('-PU')[0] # MC
            dataset	   = inDO[0].split('/')[3]

            #trial          = 'gammares_ratio_126_v2'
            #trial          = 'gammares_cc_140_v2'
            #trial          = 'gammares_ttcc_140_v5' # 24C and earlier only
            #trial          = 'gammares_llpana_pd'
            trial          = 'gammares_llpana_mc'
            #trial          = 'gammares_llpana_qcd'

            config.Data.outLFNDirBase    = "/store/user/jaking/ecalTiming/"+trial+"/"
            ##config.General.requestName   = trial+"_"+primaryDataset+"_"+runEra+"_"+runs+"_"+dataset
            ##config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_"+runs+"_dispho"
            ##config.Data.outputDatasetTag = runEra+"_"+runs+"_"+dataset
            config.General.requestName   = primaryDataset+"_"+dataset+"_"+runEra+"_t2_request"
            config.Data.outputDatasetTag = dataset+"_"+runEra

            #config.JobType.pyCfgParams   = ['globalTag=140X_dataRun3_Prompt_v2', 'outputFileName=output.root','doTwoTier=False','doDiag=True']  # 2024 tested
            #config.JobType.pyCfgParams   = ['globalTag=140X_dataRun3_Prompt_Candidate_2024_05_31_21_23_47', 'outputFileName=output.root','doTwoTier=False','doDiag=True']  # 5 GeV cali - EE cali bad
            ##config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v4','doTwoTier=False','doDiag=False']  
            #config.JobType.pyCfgParams   = ['globalTag=112X_dataRun3_Prompt_v2', 'outputFileName=output.root','doTwoTier=False','doDiag=False']]  # 2018A tested
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_PromptAnalysis_v1', 'outputFileName=output.root','doTwoTier=False','doDiag=False']] # ABCD
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v8', 'outputFileName=output.root','doTwoTier=False','doDiag=False']] # EF
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v10', 'outputFileName=output.root','doTwoTier=False','doDiag=False']] # G
#>>>>>>>>>>>>>>>>>>>     #MET 2017E
            #config.JobType.pyCfgParams   = ['globalTag=94X_dataRun2_ReReco_EOY17_v1','outputFileName=output.root','doTwoTier=False','doDiag=True'] # MET 2017E
#>>>>>>>>>>>>>>>>>>>     #GMSB
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','outputFileName=output.root','doTwoTier=False','doDiag=True'] # GMSB
#>>>>>>>>>>>>>>>>>>>     #DoubleEG/Run2017*-17Nov2017-v1/AOD good for #globalTag=94X_dataRun2_ReReco_EOY17_v2
            #config.JobType.pyCfgParams   = ['globalTag=94X_dataRun2_ReReco_EOY17_v2','outputFileName=output.root','doTwoTier=False','doDiag=True']
#>>>>>>>>>>>>>>>>>>>     #/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM
            config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','outputFileName=output.root','doTwoTier=False','doDiag=True']

#------------------------------------------------------------------------------------------------------------------
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

        subcrab( "000000-999999","",2500)# anything


########################################################

submit_run()
