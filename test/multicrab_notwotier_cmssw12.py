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

        config.Data.inputDataset = None
        #config.Data.useParent      = True
	    #config.Data.secondaryInputDataset = secInputPaths
        #config.Data.useParent      = False
        #config.Data.lumiMask     = inputDir+inputJSON
        #config.Data.splitting    = 'LumiBased'
        config.Data.splitting    = 'EventAwareLumiBased'
        #config.Data.splitting    = 'Automatic'
        #config.Data.unitsPerJob  =  300 # for auto job splitting
        config.Data.runRange  = runs #'321122-321128'
        config.Data.unitsPerJob  =  250000 # unitsPerJob = 1000 for 321122-321128 and maxMemoryMB = 4000  on EventAwareLumiBased


        #config.Data.outputDatasetTag = 'reducedRAW_EGamma_ntuple'
	     
        config.JobType.allowUndistributedCMSSW = True
        config.Data.publication      = False
        config.Site.storageSite      = 'T3_US_FNALLPC'
        config.Data.outLFNDirBase    = '/store/user/jaking/ecalTiming/'
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDataAndOpts = [

            ['/EGamma/Run2018A-17Sep2018-v2/MINIAOD'], # for ABC change to GT 102X_dataRun2_v13

            #['/EGamma/Run2022A-PromptReco-v1/MINIAOD'], # 2022 ABCD Prompt 124X_dataRun3_PromptAnalysis_v1  352499 - 355062    
            #['/EGamma/Run2022B-PromptReco-v1/MINIAOD'],	#											 		355094 - 355769											 
            #['/EGamma/Run2022C-PromptReco-v1/MINIAOD'],	#											 		355809 - 357482
            #['/EGamma/Run2022D-PromptReco-v1/MINIAOD'],	#											 		357538 - 357733
            #['/EGamma/Run2022D-PromptReco-v2/MINIAOD'],	#											 		357734 - 357902
            #['/EGamma/Run2022D-PromptReco-v3/MINIAOD'],	#											 		358381
            #['/EGamma/Run2022E-PromptReco-v1/MINIAOD'], # 2022 EF Prompt 124X_dataRun3_Prompt_v8	 		359090 - 360327
            #['/EGamma/Run2022F-PromptReco-v1/MINIAOD'], #											 		360389 - 362167
            #['/EGamma/Run2022G-PromptReco-v1/MINIAOD']  # ? 2022 G Prompt 124X_dataRun3_Prompt_v10   		362399 - 362760

	    ]
 
        for inDO in inputDataAndOpts:
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
	         #print( inDO[0] )
            primaryDataset = inDO[0].split('/')[1]
            runEra         = inDO[0].split('/')[2]
            dataset	   = inDO[0].split('/')[3]
            trial          = 'gammares_ratio_126_v2'

            config.General.requestName   = trial+"_"+primaryDataset+"_"+runEra+"_"+runs+"_"+dataset+"_dispho"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_"+runs+"_dispho"

            ##config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v4',#'nThreads='+str(config.JobType.numCores), 
            config.JobType.pyCfgParams   = ['globalTag=112X_dataRun3_Prompt_v2', 'outputFileName=output.root']  # 2018A tested
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_PromptAnalysis_v1', 'outputFileName=output.root'] # ABCD
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v8', 'outputFileName=output.root'] # EF
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v10', 'outputFileName=output.root'] # G

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

		#subcrab( "352400-358400","",2500)# Run3 ABCD
        #subcrab( "359000-362200","",2500)# Run3 EF
        #subcrab( "362300-362800","",2500)# Run3 G

        #subcrab( "316241-316245","",2500) # very small test batch 18A
        subcrab( "316000-316499","",2500) # 18A	
        #subcrab( "357101-357268","",2500)#22C?
		#subcrab( "360395-360415","",2500)#22F?


########################################################

submit_run()

