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
        #inputJSON    = 'Cert_352499-362760_13TeV_PromptReco_Collisions22_jwkgoldenjson.txt'
        #inputJSON    = 'Cert_Collisions2023_eraB_366403_367079_Golden.json'
        #inputJSON    = 'Cert_Collisions2023_eraC_367095_368823_Golden.json'
        inputJSON    = 'Cert_Collisions2023_eraD_369803_370790_Golden.json'

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        import CRABClient
        from CRABClient.UserUtilities import config
        config = config()

        config.General.workArea    = options.workArea
        #config.General.requestName = None

        config.JobType.pluginName  = 'Analysis'
        #config.JobType.psetName    = 'gammares_raw_twotier_multi.py'
        config.JobType.psetName    = 'gammares_twotier_multi.py'
        #config.JobType.numCores    = 8
        #config.JobType.maxMemoryMB = 2250 #reqmem
        #config.JobType.maxJobRuntimeMin = 1500
        config.JobType.pyCfgParams = None
        #config.JobType.inputFiles  = [ inputDir+inputPaths , inputDir+inputFilters , inputDir+inputFlags ]

        config.Data.inputDataset = None
        config.Data.useParent      = True
	#config.Data.secondaryInputDataset = secInputPaths
        #config.Data.useParent      = False
        config.Data.lumiMask     = inputDir+inputJSON
        #config.Data.splitting    = 'LumiBased'
        config.Data.splitting    = 'EventAwareLumiBased'
        #config.Data.splitting    = 'Automatic'
        #config.Data.unitsPerJob  =  300 # for auto job splitting
        config.Data.runRange  = runs #'321122-321128'
        #config.Data.unitsPerJob  =  75000 # unitsPerJob = 1000 for 321122-321128 and maxMemoryMB = 4000  on EventAwareLumiBased
        config.Data.unitsPerJob  =  25000

        #config.Data.outputDatasetTag = 'reducedRAW_EGamma_ntuple'
	     
        config.JobType.allowUndistributedCMSSW = True
        config.Data.publication      = False
        config.Site.storageSite      = 'T3_US_FNALLPC'
        #config.Data.outLFNDirBase    = '/store/user/jaking/ecalTiming/'
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDataAndOpts = [

             #['/EGamma/Run2018A-17Sep2018-v2/MINIAOD'], # for ABC change to GT globalTag=112X_dataRun3_Prompt_v2 

             ##Comm['/EGamma/Run2022A-PromptReco-v1/MINIAOD'], # 2022 ABCD Prompt 124X_dataRun3_PromptAnalysis_v1  352499 - 355062    
             ##Comm['/EGamma/Run2022B-PromptReco-v1/MINIAOD'],	#											 		355094 - 355769								
             #['/EGamma/Run2022C-PromptReco-v1/MINIAOD'],	#											 		355809 - 357482
             #['/EGamma/Run2022D-PromptReco-v1/MINIAOD'],	#											 		357538 - 357733
             #['/EGamma/Run2022D-PromptReco-v2/MINIAOD'],	#											 		357734 - 357902
             ##['/EGamma/Run2022D-PromptReco-v3/MINIAOD'],	#											 		358381

             #['/EGamma/Run2022E-PromptReco-v1/MINIAOD'], # 2022 EF Prompt 124X_dataRun3_Prompt_v10	 		359090 - 360327
             #['/EGamma/Run2022F-PromptReco-v1/MINIAOD'], #											 		360389 - 362167
             #['/EGamma/Run2022G-PromptReco-v1/MINIAOD']  # ? 2022 G Prompt 124X_dataRun3_Prompt_v10   		362399 - 362760

             #['/EGamma0/Run2023A-PromptReco-v2/MINIAOD'],    # 366323 - 366361  130X_dataRun3_Prompt_v2
             #['/EGamma0/Run2023B-PromptReco-v1/MINIAOD'],    # 366403 - 367065 130X_dataRun3_Prompt_v2
             #['/EGamma0/Run2023C-PromptReco-v1/MINIAOD'],    # 367094 - 367476 130X_dataRun3_Prompt_v3
             #['/EGamma0/Run2023C-PromptReco-v2/MINIAOD'],    # 367516 - 367758 130X_dataRun3_Prompt_v3
             ##['/EGamma0/Run2023C-PromptReco-v3/MINIAOD'],    # 367622 - 367758 130X_dataRun3_Prompt_v3
             #['/EGamma0/Run2023C-PromptReco-v4/MINIAOD'],    # 367770 - 368389 130X_dataRun3_Prompt_v3
             #['/EGamma0/Run2023D-PromptReco-v1/MINIAOD'],    # 369844 - 370496 130X_dataRun3_Prompt_v4
             #['/EGamma0/Run2023D-PromptReco-v2/MINIAOD'],    # 370666 - 370790 130X_dataRun3_Prompt_v4

             #['/EGamma1/Run2023A-PromptReco-v2/MINIAOD'],    # 366323 - 366361  130X_dataRun3_Prompt_v2
             #['/EGamma1/Run2023B-PromptReco-v1/MINIAOD'],    # 366403 - 367065 130X_dataRun3_Prompt_v2
             #['/EGamma1/Run2023C-PromptReco-v1/MINIAOD'],    # 367094 - 367515 130X_dataRun3_Prompt_v3
             #['/EGamma1/Run2023C-PromptReco-v2/MINIAOD'],    # 367516 - 367758 130X_dataRun3_Prompt_v3
             ##['/EGamma1/Run2023C-PromptReco-v3/MINIAOD'],    # 367661 - 367758 130X_dataRun3_Prompt_v3
             #['/EGamma1/Run2023C-PromptReco-v4/MINIAOD'],    # 367770 - 368412 130X_dataRun3_Prompt_v3
             ['/EGamma1/Run2023D-PromptReco-v1/MINIAOD'],    # 369844 - 370472 130X_dataRun3_Prompt_v4
             #['/EGamma1/Run2023D-PromptReco-v2/MINIAOD'],    # 370666 - 370790 130X_dataRun3_Prompt_v4

	    ]
 
        for inDO in inputDataAndOpts:
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
            #print( inDO[0] )
            primaryDataset = inDO[0].split('/')[1]
            runEra         = inDO[0].split('/')[2]
            dataset	   = inDO[0].split('/')[3]

            #trial          = 'gammares_tt_kucc_126_v4_spike'
            #trial          = 'gammares_tt_kucc_126_v5_phoclean'
            #trial          = 'gammares_tt_kucc_126_v7_diag_unclean'
            #trial          = 'gammares_tt_kucc_126_v11_diag' # added spike rechits with 2 gev cut
            #trial          = 'gammares_ttcc_140_v11_diag_mod1' # kOOT flag disable online calibration && -1.65 while setting koot flag 
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_slewc' # kOOT flag disable online calibration && -1.65 while setting koot flag + gain switch punt EB
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_exp1' # no gain and no min of zero amplitude set slew to 0
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_dslew' # bad pulse shape for double gain slew
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_exp2' # no min amplitude - no slew correction
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_gstest' # no min amplitude - has gain vs slew test
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_qfix' # no min amplitude - slew + cc > 0 set to 0
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_fix1' # use weights to zero out terms with slew issue in computeCC
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_ebsf' # added parameters to turn on and off slew correction for EB && EE seperatly - eb true ee false
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_nosf' # added parameters to turn on and off slew correction for EB && EE seperatly - eb false ee false
            #trial          = 'gammares_ttcc_140_v11_diag_ebsf_ccgt' # using 140X_dataRun3_Candidate_2024_06_26_20_41_23 GT for CC calibration in 22/23
            #trial          = 'gammares_ttcc_140_v11_diag_mod1_exp3' # changed pulse template normilization to max sample amplitude instead of sum of amplitudes
            trial          = 'gammares_ttcc_140_v11_diag_mod1_exp3' # changed pulse template normilization to max pulse template amplitude instead of sum of amplitudes

            config.Data.outLFNDirBase    = "/store/user/jaking/ecalTiming/"+trial+"/"
            #config.General.requestName   = trial+"_"+primaryDataset+"_"+runEra+"_"+runs+"_"+dataset
            #config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_"+runs+"_dispho"
            #config.Data.outputDatasetTag = runEra+"_"+runs+"_"+dataset
            config.General.requestName   = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_"+runs+"_request"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_"+runs


            ##config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v4',#'nThreads='+str(config.JobType.numCores), 
            #config.JobType.pyCfgParams   = ['globalTag=112X_dataRun3_Prompt_v2', 'outputFileName=output.root'] #'nThreads='+str(config.JobType.numCores), # run2 2018
			########## 2022 
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_PromptAnalysis_v1', 'outputFileName=output.root'] # ABCD prompt
            ##config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_v3', 'outputFileName=output.root' # ABCD rereco
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_Prompt_v10', 'outputFileName=output.root'] # EFG promp
            ##config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_v14', 'outputFileName=output.root' # E rereco
			########## 2023
            #config.JobType.pyCfgParams   = ['globalTag=130X_dataRun3_Prompt_v2', 'outputFileName=output.root'] # 23 AB prompt
            #config.JobType.pyCfgParams   = ['globalTag=130X_dataRun3_Prompt_v3', 'outputFileName=output.root'] # 23 Cv1-3 prompt
            #config.JobType.pyCfgParams   = ['globalTag=130X_dataRun3_Prompt_v3_forRun368229_v1', 'outputFileName=output.root'] # 23 Cv4 prpt
            #config.JobType.pyCfgParams   = ['globalTag=130X_dataRun3_Candidate_2023_08_08_21_30_44', 'outputFileName=output.root'] # 23 CC GT
            config.JobType.pyCfgParams   = ['globalTag=140X_dataRun3_Candidate_2024_06_26_20_41_23', 'outputFileName=output.root'] # 22/23 CC GT 140
            #config.JobType.pyCfgParams   = ['globalTag=140X_dataRun3_Prompt_v2', 'outputFileName=output.root','doTwoTier=True','doDiag=True']

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

	#subcrab( "315257-316993","",2500)#Run2 2018A

        #subcrab( "366850-366873","",2500)#

	#subcrab( "366323-367065","",2500)# Run 3 23 AB
        #subcrab( "366323-366361","","2500")# Run 3 23 A
        #subcrab( "366403-367079","","2500")# Run 3 23 B
        #subcrab( "367094-367758","",2500)# Run 3 23 C1&3
        #subcrab( "367770-368412","",2500)# Run 3 23 C4

	##subcrab( "352400-358400","",2500)# Run3 22ABCD
        ##subcrab( "359000-362200","",2500)# Run3 22EF
        ##subcrab( "362300-362800","",2500)# Run3 22G

        #subcrab( "355794-359021","",2500)# Run3 CD
        #subcrab( "359022-362760","",2500)# Run3 EFG
	#subcrab( "355890-355895","",2500)

        #subcrab( "316241-316245","",2500) # very small test batch 18A
        #subcrab( "357101-357268","",2500)#22C?
	#subcrab( "360395-360415","",2500)#22F?

        #subcrab( "369844-369999","",2500)#23D (369844-369999),(370092-370243),(370293-370580)
        #subcrab( "370293-370580","",2500)#23D 
        subcrab( "370496-370580","",2500)#23D

########################################################

submit_run()

