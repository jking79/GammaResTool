import subprocess
import sys
import os

# uses python3

def bash( bashCommand ):
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	#process = subprocess.Popen(bashCommand.split())
	output, error = process.communicate()
	return output ,error

def bashout( command ):
	#output = subprocess.check_output( command, shell=True)
	output = subprocess.run( command, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True )
	return output.stdout	

def doCommand( command ):
	output = os.system( command )
	return output

mspc = '/store/user/jaking/'
mdis = '/store/user/lpcsusylep/jaking/'
eosll = 'eos root://cmseos.fnal.gov ls '
#command = eosll+mspc+'LLPGamma/llpga_GMSB_AOD_v48/'
#command = eosll+mspc+'A/'
#command = eosll+mdis+'LLPGamma/llpga_GMSB_AOD_v58/'
#command = eosll+mdis+'/ecalTiming/tt_KUCCRes_126_Test/EGamma/'
#command = eosll+mdis+'LLPGamma/llpga_GJets_AOD_v57/'
#command = eosll+mdis+'LLPGamma/llpga_GMSB_AOD_v59/'
#command = eosll+mspc+'ecalTiming/gammares_ttcc_131_v11_diag/'
#command = eosll+mdis+'ecalTiming/EGamma/'
#command = eosll+mspc+'EGamma/'
#command = eosll+mdis+'KUCMSNtuple/GMSB_AOD_v1/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GJETS_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_WJETS_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_QCD_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GMSB_AOD_v6/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_ZJETS_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DYTT_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_JetHTi_Met50_AOD_v2/'
#command = eosll+mdis+'EcalTiming/ZeroBias/'
#command = eosll+mspc+'ecalTiming/gammares_ttcc_140_v11_diag_mod1_nosf/EGamma1/'
#command = eosll+mspc+'ecalTiming/gammares_ttcc_140_v11_diag_ebsf_ccgt/EGamma1/'
#command = eosll+mspc+'ecalTiming/gammares_ttcc_140_v11_diag_mod1_exp3/EGamma1/'
#command = eosll+mspc+'ecalTiming/gammares_ttcc_140_v11_mod1_pr_test2/EGamma1/'
#command = eosll+mspc+'ecalTiming/gammares_llpana/'
#command = eosll+mspc+'/ecalTiming/gammares_llpana_qcd/'
#command = eosll+mspc+'/ecalTiming/gammares_llpana/MET/'
command = eosll+mspc+'/ecalTiming/gammares_llpana_mc/'

version = ''
#version = '_v11_'
#version = '_noOOTAmp_'
#version = '_wthOOTAmp_'
#version = 'DYJetsToLL_M-50'

rootfile = '.root'

#dirselect = 'HTo2LongLivedTo4b'
#dirselect = '_newRtParams4_v26b_'
#dirselect = '_newRtParams3_test_v26_'
#dirselect = 'tt_kurhs_124cc5_cert'
#dirselect = '22eraC_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = '22eraC_CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'
#dirselect = 'CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'
#dirselect = 'GMSB'
#dirselect = 'AOD'
#dirselect = 'WJetsToLNu_HT-800'
#dirselect = 'QCD_HT100to200'
#dirselect = 'GMSB_L-400TeV'
dirselect = 'DY1JetsToLL_M-50'
#dirselect = 'TTJets'
#dirselect = 'v2_EGamma'
#dirselect = 'ecaltiming_dqm_132r3prompt3'
#dirselect = 'cc_140_json_v2_EGamma1_MINIAOD_Run2024C-PromptReco-v1'
#dirselect = 'MINIAOD_Run2024D-PromptReco-v1'
#dirselect = 'Run2023D-PromptReco'
#dirselect = 'QCD'
#dirselect = 'AOD_Run2017E-17Nov2017-v1'
#dirselect = 'MINIAOD'

#dirselect = ''

#debug = True
debug = False

#deep = True
deep = False

targdirs = []

dirls = bashout( command ).splitlines()
print( '************************************************')
for line in dirls:
	#print( line )
	if dirselect in line : targdirs.append( line )
    #targdirs.append( line )
print( targdirs )

for line2 in targdirs :

    subdirlist1 = []
    subdirlist2 = []
    subdirlist3 = []
    subdirlist4 = []
    filelist = []
    theFileList = ''

    #for mydir in targdirs:
    print( '-------------------------------------------------')
    print( line2 )
    print( '-------------------------------------------------')
    #subdirlist1.append( line+'/' )
    
    #if debug : print( subdirlist1 )
    #for thesubdir in subdirlist1 :
    thesubdir = line2
    command2 = command+thesubdir+'/'
    if debug : print( command2 )
    subdir2 = bashout( command2 ).rstrip().splitlines()
    for subdir in subdir2 : 
    	command3 = command+thesubdir+'/'+subdir+'/'
    	subdir3 = bashout( command3 ).rstrip().splitlines()
    	for subsubdir in subdir3 : 
    	    subdirlist2.append(thesubdir+'/'+subdir+'/'+subsubdir)
   
    if debug : print( subdirlist2 ) 
    for thesubdir2 in subdirlist2 :
        command4 = command+thesubdir2+'/'
        subdir4 = bashout( command4 ).rstrip().splitlines()
        #print( thesubdir+subdir2+'/0000/' )
        for subdir in subdir4 :
            subdirlist3.append(thesubdir2+'/'+subdir)
            #command5 = command+thesubdir+subdir+'/'
            #subdir5 = bashout( command5 ).rstrip().splitlines()
            #for subsubdir in subdir5 :
                #subdirlist3.append(thesubdir+subdir+'/'+subsubdir+'/')

    if debug : print( subdirlist3 )
    for thesubdir3 in subdirlist3 :
        command5 = command+thesubdir3+'/'
        subdir5 = bashout( command5 ).rstrip().splitlines()
        #print( thesubdir+subdir2+'/0000/' )
        for subdir in subdir5 :
            subdirlist4.append(thesubdir3+'/'+subdir)


    if debug : print( subdirlist4 )
    for subdir4 in subdirlist4:
    	lists = bashout( command+subdir4 ).rstrip().splitlines()
    	for lline in lists :
    	    if rootfile in lline : filelist.append(subdir4+lline)


    select =  line2.split("Tune")
    outfile = 'egammares_' + select[0] + 'v2.txt'
    print( outfile )
    outf = open( outfile, 'w' )
    filelist = subdirlist3
    for thefile in filelist:
    	outf.write( thefile + '\n' )
    outf.close()



