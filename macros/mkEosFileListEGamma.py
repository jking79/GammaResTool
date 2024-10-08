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
mgrp = '/store/user/lpcsusylep/jaking/'
eosll = 'eos root://cmseos.fnal.gov ls '

#command = eosll+mspc+'ecalTiming/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v2b/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v4_flipped/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v3/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v5_phoclean/EGamma/'
#command = eosll+mgrp+'ecalTiming/gammares_tt_kucc_126_v6_diag/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v7_diag/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v7_diag_unclean/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v11_diag/EGamma/'
#command = eosll+mspc+'ecalTiming/gammares_tt_kucc_126_v10_reso/EGamma/'
#command = eosll+mgrp+'LLPGamma/llpga_GMSB_AOD_v60/'
#command = eosll+mgrp+'LLPGamma/llpga_GJets_AOD_v60/'
#command = eosll+mspc+'/ecalTiming/gammares_tt_kucc_126_v11_ccEncDiag/EGamma/'
#command = eosll+mspc+'/ecalTiming/gammares_ttcc_1307_v11_diag/EGamma1/'
#command = eosll+mspc+'/ecalTiming/gammares_ttcc_140_v11_test/EGamma1/'
#command = eosll+mspc+'/ecalTiming/EGamma0/'
#command = eosll+mspc+'/ecalTiming/gammares_ttcc_140_v11_diag_mod1_exp3/EGamma1/'
#command = eosll+mspc+'/ecalTiming/gammares_llpana/'
#command = eosll+mspc+'/ecalTiming/gammares_llpana_pd/MET/'
command = eosll+mspc+'/ecalTiming/gammares_llpana_qcd/'

version = ''
#version = '_v11_'
#version = '_noOOTAmp_'
#version = '_wthOOTAmp_'

rootfile = '.root'

#dirselect = ''
#dirselect = 'gammares_tt_kucc_126_v2a_EGamma_MINIAOD_Run2022C'
#dirselect = 'gammares_ratio_126_v2_EGamma_MINIAOD_Run2018A-17Sep2018-v2_316000-316499_dispho'
#dirselect = 'gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C'
#dirselect = 'gammares_tt_kucc_126_v4_flipped_EGamma_MINIAOD_Run2022F'
#dirselect = 'gammares_tt_kucc_126_v3_EGamma_MINIAOD_Run2022D-PromptReco-v2'
#dirselect = 'gammares_tt_kucc_126_v3_EGamma_MINIAOD_Run2022D-PromptReco-v1'
#dirselect = 'gammares_tt_kucc_126_v4_flipped_EGamma_MINIAOD_Run2022G'
#dirselect = 'gammares_tt_kucc_126_v3_EGamma_MINIAOD_Run2022G'
#dirselect = 'gammares_tt_kucc_126_v5_phoclean'
#dirselect = 'gammares_tt_kucc_126_v6_diag_EGamma_MINIAOD_Run2022C-PromptReco-v1_355794-359021_dispho'
#dirselect = 'gammares_tt_kucc_126_v6_diag_EGamma_MINIAOD_Run2022E-PromptReco-v1_359022-362760_dispho'
#dirselect = 'gammares_tt_kucc_126_v5_phoclean_EGamma_MINIAOD_Run2022E-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v6_diag_EGamma_MINIAOD_Run2022G-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v7_diag_EGamma_MINIAOD_Run2022E-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v7_diag_EGamma_MINIAOD_Run2018A-17Sep2018-v2_000000-999999'
#dirselect = 'gammares_tt_kucc_126_v7_diag_unclean_EGamma_MINIAOD_Run2022E-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v10_diag_EGamma_MINIAOD_Run2022E-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v10_diag_EGamma_MINIAOD_Run2018A-17Sep2018-v2_315257-316993'
#dirselect = 'gammares_tt_kucc_126_v10_diag_EGamma_MINIAOD_Run2022G-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v10_reso_EGamma_MINIAOD_Run2022C-PromptReco-v1_355794-359021'
#dirselect = 'gammares_tt_kucc_126_v10_reso_EGamma_MINIAOD_Run2022D-PromptReco-v1_355794-359021'
#dirselect = 'gammares_tt_kucc_126_v10_reso_EGamma_MINIAOD_Run2022D-PromptReco-v2_355794-359021'
#dirselect = 'gammares_tt_kucc_126_v10_reso_EGamma_MINIAOD_Run2022E-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v10_reso_EGamma_MINIAOD_Run2022F-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v10_reso_EGamma_MINIAOD_Run2022G-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v10_diag2_EGamma_MINIAOD_Run2022G-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v11_diag_EGamma_MINIAOD_Run2022G-PromptReco-v1_359022-362760'
#dirselect = 'gammares_tt_kucc_126_v11_ccEncDiag_EGamma_MINIAOD_Run2022C-PromptReco-v1_355890-355895'
#dirselect = 'gammares_ttcc_1307_v11_diag_EGamma1_MINIAOD_Run2023B-PromptReco-v1_366323-367065'
#dirselect = 'Run2024C-PromptReco'
#dirselect = 'Run2023D-PromptReco'
dirselect = 'QCD'
#dirselect = 'AOD_Run2017E-17Nov2017-v1'

#debug = True
debug = False

#deep = True
deep = False

targdirs = []
subdirlist1 = []
subdirlist2 = []
subdirlist3 = []
filelist = []
theFileList = ''

dirls = bashout( command ).splitlines()
if debug : print( '-------------------------------------------------')
for line in dirls:
	#print( line )
	if dirselect in line : targdirs.append( line )
if debug : print( targdirs )
if deep :
	for mydir in targdirs:
		command1 = command+mydir+'/'
		subdir1 = bashout( command1 ).rstrip().splitlines()
		#print( subdir1 )
		#print( mydir+'/'+subdir1+'/' )
		for line in subdir1 : 
			#print( line )
			if version in line : 
				subdirlist1.append( mydir+'/'+line+'/' )
		#print( subdirlist1 )
else : 
	for mydir in targdirs:
		#print( mydir+'/' )
		subdirlist1.append( mydir+'/' )

if debug : print( subdirlist1 )
for thesubdir in subdirlist1 :
	command2 = command+thesubdir+'/'
	subdir2 = bashout( command2 ).rstrip().splitlines()
	#print( thesubdir+subdir2 )
	for subdir in subdir2 : 
		command3 = command+thesubdir+subdir+'/'
		subdir3 = bashout( command3 ).rstrip().splitlines()
		for subsubdir in subdir3 : 
			#print( thesubdir+subdir+'/'+subsubdir+'/' )
			subdirlist2.append(thesubdir+subdir+'/'+subsubdir+'/')


if debug : print( subdirlist2 )
for subdir2 in subdirlist2:
	lists = bashout( command+subdir2 ).rstrip().splitlines()
	for line in lists :
		if rootfile in line : 
			#print( subdir2+line )
			filelist.append(subdir2+line)

#for thefile in filelist:
#	print( thefile )

#select =  line2.split("Tune")
select = 'Met_PD_AOD_Run2017E-17Nov2017'

outfile = 'egammares_' + select + 'v2.txt'
print( outfile )
outf = open( outfile, 'w' )
#filelist = subdirlist3
for thefile in filelist:
    outf.write( thefile + '\n' )
outf.close()


	#filename = 'tmp_'+subdir2.split('/')[1]+'.root '
	#print( filename )
	#lists = bashout( "eosls "+mspc+"LLPGamma/"+subdir2 ).rstrip()
	#print( subdir2 )
	#haddcommand = "hadd -f "+filename+"`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+" | grep '\.root'`"
	#haddcommand = "`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+" | grep '\.root'`"
	#haddcommand = "`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+"`"
	#haddcommand = "hadd "+filename+lists
	#print( mspc+"LLPGamma/"+subdir2 )
	#print( haddcommand )
	#doCommand( haddcommand )
	#print( '---------------------------------------------------' )

	
#1print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	
