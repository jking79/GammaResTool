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

#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GMSB_R17_MET75_v20/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GJETS_R17_MET75_v20/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_QCD_R17_MET75_v20/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DEG_R17_v18/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R17E_MET75_v20/'
#command = eosll+mspc+'ecalTiming/gammares_llpana/'
#command = eosll+mspc+'/ecalTiming/gammares_llpana_qcd/'
#command = eosll+mspc+'/ecalTiming/gammares_llpana_v2/'
#command = eosll+mspc+'/ecalTiming/gammares_llpana_pd/'
#command = eosll+mspc+'ecalTiming/gammares_llpana_mc/'
#command = eosll+mspc+'ecalTiming/gammares_llpana_pd/'
#command = eosll+mspc+'ecalTiming/gammares_ccval/'
#command = eosll+mspc+'ecalTiming/gammares_r24fprompt/'
#command = eosll+mspc+'ecalTiming/gammares_ecaldpg_tevjets_prompt_v3/'
#command = eosll+mspc+'ecalTiming/gammares_ecaldpg_tevjets_prompt_oot3_v3/'
#command = eosll+mspc+'ecalTiming/gammares_llpana_pd/'
#command = eosll+mdis+'kuncali/gammares_cali/'
#command = eosll+mspc+'ecalTiming/gammares_r24f_cctest/'
#command = eosll+mdis+'kuncali/gammares_cali/'
#command = eosll+mspc+'ecalTiming/gammares_cali_mc/'
command = eosll+mdis+'kuncali/gammares_cali_mc/'

rootfile = '.root'

#dirselect = 'GMSB'
#dirselect = 'QCD'
#dirselect = 'TTJets'
#dirselect = 'JetHT'
#dirselect = 'JetMET1'
#dirselect = 'GJets'
#dirselect = 'DiPhotonJetsBox'
#dirselect = 'Run2017B'
#dirselect = 'MET'
#dirselect = 'DoubleEG'
#dirselect = 'MET_R17E_MET75'
#dirselect = 'EGamma1'
#dirselect = 'ZJetsToNuNu'
dirselect = 'DYJetsToLL'


#dirselect = ''

#version = 'Run2018D'
#version = 'GJets'
#version = 'QCD'
#version = 'GMSB'
#version = 'JetHT' 
#version = 'DiPhotonJetsBox'
#version = 'MET'
#version = 'DoubleEG'
#version = 'MINIAOD'
#version = 'EGamma1'
#version = 'Run2017E'
#version = 'JetMET1'
#version = 'ZJetsToNuNu'
version = 'RunIIFall17DRPremix'

debug = True
#debug = False

targdirs = []

dirls = bashout( command ).splitlines()
print( dirls )
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
    if 'none' not in version :
        for subdir in subdir2 : 
            command3 = command+thesubdir+'/'+subdir+'/'
            if debug : print( command3 )
            subdir3 = bashout( command3 ).rstrip().splitlines()
            for subsubdir in subdir3 : 
                if version in subdir : subdirlist2.append(thesubdir+'/'+subdir+'/'+subsubdir)
   
    if debug : print( subdirlist2 ) 
    if 'none' in version : subdirlist2 = subdir2
    for thesubdir2 in subdirlist2 :
        command4 = command+thesubdir2+'/'
        subdir4 = bashout( command4 ).rstrip().splitlines()
        #print( thesubdir+subdir2+'/0000/' )
        for subdir in subdir4 :
            subdirlist3.append(thesubdir2+'/'+subdir+'/')
            #command5 = command+thesubdir+subdir+'/'
            #subdir5 = bashout( command5 ).rstrip().splitlines()
            #for subsubdir in subdir5 :
                #subdirlist3.append(thesubdir+subdir+'/'+subsubdir+'/')
    
    
    if debug : print( subdirlist3 )
    for subdir2 in subdirlist3:
    	lists = bashout( command+subdir2 ).rstrip().splitlines()
    	for lline in lists :
    		if rootfile in lline : filelist.append(subdir2+lline)
   
    select =  line2.split("Tune")
    #outfile = 'kuntuple_' + select[0] + '_p3_v21.txt'
    #outfile = 'egammares_' + select[0] + '_v21.txt'
    outfile = 'egammares_' + select[0] + 'MINIAOD_RunIIFall17DRPremix_cali_v21.txt'
    #print( outfile )
    outf = open( outfile, 'w' )
    for thefile in filelist:
    	outf.write( thefile + '\n' )
    outf.close()



