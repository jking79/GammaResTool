//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSRootHelperBaseClass.hh"

#ifndef KUCMSTimeCalibrationClass
#define KUCMSTimeCalibrationClass

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////////////////////
//
//  this code assumes that the ecal_config directory is included beside this class
//  ecal_config should contain :
//
//  	fullinfo_detids_EB.txt
//  	fullinfo_detids_EE.txt
//  	fullinfo_v2_detids_EB.txt
//  	fullinfo_v2_detids_EE.txt
//  	howto.txt
//  	reducedinfo_detids.txt
//  	rhid_i12_list.txt
//  	rhid_info_list.txt 
//
/////////////////////////////////////////////////////////////////////////////////////////

class KUCMSTimeCalibration : KUCMSRootHelperBaseClass {

	public:

    KUCMSTimeCalibration();
	~KUCMSTimeCalibration();


	private:

	enum ECAL {EB, EM, EP, NONE};

    std::string detIDConfigEB;
    std::string detIDConfigEE;

    std::string caliFileDir;
    std::string caliRunConfig;
    std::string smearConfig;
    std::string caliTFile;

    std::string xtalHistMapName;
    std::string ttHistMapName;
    TFile* caliRootFile;

	struct DetIDStruct {

		DetIDStruct() {}
		DetIDStruct( const int ni1, const int ni2, const int nTT, const int & necal, const float nphi, const float neta )
			: i1(ni1), i2(ni2), TT(nTT), ecal(necal), phi(nphi), eta(neta) {}

		int i1; // EB: iphi, EE: ix
		int i2; // EB: ieta, EE: iy
		int TT; // trigger tower
		int ecal; // EB, EM, EP
		float phi; // xtal phi
		float eta; // xtal eta

	};//<<>>struct DetIDStruct

	struct caliMapStruct {

		caliMapStruct(){}
		caliMapStruct( TH2F* thistfile, std::string histname, bool isnew, bool isopen, int last ) 
			: histFile(thistfile), histName(histname), isNew(isnew), isOpen(isopen), lastRun(last) {}

		TH2F* histFile;
		std::string histName;
		bool isNew; // has not been saved to TFile yet
		bool isOpen; // is open for use, false -> read to use, true -> not all runs processed
		int lastRun;

	};//<<>>struct caliMapStruct

	struct CaliRunStruct {

		CaliRunStruct() {}
		CaliRunStruct( std::string tmpxtalmap, int tstart, int tend, float tlumi, int last ) 
			: histMaps(tmpxtalmap), startRun(tstart), endRun(tend), lumi(tlumi), lastRun(last) {}	

		int startRun;
		int endRun;
		std::string histMapName;
		int lastRun;
		float lumi;

	};//<<>>TimeCaliTagStruct

    struct lumiRunStruct {

		lumiRunStruct() {}
		lumiRunStruct( int trun, int tfill, int tlumi )
			: run(trun), fill(tfill), lumi(tlumi) {}

		int run;
		int fill;
		int lumi

	};//<<>>lumiRunStruct

	struct smearParameters {

		float noise;
		float stochastic;
		float constant;

	};//smearParameters

	struct smearRunStruct {

		smearRunMap() {};
		smearRunMap( std::string source, std::string target, float noise, float stoc, float cons, int strt, int end )
			: sourceTag(source), targetTag(target), smear.noise(noise), smear.stochastic(stoc), 
			  smear.constant(cons), startRun(strt), endRun(end) {}
        smearRunMap( std::string source, std::string target, smearParameters spars, int strt, int end )
			: sourceTag(source), targetTag(target), smear(spars), startRun(strt), endRun(end) {}

		smearParameters smear;
		std::string sourceTag;
		std::string targetTag;
		int startRun;
		int endRun;

	};//<<>>smearRunMap

    std::string curTTIov;
    std::string curXIov;
	std::string curTag; // used to change rereco version or campaian calibrations are for

    std::map<UInt_t,DetIDStruct> DetIDMap; // map of information by detid for each crystal

	std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet; // str is tag, int is internal iov/era number  - final cali map
    std::map<std::string,std::map<int,CaliRunStruct>> TTCaliRunMapSet; // str is tag, int is internal iov/era number  - final cali map
	std::map<std::string,caliMapStruct> CaliMaps; // str is label for specific xtalcalimap - stored as ref in calirunmap
    std::map<std::string,caliMapStruct> TTCaliMaps; // intenal for use to create calirunmaps,  int is run, str is tag

    std::map<std::string,std::map<int,lumiRunStruct>> lumiRunMaps;// str is settag, int is run - ref for maps - give lumi by run
    std::map<std::string,smearRunStruct> smearRunMap; // map of smear paramerts

    std::map<std::string,std::map<int,int>> iovMaps; // < tag, < start run, end rn >>

	public:

	void SetupDetIDsEB();
	void SetupDetIDsEE();
	void SetupIovMaps();

    void SetupTTIovMap( std::string tag );

	void ReadCaliRunFile();
    void SaveCaliRunFile();
    void ReadTTRunFile();
    void SaveTTRunFile();
	void ReadSmearFile();
    void SaveSmearFile();
	void ReadLumiFile( std::string lumifile, std::string tag );// only need to create caliRunMap entries - made with brilcalc

    void SetupCaliMaps();// loads up all starting info
    void SaveCaliMaps();// save existing calimaps to root TFile

	// use to create calibration files and add to DB
	// work in progress - still thinking this area trhough
	void makeCaliMaps();
	void makeSmearMaps();

	// accessors to utilize calibration information
	float getCalibration( uInt rhid, int run, std::string tag  ); // tag indicates which calibration set to use
    float getCalibration( uInt rhid, int run )
			{ getCalibration( rhid, run, curTag ); };
	float getSmearedTime( float rhtime, float rhamp, uInt rhid, int run, std::string stag );//tag indicates which smear to use
    float getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int run, std::string ctag, std::string stag )
    float getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int run, std::string stag )
			{ getSmrdCalibTime( rhtime, rhamp, rhid, run, curTag, stag ); }; 
	DetIDStruct& getDetIdInfo( uInt rhid );

	void setTag( std::string tag ){ curTag = tag; }; // used to change rereco version or campaian calibrations are for
    void setXIov( std::string tag ){ curXIov = tag; };
    void setTTIov( std::string tag ){ curTTIov = tag; };
	// function to create calibration and smear information 
	void makeTTCaliMapEGR( std::string inputFileName ); // make TT maps from egamma res ntuple from lpc/jwk space

};//<<>>class KUCMSTimeCalibration : KUCMSRootHelperBaseClass

//////////////////////////////////////////////////////////////////////////////////////////
//  Class Object code
///////////////////////////////////////////////////////////////////////////////////////////

KUCMSTimeCalibration::KUCMSTimeCalibration(){

	// parameter and setup hardcoded in constructer - some of this needs moved to "run" code

    caliFileDir = "ecal_config/";// move
	detIDConfigEB = "ecal_config/fullinfo_v2_detids_EB.txt";
    detIDConfigEE = "ecal_config/fullinfo_v2_detids_EE.txt";
	caliRunConfig = "ecal_config/caliRunConfig.txt";
    caliSmearConfig = "ecal_config/caliSmearConfing.txt";
    caliTTConfig = "ecal_config/caliTTConfig.txt";

    xtalHistMapName = "AveTTRecTimeMap"
    ttHistMapName = "AveXtalRecTimeMap"

    caliTFile = "ecal_config/caliHistsTFile.root";// name configurable ?
	caliRootFile = TFile::Open( caliTFile.c_str(), "update" );

    curTag = "default";
	curTTIov = "default";
    curXIov = "default";

    SetupDetIDsEB();
    SetupDetIDsEE();
    SetupIovMaps(); // this info hardcoded from DB

	///////////////////////////////////////////////////////////////////
	// this section should be moved to run file -----------------------
	// no need to "recreate" cali maps every time, they will be saved -
	///////////////////////////////////////////////////////////////////

	// for mc use "mc" for TTIov and XIov

    //-----//////////  making tt cali file :
    std::string r2ulTag( "r2ul" );
    setTTIov( r2ulTag );
	std::string inputfilename( "someEGRfile.txt" );
	makeTTCaliMapEGR( inputfilename );
	//-----//////////  making xtal cali file :
	std::string xiovtag( "prompt" );
    setXIov( xiovtag );
	// same as above std::string inputfilename( "someEGRfile.txt" );
	makeXCaliMapEGR( inputfilename );

	// make smear maps

	/////////////////////////////////////////////////////////////////////
	//-------------------------------------------------------------------

	//ReadCaliRunFile();
	//ReadTTRunFile();
	//ReadSmearFile();

	//SetupCaliMaps();

};//<<>>KUCMSTimeCalibration()   

KUCMSTimeCalibration::~KUCMSTimeCalibration(){

	//SaveCaliRunFile();
	SaveTTRunFile();
	//SaveSmearFile();

	caliRootFile.Close();

}//<<>>KUCMSTimeCalibration::~KUCMSTimeCalibration()

void KUCMSTimeCalibration::SetupDetIDsEB(){

    std::ifstream infile( detIDConfigEB, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    float phi, eta;
    std::string pos;

    while( infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos
                >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM >> phi >> eta ){

        DetIDMap[cmsswId] = {iphi,ieta,TT25,0,phi,eta};

    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )

void KUCMSTimeCalibration::SetupDetIDsEE(){

    std::ifstream infile( detIDConfigEE, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    float phi, eta;
    std::string EE;

    while( infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC
                >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant >> phi >> eta ){

        int ec = 1;
        if( side > 0 ) ec = 2;
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec,phi,eta};

    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )

void KUCMSTimeCalibration::SetupIovMaps(){

	//Prompt 
	std::map<int,int> promptIovMap;
	promptIovMap[1] = 189158;//03/02/2012
	promptIovMap[189159] = 203827;//26/03/2012
    promptIovMap[203828] = 208838;//28/09/2012
    promptIovMap[208839] = 253983;//11/12/2012
    promptIovMap[253984] = 276811;//11/08/2015
    promptIovMap[276812] = 293998;//14/07/2016
    promptIovMap[293999] = 296398;
    promptIovMap[296399] = 297726;
    promptIovMap[297727] = 300201;
    promptIovMap[300202] = 301486;
    promptIovMap[301487] = 304475;
    promptIovMap[304476] = 307554;//05/10/2017
    promptIovMap[307555] = 314766;
    promptIovMap[314767] = 315343;
    promptIovMap[315344] = 316360;
    promptIovMap[316361] = 316569;
    promptIovMap[316570] = 318744;
    promptIovMap[318745] = 319329;
    promptIovMap[319330] = 321506;
    promptIovMap[321507] = 322717;
    promptIovMap[322718] = 323412;
    promptIovMap[323413] = 324305;
    promptIovMap[324306] = 327238;
    promptIovMap[327239] = 356513;//25/11/2018
    promptIovMap[356514] = 357289;
    promptIovMap[357290] = 358883;
    promptIovMap[358884] = 359420;
    promptIovMap[359421] = 360089;
    promptIovMap[360090] = 360981;
    promptIovMap[360982] = 361416;
    promptIovMap[361417] = 362522;
    promptIovMap[362523] = 367094;//24/11/2022
    promptIovMap[367095] = 367515;
    promptIovMap[367516] = 367882;
    promptIovMap[367884] = 368825;
    promptIovMap[368826] = 369801;
    promptIovMap[369802] = 369912;
    promptIovMap[369913] = 370496;
    promptIovMap[370497] = 372308;
    promptIovMap[372309] = 372750;//"29/08/2023";//cc start
    promptIovMap[372751] = 373539;
    promptIovMap[373540] = 374564;
    promptIovMap[374565] = 375351;
    promptIovMap[375352] = 378673;//"25/10/2023";
    promptIovMap[378674] = 378793;
    promptIovMap[378794] = 379190;
    promptIovMap[379191] = 379323;
    promptIovMap[379324] = 379366;
    promptIovMap[379367] = 379433;
    promptIovMap[379434] = 380065;
    promptIovMap[380066] = 381383;
    promptIovMap[381384] = 382007;
    promptIovMap[382008] = 382208;
    promptIovMap[382209] = 382959;
    promptIovMap[382960] = 385285;
    promptIovMap[385286] = 999999;//"06/09/2024";
	iovMaps["prompt"] = promptIovMap;

	std::map<int,int> mcIovMap;
	mvIovMap[0] = 999999;
	iovMaps["mc"] = mcIovMap;

    float minttlumi = 10000000;// in /ub
    std::string r2ulTag( "r2ul" );
    std::string lumiUL16Config = "UL2016_runlumi.txt";
    std::string lumiUL17Config = "UL2017_runlumi.txt";
    std::string lumiUL18Config = "UL2018_runlumi.txt";

    ReadLumiFile( caliFileDir+lumi16Config, r2ulTag );
    ReadLumiFile( caliFileDir+lumi17Config, r2ulTag );
    ReadLumiFile( caliFileDir+lumi18Config, r2ulTag );
    SetupTTIovMap( r2ulTag, minttlumi );

}//<<>>void KUCMSTimeCalibration::SetupEraIovMap()

void KUCMSTimeCalibration::SetupTTIovMap( std::string tag, float minTTLumi ){

	bool newrange( true );
	int start = 0;
	int end = 0;
	float lumisum = 0;
	std::map<int,int> ttIovMap;
	for( auto& lumirun : lumiRunMaps[tag] ){
	
		if( newrange ){ start = lumirun.second.run; newrange = false; }
		lumisum += lumirun.second.lumi;
		if( lumisum > minTTLumi ){ 
			end = lumirun.second.run - 1;
			if( end <= start ) end = start;
			ttIovMap[start] = end; 
			newrange = true;
			lumisum = 0;
		}//<<>>if( lumisum > minTTLumi )

	}//<<>>for( auto& lumirun : lumiRunMap )
	if( not newrange ) ttIovMap[start] = 999999;
	iovMaps[tag] = ttIovMap;		

}//<<>>void KUCMSTimeCalibration::makeTTIovMap()

void KUCMSTimeCalibration::ReadLumiFile( std::string lumifile, std::string tag ){

    std::ifstream infile( lumifile, std::ios::in);
    std::string infilestr;
	while( std::getline( iflstream, infilestr ) ){

		if( infilestr[0] != '#' ){

			std::stringstream ss(infilestr);
    		std::string runfill, date, time;
    		int run, fill, nls, ncms;
    		float delivered, recorded;

			ss >> runfill >> date >> time >> nls >> ncms >> delivered >> recorded;
			std::stringstream srun( runfill.substr( 0, 5 ) );
			srun >> run;
			std::stringstream sfil( runfill.substr( 7, 10 ) );
			sfil >> fill;		

 			lumiRunMaps[tag][run] = { run, fill, lumi };

		}//<<>>if( infilestr[0] != '#' )

	}//<<>>while(std::getline(iflstream,infilestr))	

}//<<>>void KUCMSTimeCalibration::upLoadLumiFile( std::string lumifile )

void KUCMSTimeCalibration::ReadCaliRunFile(){

    std::ifstream infile( caliRunConfig, std::ios::in);
    std::string caliMapName, tag;
	int srun, erun, lrun;
    float lumi;
	bool first( true );
    while( infile >> tag >> srun >> erun >> caliMapName >> lumi >> lrun ){
		CaliRunMapSet[tag][srun] = { caliMapName, srun, erun, lumi, lrun };
		if( first ){ first = false; curTag = tag; }
    }//<<>>while (infile >>
	infile.close();

}//<<>>void ReadTimeCaliTagFile()

void KUCMSTimeCalibration::SaveCaliRunFile(){

    std::iostream outfile( caliRunConfig, std::ios::out | std::ios::trunc );
	for( auto& calirunmap : CaliRunMapSet ){

		std::string tag = calirunmap.first;
		for( auto& calirunsct : calirunmap.second ){
		
			int srun = calirunsct.second.startRun;
        	int erun = calirunsct.second.endRun;
        	std::string mapName = calirunsct.second.histMapName;
        	float lumi = calirunsct.second.lumi;
			int lrun = calirunsct.second.lastRun;

			outfile << tag << srun << erun << mapName << lumi << lrun;

		}//<<>>for( auto& calirunsct : calirunmap[tag] )

	}//<<>>for( auto& calirunmap : CaliRunMapSet )
	outfile.close();

}//<<>>void ReadTimeCaliTagFile()

void KUCMSTimeCalibration::ReadTTRunFile(){

    std::ifstream infile( caliTTConfig, std::ios::in);
    std::string caliMapName, tag;
    int srun, erun, lrun;
    float lumi;
    while( infile >> tag >> srun >> erun >> caliMapName >> lumi >> lrun ){
        TTCaliRunMapSet[tag][srun] = { caliMapName, srun, erun, lumi, lrun };
    }//<<>>while (infile >>
	infile.close();

}//<<>>void KUCMSTimeCalibration::ReadTTRunFile()

//std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet;

void KUCMSTimeCalibration::SaveTTRunFile(){
    
    std::iostream outfile( caliTTConfig, std::ios::out | std::ios::trunc );
    for( auto& calirunmap : CaliRunMapSet ){
        
        std::string tag = calirunmap.first;
		for( auto& calirunsct : calirunmap.second ){

        	int srun = calirunsct.second.startRun;
        	int erun = calirunsct.second.endRun;
        	std::string mapName = calirunsct.second.histMapName;
        	float lumi = calirunsct.second.lumi;
			int lrun = calirunsct.second.lastRun;        

        	outfile << tag << srun << erun << mapName << lumi << lrun;
    
		}//<<>>for( auto& calirunsct : calirunmap )

    }//<<>>for( auto& calirunmap : CaliRunMapSet )
	outfile.close();

}//<<>>void ReadTimeCaliTagFile()

void KUCMSTimeCalibration::ReadSmearFile(){

    std::ifstream infile( caliSmearConfig, std::ios::in);
	float noise, stochastic, constant;
    std::string sourceTag, targetTag;
    int startRun, endRun;
	while( infile >> sourceTag >> targetTag >> startRun >> endRun >> noise >> stochastic >> constant ){
		smearRunMap[sourceTag] = { sourceTag, targetTag, noise, stochastic, constant, startRun, endRun };
	}//<<>>while( infile >> sourceTag >> t
    infile.close();

}//<<>>void KUCMSTimeCalibration::ReadSmearFile()

void KUCMSTimeCalibration::SaveSmearFile(){

    std::iostream outfile( caliSmearConfig, std::ios::out | std::ios::trunc );
    for( auto& runmap : smearRunMap ){

		std::string smrtag = runmap.first;
        std::string stag = runmap.second.sourceTag;
        std::string ttag = runmap.second.targetTag;
        int startRun = runmap.second.startRun;
        int endRun = runmap.second.endRun;
        float noise = runmap.second.smear.noise;
        float stochastic = runmap.second.smear.stochastic;
		float constant = runmap.second.smear.constant;

        outfile << smrtag << sourceTag << targetTag << startRun << endRun << noise << stochastic << constant;

    }//<<>>for( auto& calirunmap : CaliRunMapSet )
    outfile.close();

}//<<>>void KUCMSTimeCalibration::SaveSmearFile(){

//std::string ttfilename = ttHistMapName+name+std::to_string(run);
//		AveTTRecTimeMap + pd/camp/tag name + start run
//		AveXtalRecTimeMap + pd/camp/tag name + start run

void KUCMSTimeCalibration::SetupCaliMaps(){

	caliRootFile->cd();
	for( auto& calirunmap : TTCaliRunMapSet ){
		for( auto& calirunsct : calirunmap.second ){
			std::string ttfilename( calirunsct.second.histMapName );
			int erun( calirunsct.second.endRun );
			int lrun( calirunsct.second.lastRun );
			bool isopen( lrun < erun );
			if( ttfilename != "none" && TTCaliMaps.find(ttfilename) == TTCaliMaps.end() ){ 
				TH2F* hist = (TH2F*)caliFile->Get(ttfilename);
				TTCaliMaps[ttfilename] = { hist, ttfilename, false, isopen, lrun }; // histfile histname isnew isopen lastrun
			}//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )
		}//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

    for( auto& calirunmap : CaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
        	std::string xtfilename( calirunsct.second.histMapName );
            int erun( calirunsct.second.endRun );
            int lrun( calirunsct.second.lastRun );
            bool isopen( lrun < erun );
        	if( xtfilename != "none" && XtalCaliMaps.find(xtfilename) == XtalCaliMaps.end() ){
            	TH2F* hist = (TH2F*)caliFile->Get(xtfilename);
            	XtalCaliMaps[xtfilename] = { hist, xtfilename, false, isopen, lrun };
        	}//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )
		}//<<>>for( auto& calirunsct : calirunmap )
	}//<<>>for( auto& calirunmap : CaliRunMapSet )

}//<<>>void KUCMSTimeCalibration::SetupCaliMaps()

void KUCMSTimeCalibration::SaveCaliMaps(){

	caliRootFile->cd();
	for( auto& calimapsct : TTCaliMaps ){

		if( calimapsct.second.isNew || calimapsct.second.isOpen ) calimapsct.second.histFile->Write();	
		calimapsct.second.histFile->Close();

	}//<<>>for( auto& calimapsct : TTCaliMaps )
    for( auto& calimapsct : XtalCaliMaps ){

        if( calimapsct.second.isNew || calimapsct.second.isOpen  ) calimapsct.second.histFile->Write();
        calimapsct.second.histFile->Close();

    }//<<>>for( auto& calimapsct : XtalCaliMaps )

}///<<>>void KUCMSTimeCalibration::SaveCaliMaps()

DetIDStruct& KUCMSTimeCalibration::getDetIdInfo( uInt rhid ){

	return DetIDMap[rhid];

}//<<>>DetIDStruct KUCMSTimeCalibration::getDetIdInfo( uInt rhid )

//std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet;

float KUCMSTimeCalibration::getCalibration( uInt rhid, int run, std::string tag ){

	auto& idinfo = DetIDMap[rhid];
	int iEta = fill_idinfo.i2;
    int iPhi = fill_idinfo.i1;
    bool isEB = (fill_idinfo.ecal == ECAL::EB);
	float xtaltime = 0.f;	
	//if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
	if( not isEB ){ std::cout << "Calibration for EE is not supported." << std::endl; return 0.f; }
    for( auto& calirunmap : CaliRunMapSet[tag] ){	
		if( run >= calirunmap.second.startRun && run <= calirunmap.second.endRun ){
			XMapName = calirunmap.second.XtalMap;
			xtaltime = XtalCaliMaps[XMapName].hist->GetBinContent( iEta + 86, iPhi );
		}//<<>>if( run >= calirunmap.second.startRun
	}//<<>>for( auto& calirunmap : CaliRunMapSet )
    return xtaltime

}//<<>>float KUCMSTimeCalibration::getCalibration( std::string tag )

//std::map<std::string,smearRunStruct> smearRunMap; // map of smear paramerts
float KUCMSTimeCalibration::getSmearedTime( float rhtime, float rhamp, uInt rhid, int run, std::string stag ){

    double noise = smearRunMap[stag].smear_n;
    double stachastic = smearRunMap[stag].smear_s;
    double constant = smearRunMap[stag].smear_c;
    double resolution = std::sqrt( ( pow(noise/rhamp,2) + POW(stachastic,2)/rhamp  + 2*pow(constant,2) )/2 );
    if( resolution <= 0 ){ std::cout << "No smearing values set for this tag !!!" << std::endl; return time; }
    float smearedtime = getRandom->Gaus( rhtime, resolution );
    return smearedtime;

}//<<>>float KUCMSTimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

float KUCMSTimeCalibration::getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int run, std::string ctag, std::string stag ){

    float calibration = getCalibration( rhid, run, ctag );
    float smrdCalibTime = getSmearedTime( rhtime+calibration, rhamp, rhid, stag );
    return smrdCalibTime;

}//<<>>float KUCMSTimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

// assume full runs are present in input files
// input file will give the name of a text file "list" of root files, the path to the root files, and the tag to use for those files 
// check if run range started
// check if run range complete 
//     std::map<std::string,std::map<int,int>> iovMaps; // < tag, < start run, end rn >>
//	- if started but not complete - is run needed ( lastrun )
//	- if not started - start new cali map
void KUCMSTimeCalibration::makeTTCaliMapEGR( std::string inputFileName ){

	makeCaliMapsEGR( inputFileName, true );

}//<<>>void KUCMSTimeCalibration::makeTTCaliMapEGR( std::string inputFileName )

void KUCMSTimeCalibration::makeXCaliMapEGR( std::string inputFileName ){

    makeCaliMapsEGR( inputFileName, false );

}//<<>>void KUCMSTimeCalibration::makeXCaliMapEGR( std::string inputFileName )


void KUCMSTimeCalibration::makeCaliMapsEGR( std::string inputFileName, bool doTT ){

// make generic fill with specialized calling functions for input and TT v X
// assume bool doTT to inicate TT of X Cali set

  	const double offset = 0.0;
    const int bin_offset = 86;
    const float minRhEnergy = 5.0;   

    // Declaration of leaf types
    UInt_t run;
    std::vector<uInt> *rhID = 0;
    std::vector<float> *rhRtTime = 0;
    std::vector<float> *rhEnergy = 0;
	std::vector<float> *rhAmp = 0;

    // List of branches
    TBranch *b_run;   //!
    TBranch *b_rhID;   //!
    TBranch *b_rhRtTime;   //!
    TBranch *b_rhEnergy;   //!
    TBranch *b_rhAmp;

	if( doTT && iovMaps.find(curTTIov) == iovMaps.end() ){ std::cout << " TT Iov not valid !!" << std::endl; return; }
    if( not doTT && iovMaps.find(curXIov) == iovMaps.end() ){ std::cout << " X Iov not valid !!" << std::endl; return; }

    std::ifstream infilelist(inputFileName);
    std::string infilestr;
    while (std::getline(infilelist,infilestr)){

        std::stringstream ss(infilestr);
        std::string infilename, tag;
        int srun, erun;
		
        ss >> infilename >> srunstr >> erunstr >> type >> tag;
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "For Run " << srunstr << " to Run " << erunstr << std::endl;
        std::cout << "Producing maps for " << tag << std::endl;

		if( not doTT && TTCaliRunMapSet.find(tag) == iovMaps.end() ){ std::cout << " No TT maps for this tag !!" << std::endl; return; }

        std::ifstream infile(infilename);
        std::string instr;
        auto fInTree = new TChain(treename.c_str());
        std::cout << "Adding files to TChain." << std::endl;
        const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");// nput parameter!
        //const std::string eosdir("root://cmseos.fnal.gov//store/user/");  
		const std::string indir("/ecalTiming/gammares_llpana_pd/");// input paramter !  
        while (std::getline(infile,instr)){
            auto tfilename = eosdir + indir + instr;
            std::cout << "-";
            fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
        std::cout << std::endl;

       	run = 0;
        rhID = 0;
        rhRtTime = 0;
        rhEnergy = 0;

        fInTree->SetBranchAddress("run", &run, &b_run);
        fInTree->SetBranchAddress("rhID", &rhID, &b_rhID);
        fInTree->SetBranchAddress("rhRtTime", &rhRtTime, &b_rhRtTime);
        fInTree->SetBranchAddress("rhEnergy", &rhEnergy, &b_rhEnergy);

		int prevrun = 999999;
        auto nEntries = fInTree->GetEntries();
        std::cout << "Starting entry loops for " << nEntries << " Entries" <<  std::endl;
        if( debug ) nEntries = 100;
        for (Long64_t centry = 0; centry < nEntries; centry++){

            if( centry%1000000 == 0 or centry == 0){
                std::cout << "Proccessed " << centry << " of " << nEntries;
                std::cout << " (" << static_cast<float>((10000*centry)/nEntries)/(100) << "%)" << std::endl;
            }//<<>>if( centry%1000000 == 0 or centry == 0)

            auto entry = fInTree->LoadTree(centry);

            b_run->GetEntry(entry);   //!
            b_rhID->GetEntry(entry);   //!
            b_rhRtTime->GetEntry(entry);   //!
            b_rhEnergy->GetEntry(entry);

            if( run < srun || run > erun ) continue;

        	// check if cali map for this run exists
        	// if exists check to see if complete  
        	// 		if complete continue
        	// 		if not complete - fill and check for completion
        	// if not exist find iov range for new map

//std::map<std::string,std::map<int,CaliRunStruct>> TTCaliRunMapSet; // str is tag, int is internal iov/era number  - final cali map
//    struct CaliRunStruct {
//        int startRun;
//        int endRun;
//        std::string histMapName;
//        int lastRun;
//        float lumi;
//    };//<<>>TimeCaliTagStruct

			auto& calirunset = doTT ? TTCaliRunMapSet : CaliRunMapSet;
			if( calirunset.find(tag) != calirunset.end() ){ // tag exists
				if( calirunset[tag].find(run) != calirunset[tag].end() ){ // run exists
					auto& therun = calirunset[tag][run];
					if( thecai.endRun == thecali.lastRun ) continue; // do nothing - this map is complete
					if( run > thecali.lastRun ){
						// make info to do fill map and update lastRun of entry -------------------------------------
					}//<<>>if( run > thecali.lastRun )					
				} else { // run does not exist
					// make info to do fill map and make new run entry --------------------------------------------------
				}//<<>> if( TTCaliRunMapSet[tag].find(run) != TTCaliRunMapSet[tag].end() )
			} else { //tag does not exit
				// make info to do fill map and creat new tag + run entry
			}//<<>>if( TTCaliRunMapSet.find(tag) != TTCaliRunMapSet.end() )






			// finding iov range for new run entry
			auto& iovset = doTT ? iovMaps[curTTIov] : iovMaps[curXIov];
        	for( auto& iovmap : iovset ){

            	if( run == iovmap.first ){}// first run in range - start new map  endrun = iovmap.second

            	if( run >= iovmap.first && run < iovmap.second ){}// not first run in range no current partial mapping ?

        	}//<<>>for( auto& iovmap : iovMaps[tag] )




		}//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)

	}//<<>>while (std::getline(infilelist,infilestr))

	// fill map hist -----------------------------------------------


}//<<>>void KUCMSTimeCalibration::makeTTCaliMap( std::string inputFileName )

#endif
//-------------------------------------------------------------------------------------------------------------------

