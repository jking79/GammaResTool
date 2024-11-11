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

    const std::string detIDConfigEB;
    const std::string detIDConfigEE;
    const std::string caliFileDir;
    const std::string caliRunConfig;


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
		caliMapStruct( TH2F* thistfile, TFile* ttreefile ) 
			: histFile(thistfile), treeFile(ttreefile) {}

		TH2F* histFile;
		TFile* treeFile;

	};//<<>>struct caliMapStruct

	struct CaliRunStruct {

		TimeCaliTagStruct() {}
		TimeCaliTagStruct( std::string tmpxtalmap, std::string tmpttmap, float tn, float ts, float tc, float tlumi ) 
			: XtalMaps(tmpxtalmap), TTMaps(tmpttmap), smrn(tn), smrs(ts), smrc(tc), lumi(tlumi) {}	

		std::string XtalMap;
		std::string TTMap;
		float smearn;
		float smears;
		float seamrc;
		float lumi;

	};//<<>>TimeCaliTagStruct

    std::map<UInt_t,DetIDStruct> DetIDMap;
	std::map<int,CaliRunStruct> CaliRunMap;
	std::map<std::string,caliMapStruct> XtalCaliMaps;
    std::map<std::string,caliMapStruct> TTCaliMaps;

	public:

	void SetupDetIDsEB();
	void SetupDetIDsEE();

	void ReadCaliRunFile();
	void SaveCaliRunFile();
	void SetupCaliMaps();
	//void ReadLumiFile();
	//void SaveLumifile();

	void upLoadCaliFile( std::string infilelist );
	void upLoadLumiMap( std::string lumifile );

	float getCalibration( uInt rhid, int run );
    float getSmearedTime( float rhtime, float rhamplitude, uInt rhid, int run );
    float getSmearedCalibratedTime( float rhtime, float rhamplitude, uInt rhid, int run );
	DetIDStruct& getDetIdInfo( uInt rhid );

};//<<>>class KUCMSTimeCalibration : KUCMSRootHelperBaseClass

//////////////////////////////////////////////////////////////////////////////////////////
//  Class Object code
///////////////////////////////////////////////////////////////////////////////////////////

KUCMSTimeCalibration::KUCMSTimeCalibration(){


	detIDConfigEB = "ecal_config/fullinfo_v2_detids_EB.txt";
    detIDConfigEE = "ecal_config/fullinfo_v2_detids_EE.txt";
	caliRunConfig = "ecal_config/calibrationRunConfig.txt";
	caliFileDir = "cali_root_files/";

    SetupDetIDsEB();
	SetupDetIDsEE();
	ReadCaliRunFile();
	SetupXtalMaps();
	SetupTTMaps();

};//<<>>KUCMSTimeCalibration()   

KUCMSTimeCalibration::~KUCMSTimeCalibration(){

	for( auto& calimap : XtalMaps ){

		if( calimap.second.histFile ) calimap.second.histFile.Close();
        if( calimap.second.treeFile ) calimap.second.treeFile.Close();		

	}//<<>>for( auto calitag : TimeCaliTagMap )

    for( auto& calimap : TTMaps ){

        if( calimap.second.histFile ) calimap.second.histFile.Close();
        if( calimap.second.treeFile ) calimap.second.treeFile.Close();

    }//<<>>for( auto calitag : TimeCaliTagMap )

	SaveCaliRunFile();

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

void KUCMSTimeCalibration::upLoadLumiFile( std::string lumifile ){

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

			if( CaliRunMap.find(run) == CaliRunMap.end() ){
				
				CaliRunStruct newrun( NULL, NULL, "none", "none", 0, 0, 0, lumi );
 				CaliRunMap[run] = newrun;

			} else CaliRunMap[run].lumi = recorded;
		
		}//<<>>if( infilestr[0] != '#' )

	}//<<>>while(std::getline(iflstream,infilestr))	


}//<<>>void KUCMSTimeCalibration::upLoadLumiFile( std::string lumifile )

void KUCMSTimeCalibration::ReadCaliRunFile(){

    std::ifstream infile( timeTagConfig, std::ios::in);
    std::string TTCaliMapName, XtalCaliMapName;
	int run;
    float n, s, c, lumi;
    while( infile >> run >> TTCaliMapName >> XtalCaliMapName >> n >> s >> c >> lumi ){
		CaliRunMap[run] = { TTCaliMapName, XtalCaliMapName, n, s, c, lumi };
    }//<<>>while (infile >>

}//<<>>void ReadTimeCaliTagFile()

void KUCMSTimeCalibration::SetupCaliMaps(){

	for( auto& calirunmap : CaliRunMap ){

		std::string ttfilename( calirunmap.second.TTCaliMapName );
		if( ttfilename != "none" && TTMaps.find(ttfilename) == TTMaps.end() ){ 
			TFile* caliFile = TFile::Open( ttfilename.c_str(), "read" );
			TH2F* hist = (TH2F*)caliFile->Get("AveTTRecTimeMap");
			caliMapStruct newmap( hist, caliFile );
			TTMaps[ttfilename] = newmap;
		}//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )

        std::string xtfilename( calirunmap.second.XtalCaliMapName );
        if( xtfilename != "none" && TTMaps.find(xtfilename) == XtalMaps.end() ){
            TFile* caliFile = TFile::Open( xtfilename.c_str(), "read" );
            TH2F* hist = (TH2F*)caliFile->Get("AveXtalRecTimeMap");
            caliMapStruct newmap( hist, caliFile );
            XtalMaps[xtfilename] = newmap;
        }//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )

	}//<<>>for( auto& calirunmap : CaliRunMap )

}//<<>>void KUCMSTimeCalibration::SetupCaliMaps()

DetIDStruct& KUCMSTimeCalibration::getDetIdInfo( uInt rhid ){

	return DetIDMap[rhid];

}//<<>>DetIDStruct KUCMSTimeCalibration::getDetIdInfo( uInt rhid )

float KUCMSTimeCalibration::getCalibration( uInt rhid ){

	auto& idinfo = DetIDMap[rhid];
	int iEta = fill_idinfo.i2;
    int iPhi = fill_idinfo.i1;
    bool isEB = (fill_idinfo.ecal == ECAL::EB);
	
	if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
	if( not isEB ){ std::cout << "Calibration for EE is not supported." << std::endl; return 0.f; }

	int ttphi = (iPhi-1)/5;
    int ttabseta = (std::abs(iEta)-1)/5;
    int tteta = ( iEta < 0 ) ? ttabseta : ttabseta + 17;
    float xtaltime = CurrentTag.XtalCaliMap->GetBinContent( iEta + 86, iPhi );
    float tttime = CurrentTag.TTCaliMap->GetBinContent( tteta, ttphi );

	return xtaltime + tttime;

}//<<>>float KUCMSTimeCalibration::getCalibration( std::string tag )

float KUCMSTimeCalibration::getSmearedCalibratedTime( float rhtime, float rhamplitude, uInt rhid ){

    if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
	float calibration = getCalibration( rhid );
	float smearedCalibratedTime = getSmearedTime( rhtime+calibration, rhamplitude, rhid );	
	return smearedCalibratedTime;

}//<<>>float KUCMSTimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

float KUCMSTimeCalibration::getSmearedTime( float rhtime, float rhamplitude, uInt rhid ){

    if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
    double noise = CurrentTag.smear_n;
    double stachastic = CurrentTag.smear_s;
    double constant = CurrentTag.smear_c;
    double resolution = std::sqrt( ( pow(noise/rhamplitude,2) + POW(stachastic,2)/rhamplitude  + 2*pow(constant,2) )/2 );
    if( resolution <= 0 ){ std::cout << "No smearing values set for this tag !!!" << std::endl; return time; }
    float smearedtime = getRandom->Gaus( rhtime, resolution );
    return smearedtime;

}//<<>>float KUCMSTimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

#endif
//-------------------------------------------------------------------------------------------------------------------

