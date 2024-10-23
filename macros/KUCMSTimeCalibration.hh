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

	std::string CurrentTag;
    const std::string detIDConfigEB;
    const std::string detIDConfigEE;
    const std::string timeTagConfig;

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

	struct TimeCaliTagStruct {

		TimeCaliTagStruct() {}
		TimeCaliTagStruct( const TH2F* xtal, const TH2F* tt, const float n, const float s, const float c, const TFile* xf, const TFile* ttf )
			: XtalCaliMap(xtal), TTCaliMap(tt), smear_n(n), smear_s(s), seamr_c(c), XtalCaliMapFile(xf), TTCaliMapFile(ttf) {}	

        TFile* XtalCaliMapFile;
		TH2F* XtalCaliMap;
        TFile* TTCaliMapFile;
		TH2F* TTCaliMap;
		float smear_n;
		float smear_s;
		float seamr_c;

	};//<<>>TimeCaliTagStruct

	bool validCurrentTag;
	TimeCaliTagStruct CurrentTag;
    std::map<UInt_t,DetIDStruct> DetIDMap;
	std::map<std::string,TimeCaliTagStruct> TimeCaliTagMap;

	public:

	void SetupDetIDsEB();
	void SetupDetIDsEE();
	void ReadTimeCaliTagFile();

	void setTag( std::string tag );
	float getCalibration( uInt rhid );
	float getSmearedTime( float rhtime, float rhamplitude, uInt rhid );
    float getSmearedCalibratedTime( float rhtime, float rhamplitude, uInt rhid );
	DetIDStruct& getDetIdInfo( uInt rhid );

};//<<>>class KUCMSTimeCalibration : KUCMSRootHelperBaseClass

//////////////////////////////////////////////////////////////////////////////////////////
//  Class Object code
///////////////////////////////////////////////////////////////////////////////////////////

KUCMSTimeCalibration::KUCMSTimeCalibration(){


	CurrentTag = { NULL, NULL, 0.f, 0.f, 0.f, NULL, NULL };
	detIDConfigEB = "ecal_config/fullinfo_v2_detids_EB.txt";
    detIDConfigEE = "ecal_config/fullinfo_v2_detids_EE.txt";
	timeTagConfig = "ecal_config/timeCalibrationTagConfig.txt";
	caliFileDir = "cali_root_files/";
	validCurrentTag = false;

    SetupDetIDsEB();
	SetupDetIDsEE();
	ReadTagFile();


};//<<>>KUCMSTimeCalibration()   

KUCMSTimeCalibration::~KUCMSTimeCalibration(){

	for( auto& calitag : TimeCaliTagMap ){

		if( calitag.second.XtalCaliMapFile ) calitag.second.XtalCaliMapFile.Close();
        if( calitag.second.TTCaliMapFile ) calitag.second.TTCaliMapFile.Close();		

	}//<<>>for( auto calitag : TimeCaliTagMap )

}//<<>>KUCMSTimeCalibration::~KUCMSTimeCalibration()

void KUCMSTimeCalibration::setTag( std::string tag ){

	std::map<std::string,TimeCaliTagStruct>::iterator taginfo = TimeCaliTagMap.find(tag);
	if( taginfo != TimeCaliTagMap.end() ){ CurrentTag = taginfo->second; validCurrentTag = true; }
	else std::cout << "No such tag in calibration tag map." << std::endl; 

}//<<>>void KUCMSTimeCalibration::setTag( std::string tag )

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

void KUCMSTimeCalibration::ReadTimeCaliTagFile(){

    std::ifstream infile( timeTagConfig, std::ios::in);
    std::string tag, TTCaliMapName, XtalCaliMapName;
    float n, s, c;

    while( infile >> tag >> TTCaliMapName >> XtalCaliMapName >> n >> s >> c ){

		TFile* xtalCaliFile(NULL);
        TH2F* XtalCaliMap(NULL);
        TFile* ttCaliFile(NULL);
        TH2F* TTCaliMap(NULL);
		if( TTCaliMapName != "none" ){  
			ttCaliFile = TFile::Open( TTCaliMapName.c_str(), "read" ); }
			TTCaliMap = (TH2F*)ttCaliFile->Get("AveTTRecTimeMap");
		}//<<>>if( TTCaliMapName != "none" )
        if( XtalCaliMapName != "none" ){  
			xtalCaliFile = TFile::Open( TTCaliMapName.c_str(), "read" ); }
			XtalCaliMap = (TH2F*)xtalCaliFile->Get("AveXtalRecTimeEBMap");
		}//<<>>if( TTCaliMapName != "none" )
		TimeCaliTagMap[tag] = { TTCaliMap, XtalCaliMap, n, s, c, xtalCaliFile, ttCaliFile };

    }//<<>>while (infile >>

}//<<>>void ReadTimeCaliTagFile()

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

