// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TText.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFormula.h"
#include "Math/PositionVector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TGraph.h"
#include "TMathBase.h"

// STL includes
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <sys/stat.h>

#include "list_files/egammares_hist_base_v4.hh"

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
//#define nEBEEMaps 36
#define nEBEEMaps 0
#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.14159265358979323846 // pie ... 
#define EVPHR 36000000.0
#define nCrys 4

#define CFlt  const float
#define CDbl  const double
#define CVFlt const std::vector<float>
#define CVDbl const std::vector<double>

#define DEBUG false
//#define DEBUG true  // ? do both debug flags?  debug & DEBUG

typedef unsigned int uInt;

enum ECAL {EB, EM, EP, NONE};

struct DetIDStruct { 
  	DetIDStruct() {}
  	DetIDStruct(const Int_t rhid, const Int_t ni1, const Int_t ni2, float x, float y, float z, int necal) : id(rhid), i1(ni1), i2(ni2), rhx(x), rhy(y), rhz(z), ecal(necal){}
	Int_t id;
  	Int_t i1; // EB: iphi, EE: ix
  	Int_t i2; // EB: ieta, EE: iy
  	float rhx; // trigger tower
    float rhy;
    float rhz;
  	Int_t ecal; // EB, EM, EP
};

struct CaliStruct {
    CaliStruct() {}
    CaliStruct(const Int_t ni1, const Int_t ni2, int necal, float time, int detid) : i1(ni1), i2(ni2), ecal(necal), ctime(time), id(detid) {}
    Int_t i1; // EB: iphi, EE: ix
    Int_t i2; // EB: ieta, EE: iy
    Int_t ecal; // EB, EM, EP
	float ctime;
	Int_t id;
};

// helper functions ------------------------------------------------------------

void SetupDetIDs( std::map<UInt_t,DetIDStruct> &DetIDMap )
{   
    const std::string detIDConfigEB("ecal_config/rhid_info_list.txt");
    std::ifstream infile( detIDConfigEB, std::ios::in);
    
    UInt_t cmsswId; 
    Int_t iphi, ieta, side;
	float rhx, rhy, rhz;
    
	//std::map<UInt_t,DetIDStruct> DetIDMap;
    while( infile >> cmsswId >> iphi >> ieta >> side >> rhx >> rhy >> rhz ){   
        //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << 0 << std::endl;
        DetIDStruct info(cmsswId, iphi, ieta, rhx, rhy, rhz, side);
		uInt loc = (side+2)*1000000 + (ieta+100)*1000 + iphi;
		DetIDMap[loc] = info;
        //auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

// math functions -------------------------

const auto sq2      (CFlt x){return x*x;}
const auto sq2      (CDbl x){return x*x;}
const auto rad2     (CFlt x, CFlt y, CFlt z = 0.f){return x*x+y*y+z*z;}
const auto hypo     (CFlt x, CFlt y, CFlt z = 0.f){return std::sqrt(rad2(x,y,z));}
const auto phi      (CFlt x, CFlt y){return std::atan2(y,x);}
const auto theta    (CFlt r, CFlt z){return std::atan2(r,z);}
const auto eta      (CFlt x, CFlt y, CFlt z){return -1.0f*std::log(std::tan(theta(hypo(x,y),z)/2.f));}
const auto effMean  (CFlt x, CFlt y){return (x*y)/sqrt(x*x+y*y);}
const auto dltIPhi  (CFlt x, CFlt y){auto dp(x-y); if( dp > 180 ){dp-=360.0;} else if( dp < -180 ){ dp+=360.0;} return dp;}
const auto dltPhi   (CFlt x, CFlt y){auto dp(x-y);if(dp>PI) dp-=2*PI; else if(dp<=-PI) dp+=2*PI; return dp;}
const auto dltAngle (CFlt x, CFlt y){auto dp(x-y);if(dp>=2*PI) dp-=2*PI; else if(dp<=-2*PI) dp+=2*PI; return dp;}
const auto max      (CVFlt x){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}

const auto deltaR2  (CDbl e0, CDbl e1, CDbl p0, CDbl p1 ){ auto dp(p1-p0); if(dp>PI) dp-=2*PI; else if(dp<=-PI) dp+=2*PI; return sq2(dp)+sq2(e1-e0);}
const auto deltaR   (CDbl e0, CDbl e1, CDbl p0, CDbl p1 ){ return std::sqrt(deltaR2(e0,e1,p0,p1));}

// histogram functions -------------------------------------------------

std::string addstr( std::string current, std::string input ){ return (current+input); }//<<>>std::string addstr( std::string current, char* input )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////       Class Declaration 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class makehists : egammares_hist_base {

	public:

	void llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename, std::string califilename, 
								int brun, int erun, std::string fhtitle  );	
	void initHists( std::string fHTitle );
	void getBranches( Long64_t entry );
	void eventLoop( Long64_t entry );
 	void endJobs();	
	void initCali( std::string califilename ); 
	void makeEBEEMaps( std::vector<unsigned int> rhcol );
	void makeEBEEMaps( int phoit );
	int getRhIdx( uInt rhDetID );

	std::map<UInt_t,DetIDStruct> DetIDMap;

	TH2F *icmap[6];
	int startRun, endRun;
    float totrhs, totrhs0, totrhs05, totrhs1, totrhs2, totrhs5, totrhs10, encrhs;

    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];

    int nMaps;
    bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
	
	//std::vector<std::vector<float>> cryscctimes, crysrttimes;
    //TH1D *hist1d136[nCrys], *hist1d140[nCrys], *hist1d137[nCrys], *hist1d138[nCrys], *hist1d139[nCrys], *hist1d143[nCrys], *hist1d144[nCrys];
    TH1D *hist1d140[nCrys], *hist1d137[nCrys], *hist1d138[nCrys], *hist1d144[nCrys];

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class Functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int makehists::getRhIdx( uInt rhDetID ){

    for( int idx = 0; idx < rhCaliID->size(); idx++ ){ if( rhDetID == (*rhCaliID)[idx] ) return idx; }
    //std::cout << " -- !! no rhDetID to (*rhCaliID)[idx] match !! ---------------------- " << std::endl;
    return -1;

}//<<>>int makehists::getRhIdx( int rhDetID )
////----------------------------------------  Make Hists function call  ---------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

void makehists::llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename, std::string califilename, 
										int brun, int erun, std::string fhtitle ){

    //bool debug = true; // this is for main loop only
	bool debug = false;

    SetupDetIDs(DetIDMap);
    auto fCaliFile = TFile::Open(califilename.c_str(), "read");
    std::string calistring[2] = { "AveXtalRatioRecTime", "AveXtalKuccRecTime" };
    std::string cmbs = calistring[0];
    std::string ebmapstring(cmbs+"EBMap");
    std::string epmapstring(cmbs+"EPMap");
    std::string emmapstring(cmbs+"EMMap");
    icmap[0] = (TH2F*)fCaliFile->Get(ebmapstring.c_str());
    icmap[1] = (TH2F*)fCaliFile->Get(epmapstring.c_str());
    icmap[2] = (TH2F*)fCaliFile->Get(emmapstring.c_str());

	std::map<UInt_t,CaliStruct>  calimap;
	for( int iphi(0); iphi < 361; iphi++ ){ 
		for( int ieta(-85); ieta < 86; ieta++ ){
			float calitime = icmap[0]->GetBinContent( ieta + 86, iphi ); 
			uInt loc = (2)*1000000 + (ieta+100)*1000 + iphi;
			if( DetIDMap.find(loc) != DetIDMap.end() ){	
			//	std::cout << ieta << " " << iphi << " 0 " << calitime << std::endl;
				auto info = DetIDMap[loc];
				CaliStruct payload(info.i2, info.i1, info.ecal, calitime, info.id);
				int index = ( info.i2 + 100 )*1000 + info.i1;
				calimap[index] = payload;
				//std::cout << info.i2 << " " << info.i1 << " " << ieta << " " << iphi << " " << calitime << std::endl;
			}//<<>>if( DetIDMap.find(loc) != DetIDMap.end() )
		}//<<>>for( int eta(-85); eta < 86; rta++ )
	}//<<>>for( int phi(0); phi < 361; phi++ )
    for( int iphi(0); iphi < 101; iphi++ ){
        for( int ieta(0); ieta < 101; ieta++ ){
            float calitime = icmap[1]->GetBinContent( ieta, iphi );
            uInt loc = (3)*1000000 + (ieta+100)*1000 + iphi;
            if( DetIDMap.find(loc) != DetIDMap.end() ){
			//	std::cout << ieta << " " << iphi << " 1 " << calitime << std::endl;
                auto info = DetIDMap[loc];
                CaliStruct payload(info.i1, info.i2, info.ecal, calitime, info.id);
                int index = info.i2*1000 + info.i1 + 400000;
                calimap[index] = payload;
                //std::cout << info.i1 << " " << info.i2 << " " << iphi << " " << ieta << " " << calitime << std::endl;
            //    std::cout << info.rhx << " " << info.rhy << " " << info.rhz << " " << calitime << std::endl;
            }//<<>>if( DetIDMap.find(loc) != DetIDMap.end() )
        }//<<>>for( int eta(-85); eta < 86; rta++ )
    }//<<>>for( int phi(0); phi < 361; phi++ )
    for( int iphi(0); iphi < 101; iphi++ ){
        for( int ieta(0); ieta < 101; ieta++ ){
            float calitime = icmap[2]->GetBinContent( ieta, iphi );
            uInt loc = (1)*1000000 + (ieta+100)*1000 + iphi;
            if( DetIDMap.find(loc) != DetIDMap.end() ){
			//	std::cout << ieta << " " << iphi << " -1 " << calitime << std::endl;
                auto info = DetIDMap[loc];
                CaliStruct payload(info.i1, info.i2, info.ecal, calitime, info.id);
                int index = info.i2*1000 + info.i1 + 200000;
                calimap[index] = payload;
                //std::cout << info.i1 << " " << info.i2 << " " << iphi << " " << ieta << " " << calitime << std::endl;
            //    std::cout << info.rhx << " " << info.rhy << " " << info.rhz << " " << calitime << std::endl;
            }//<<>>if( DetIDMap.find(loc) != DetIDMap.end() )
        }//<<>>for( int eta(-85); eta < 86; rta++ )
    }//<<>>for( int phi(0); phi < 361; phi++ )

	for( auto & pl : calimap ) { float ct = pl.second.ctime; if( ct != 0  ) ct *= -1; std::cout << pl.second.i1 << " " << pl.second.i2 << " " << pl.second.ecal << " " << ct << " " << pl.second.id << std::endl; }

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

    std::cout << "Thats all Folks!!" << std::endl;
}//<<>>void llpgana_hist_maker

//------------------------------------------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------------------------------------------
//
//// ----------------------------------------------- event loop -------------------------------------------------------------------------

void makehists::eventLoop( Long64_t entry ){

        //------------------------------------------------------------------------------------------------

}//<<>>void makehists::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void makehists::initCali( std::string califilename ){

    int nAlgos(2);
	if( califilename == "none" ) return;
	auto fCaliFile = TFile::Open(califilename.c_str(), "read");
    std::string calistring[nAlgos] = { "AveXtalRatioRecTime", "AveXtalKuccRecTime" };	
	for( int i = 0; i < nAlgos; i++ ){
		std::string cmbs = calistring[i];
    	std::string ebmapstring(cmbs+"EBMap");
    	std::string epmapstring(cmbs+"EPMap");
    	std::string emmapstring(cmbs+"EMMap");
    	icmap[i*3] = (TH2F*)fCaliFile->Get(ebmapstring.c_str());
    	icmap[i*3+1] = (TH2F*)fCaliFile->Get(epmapstring.c_str());
    	icmap[i*3+2] = (TH2F*)fCaliFile->Get(emmapstring.c_str());
	}//<<>>for( int i = 0; i < 2; i++ )
	std::cout << "Getting Cali maps : " << icmap[0] << " && " << icmap[3] << std::endl;

}//<<>>void makehists::initCali( TFile* fCaliFile )

void makehists::getBranches( Long64_t entry ){

}//<<>>void makehists::getBranches( Long64_t entry )

void makehists::endJobs(){


}//<<>>void makehists::endJobs()

void makehists::initHists( std::string fHTitle ){

}//<<>>void makehists::initHists()

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                auto indir = "ecalTiming/EGamma0/";

                //auto infilename = "egammares_gammares_cc_140_v2_EGamma0_MINIAOD_Run2024C-PromptReco-v1_378981-379349_dispho_v2.txt";
                auto infilename = "list_files/egammares_gammares_cc_140_json_v2_EGamma0_MINIAOD_Run2024C-PromptReco-v1_379415-380238_dispho_v2.txt";
                //auto infilename = "list_files/egammares_gammares_cc_140_v2_EGamma1_MINIAOD_Run2024D-PromptReco-v1_380306-380947_dispho_v2.txt";

				int brun = 0;
            	int erun = 999999;
                auto califilename = "cali_root_files/tt_KUCCRes_1404_v12_run3_2024C_eg0_Prompt_379415_380238_Cali.root";

                //auto outfilename = "egammares_diag_2022G_Prompt_359022-362760_v11_IDCT_EB_E10.root";
                auto outfilename = "egammares_diag_24C_EGamma0_v12_gold_iov3_EB.root";

                //auto fhtitle = "Run2023D 13.2.3 EB ";
                auto fhtitle = "Run2024C 14_0_4 EB ";

				makehists base;				
                base.llpgana_hist_maker( indir, infilename, outfilename, califilename, brun, erun, fhtitle );
    //}
    return 1;
}//<<>>int main ( int argc, char *argv[] )



