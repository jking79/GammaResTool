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

using namespace std;

typedef unsigned int uInt;

template <typename T> std::string to_string(T value)
{  
   //create an output string stream
   std::ostringstream os ;
   
   //throw the value into the string stream
   os << value ;
   
   //convert the string stream into a string and return
   return os.str() ;

//you can also do this
//std::string output;
//os >> output;  //throw whats in the string stream into the string
}

enum ECAL {EB, EM, EP, NONE};
std::string ecal_config_path("ecal_config/");

struct DetIDStruct
{
  DetIDStruct() {}
  DetIDStruct(const Int_t inI1, const Int_t inI2, const Int_t inTT, const Int_t & inEcal) : i1(inI1), i2(inI2), TT(inTT), ecal(inEcal)  {}

  Int_t i1; // EB: iphi, EE: ix
  Int_t i2; // EB: ieta, EE: iy
  Int_t TT; // trigger tower
  Int_t ecal; // EB, EM, EP
};

void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )
{
    const std::string detIDConfigEB(ecal_config_path+"fullinfo_detids_EB.txt");
    std::ifstream infile( detIDConfigEB, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    TString pos;

    while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM)
    {
        //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << EB << std::endl;
        DetIDMap[cmsswId] = {iphi,ieta,TT25,EB};
        //auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )
{
    const std::string detIDConfigEE(ecal_config_path+"fullinfo_detids_EE.txt");
    std::ifstream infile( detIDConfigEE, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    TString EE;

    while (infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant)
    {
        ECAL ec = EM;
        if( side > 0 ) ec = EP;
        //std::cout << "DetID Input Line: " << cmsswId << " " << ix << " "  << iy << " " << ec << std::endl; 
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec};
        //auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

std::string addstr( std::string current, std::string input ){ return (current+input); }

void SetupDetIDEBHists( std::map<UInt_t,TH1F *> &DetIdHists )
{
    const std::string detIDConfigEB(ecal_config_path+"fullinfo_detids_EB.txt");
    std::ifstream infile( detIDConfigEB, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    TString pos;

    while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM)
    {
        //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << EB << std::endl;
        //DetIDMap[cmsswId] = {iphi,ieta,TT25,EB};
        std::string histname = addstr( "XtalTimeHist", to_string(cmsswId) );		
		DetIdHists[cmsswId] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",500,-5,5);
        //auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

const auto sq2(const float x){return x*x;}
const auto sq2(const double x){return x*x;}

void fillRatioHist(TH1D* numi, TH1D* denom, TH1D* result ){

    const auto nbins = numi->GetNbinsX();
    for (auto ibin = 0; ibin <= nbins; ibin++){
        auto nc = numi->GetBinContent(ibin);
        auto ncer = numi->GetBinError(ibin);
        auto dc = denom->GetBinContent(ibin);
        auto dcer = denom->GetBinError(ibin);
        auto ratio(0.0);
        auto rerr(0.0);
        if( dc > 20 ){
            ratio = nc/dc;
			rerr = std::sqrt( sq2(ncer/dc)-(sq2(ratio)/dc) );
			//rerr = std::sqrt((sq2(ncer/dc)+sq2((nc/sq2(dc))*dcer))/dc);
        }//<<>>if( dc > 0 )
        result->SetBinContent(ibin,ratio);
        result->SetBinError(ibin,rerr);
    }//<<>>for (auto ibinX = 1; ibinX <= nXbins; ibinX++)

}//<<>>fillRatioHist(TH1F* numi, TH1F* denom, TH1F* result )

void wc_ku_InterCali_aveRecHit_mini( string indir, string infilelistname, string outfilename, bool filterRecHits ){

    const int  nAlgos = 1; // Mini, MfootCCStc
    const double offset = 0.0;
    const int bin_offset = 86;
	const float minRhEnergy = 5.0;

	//const bool debug = true;
    const bool debug = false;
	const bool useEnergy = true;

    const string treename("tree/llpgtree");

    std::cout << "Opening Outfile : " << outfilename << std::endl;
    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

    TH2F * IcMapEB[nAlgos];
    TH2F * IcMapErrEB[nAlgos];
    TH2F * IcMapOccEB[nAlgos];
    TH1F * IcDistEB[nAlgos];
    TH1F * IcDistErrEB[nAlgos];
    TH2F * IcMapEP[nAlgos];
    TH2F * IcMapEM[nAlgos];

    TH2F * IcDistMeanEBTTMap;
    TH2F * IcDistMeanEBTTMapOcc;
    TH2F * IcDistMeanEBMDMap;
    TH2F * IcDistMeanEBMDMapOcc;
    TH2F * IcDistMeanEBSMMap;
    TH2F * IcDistMeanEBSMMapOcc;

    TH1F * IcDistMeanEBEta;
    TH1F * IcDistMeanErrEBEta;
    TH1F * IcDistMeanEBPhi;
    TH1F * IcDistMeanErrEBPhi;
    TH1F * IcMapETEta;
    TH1F * IcMapETPhi;
    TH1F * IcMapETErrEta;
    TH1F * IcMapETErrPhi;
    TH1F * IcDistMeanEBTT;
    TH1F * IcDistMeanErrEBTT;
    TH1F * IcDistMeanEBMD;
    TH1F * IcDistMeanErrEBMD;
    TH1F * IcDistMeanEBSM;
    TH1F * IcDistMeanErrEBSM;

    string algostring[1] = { "" };
    for( auto i = 0; i < nAlgos; i++){
        string hnameEB( "AveXtal"+algostring[i]+"RecTimeEBMap");
        string htitleEB( "AveXtal "+algostring[i]+" RecTimeEBMap EB ");
        IcMapEB[i] = new TH2F(hnameEB.c_str(),htitleEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
        IcMapEB[i]->Sumw2();
        string hnameEP( "AveXtal"+algostring[i]+"RecTimeEPMap");
        string htitleEP( "AveXtal "+algostring[i]+" RecTimeEPMap EP ");
        IcMapEP[i] = new TH2F(hnameEP.c_str(),htitleEP.c_str(),100,0.5,100.5,100,0.5,100.5);
        IcMapEP[i]->Sumw2();
        string hnameEM( "AveXtal"+algostring[i]+"RecTimeEMMap");
        string htitleEM( "AveXtal "+algostring[i]+" RecTimeEBMap EM ");
        IcMapEM[i] = new TH2F(hnameEM.c_str(),htitleEM.c_str(),100,0.5,100.5,100,0.5,100.5);
        IcMapEM[i]->Sumw2();
        string hnameOccEB( "AveXtal"+algostring[i]+"OccEBMap");
        string htitleOccEB( "AveXtal "+algostring[i]+" OccEBMap EB ");
        IcMapOccEB[i] = new TH2F(hnameOccEB.c_str(),htitleOccEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
        IcMapOccEB[i]->Sumw2();
        string hnameErrEB( "AveXtal"+algostring[i]+"RecTimeErrEBMap");
        string htitleErrEB( "AveXtal "+algostring[i]+" RecTimeErrEBMap EB ");
        IcMapErrEB[i] = new TH2F(hnameErrEB.c_str(),htitleErrEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
        IcMapErrEB[i]->Sumw2();
        string hnameDistEB( "AveXtal"+algostring[i]+"EBdist");
        string htitleDistEB( "AveXtal "+algostring[i]+" EBDist EB ");
        IcDistEB[i] = new TH1F(hnameDistEB.c_str(),htitleDistEB.c_str(),320,-4,4);
        IcDistEB[i]->Sumw2();
        string hnameErrDistEB( "AveXtal"+algostring[i]+"EBErrdist");
        string htitleErrDistEB( "AveXtal "+algostring[i]+" EBErrDist EB ");
        IcDistErrEB[i] = new TH1F(hnameErrDistEB.c_str(),htitleErrDistEB.c_str(),200,0,0.2);
        IcDistErrEB[i]->Sumw2();
    }//<<>>for( auto i = 0; i < nAlgos; i++)

    IcDistMeanEBEta = new TH1F("AveEtaRecTimeEta","AveEtaRecTimeEta;iEta;MeanTime [ns]",171,-85.5,85.5);
    IcDistMeanErrEBEta = new TH1F("AvePhiRecTimeEtaDist","AvePhiRecTimeEtaDist;MeanTime [ns]",160,-1,1);
    IcDistMeanEBEta->Sumw2();
    IcDistMeanErrEBEta->Sumw2();

    IcDistMeanEBPhi = new TH1F("AveEtaRecTimePhi","AveEtaRecTimePhi;iPhi;MeanTime [ns]",360,0.5,360.5);
    IcDistMeanErrEBPhi = new TH1F("AveEtaRecTimePhiDist","AveEtaRecTimePhiDist;MeanTime [ns]",160,-1,1);
    IcDistMeanEBPhi->Sumw2();
    IcDistMeanErrEBPhi->Sumw2();

    IcMapETEta = new TH1F("AveEtaRecMMTimeEta","AveEtaRecMeanMeanTimeEta;iEta;Mean MeanTime [ns]",171,-85.5,85.5);
    IcMapETPhi = new TH1F("AvePhiRecMMTimePhi","AvePhiRecMeanMeanTimePhi;iPhi;Mean MeanTime [ns]",360,0.5,360.5);
    IcMapETPhi->Sumw2();
    IcMapETEta->Sumw2();

    IcMapETErrEta = new TH1F("AveEtaRecMMTimeEtaErr","AveEtaRecMMTimeEtaErrt;Mean MeanTime [ns]",160,-1,1);
    IcMapETErrPhi = new TH1F("AvePhiRecMMTimePhiErr","AvePhiRecMMTimePhiErr;Mean MeanTime [ns]",160,-1,1);
    IcMapETErrEta->Sumw2();
    IcMapETErrPhi->Sumw2();

    IcDistMeanEBTTMap = new TH2F("AveTTRecTimeMap","AveTTRecTimeMap",34,0,34,72,0,72);
    IcDistMeanEBTT = new TH1F("AveTTRecTimeDist","AveTTRecTimeDist;MeanTime [ns]",160,-4,4);
    IcDistMeanErrEBTT = new TH1F("AveTTRecTimeErr","AveTTRecTimeErr;MeanTime Error [ns]",200,0,0.2);
    IcDistMeanEBTTMapOcc = new TH2F("AveTTRecTimeOccMap","AveTTRecTimeOccMap",34,0,34,72,0,72);
    IcDistMeanEBTTMap->Sumw2();
    IcDistMeanEBTT->Sumw2();
    IcDistMeanErrEBTT->Sumw2();
    IcDistMeanEBTTMapOcc->Sumw2();

	int nTTCellMaps = 49;
	TH1F * IcDistTTTimes[nTTCellMaps];
	for( int i = 0; i < nTTCellMaps; i++ ){ 
		IcDistTTTimes[i] = new TH1F(addstr("AveMDRecTimeDist_",to_string(i)).c_str(),"AveMDRecTimeDist;MeanTime [ns]",480,-24,24); }

    IcDistMeanEBMDMap = new TH2F("AveMDRecTimeMap","AveMDRecTimeMap",8,0,8,18,0,18);
    IcDistMeanEBMD = new TH1F("AveMDRecTimeDist","AveMDRecTimeDist;MeanTime [ns]",120,-3,3);
    IcDistMeanErrEBMD = new TH1F("AveMDRecTimeErr","AveMDRecTimeErr;MeanTime Error [ns]",400,0,0.04);
    IcDistMeanEBMDMapOcc = new TH2F("AveMDRecTimeOccMap","AveMDRecTimeOccMap",8,0,8,18,0,18);
    IcDistMeanEBMDMap->Sumw2();
    IcDistMeanEBMD->Sumw2();
    IcDistMeanErrEBMD->Sumw2();
    IcDistMeanEBMDMapOcc->Sumw2();

    IcDistMeanEBSMMap = new TH2F("AveSMRecTimeMap","AveSMRecTimeMap",2,0,2,18,0,18);
    IcDistMeanEBSM = new TH1F("AveSMRecTimeDist","AveSMRecTimeDist;MeanTime [ns]",80,-2,2);
    IcDistMeanErrEBSM = new TH1F("AveSMRecTimeErr","AveSMRecTimeErr;MeanTime Error [ns]",200,0,0.02);
    IcDistMeanEBSMMapOcc = new TH2F("AveSMRecTimeOccMap","AveSMRecTimeOccMap",2,0,2,18,0,18);
    IcDistMeanEBSMMap->Sumw2();
    IcDistMeanEBSM->Sumw2();
    IcDistMeanErrEBSM->Sumw2();
    IcDistMeanEBSMMapOcc->Sumw2();

    auto distRhTime = new TH1F("distrhtime","rh Time",2000,-50,50);
    auto distRhTimeTof = new TH1F("distrhtimetof","rh |Time|<25 + TOF",1000,-25,25);

    std::cout << "Setting up DetIDs." << std::endl;
    std::map< UInt_t,DetIDStruct > DetIDMap;
    SetupDetIDsEB( DetIDMap );
    SetupDetIDsEE( DetIDMap );
	//for( auto iter : DetIDMap ){ std::cout << iter.first << " " << iter.second.i1 << " " << iter.second.i2 << std::endl; }
    std::cout << "Setting up DetIDHists." << std::endl;
	std::map< UInt_t,TH1F * > DetIdHists;
	SetupDetIDEBHists( DetIdHists );
    std::cout << "DetIDHists Filled." << std::endl;

    std::map<UInt_t,Float_t> sumXtalMiniRecTime;
    std::map<UInt_t,Float_t> sumXtal2MiniRecTime;
    std::map<UInt_t,UInt_t> numXtalMiniRecTime;
	std::map<UInt_t,Float_t> xtalMeanRtRecTime;

    std::map<int,Float_t> sumXtalEtaRecTime;
    std::map<int,Float_t> sumXtalPhiRecTime;
    std::map<int,Float_t> sumXtal2EtaRecTime;
    std::map<int,Float_t> sumXtal2PhiRecTime;
    std::map<int,UInt_t> numXtalEtaRecTime;
    std::map<int,UInt_t> numXtalPhiRecTime;

    std::map<int,Float_t> sumXtalEtaRecMMTime;
    std::map<int,Float_t> sumXtalPhiRecMMTime;
    std::map<int,Float_t> sumXtal2EtaRecMMTime;
    std::map<int,Float_t> sumXtal2PhiRecMMTime;
    std::map<int,UInt_t> numXtalEtaRecMMTime;
    std::map<int,UInt_t> numXtalPhiRecMMTime;

    std::map<int,Float_t> sumXtalTTRecTime;
    std::map<int,Float_t> sumXtal2TTRecTime;
    std::map<int,UInt_t> numXtalTTRecTime;
    std::map<int,Float_t> sumXtalMDRecTime;
    std::map<int,Float_t> sumXtal2MDRecTime;
    std::map<int,UInt_t> numXtalMDRecTime;
    std::map<int,Float_t> sumXtalSMRecTime;
    std::map<int,Float_t> sumXtal2SMRecTime;
    std::map<int,UInt_t> numXtalSMRecTime;

    // Declaration of leaf types
    UInt_t run;
	ULong64_t event;
    vector<uInt> *rhID;
    vector<float> *rhRtTime;
    vector<float> *rhEnergy;
    vector<bool> *rhRtisOOT = 0;
    vector<bool> *rhisWeird = 0;
    vector<bool> *rhisDiWeird = 0;
    vector<float> *rhSwCross = 0;

    // List of branches
    TBranch *b_run;   //!
    TBranch *b_event;   //!
    TBranch *b_rhID;   //!
    TBranch *b_rhRtTime;   //!
    TBranch *b_rhEnergy;   //!
    TBranch *b_rhRtisOOT;
    TBranch *b_rhisWeird;
    TBranch *b_rhisDiWeird;
    TBranch *b_rhSwCross;

    std::ifstream infilelist(infilelistname);
    std::string infilestr;
    while (std::getline(infilelist,infilestr)){

    	std::stringstream ss(infilestr);
        std::string infilename;
        std::string srunstr;
        std::string erunstr;
        ss >> infilename >> srunstr >> erunstr;
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "For Run " << srunstr << " to Run " << erunstr << std::endl;
        auto srun = std::stoi(srunstr);
        auto erun = std::stoi(erunstr);

    	std::ifstream infile(infilename);
    	std::string instr;
        auto fInTree = new TChain(treename.c_str());
        std::cout << "Adding files to TChain." << std::endl;
		const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");
		//const std::string eosdir("root://cmseos.fnal.gov//store/user/");	
        while (std::getline(infile,instr)){
			auto tfilename = eosdir + indir + instr;
         	//auto tfilename = indir + "/" + str;
         	//std::cout << "--  adding file: " << tfilename << std::endl;
			std::cout << "-";
         	fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
		std::cout << std::endl;

        run = 0;
		event = 0;

        rhID = 0;
    	rhRtTime = 0;
   		rhEnergy = 0;

   		rhRtisOOT = 0;
   		rhisWeird = 0;
   		rhisDiWeird = 0;
   		rhSwCross = 0;

        fInTree->SetBranchAddress("run", &run, &b_run);
        fInTree->SetBranchAddress("event", &event, &b_event);

        fInTree->SetBranchAddress("rhID", &rhID, &b_rhID);
        fInTree->SetBranchAddress("rhRtTime", &rhRtTime, &b_rhRtTime);
		fInTree->SetBranchAddress("rhEnergy", &rhEnergy, &b_rhEnergy);

   		fInTree->SetBranchAddress("rhRtisOOT", &rhRtisOOT, &b_rhRtisOOT);
   		fInTree->SetBranchAddress("rhisWeird", &rhisWeird, &b_rhisWeird);
   		fInTree->SetBranchAddress("rhisDiWeird", &rhisDiWeird, &b_rhisDiWeird);
   		fInTree->SetBranchAddress("rhSwCross", &rhSwCross, &b_rhSwCross);
 
         // >> calcs  <<
     
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
            b_event->GetEntry(entry);   //!

            b_rhID->GetEntry(entry);   //!
            b_rhRtTime->GetEntry(entry);   //!
			b_rhEnergy->GetEntry(entry);

         	b_rhRtisOOT->GetEntry(entry);
         	b_rhisWeird->GetEntry(entry);
         	b_rhisDiWeird->GetEntry(entry);
         	b_rhSwCross->GetEntry(entry);

            if( debug ) std::cout << " - GetEntry : Energy "  << std::endl;
			//if( useEnergy ) b_rhEnergy->GetEntry(entry);   //!
			// time cuts ? non 0 and within -25 to 25 ? for AOD
			if( run < srun || run > erun ) continue;

            const auto nRecHits1 = rhID->size(); //(cluster[ipho0])->size();
            if( debug ) std::cout << "Looping over first recHits"  << std::endl;
            for ( auto rh_i = 0U; rh_i < nRecHits1; rh_i++ ){
				//auto rhe = useEnergy ? (*rhEnergy)[rh_i] : 999;	
				//if( rhe < minRhEnergy ) continue; 
                bool underMinEnergy = (*rhEnergy)[rh_i] < 5.0;
                bool rhTimeZero = (*rhRtTime)[rh_i] == 0.0;
                bool timeOutOfRange = std::abs((*rhRtTime)[rh_i]) > 25.0;
				bool isOOT = (*rhRtisOOT)[rh_i];
				bool isWeird = (*rhisWeird)[rh_i];
                bool isDiWeird = (*rhisDiWeird)[rh_i];
                bool hasHighSwissCross = (*rhSwCross)[rh_i] > 0.96;
				bool badRecHit = filterRecHits && ( isOOT or isWeird or isDiWeird or hasHighSwissCross );
                if( underMinEnergy or rhTimeZero or badRecHit ) continue;
                //if( timeOutOfRange ) continue;
                auto id_i = (*rhID)[rh_i];
                const auto & fill_idinfo = DetIDMap[id_i];
                auto iEta = fill_idinfo.i2;
                auto iPhi = fill_idinfo.i1;
                bool isEB = (fill_idinfo.ecal == ECAL::EB);
                auto Mini_t_i = (*rhRtTime)[rh_i];
                if( isEB ) distRhTime->Fill(Mini_t_i);
                //if( isEB ) distRhTimeTof->Fill(Mini_t_i);
     	     	if( debug ) std::cout << "Getting maps " << std::endl;
				//DetIdHists[id_i]->Fill(Mini_t_i);
                sumXtalMiniRecTime[id_i] += Mini_t_i; 
				numXtalMiniRecTime[id_i] += 1;
                sumXtal2MiniRecTime[id_i] += Mini_t_i*Mini_t_i;
                int ttphi = (iPhi-1)/5;
                int tteta = (std::abs(iEta)-1)/5;
                int ttCellPhi = ((ttphi+7)%10 == 0 ) ? (ttphi+7)/10 : -100;
				int ttacteta = ( iEta < 0 ) ? tteta : tteta + 17;
				int ttCellEta = ((ttacteta+2)%5 == 0 ) ? (ttacteta+2)/5 : -100;
                int ttCellIdx = ((ttCellPhi-1)*7)+(ttCellEta-1);
                int ttidx = ( iEta < 0 ) ? (ttphi+(tteta+1)*100) * -1 : (ttphi+tteta*100);
                int mdphi = (iPhi-1)/20;
                int mdeta = ( std::abs(iEta) < 26 ) ? 0 : 1 + (std::abs(iEta)-26)/20;
                int mdidx = ( iEta < 0 ) ? (mdphi+(mdeta+1)*100) * -1 : (mdphi+mdeta*100);
                int smphi = (iPhi-1)/20;
                int smidx = ( iEta < 0 ) ? (smphi+1) * -1 : smphi;
                if( isEB ){

					DetIdHists[id_i]->Fill(Mini_t_i);
                    sumXtalEtaRecTime[iEta] += Mini_t_i;
                    numXtalEtaRecTime[iEta] += 1;
                    sumXtal2EtaRecTime[iEta] += Mini_t_i*Mini_t_i;

                    sumXtalPhiRecTime[iPhi] += Mini_t_i;
                    numXtalPhiRecTime[iPhi] += 1;
                    sumXtal2PhiRecTime[iPhi] += Mini_t_i*Mini_t_i;

                    sumXtalTTRecTime[ttidx] += Mini_t_i;
                    sumXtal2TTRecTime[ttidx] += Mini_t_i*Mini_t_i;
                    numXtalTTRecTime[ttidx] += 1;
					if( ttCellIdx > -1 && ttCellIdx < 100 ){ 
						if( ttCellIdx < nTTCellMaps ) IcDistTTTimes[ttCellIdx]->Fill(Mini_t_i);
						//std::cout << " TTCell idx : " << ttCellIdx << " TTCell E-P : " << ttCellEta << " - ";
						//std::cout << ttCellPhi << " Act E-P : " << iEta << " - " << iPhi << " TT : " << ttidx << std::endl;
					}//<<>>if( ttCellIdx > -1 && ttCellIdx < 100 )
                    sumXtalMDRecTime[mdidx] += Mini_t_i;
                    sumXtal2MDRecTime[mdidx] += Mini_t_i*Mini_t_i;
                    numXtalMDRecTime[mdidx] += 1;

                    sumXtalSMRecTime[smidx] += Mini_t_i;
                    sumXtal2SMRecTime[smidx] += Mini_t_i*Mini_t_i;
                    numXtalSMRecTime[smidx] += 1;

                }//<<>>if( fill_idinfo.ecal == ECAL::EB )
             }//<<>>for (auto i = 0U; i < nRecHits1; i++) // end loop over rechits
             if( debug ) std::cout << "RecHits Loop done "<< std::endl;
         }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)
		delete fInTree;

    } // while (std::getline(infilelist,infiles))

	std::cout << "Filling Calibration Maps" <<  std::endl;
	//auto fitForm = new TFormula("fitFormula","[0]*exp(-0.5*((x-[1])/[2])**2)");
    //auto fitFunc  = new TF1("gfit","fitFormula",-3.0,3.0);
    std::map<UInt_t,Float_t> *  icmaps[nAlgos] = {&sumXtalMiniRecTime};
    std::map<UInt_t,UInt_t> *  nicmaps[nAlgos] = {&numXtalMiniRecTime};
    std::map<UInt_t,Float_t> *  ic2maps[nAlgos] = {&sumXtal2MiniRecTime};
	std::map<UInt_t,Float_t> *  meanMaps[nAlgos] = {&xtalMeanRtRecTime};
    for( auto ai = 0; ai < nAlgos; ai++ ){
         for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it){
            const auto & fill_idinfo = DetIDMap[it->first];
            float map_time = (((*icmaps[ai])[it->first])/((*nicmaps[ai])[it->first])) + offset; 
				// - (drift/(icmaps[ai]->size()))) + offset;
            float map_occ = (*nicmaps[ai])[it->first];
            float map_err = sqrt((((*ic2maps[ai])[it->first])/map_occ - map_time*map_time)/map_occ);
			//if( debug ) std::cout << "Fill hist for Algo " << ai << " at " << fill_idinfo.i2; 
			//if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << " for iter " << std::endl;
            if( fill_idinfo.ecal == ECAL::EB ){
		   		if( debug ) std::cout << "Fill EB hist for Algo " << ai << " at " << fill_idinfo.i2 << " "; 
                if( debug ) std::cout << fill_idinfo.i1 << " with " << map_time << std::endl;
				if( filterRecHits ){
					auto fitFunc  = new TF1("gfit","gaus",-3.0,3.0);
					TH1F* thisHist = DetIdHists[it->first]; 
					auto thefit = thisHist->Fit("gfit","QNS");
					float ftime = 1.f*fitFunc->GetParameter(1);
					//map_time = ( std::abs(ftime) < 5.0 ) ? time : map_time;
					//map_err = (fitFunc->GetParameter(1))/std::sqrt(map_occ);
					float ferr = 1.f*fitFunc->GetParError(1)/std::sqrt(map_occ);
					if( ferr > 1.f ){ map_time = ftime; map_err = ferr; }
					//std::cout << " FitFill : " << map_time << " : " << map_err << std::endl;
					delete fitFunc;					
				}//<<>>if( filterRecHits )
            	(IcMapEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
                (IcMapOccEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_occ );
                (IcMapErrEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_err );
                (IcDistEB[ai])->Fill(map_time);
                (IcDistErrEB[ai])->Fill(map_err);
				(*meanMaps[ai])[it->first] = map_time;
            } else if( fill_idinfo.ecal == ECAL::EP ){
                if( debug ) std::cout << "Fill EP hist for Algo " << ai << " at " << fill_idinfo.i2; 
                if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                (IcMapEP[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            } else if( fill_idinfo.ecal == ECAL::EM ){
                if( debug ) std::cout << "Fill EM hist for Algo " << ai << " at " << fill_idinfo.i2;
                if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                (IcMapEM[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            }//<<>>if( fill_idinfo.ecal == ECAL::EB )
         }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)
    }//<<>>for( auto ai = 0; ai < nAlgos; ai++ )

    for( std::map<int,Float_t>::iterator it=sumXtalEtaRecTime.begin(); it!=sumXtalEtaRecTime.end(); ++it){
       const auto iEta = it->first;
       const auto & map_time = sumXtalEtaRecTime[iEta]/numXtalEtaRecTime[iEta]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalEtaRecTime[iEta];
       const auto & map_err = sqrt((sumXtal2EtaRecTime[iEta]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcDistMeanEBEta->SetBinContent(iEta+86,map_time);
           IcDistMeanEBEta->SetBinError(iEta+86,map_err);
           IcDistMeanErrEBEta->Fill( map_time );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalPhiRecTime.begin(); it!=sumXtalPhiRecTime.end(); ++it){
       const auto iPhi = it->first;
       const auto & map_time = sumXtalPhiRecTime[iPhi]/numXtalPhiRecTime[iPhi]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalPhiRecTime[iPhi];
       const auto & map_err = sqrt((sumXtal2PhiRecTime[iPhi]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcDistMeanEBPhi->SetBinContent(iPhi,map_time);
           IcDistMeanEBPhi->SetBinError(iPhi,map_err);
           IcDistMeanErrEBPhi->Fill( map_time );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalEtaRecMMTime.begin(); it!=sumXtalEtaRecMMTime.end(); ++it){
       const auto iEta = it->first;
       const auto & map_time = sumXtalEtaRecMMTime[iEta]/numXtalEtaRecMMTime[iEta]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalEtaRecMMTime[iEta];
       const auto & map_err = sqrt((sumXtal2EtaRecMMTime[iEta]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           //IcMapETEta->Fill( iEta, map_time );
           IcMapETErrEta->Fill( map_time );
           IcMapETEta->SetBinContent(iEta+86,map_time);
           IcMapETEta->SetBinError(iEta+86,map_err);
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalEtaRecMMTime.begin(); it!=sumXtalEtaRecMMTime.end(); ++it){
       const auto iEta = it->first;
       const auto & map_time = sumXtalEtaRecMMTime[iEta]/numXtalEtaRecMMTime[iEta]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalEtaRecMMTime[iEta];
       const auto & map_err = sqrt((sumXtal2EtaRecMMTime[iEta]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           //IcMapETEta->Fill( iEta, map_time );
           IcMapETErrEta->Fill( map_time );
           IcMapETEta->SetBinContent(iEta+86,map_time);
           IcMapETEta->SetBinError(iEta+86,map_err);
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalPhiRecMMTime.begin(); it!=sumXtalPhiRecMMTime.end(); ++it){
       const auto iPhi = it->first;
       const auto & map_time = sumXtalPhiRecMMTime[iPhi]/numXtalPhiRecMMTime[iPhi]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalPhiRecMMTime[iPhi];
       const auto & map_err = sqrt((sumXtal2PhiRecMMTime[iPhi]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcMapETErrPhi->Fill( map_time );
           IcMapETPhi->SetBinContent(iPhi,map_time);
           IcMapETPhi->SetBinError(iPhi,map_err);
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    std::cout << "TT : map fills " << std::endl;
    for( std::map<int,Float_t>::iterator it=sumXtalTTRecTime.begin(); it!=sumXtalTTRecTime.end(); ++it){
       const int iTT = it->first;
       const int eta = iTT/100;
       const int iPhiTT = std::abs(iTT) - std::abs(eta)*100;
       const int iEtaTT = eta + 17;
       //std::cout << "TT : iTT = " << iTT << " phi: " << iPhiTT << " eta: " << iEtaTT << std::endl;
       const float & map_time = sumXtalTTRecTime[iTT]/numXtalTTRecTime[iTT]; // - (drift/(icmaps[ai]->size()))) + offset;
       const int & map_occ = numXtalTTRecTime[iTT];
       const float & map_err = sqrt((sumXtal2TTRecTime[iTT]/map_occ - map_time*map_time)/map_occ);
       //std::cout << "TT vals : time: " << map_time << " occ: " << map_occ << " err: " << map_err << std::endl;
       //if( fill_idinfo.ecal == ECAL::EB ){
           (IcDistMeanEBTTMap)->Fill(iEtaTT,iPhiTT,map_time);
           (IcDistMeanEBTTMapOcc)->Fill(iEtaTT,iPhiTT,map_occ);
           IcDistMeanEBTT->Fill( map_time );
           IcDistMeanErrEBTT->Fill( map_err );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    std::cout << "MD : map fills " << std::endl;
    for( std::map<int,Float_t>::iterator it=sumXtalMDRecTime.begin(); it!=sumXtalMDRecTime.end(); ++it){
       const int iMD = it->first;
       const int eta = iMD/100;
       const int iPhiMD = std::abs(iMD) - std::abs(eta)*100;
       const int iEtaMD = eta + 4;
       //std::cout << "MD : iMD = " << iMD << " phi: " << iPhiMD << " eta: " << iEtaMD << std::endl;
       const float & map_time = sumXtalMDRecTime[iMD]/numXtalMDRecTime[iMD]; // - (drift/(icmaps[ai]->size()))) + offset;
       const int & map_occ = numXtalMDRecTime[iMD];
       const float & map_err = sqrt((sumXtal2MDRecTime[iMD]/map_occ - map_time*map_time)/map_occ);
       //std::cout << "MD vals : time: " << map_time << " occ: " << map_occ << " err: " << map_err << std::endl;
       //if( fill_idinfo.ecal == ECAL::EB ){
           (IcDistMeanEBMDMap)->Fill(iEtaMD,iPhiMD,map_time);
           (IcDistMeanEBMDMapOcc)->Fill(iEtaMD,iPhiMD,map_occ);
           IcDistMeanEBMD->Fill( map_time );
           IcDistMeanErrEBMD->Fill( map_err );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    std::cout << "SM : map fills " << std::endl;
    for( std::map<int,Float_t>::iterator it=sumXtalSMRecTime.begin(); it!=sumXtalSMRecTime.end(); ++it){
       const int iSM = it->first;
       //const int eta = iSM/100;
       const int iPhiSM = ( iSM < 0 ) ? std::abs(iSM) - 1 : std::abs(iSM);
       const int iEtaSM = ( iSM < 0 ) ? 0 : 1;
       //std::cout << "SM : iSM = " << iSM << " phi: " << iPhiSM << " eta: " << iEtaSM << std::endl;
       const float & map_time = sumXtalSMRecTime[iSM]/numXtalSMRecTime[iSM]; // - (drift/(icmaps[ai]->size()))) + offset;
       const int & map_occ = numXtalSMRecTime[iSM];
       const float & map_err = sqrt((sumXtal2SMRecTime[iSM]/map_occ - map_time*map_time)/map_occ);
       //std::cout << "SM vals : time: " << map_time << " occ: " << map_occ << " err: " << map_err << std::endl;
       //if( fill_idinfo.ecal == ECAL::EB ){
           (IcDistMeanEBSMMap)->Fill(iEtaSM,iPhiSM,map_time);
           (IcDistMeanEBSMMapOcc)->Fill(iEtaSM,iPhiSM,map_occ);
           IcDistMeanEBSM->Fill( map_time );
           IcDistMeanErrEBSM->Fill( map_err );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)



    for( auto ai = 0; ai < nAlgos; ai++ ){ (*icmaps[ai]).clear(); (*nicmaps[ai]).clear(); (*meanMaps[ai]).clear(); }

    fOutFile->cd();

    std::cout << "Write AveXtal Rechit Time Maps" << std::endl;

    distRhTime->Write();
    distRhTimeTof->Write();
    delete distRhTime;
    delete distRhTimeTof;

	for( auto detidHist : DetIdHists ){ delete detidHist.second; }
    for( int i = 0; i < nTTCellMaps; i++ ){ IcDistTTTimes[i]->Write(); delete IcDistTTTimes[i]; }

    for( auto i = 0; i < nAlgos; i++){

         IcMapEB[i]->Write();
         IcMapErrEB[i]->Write();
         IcMapOccEB[i]->Write();
         IcDistEB[i]->Write();
         IcDistErrEB[i]->Write();
         IcMapEP[i]->Write();
         IcMapEM[i]->Write();

         delete IcMapEB[i];
         delete IcMapErrEB[i];
         delete IcMapOccEB[i];
         delete IcDistEB[i];
         delete IcDistErrEB[i];
         delete IcMapEP[i];
         delete IcMapEM[i];

    }//<<>>for( auto i = 0; i < nAlgos; i++)

    IcMapETErrPhi->Write();
    IcMapETPhi->Write();
    IcMapETErrEta->Write();
    IcMapETEta->Write();
    IcDistMeanErrEBPhi->Write();
    IcDistMeanEBPhi->Write();
    IcDistMeanErrEBEta->Write();
    IcDistMeanEBEta->Write();

    IcDistMeanEBTTMap->Write();
    IcDistMeanEBTTMapOcc->Write();
    IcDistMeanEBTT->Write();
    IcDistMeanErrEBTT->Write();

    IcDistMeanEBMDMap->Write();
    IcDistMeanEBMDMapOcc->Write();
    IcDistMeanEBMD->Write();
    IcDistMeanErrEBMD->Write();

    IcDistMeanEBSMMap->Write();
    IcDistMeanEBSMMapOcc->Write();
    IcDistMeanEBSM->Write();
    IcDistMeanErrEBSM->Write();

    delete IcMapETErrPhi;
    delete IcMapETPhi;
    delete IcMapETErrEta;
    delete IcMapETEta;
    delete IcDistMeanErrEBPhi;
    delete IcDistMeanEBPhi;
    delete IcDistMeanErrEBEta;
    delete IcDistMeanEBEta;

    delete IcDistMeanEBTTMap;
    delete IcDistMeanEBTTMapOcc;
    delete IcDistMeanEBTT;
    delete IcDistMeanErrEBTT;

    delete IcDistMeanEBMDMap;
    delete IcDistMeanEBMDMapOcc;
    delete IcDistMeanEBMD;
    delete IcDistMeanErrEBMD;

    delete IcDistMeanEBSMMap;
    delete IcDistMeanEBSMMapOcc;
    delete IcDistMeanEBSM;
    delete IcDistMeanErrEBSM;

    //delete fInTree;
    delete fOutFile;

}//<<>>void wc_ku_InterCali_aveRecHit_mini( string indir, string infilelistname, string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {

        //auto indir = "/ecalTiming/gammares_llpana_pd/MET/";
		//auto indir = "ecalTiming/gammares_llpana/MET/";
        //auto indir = "/ecalTiming/gammares_llpana_qcd/";
        //auto indir = "/ecalTiming/gammares_llpana_v2/MET/";
		//auto indir = "/ecalTiming/gammares_llpana_pd/";
        //auto indir = "ecalTiming/gammares_llpana/";
        //auto indir = "/ecalTiming/gammares_llpana_mc/"; 
        //auto indir = "/ecalTiming/gammares_r24f_prompt/";
        auto indir = "/ecalTiming/gammares_ECAL_CC_HCAL_DI-v3/";

        //auto infilelist = "master_list_files/egammares_Met_PD_AOD_Run2017E-17Nov2017v2_reso_califilelist.txt";
		//auto infilelist = "master_list_files/egammares_MetPD_MINIAOD_Run2017E_304475_califilelist.txt";
        //auto infilelist = "master_list_files/egammares_Met_PD_AOD_Run2017E-17Nov2017v2_reso_califilelist.txt";
        //auto infilelist = "master_list_files/kuntuple_QCDHT_Met75_R17_v20_califilelist.txt";
        //auto infilelist = "master_list_files/egammares_MetPD_MINIAOD_Run2017E_Single_califilelist.txt";
        //auto infilelist = "master_list_files/egammares_MetPD_MINIAOD_Run2017E_v2_Full_califilelist.txt";
        //auto infilelist = "master_list_files/egammares_DEGPD_AOD_Run2017F_califilelist.txt";
        //auto infilelist = "master_list_files/egammares_EGMPD_MINIAOD_Run2018D_v2_327238_califilelist.txt";
        //auto infilelist = "master_list_files/kuntuple_QCDHT_Met75_R17_v20_califilelist.txt";
        //auto infilelist = "master_list_files/kuntuple_QCDHT100200_Met75_R17_v20_infileslist.txt";
        //auto infilelist = "master_list_files/egammares_DYJetsToLL_Met75_R17_v20_califilelist.txt";
        //auto infilelist = "master_list_files/egammares_DEGPD_AOD_Run2017E_califilelist.txt";
        auto infilelist = "kucmsTimeCaliR24FCCvRt_Cali_TFile.txt";

        //auto outfilename = "KURes_14011_v12_MetPD_AOD_Run2017E_Cali_FilterFit_Full.root";
        //auto outfilename = "KURes_14011_v12_MetPD_MINIAOD_Run2017E_Cali_Filtered_304475.root";
        //auto outfilename = "KURes_14011_v12_MetPD_MINIAOD_Run2017E_Cali_Single.root";
        //auto outfilename = "KURes_14011_v12_QCD_AOD_Met75_R17_Cali_Filtered.root";
        //auto outfilename = "KURes_14011_v12_MetPD_MINIAOD_Run2017E_v2_Cali_FilterFit_Full.root";
        //auto outfilename = "KURes_14011_v12_DEGPD_AOD_Run2017F_Cali_FilterFit_305040_305365.root";
        //auto outfilename = "KURes_14011_v12_EGMPD_MINIAOD_Run2018D_Cali_FilterFit_324305_327239.root";
        //auto outfilename = "KURes_14011_v12_QCD_AODSIM_Met75_R17_Cali_FilterFit.root";
        //auto outfilename = "KURes_14011_v12_QCDHT100200_AODSIM_Met75_R17_Cali_FilterFit.root";
        //auto outfilename = "KURes_14011_v12_DY1JetsToLL_AODSIM_Met75_R17_Cali_FilterFit.root";
		//auto outfilename = "KURes_14011_v12_DEGPD_AOD_Run2017E_Cali_FilterFit_301487_304475.root";
        auto outfilename = "KURes_14011_v12_DEGPD_MiniAOD_ECAL_CC_HCAL_DI-v3_Cali.root";

		//bool filterRecHits = true;
        bool filterRecHits = false;

        wc_ku_InterCali_aveRecHit_mini( indir, infilelist, outfilename, filterRecHits );

    //}//<<>>if( argc != 4 )
    return 1;

}//<<>>int main ( int argc, char *argv[] )

