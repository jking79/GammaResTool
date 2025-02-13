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
  	DetIDStruct(const Int_t ni1, const Int_t ni2, const Int_t nTT, const Int_t & necal) : i1(ni1), i2(ni2), TT(nTT), ecal(necal){}
  	Int_t i1; // EB: iphi, EE: ix
  	Int_t i2; // EB: ieta, EE: iy
  	Int_t TT; // trigger tower
  	Int_t ecal; // EB, EM, EP
};

// helper functions ------------------------------------------------------------

void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )
{   
    const std::string detIDConfigEB("ecal_config/fullinfo_detids_EB.txt");
    std::ifstream infile( detIDConfigEB, std::ios::in);
    
    UInt_t cmsswId, dbID; 
    Int_t hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    TString pos;
    
    while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM)
    {   
        //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << 0 << std::endl;
        DetIDMap[cmsswId] = {iphi,ieta,TT25,0};
        //auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )
{
    const std::string detIDConfigEE("ecal_config/fullinfo_detids_EE.txt");
    std::ifstream infile( detIDConfigEE, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    TString EE;

    while (infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant)
    {
        int ec = 1;
        if( side > 0 ) ec = 2;
        //std::cout << "DetID Input Line: " << cmsswId << " " << ix << " "  << iy << " " << ec << std::endl; 
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec};
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
const char* addstrc( std::string current, std::string input ){ return (current+input).c_str(); }//<<>>std::string addstr( std::string current, char* input )

void fillOUHist1F( float val, float low, float high, float div, TH1F * hist ){

    auto step = ((high-low)/div)/2;
    if( val < low ) hist->Fill( low+step );
    else if ( val > high ) hist->Fill( high-step );
    else hist->Fill( val );

}//<<>>void fillOUHist1F( float val, float low, float high, TH1F & hist )

void fillTH1( float val, TH1F *& hist ){

   auto nBins = hist->GetNbinsX();
   auto low = hist->GetBinCenter(1);
   auto high = hist->GetBinCenter(nBins);
   if( val < low ) hist->Fill( low );
   else if ( val > high ) hist->Fill( high );
   else hist->Fill( val );

}//<<>>void fillTH1F( float val, TH1F *& hist )

void fillTH1( float val, TH1D* hist ){

    auto nBins = hist->GetNbinsX();
    auto low = hist->GetBinCenter(1);
    auto high = hist->GetBinCenter(nBins);
    if( val < low ) hist->Fill( low );
    else if ( val > high ) hist->Fill( high );
    else hist->Fill( val );

}//<<>>void fillTH1( float val, TH1D* hist )

void normTH2D(TH2D* hist){

    std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;

    const auto nXbins = hist->GetNbinsX();
    const auto nYbins = hist->GetNbinsY();

    for (auto ibinX = 1; ibinX <= nXbins; ibinX++){

        const auto norm = hist->Integral(ibinX,ibinX,1,nYbins);
        if( norm == 0.0 ) continue;
        for (auto ibinY = 1; ibinY <= nYbins; ibinY++){

            // get content/error
            auto content = hist->GetBinContent(ibinX,ibinY);
            auto error   = hist->GetBinError  (ibinX,ibinY);
            // set new contents
            content /= norm;
            error /= norm;
            hist->SetBinContent(ibinX,ibinY,content);
            hist->SetBinError  (ibinX,ibinY,error);

        }//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
    }//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){

}//<<>>void NormTH2D(TH2D* hist){

void normTH1D(TH1D* hist){

    std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;

    const auto nBins = hist->GetNbinsX();
    const auto norm = hist->Integral();
    for (auto ibinX = 1; ibinX <= nBins; ibinX++){

        if( norm == 0.0 ) continue;
        // get content/error
        auto content = hist->GetBinContent(ibinX);
        auto error   = hist->GetBinError(ibinX);
        // set new contents
        content /= norm;
        error /= norm;
        hist->SetBinContent(ibinX,content);
        hist->SetBinError  (ibinX,error);

    }//<<>>for (auto ibinX = 1; ibinX <= nBins; ibinX++)

}//<<>>void NormTH1D(TH1D* hist)

void profileTH2D(TH2D* nhist, TH1D* prof, TH1D* fithist){

    std::cout << "Profile " << " hist: " << nhist->GetName() << std::endl;

    const auto nXBins = nhist->GetNbinsX();
    //const auto nYBins = nhist->GetNbinsY();
    for (auto ibinX = 1; ibinX <= nXBins; ibinX++){

        auto phist = (TH1F*)nhist->ProjectionY("temp",ibinX,ibinX);

        auto mean = phist->GetMean();
        auto stdv = phist->GetStdDev();
        auto norm = phist->GetBinContent(phist->GetMaximumBin());
        auto high = mean + 0.2*stdv;
        auto low = mean - 0.2*stdv;
        //std::cout << " - Profile: m " << mean << " s " << stdv << " h " << high << " l " << low << " n " << norm << std::endl;
        if( abs(stdv) > 0.01 && abs(norm) > 1 ){
            auto tmp_form = new TFormula("tmp_formula","[0]*exp(-0.5*((x-[1])/[2])**2)");
            auto tmp_fit  = new TF1("tmp_fit",tmp_form->GetName(),low,high);
            tmp_fit->SetParameter(0,norm); //tmp_fit->SetParLimits(0,norm/2,norm*2);
            tmp_fit->SetParameter(1,mean); //tmp_fit->SetParLimits(1,-2,2);
            tmp_fit->SetParameter(2,stdv); //tmp_fit->SetParLimits(2,0,10);
            phist->Fit(tmp_fit->GetName(),"RBQ0");
            auto fmean = tmp_fit->GetParameter(1);
            auto error = tmp_fit->GetParError(1);
            auto fNdf = tmp_fit->GetNDF();
            auto fProb = tmp_fit->GetProb();

            // set new contents
            if( fNdf > 0 && fProb > 0.05 && error < 1.0 ){
                fithist->SetBinContent( ibinX, fProb );
                fithist->SetBinError( ibinX, 0 );
                prof->SetBinContent( ibinX, fmean );
                prof->SetBinError( ibinX, error );
            }//<<>>if( fmean < 1 && error < 0.1 )

            //delete tmp_form;
            delete tmp_fit;
        }//<<>>if( stdv > 0.01 )

    }//<<>>for (auto ibinX = 1; ibinX <= nBins; ibinX++)

}//<<>>void profileTH2D(TH2D* hist, TH1D* prof)

void thresDivTH2D(TH2D* numi, TH2D* denom, float thres){

    std::cout << "Threshold Division - " << " hist: " << numi->GetName() << std::endl;
    const auto nXbins = numi->GetNbinsX();
    const auto nYbins = numi->GetNbinsY();
    for (auto ibinX = 1; ibinX <= nXbins; ibinX++){
        for (auto ibinY = 1; ibinY <= nYbins; ibinY++){
            // get content/error
            auto ncontent = numi->GetBinContent(ibinX,ibinY);
            auto nerror   = numi->GetBinError  (ibinX,ibinY);
            auto dcontent = denom->GetBinContent(ibinX,ibinY);
            auto derror   = denom->GetBinError  (ibinX,ibinY);
            // set new contents
            auto content(0.0);
            auto error(0.0);
            if( dcontent > thres ){ content = ncontent/dcontent; error = nerror/derror; }
            numi->SetBinContent(ibinX,ibinY,content);
            numi->SetBinError  (ibinX,ibinY,error);
        }//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
    }//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){

}//<<>>void thresDivTH2D(TH2D* numi, TH2D* denom, float thres){

void fillMeanHist(TH1D* numi, TH1D* denom, TH1D* result ){

    const auto nbins = numi->GetNbinsX();
    for (auto ibin = 0; ibin <= nbins; ibin++){
        auto nc = numi->GetBinContent(ibin);
        auto ncer = numi->GetBinError(ibin);
        auto dc = denom->GetBinContent(ibin);
        auto dcer = denom->GetBinError(ibin);
        float ratio(0.0);
        float rerr(0.0);
		float mnerr(0.0);
        if( dc > 20 ){
            ratio = nc/dc;
			mnerr = std::sqrt( sq2(ncer/dc)-(sq2(ratio)/dc) ); 
        }//<<>>if( dc > 0 )
        result->SetBinContent(ibin,ratio);
        result->SetBinError(ibin,mnerr);
    }//<<>>for (auto ibinX = 1; ibinX <= nXbins; ibinX++)

}//<<>>fillRatioHist(TH1F* numi, TH1F* denom, TH1F* result )

void fillRatioHist(TH1D* numi, TH1D* denom, TH1D* result ){

    const auto nbins = denom->GetNbinsX();
    for (auto ibin = 1; ibin < nbins; ibin++){
        auto nc = numi->GetBinContent(ibin);
        auto ncer = numi->GetBinError(ibin);
        auto dc = denom->GetBinContent(ibin);
        auto dcer = denom->GetBinError(ibin);
        if( dc > 0 ){
            auto ratio = nc/dc;
			auto rerr = std::sqrt( ((nc+1)*(dc-nc+1)) / (sq2(dc+2)*(dc+3)) );	
        	result->SetBinContent(ibin,ratio);
        	result->SetBinError(ibin,rerr);
        }//<<>>if( dc > 0 )
    }//<<>>for (auto ibinX = 1; ibinX <= nXbins; ibinX++)

}//<<>>fillRatioHist(TH1F* numi, TH1F* denom, TH1F* result )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////       Class Declaration 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class makehists : egammares_hist_base {

	public:

	void llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename, std::string califilename, 
								int brun, int erun, std::string fhtitle, bool useAllData  );	
	void initHists( std::string fHTitle );
	void getBranches( Long64_t entry );
	void eventLoop( Long64_t entry );
 	void endJobs();	
	void initCali( std::string califilename ); 
	void makeEBEEMaps( std::vector<unsigned int> rhcol );
	void makeEBEEMaps( int phoit );
	int getRhIdx( uInt rhDetID );

	std::map<UInt_t,DetIDStruct> DetIDMap;

	TH2F *icmap[3];
    bool useAll;
	int startRun, endRun;
    float totrhs, totrhs0, totrhs05, totrhs1, totrhs2, totrhs5, totrhs10, encrhs;

    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];

    int nMaps;
    bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
	
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class Functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int makehists::getRhIdx( uInt rhDetID ){

    for( int idx = 0; idx < rhCaliID->size(); idx++ ){ if( rhDetID == (*rhCaliID)[idx] ) return idx; }
    //std::cout << " -- !! no rhDetID to (*rhCaliID)[idx] match !! ---------------------- " << std::endl;
    return -1;

}//<<>>int makehists::getRhIdx( int rhDetID )

void makehists::makeEBEEMaps( std::vector<unsigned int> rhcol ){

    if( DEBUG ) std::cout << " -- looping over ebee maps:" << std::endl;
    if( nMaps >= nEBEEMaps ) return;
    bool fill(false);
    auto nrh = rhcol.size();
    if( nMaps < 6 ){ if( nrh < 25 ) fill = true; }
    else if( nMaps < 12 ){ if( nrh < 25 ) fill = true; }
    else if( nMaps < 18 ){ if( nrh < 45 ) fill = true; }
    else if( nMaps < 24 ){ if( nrh < 65 ) fill = true; }
    else if( nMaps < 30 ){ if( nrh < 85 ) fill = true; }
    else { if( nrh >= 85 ) fill = true; }
    if( DEBUG ) std::cout << " -- start fill of ebee maps:" << std::endl;
    if(fill){
        if( DEBUG ) std::cout << " -- Filling ebeeMapT : " << nMaps << " with nRH : " << nrh << std::endl;
        for( int idx = 0; idx < nrh; idx++ ){
            const auto rhIDX = getRhIdx(rhcol[idx]);
            auto idinfo = DetIDMap[rhcol[idx]];
            const auto rhEtaPos = idinfo.i2;//recHitPos.eta();
            const auto rhPhiPos = idinfo.i1;//recHitPos.phi();
            //auto res = 1/(sq2(3.64/(*rhEnergy)[rhIDX])+0.18);
			//auto res = (*rhEnergy)[rhIDX];
            ebeeMapP[nMaps]->Fill( rhEtaPos, rhPhiPos,1);
            ebeeMapP[nMaps]->Fill( nMaps, 1,1);
            ebeeMapT[nMaps]->Fill( rhEtaPos, rhPhiPos,(*rhCaliCCTime)[rhIDX]);
            ebeeMapR[nMaps]->Fill( rhEtaPos, rhPhiPos,(*rhEnergy)[rhIDX]);
        }//<<>>for( idx = 0; idx < nRecHits; idx++ )
        nMaps++;
    }//<<>>if(fill) 

}//<<>>void makehists::makeEBEEMaps( vector<unsigned int> )


void makehists::makeEBEEMaps( int phoit ){

    if( DEBUG ) std::cout << " - looping photons : making ebee maps" << std::endl;
    auto isEB = true; //(*phoIsEB)[phoit];
    auto rhcol = (*phoRhIds)[phoit];
    if( isEB ) makeEBEEMaps(rhcol);
    if( DEBUG ) std::cout << " - looping photons : Finished making ebee maps" << std::endl;
    return;

}//int makehists::makeEBEEMaps( int phoit )

//-----------------------------------------------------------------------------------------------------------------------------
////----------------------------------------  Make Hists function call  ---------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

void makehists::llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename, std::string califilename, 
										int brun, int erun, std::string fhtitle, bool useAllData ){

    //bool debug = true; // this is for main loop only
	bool debug = false;

    // these need to be changed to input paramters for public verions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    const std::string disphotreename("tree/llpgtree");
    const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");
    //const std::string listdir("llpgana_list_files/");

    std::cout << "Producing Histograms for : " << outfilename << std::endl;
    std::ifstream infile(infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infilelist << " >> " << fInTree << std::endl;
    std::string str;
    while (std::getline(infile,str)){
        auto tfilename = eosdir + indir + str;
        if(debug) std::cout << "--  adding file: " << tfilename << std::endl;
        fInTree->Add(tfilename.c_str());
    }

	Init(fInTree);
	initHists(fhtitle);

	useAll = useAllData;

	startRun = brun;
	endRun = erun;

    SetupDetIDsEB(DetIDMap);
    SetupDetIDsEE(DetIDMap);

	initCali(califilename);

    std::cout << "Setting up For Main Loop." << std::endl;
    auto nEntries = fInTree->GetEntries();
    if(debug) nEntries = 1000;
    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
    int npace = nEntries/20;
    for (Long64_t centry = 0; centry < nEntries; centry++){
        if( centry%npace == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " (" << static_cast<float>((10000*centry)/nEntries)/(100) << "%)" << std::endl;
        if(debug) std::cout << "*****************************************************************************" << std::endl;
        auto entry = fInTree->LoadTree(centry);
		if(debug) std::cout << " - getBranches " << std::endl;
		getBranches(entry);
		if(debug) std::cout << " - eventLoop " << std::endl;
		eventLoop(entry);
    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop

    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

	endJobs();
    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }

    for( int it = 0; it < nEBEEMaps; it++ ){
        ebeeMapP[it]->Write(); delete ebeeMapP[it];
        ebeeMapT[it]->Write(); delete ebeeMapT[it];
        ebeeMapR[it]->Write(); delete ebeeMapR[it];
    }//<<>>for( int it = 0; it < nEBEEMaps; it++ )

    fOutFile->Close();
    std::cout << "llpgana_hist_maker : Thats all Folks!!" << std::endl;
}//<<>>void llpgana_hist_maker

//------------------------------------------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------------------------------------------
//
//// ----------------------------------------------- event loop -------------------------------------------------------------------------

void makehists::eventLoop( Long64_t entry ){

	if( DEBUG ) std::cout << " Event Loop Start -- Checking Validations " << std::endl;
	auto goodRunRange = (run >= startRun) && (run <= endRun);
    if( DEBUG ) std::cout << "Run : " << run << " in range " << startRun << " to " << endRun << " Good? " << goodRunRange << std::endl;
    auto goodLocRHs = (*resRhID)[0] != 0 && (*resRhID)[1] != 0;
    auto goodGloRHs = (*resRhID)[2] != 0 && (*resRhID)[3] != 0;
    auto goodLocTimes = (*resRtTime)[0] > -99 && (*resRtTime)[1] > -99;
    auto goodGloTimes = (*resRtTime)[2] > -99 && (*resRtTime)[3] > -99;

    if( DEBUG ) std::cout << " -- Setting Selections " << std::endl;
	auto minLocEnergy = goodLocRHs ? (*resE)[0] > 10.0 && (*resE)[1] > 10.0 : false;
    auto minGloEnergy = goodGloRHs ? (*resE)[2] > 10.0 && (*resE)[3] > 10.0 : false;
	// ""CutSelct"" are pass filters -> if "True" then events will have the stated condition
    auto sieieLocCut = false; //goodLocRHs ? (*phoSigmaIEtaIEta)[0] < 0.013 : false;//less then 
	auto sieieGloCut = false; //goodLocRHs ? (*phoSigmaIEtaIEta)[1] < 0.013 && (*phoSigmaIEtaIEta)[2] < 0.013 : false;//less then 

    auto locCutSelct = true;
	//auto locCutSelct = minLocEnergy;
    //auto locCutSelct = sieieLocCut;

    auto gloCutSelct = true;
	//auto gloCutSelct = minGloEnergy;
    //auto gloCutSelct = sieieGloCut;

	//auto mapLocal = true;
    auto mapLocal = false;
	//auto mapGlobal = not mapLocal;
    auto mapGlobal = false;
	//auto useCali = true;
    //auto useCali = false;

	if( DEBUG ) std::cout << " -- Seting Eta range " << std::endl;
	// --------- eta cuts ------------------------
	auto doEBOnly = true;
    //auto doEBOnly = false;
	auto isLocEB = (DetIDMap[(*resRhID)[0]].ecal == EB) && (DetIDMap[(*resRhID)[1]].ecal == EB);
	auto isGloEB = (DetIDMap[(*resRhID)[2]].ecal == EB) && (DetIDMap[(*resRhID)[3]].ecal == EB);
	auto locEBSelct = doEBOnly ? isLocEB : true;
	auto gloEBSelct = doEBOnly ? isGloEB : true;

	if( goodRunRange ){
		if( DEBUG ) std::cout << "Run : " << run << " in range " << startRun << " to " << endRun << std::endl;

		hist1d[0]->Fill(run);
        hist1d[69]->Fill((event/EVPHR)*100);

		if( DEBUG ) std::cout << " -- Looping rh info " << std::endl;
		for( int it=0; it < (*rhEnergy).size(); it++){

			if( DEBUG ) std::cout << " - Rechit Cali info" << std::endl;
			auto rhIdInfo = DetIDMap[(*rhID)[it]];
			if( doEBOnly && rhIdInfo.ecal != EB ) continue;
            if( DEBUG ) std::cout << " - Rechit Cali map info" << std::endl;
			// icmap lookup only works with EB rhs !!!!!!
			float rhRtCali = ( icmap[0] != NULL ) ? icmap[0]->GetBinContent(rhIdInfo.i2 + 86, rhIdInfo.i1) : 0.f;
			float caliRtTime = (*rhRtTime)[it]-rhRtCali;

			if( DEBUG ) std::cout << " - Rechit loop  1" << std::endl;
        	hist1d[67]->Fill((*rhRtTime)[it]);
            hist1d[70]->Fill(caliRtTime);
            if( DEBUG ) std::cout << " - Rechit loop  1a" << std::endl;
            hist1d[72]->Fill(rhRtCali);
        	hist2d[25]->Fill((*rhEnergy)[it],caliRtTime);
			if( DEBUG ) std::cout << " - Rechit loop  1b" << std::endl;
			if( (*rhEnergy)[it] > 4.0 ){// plot greater then 4 GeV
				(*rhRtisOOT)[it] ? hist1d[127]->Fill(caliRtTime) : hist1d[129]->Fill(caliRtTime);
			}//<<>>if( (*rhEnergy)[it] > 4.0 )

            if( DEBUG ) std::cout << " - Rechit loop  2" << std::endl;
            hist1d[83]->Fill((*rhEnergy)[it]);
            hist1d[84]->Fill((*rhRtisOOT)[it]);
            if( DEBUG ) std::cout << " - Rechit loop  2a" << std::endl;
            hist1d[85]->Fill((*rhisWeird)[it]);
            hist1d[86]->Fill((*rhisDiWeird)[it]);
            hist1d[87]->Fill((*rhSwCross)[it]);
            if( DEBUG ) std::cout << " - Rechit loop  2b" << std::endl;
            hist1d[88]->Fill((*rhisGS6)[it]);
            hist1d[89]->Fill((*rhisGS1)[it]);
            hist1d[90]->Fill((*rhadcToGeV)[it]);
            hist1d[91]->Fill((*rhpedrms12)[it]);

			if( DEBUG ) std::cout << " - Rechit loop  3" << std::endl;

            if( (*rhisGS6)[it] == true && (*rhisGS1)[it] == true ){ 
				hist2d[119]->Fill((*rhEnergy)[it],caliRtTime);
				if( (*rhRtisOOT)[it] == false ) hist2d[118]->Fill((*rhEnergy)[it],caliRtTime);
			}//<<>>if( (*rhisGS6)[it] == true && (*rhisGS1)[it] == true )
            if( (*rhisGS6)[it] == true && (*rhisGS1)[it] == false ){
                hist2d[121]->Fill((*rhEnergy)[it],caliRtTime);
                if( (*rhRtisOOT)[it] == false ) hist2d[120]->Fill((*rhEnergy)[it],caliRtTime);
            }//<<>>if( (*rhisGS6)[it] == true && (*rhisGS1)[it] == true )
            if( (*rhisGS6)[it] == false && (*rhisGS1)[it] ==true ){
                hist2d[123]->Fill((*rhEnergy)[it],caliRtTime);
                if( (*rhRtisOOT)[it] == false ) hist2d[122]->Fill((*rhEnergy)[it],caliRtTime);
            }//<<>>if( (*rhisGS6)[it] == true && (*rhisGS1)[it] == true )
            if( (*rhisGS6)[it] == false && (*rhisGS1)[it] ==false ){
                hist2d[125]->Fill((*rhEnergy)[it],caliRtTime);
                if( (*rhRtisOOT)[it] == false ) hist2d[124]->Fill((*rhEnergy)[it],caliRtTime);
            }//<<>>if( (*rhisGS6)[it] == true && (*rhisGS1)[it] == true )

			if( (*rhEnergy)[it] > 10.0 ){// plot greater then 10 GeV 
				hist2d[108]->Fill((*rhSwCross)[it],caliRtTime); 
				hist2d[117]->Fill((*rhSwCross)[it],(*rhEnergy)[it]);
				auto toposel = not ( (*rhisWeird)[it] || (*rhisDiWeird)[it] );
				if( toposel ){ hist2d[110]->Fill((*rhSwCross)[it],caliRtTime); hist2d[116]->Fill((*rhSwCross)[it],(*rhEnergy)[it]); }
				auto rttimesel = not ( (*rhisWeird)[it] || (*rhisDiWeird)[it] || (*rhRtisOOT)[it] );
            	if( rttimesel ){ hist2d[112]->Fill((*rhSwCross)[it],caliRtTime); hist2d[115]->Fill((*rhSwCross)[it],(*rhEnergy)[it]); }
            }//<<>>if( (*rhEnergy)[it] > 10.0 )

			auto isSpike = (*rhSwCross)[it] > 0.95;
			auto isSpikeRtSel = (*rhisWeird)[it] || (*rhisDiWeird)[it] || (*rhRtisOOT)[it];						
			if( isSpike && (*rhEnergy)[it] > 10.0 ){
				hist1d[131]->Fill((*rhEnergy)[it]);
				if( isSpikeRtSel ) hist1d[132]->Fill((*rhEnergy)[it]);
			}//<<>>if( isSpike )
			if( DEBUG ) std::cout << " - Rechit loop  4" << std::endl;

		}//<<>>for( int it=0; it < (*rhEnergy).size(); it++)

        //------------------------------------------------------------------------------------------------

        if( goodLocRHs && goodLocTimes && locCutSelct && locEBSelct ){
			if( DEBUG ) std::cout << " - Local Phos/seeds Fill  " << std::endl;

			auto dooutl = false; //(*resCCTime)[0] == 0 || (*resCCTime)[1] == 0 ;
            //auto dooutl = true; 
            auto idinfoL0 = DetIDMap[(*resRhID)[0]];
            auto idinfoL1 = DetIDMap[(*resRhID)[1]];
			auto isSRU = (idinfoL0.TT == idinfoL1.TT); // true = same, fasle = different			
            if( dooutl ) std::cout << " Fetching Local Cali values with Rt : " << (*resRtTime)[0] << " & " << (*resRtTime)[1] << std::endl;
            if( dooutl ) std::cout << "  - for : " << (*resRhID)[0] << " && " << (*resRhID)[1] << std::endl;
            if( dooutl ) std::cout << "  - with : " << idinfoL0.i2 << " && " << idinfoL0.i1 << std::endl;
            if( dooutl ) std::cout << "  - with : " << idinfoL1.i2 << " && " << idinfoL1.i1 << std::endl;
            auto caliRtL0 = ( icmap[0] != NULL ) ? icmap[0]->GetBinContent(idinfoL0.i2 + 86, idinfoL0.i1) : 0.f;
            auto caliRtL1 = ( icmap[0] != NULL ) ? icmap[0]->GetBinContent(idinfoL1.i2 + 86, idinfoL1.i1) : 0.f;

			//auto caliSel = caliRtL0 > 0.7 && caliRtL0 < 1.25 && caliRtL1 > 0.7 && caliRtL1 < 1.25;
            auto caliSel = true;
			if( caliSel ){ //------------- Cali Cut ( if/then stament )

	            auto ampl0a = (*resAmp)[0];//  same as resAmp->index(0);
	            auto ampl1a = (*resAmp)[1];
	            auto ampl0 = (*resE)[0];//  same as resAmp->index(0);
	            auto ampl1 = (*resE)[1];
	            auto difAmpL = ampl0 - ampl1;
	            //auto effAmpL = (ampl0*ampl1)/sqrt(sq2(ampl0)+sq2(ampl1));
	            auto effAmpL = (ampl0 + ampl1)/2;
	            auto doeAmpL = difAmpL/effAmpL;
	
	            auto effampl = (ampl0*ampl1)/sqrt(sq2(ampl0)+sq2(ampl1));
	            auto dtrtloc = ((*resRtTime)[0] - caliRtL0 ) - ((*resRtTime)[1] - caliRtL1 );
	
	            auto nlphrh = 0;//((*phoRhIds)[0]).size();
	
				//hist2d[46]->Fill(event/EVPHR,(*phoSigmaIEtaIEta)[0]);
				if( mapLocal ) makeEBEEMaps(0);
	
				//int locIdx( -1 );
				//for( int idx = 0; idx < phoSelType->size(); idx++ ){ if( (*phoSelType)[idx] == 0 ){ locIdx = idx; break; }} 
				//if( locIdx == -1 ) std::cout << " BAD LOCAL PHO INDEX " << std::endl;
	
				if( DEBUG ) std::cout << " - TT0 : " << idinfoL0.TT << " TT1 : " << idinfoL1.TT << std::endl;
				if( isSRU ){
					if( DEBUG ) std::cout << " -- SRU " << std::endl;
	
		            hist1d[47]->Fill(caliRtL0);
		            hist1d[48]->Fill(caliRtL1);
		
		        	hist1d[1]->Fill(doeAmpL);
		        	hist2d[0]->Fill(difAmpL,effAmpL);
		        	hist2d[2]->Fill(ampl0,ampl1);
	
	                hist2d[56]->Fill(effampl,dtrtloc);
	
		            hist1d[39]->Fill((*resRtTime)[0]-caliRtL0);
		            hist1d[40]->Fill((*resRtTime)[1]-caliRtL1);
	                hist1d[45]->Fill((*resRtTime)[0]-caliRtL0-((*resRtTime)[1]-caliRtL1));
		
		            hist1d[55]->Fill((*resTOF)[0]);
		            hist1d[56]->Fill((*resTOF)[1]);
		            hist1d[59]->Fill((*resE)[0]);
		            hist1d[60]->Fill((*resE)[1]);
		            hist1d[63]->Fill(ampl0a);
		            hist1d[64]->Fill(ampl1a);
		
					hist2d[13]->Fill((*resE)[0],ampl0a);
		            hist2d[14]->Fill((*resE)[0],(*resRtTime)[0]-caliRtL0);
		
		            hist2d[16]->Fill((*resE)[1],ampl1a);
		            hist2d[17]->Fill((*resE)[1],(*resRtTime)[1]-caliRtL1);
	
				} else { //<<>>if( isSRU )
					if( DEBUG ) std::cout << " -- DRU " << std::endl;
				
	                hist1d[97]->Fill(caliRtL0);
	                hist1d[98]->Fill(caliRtL1);
	    
	                hist1d[101]->Fill(doeAmpL);
	                hist2d[35]->Fill(difAmpL,effAmpL);
	                hist2d[36]->Fill(ampl0,ampl1);
	    
	                hist2d[57]->Fill(effampl,dtrtloc);
	
	                hist1d[117]->Fill((*resRtTime)[0]-caliRtL0);
	                hist1d[118]->Fill((*resRtTime)[1]-caliRtL1);
	    
	                hist1d[121]->Fill((*resTOF)[0]);
	                hist1d[122]->Fill((*resTOF)[1]);
	                hist1d[123]->Fill((*resE)[0]);
	                hist1d[124]->Fill((*resE)[1]);
	                hist1d[125]->Fill(ampl0a);
	                hist1d[126]->Fill(ampl1a);
	    
	                hist2d[40]->Fill((*resE)[0],ampl0a);
	                hist2d[41]->Fill((*resE)[0],(*resRtTime)[0]-caliRtL0);
	    
	                hist2d[43]->Fill((*resE)[1],ampl1a);
	                hist2d[44]->Fill((*resE)[1],(*resRtTime)[1]-caliRtL1);
	
	            }//<<>>if( isSRU )

			}//<<>>if( caliselect cut //-------------  Cali Cut ( if/then stament )
			//std::cout << " Finished Local" << std::endl;
		}//<<>>if( goodLocRHs )

        if( goodGloRHs && goodGloTimes && gloCutSelct && gloEBSelct ){
			if( DEBUG ) std::cout << " - Global Phos/seeds Fill  " << std::endl;

            int gloIdx0( -1 );
            int gloIdx1( -1 );
            for( int idx = 0; idx < phoSelType->size(); idx++ ){ 
				if( (*phoSelType)[idx] == 1 ) gloIdx0 = idx; if( (*phoSelType)[idx] == 2 ) gloIdx1 = idx; }

            auto dooutg = false; //(*resCCTime)[2] == 0 || (*resCCTime)[3] == 0 ;

			if( DEBUG ) std::cout << " --- glo a0 " << std::endl;
            if( dooutg ) std::cout << " Fetching Global Cali values for CC : " << (*resCCTime)[2] << " & " << (*resCCTime)[3] << std::endl;
            if( dooutg ) std::cout << " Fetching Global Cali values with Rt : " << (*resRtTime)[2] << " & " << (*resRtTime)[3] << std::endl;
            if( dooutg ) std::cout << "  - for : " << (*resRhID)[2] << " && " << (*resRhID)[3] << std::endl;
            auto idinfoG0 = DetIDMap[(*resRhID)[2]];
            auto idinfoG1 = DetIDMap[(*resRhID)[3]];
			if( DEBUG ) std::cout << " --- glo a1 " << std::endl;
            if( dooutg ) std::cout << "  - with : " << idinfoG0.i2 << " && " << idinfoG0.i1 << std::endl;
            if( dooutg ) std::cout << "  - with : " << idinfoG1.i2 << " && " << idinfoG1.i1 << std::endl;
            auto caliRtG0 = ( icmap[0] != NULL ) ? icmap[0]->GetBinContent(idinfoG0.i2 + 86, idinfoG0.i1) : 0.f;// )
            auto caliRtG1 = ( icmap[0] != NULL ) ? icmap[0]->GetBinContent(idinfoG1.i2 + 86, idinfoG1.i1) : 0.f;// )

			if( DEBUG ) std::cout << " --- glo a " << std::endl;

			auto caliSel = true;
            //auto caliSel = false;
            if( caliSel ){ //-------------  CC Cali Cut ( if/then stament )

	        	auto ampg0a = (*resAmp)[2];
	        	auto ampg1a = (*resAmp)[3];
	            auto ampg0 = (*resE)[2];
	            auto ampg1 = (*resE)[3];
	        	auto difAmpG = ampg0 - ampg1;
	        	//auto effAmpG = (ampg0*ampg1)/sqrt(sq2(ampg0)+sq2(ampg1));
	            auto effAmpG = (ampg0 + ampg1)/2;
				auto doeAmpG = difAmpG/effAmpG;
				auto effampg = (ampg0*ampg1)/sqrt(sq2(ampg0)+sq2(ampg1)); 
	            auto dtrtglo = ((*resRtTime)[2] - caliRtG0 ) - ((*resRtTime)[3] - caliRtG1 );
	
	            auto ngphrh0 = ((*phoRhIds)[gloIdx0]).size();
	            auto ngphrh1 = ((*phoRhIds)[gloIdx1]).size();
	
	            if( mapGlobal ) makeEBEEMaps(1);
	            if( mapGlobal ) makeEBEEMaps(2);
	
				if( DEBUG ) std::cout << " --- glo b " << std::endl;
	
				hist1d[2]->Fill(doeAmpG);
	        	hist2d[1]->Fill(difAmpG,effAmpG);
	        	hist2d[3]->Fill(ampg0,ampg1);
	
	            hist2d[55]->Fill(effampg,dtrtglo);
	
	        	hist1d[3]->Fill(phoDiMass);
	        	hist1d[4]->Fill(phoDiAngle);
	        	hist1d[5]->Fill(phoDiDr);
	            hist1d[92]->Fill(phoDiEta);
	            hist1d[93]->Fill(phoDiPhi);
	
				hist1d[32]->Fill((*phoPt)[gloIdx0],(*phoPt)[gloIdx1]);
	            hist1d[33]->Fill((*phoEnergy)[gloIdx0],(*phoEnergy)[gloIdx1]);
	
	            hist1d[17]->Fill((*phoEnergy)[gloIdx0]);
	            hist1d[18]->Fill((*phoPt)[gloIdx0]);
	            hist1d[19]->Fill((*phoEta)[gloIdx0]);
	            hist1d[20]->Fill((*phoPhi)[gloIdx0]);
	            hist1d[21]->Fill((*phoHadOverEM)[gloIdx0]);
	            hist1d[22]->Fill((*phoSigmaIEtaIEta)[gloIdx0]);
	            hist1d[75]->Fill((*phoCov2IEtaIEta)[gloIdx0]);
	            hist1d[78]->Fill((*phoCov2IEtaIPhi)[gloIdx0]);
	            hist1d[81]->Fill((*phoCov2IPhiIPhi)[gloIdx0]);
	            hist1d[23]->Fill((*phoEcalRHSumEtConeDR04)[gloIdx0]);
	            hist1d[24]->Fill((*phoHcalTwrSumEtConeDR04)[gloIdx0]);
	            hist1d[25]->Fill((*phoTrkSumPtSolidConeDR04)[gloIdx0]);
	            hist1d[26]->Fill((*phoTrkSumPtHollowConeDR04)[gloIdx0]);
	            hist1d[27]->Fill((*phoR9)[gloIdx0]);
	            hist1d[95]->Fill(ngphrh0);
	    		hist2d[28]->Fill(phoDiMass,(*phoEnergy)[gloIdx0]);
	            hist2d[30]->Fill(phoDiMass,(*phoPt)[gloIdx0]);
	
				if( DEBUG ) std::cout << " --- glo c " << std::endl;
	
				hist2d[4]->Fill((*phoEnergy)[gloIdx0],(*phoHadOverEM)[gloIdx0]);
	            hist2d[7]->Fill((*phoEnergy)[gloIdx0],(*phoSigmaIEtaIEta)[gloIdx0]);
	            hist2d[10]->Fill((*phoEnergy)[gloIdx0],(*phoR9)[gloIdx0]);
	            hist2d[49]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoR9)[gloIdx0]);
	
				hist2d[51]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoSigmaIEtaIEta)[gloIdx1]);
	
	            hist1d[28]->Fill((*phoEnergy)[gloIdx1]);
	            hist1d[29]->Fill((*phoPt)[gloIdx1]);
	            hist1d[30]->Fill((*phoEta)[gloIdx1]);
	            hist1d[31]->Fill((*phoPhi)[gloIdx1]);
	            hist1d[32]->Fill((*phoHadOverEM)[gloIdx1]);
	            hist1d[33]->Fill((*phoSigmaIEtaIEta)[gloIdx1]);
	            hist1d[76]->Fill((*phoCov2IEtaIEta)[gloIdx1]);
	            hist1d[79]->Fill((*phoCov2IEtaIPhi)[gloIdx1]);
	            hist1d[82]->Fill((*phoCov2IPhiIPhi)[gloIdx1]);
	            hist1d[34]->Fill((*phoEcalRHSumEtConeDR04)[gloIdx1]);
	            hist1d[35]->Fill((*phoHcalTwrSumEtConeDR04)[gloIdx1]);
	            hist1d[36]->Fill((*phoTrkSumPtSolidConeDR04)[gloIdx1]);
	            hist1d[37]->Fill((*phoTrkSumPtHollowConeDR04)[gloIdx1]);
	            hist1d[38]->Fill((*phoR9)[gloIdx1]);
	            hist1d[96]->Fill(ngphrh1);
	            hist2d[29]->Fill(phoDiMass,(*phoEnergy)[gloIdx1]);
	            hist2d[31]->Fill(phoDiMass,(*phoPt)[gloIdx1]);
	
	            hist2d[32]->Fill((*phoPt)[gloIdx0],(*phoPt)[gloIdx1]);
	            hist2d[33]->Fill((*phoEnergy)[gloIdx0],(*phoEnergy)[gloIdx1]);
	            hist2d[34]->Fill(phoDiMass,phoDiDr);
	
	            hist2d[5]->Fill((*phoEnergy)[gloIdx1],(*phoHadOverEM)[gloIdx1]);
	            hist2d[8]->Fill((*phoEnergy)[gloIdx1],(*phoSigmaIEtaIEta)[gloIdx1]);
	            hist2d[11]->Fill((*phoEnergy)[gloIdx1],(*phoR9)[gloIdx1]);
	            hist2d[50]->Fill((*phoSigmaIEtaIEta)[gloIdx1],(*phoR9)[gloIdx1]);
	
				if( DEBUG ) std::cout << " --- glo d " << std::endl;
	
	            hist1d[43]->Fill((*resRtTime)[2]-caliRtG0);
	            hist1d[44]->Fill((*resRtTime)[3]-caliRtG1);
	            hist1d[46]->Fill((*resRtTime)[2]-caliRtG0-((*resRtTime)[3]-caliRtG1));
	
	            hist1d[51]->Fill(caliRtG0);
	            hist1d[52]->Fill(caliRtG1);
	
	            hist1d[57]->Fill((*resTOF)[2]);
	            hist1d[58]->Fill((*resTOF)[3]);
	            hist1d[61]->Fill((*resE)[2]);
	            hist1d[62]->Fill((*resE)[3]);
	            hist1d[65]->Fill(ampg0a);
	            hist1d[66]->Fill(ampg1a);
	
	            hist2d[19]->Fill((*resE)[2],ampg0a);
	            hist2d[20]->Fill((*resE)[2],(*resRtTime)[2]-caliRtG0);
	
	            hist2d[22]->Fill((*resE)[3],ampg1a);
	            hist2d[23]->Fill((*resE)[3],(*resRtTime)[3]-caliRtG1);
	
	            hist2d[90]->Fill((*phoEnergy)[gloIdx0],(*phoHadOverEM)[gloIdx0]);
	            hist2d[91]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoHadOverEM)[gloIdx0]);
	            hist2d[92]->Fill((*phoEnergy)[gloIdx0],(*phoCov2IEtaIEta)[gloIdx0]);
	            hist2d[93]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoCov2IEtaIEta)[gloIdx0]);
	            hist2d[94]->Fill((*phoEnergy)[gloIdx0],(*phoCov2IEtaIPhi)[gloIdx0]);
	            hist2d[95]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoCov2IEtaIPhi)[gloIdx0]);
	            hist2d[96]->Fill((*phoEnergy)[gloIdx0],(*phoCov2IPhiIPhi)[gloIdx0]);
	            hist2d[97]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoCov2IPhiIPhi)[gloIdx0]);
	            hist2d[98]->Fill((*phoEnergy)[gloIdx0],(*phoEcalRHSumEtConeDR04)[gloIdx0]);
	            hist2d[99]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoEcalRHSumEtConeDR04)[gloIdx0]);
	            hist2d[100]->Fill((*phoEnergy)[gloIdx0],(*phoHcalTwrSumEtConeDR04)[gloIdx0]);
	            hist2d[101]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoHcalTwrSumEtConeDR04)[gloIdx0]);
	            hist2d[102]->Fill((*phoEnergy)[gloIdx0],(*phoTrkSumPtSolidConeDR04)[gloIdx0]);
	            hist2d[103]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoTrkSumPtSolidConeDR04)[gloIdx0]);
	            hist2d[104]->Fill((*phoEnergy)[gloIdx0],(*phoTrkSumPtHollowConeDR04)[gloIdx0]);
	            hist2d[105]->Fill((*phoSigmaIEtaIEta)[gloIdx0],(*phoTrkSumPtHollowConeDR04)[gloIdx0]);

			}//-------------  Cali Cut ( if/then stament )
            //std::cout << " Finished Global" << std::endl;
		}//<<>>if( goodGloRHs )
	}//<<>>if( run > startRun && run < endRun )
}//<<>>void makehists::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void makehists::initCali( std::string califilename ){

    int nAlgos(1);
	bool setnull = false;
	if( califilename == "none" ) setnull = true;
	auto fCaliFile = ( not setnull ) ? TFile::Open(califilename.c_str(), "read") : NULL;
    std::string calistring[nAlgos] = { "AveXtalRecHitTime"};	
	for( int i = 0; i < nAlgos; i++ ){
		std::string cmbs = calistring[i];
    	std::string ebmapstring(cmbs+"EBMap");
    	std::string epmapstring(cmbs+"EPMap");
    	std::string emmapstring(cmbs+"EMMap");
    	icmap[i*3] = setnull ? NULL : (TH2F*)fCaliFile->Get(ebmapstring.c_str());
    	icmap[i*3+1] = setnull ? NULL : (TH2F*)fCaliFile->Get(epmapstring.c_str());
    	icmap[i*3+2] = setnull ? NULL : (TH2F*)fCaliFile->Get(emmapstring.c_str());
	}//<<>>for( int i = 0; i < 2; i++ )
	std::cout << "Getting Cali maps : " << icmap[0] << " && " << icmap[1] << " && " << icmap[2] << std::endl;

}//<<>>void makehists::initCali( TFile* fCaliFile )

void makehists::getBranches( Long64_t entry ){

    b_run->GetEntry(entry);   //!
    b_lumi->GetEntry(entry);   //!
    b_event->GetEntry(entry);   //!

	if( DEBUG ) std::cout << " -- getting Banches rhCali " << std::endl;
    b_rhCaliID->GetEntry(entry);   //!
    b_rhCaliEnergy->GetEntry(entry);   //!
    b_rhCaliRtTime->GetEntry(entry);   //!

    if( DEBUG ) std::cout << " -- getting Banches res " << std::endl;
    b_resRhID->GetEntry(entry);   //!
    b_resAmp->GetEntry(entry);   //!
    b_resE->GetEntry(entry);   //!
    b_resRtTime->GetEntry(entry);   //!
    b_resTOF->GetEntry(entry);   //!

    if( DEBUG ) std::cout << " -- getting Banches rh " << std::endl;
    b_rhEnergy->GetEntry(entry);   //!
    b_rhID->GetEntry(entry);   //!
    if( DEBUG ) std::cout << " -- getting Banches rh 1" << std::endl;
    b_rhRtTime->GetEntry(entry);   //!
    b_rhCCTime->GetEntry(entry);   //!
    if( DEBUG ) std::cout << " -- getting Banches rh koot stuff" << std::endl;
    b_rhRtisOOT->GetEntry(entry);   //!
    b_rhisWeird->GetEntry(entry);   //!
    b_rhisDiWeird->GetEntry(entry);   //!
    b_rhSwCross->GetEntry(entry);   //!
    if( DEBUG ) std::cout << " -- getting Banches rh amp stuff" << std::endl;
    b_rhisGS6->GetEntry(entry);
    b_rhisGS1->GetEntry(entry);
    b_rhadcToGeV->GetEntry(entry);
    b_rhpedrms12->GetEntry(entry);


    if( DEBUG ) std::cout << " --- getting Banches pho " << std::endl;
    b_phoEnergy->GetEntry(entry);   //!
    b_phoRhIds->GetEntry(entry);   //!
    b_phoPt->GetEntry(entry);   //!
    b_phoEta->GetEntry(entry);   //!
    b_phoPhi->GetEntry(entry);   //!
    b_phoHadOverEM->GetEntry(entry);   //!
    b_phoSigmaIEtaIEta->GetEntry(entry);   //!

    if( DEBUG ) std::cout << " --- getting Banches pho 1" << std::endl;
    b_phoCov2IEtaIEta->GetEntry(entry);   //!
    b_phoCov2IEtaIPhi->GetEntry(entry);   //!
    b_phoCov2IPhiIPhi->GetEntry(entry);   //!
    if( DEBUG ) std::cout << " --- getting Banches pho 2" << std::endl;
    b_phoEcalRHSumEtConeDR04->GetEntry(entry);   //!
    b_phoHcalTwrSumEtConeDR04->GetEntry(entry);   //!
    b_phoTrkSumPtSolidConeDR04->GetEntry(entry);   //!
    b_phoTrkSumPtHollowConeDR04->GetEntry(entry);   //!
    b_phoR9->GetEntry(entry);   //!
    b_phoSelType->GetEntry(entry);   //!

    if( DEBUG ) std::cout << " --- getting Banches pho 3" << std::endl;
    b_phoDiMass->GetEntry(entry);   //!
    b_phoDiAngle->GetEntry(entry);   //!
    b_phoDiDr->GetEntry(entry);   //!
    b_phoDiPhi->GetEntry(entry);   //!
    b_phoDiEta->GetEntry(entry);   //!

}//<<>>void makehists::getBranches( Long64_t entry )

void makehists::endJobs(){

    if( DEBUG ) std::cout << " Starting End jobs " << std::endl;

	fillRatioHist(hist1d[132],hist1d[131],hist1d[133]);

}//<<>>void makehists::endJobs()

void makehists::initHists( std::string fHTitle ){

    totrhs = 0;
    totrhs0 = 0;
    totrhs05 = 0;
    totrhs1 = 0;
    totrhs2 = 0;
    totrhs5 = 0;
    totrhs10 = 0;
    encrhs = 0;

    nMaps = 0;
    for(int it=0; it<nEBEEMaps; it++){
        fMap[it] = false;
        std::string label(";iEta;iPhi");
        std::string stt1("ebeeMapPhoCluster_"+std::to_string(it));
        ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
        std::string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
        ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
        std::string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
        ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
    }//<<>>for(int it=0; it<nEBEEMaps; it++)

    for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < nHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

    //------ 1D Hists --------------------------------------------------------------------------

	hist1d[0] = new TH1D("run",addstr(fHTitle,"Run;Run").c_str(),88000,300000,400000); 
	hist1d[69] = new TH1D("event",addstr(fHTitle,"Event;Run-Evt??").c_str(),10000,0,10000);

    hist1d[1] = new TH1D("difampoeffamp_sruloc",addstr(fHTitle,"SRU Local Diff E / Ave E; DiffE/AveE").c_str(),300,-0.3,0.3);
    hist1d[2] = new TH1D("difampoeffamp_glo",addstr(fHTitle,"Global Diff E / Ave E; DiffE/AveE").c_str(),300,-3,3);

    hist1d[3] = new TH1D("phoDiMass",addstr(fHTitle,"phoDiMass").c_str(),350,55,125);
    hist1d[4] = new TH1D("phoDiAngle",addstr(fHTitle,"phoDiAngle").c_str(),335,-0.1,3.25);
    hist1d[5] = new TH1D("phoDiDr",addstr(fHTitle,"phoDiDr").c_str(),300,0,6);

    hist1d[17] = new TH1D("phoEnergy_Glo0",addstr(fHTitle,"phoEnergy_Glo0").c_str(),1000,0,1000);
    hist1d[28] = new TH1D("phoEnergy_Glo1",addstr(fHTitle,"phoEnergy_Glo1").c_str(),1000,0,1000);

    hist1d[18] = new TH1D("phoPt_Glo0",addstr(fHTitle,"phoPt_Glo0").c_str(),500,0,500);
    hist1d[29] = new TH1D("phoPt_Glo1",addstr(fHTitle,"phoPt_Glo1").c_str(),500,0,500);

    hist1d[19] = new TH1D("phoEta_Glo0",addstr(fHTitle,"phoEta_Glo0").c_str(),700,-3.5,3.5);
    hist1d[30] = new TH1D("phoEta_Glo1",addstr(fHTitle,"phoEta_Glo1").c_str(),700,-3.5,3.5);

    hist1d[20] = new TH1D("phoPhi_Glo0",addstr(fHTitle,"phoPhi_Glo0").c_str(),700,-3.5,3.5);
    hist1d[31] = new TH1D("phoPhi_Glo1",addstr(fHTitle,"phoPhi_Glo1").c_str(),700,-3.5,3.5);

    hist1d[21] = new TH1D("phoHadOverEM_Glo0",addstr(fHTitle,"phoHadOverEM_Glo0").c_str(),100,0,1.0);
    hist1d[32] = new TH1D("phoHadOverEM_Glo1",addstr(fHTitle,"phoHadOverEM_Glo1").c_str(),100,0,1.0);

    hist1d[22] = new TH1D("phoSigmaIEtaIEta_Glo0",addstr(fHTitle,"phoSigmaIEtaIEta_Glo0").c_str(),1000,0,0.1);
    hist1d[33] = new TH1D("phoSigmaIEtaIEta_Glo1",addstr(fHTitle,"phoSigmaIEtaIEta_Glo1").c_str(),1000,0,0.1);

    hist1d[23] = new TH1D("phoEcalRHSumEtConeDR04_Glo0",addstr(fHTitle,"phoEcalRHSumEtConeDR04_Glo0").c_str(),350,0,350);
    hist1d[34] = new TH1D("phoEcalRHSumEtConeDR04_Glo1",addstr(fHTitle,"phoEcalRHSumEtConeDR04_Glo1").c_str(),350,0,350);

    hist1d[24] = new TH1D("phoHcalTwrSumEtConeDR04_Glo0",addstr(fHTitle,"phoHcalTwrSumEtConeDR04_Glo0").c_str(),350,0,350);
    hist1d[35] = new TH1D("phoHcalTwrSumEtConeDR04_Glo1",addstr(fHTitle,"phoHcalTwrSumEtConeDR04_Glo1").c_str(),350,0,350);

    hist1d[25] = new TH1D("phoTrkSumPtSolidConeDR04_Glo0",addstr(fHTitle,"phoTrkSumPtSolidConeDR04_Glo0").c_str(),350,0,350);
    hist1d[36] = new TH1D("phoTrkSumPtSolidConeDR04_Glo1",addstr(fHTitle,"phoTrkSumPtSolidConeDR04_Glo1").c_str(),350,0,350);

    hist1d[26] = new TH1D("phoTrkSumPtHollowConeDR04_Glo0",addstr(fHTitle,"phoTrkSumPtHollowConeDR04_Glo0").c_str(),350,0,350);
    hist1d[37] = new TH1D("phoTrkSumPtHollowConeDR04_Glo1",addstr(fHTitle,"phoTrkSumPtHollowConeDR04_Glo1").c_str(),350,0,350);

    hist1d[27] = new TH1D("phoR9_Glo0",addstr(fHTitle,"phoR9_Glo0").c_str(),100,0,1);
    hist1d[38] = new TH1D("phoR9_Glo1",addstr(fHTitle,"phoR9_Glo1").c_str(),100,0,1);

	hist1d[39] = new TH1D("seedRtTime_sruLoc0",addstr(fHTitle,"seedRtTime SRU Loc0").c_str(),500,-25,25);
    hist1d[40] = new TH1D("seedRtTime_sruLoc1",addstr(fHTitle,"seedRtTime SRU Loc1").c_str(),500,-25,25);
    hist1d[43] = new TH1D("seedRtTime_Glo0",addstr(fHTitle,"seedRtTime_Glo0").c_str(),500,-25,25);
    hist1d[44] = new TH1D("seedRtTime_Glo1",addstr(fHTitle,"seedRtTime_Glo1").c_str(),500,-25,25);

    hist1d[45] = new TH1D("seedRtTimeDiff_sruLoc",addstr(fHTitle,"seedRtTimeDiff SRU Loc").c_str(),500,-25,25);
    hist1d[46] = new TH1D("seedRtTimeDiff_Glo",addstr(fHTitle,"seedRtTimeDiff_Glo").c_str(),500,-25,25);

    hist1d[47] = new TH1D("seedRtCali_sruLoc0",addstr(fHTitle,"seedRtCali SRU Loc0").c_str(),1000,-5,5);
    hist1d[48] = new TH1D("seedRtCali_sruLoc1",addstr(fHTitle,"seedRtCali SRU Loc1").c_str(),1000,-5,5);
    hist1d[51] = new TH1D("seedRtCali_Glo0",addstr(fHTitle,"seedRtCali_Glo0").c_str(),1000,-5,5);
    hist1d[52] = new TH1D("seedRtCali_Glo1",addstr(fHTitle,"seedRtCali_Glo1").c_str(),1000,-5,5);

    hist1d[55] = new TH1D("seedTOF_sruLoc0",addstr(fHTitle,"seedTOF SRU Loc0").c_str(),250,-1.25,1.25);
    hist1d[56] = new TH1D("seedTOF_sruLoc1",addstr(fHTitle,"seedTOF SRU Loc1").c_str(),250,-1.25,1.25);
    hist1d[57] = new TH1D("seedTOF_Glo0",addstr(fHTitle,"seedTOF_Glo0").c_str(),250,-1.25,1.25);
    hist1d[58] = new TH1D("seedTOF_Glo1",addstr(fHTitle,"seedTOF_Glo1").c_str(),250,-1.25,1.25);

    hist1d[59] = new TH1D("seedEnergy_sruLoc0",addstr(fHTitle,"seedEnergy SRU Loc0").c_str(),500,0,500);
    hist1d[60] = new TH1D("seedEnergy_sruLoc1",addstr(fHTitle,"seedEnergy SRU Loc1").c_str(),500,0,500);
    hist1d[61] = new TH1D("seedEnergy_Glo0",addstr(fHTitle,"seedEnergy_Glo0").c_str(),500,0,500);
    hist1d[62] = new TH1D("seedEnergy_Glo1",addstr(fHTitle,"seedEnergy_Glo1").c_str(),500,0,500);

    hist1d[63] = new TH1D("seedAmplitude_sruLoc0",addstr(fHTitle,"seedAmplitude SRU Loc0").c_str(),1000,0,1000);
    hist1d[64] = new TH1D("seedAmplitude_sruLoc1",addstr(fHTitle,"seedAmplitude SRU Loc1").c_str(),1000,0,1000);
    hist1d[65] = new TH1D("seedAmplitude_Glo0",addstr(fHTitle,"seedAmplitude_Glo0").c_str(),1000,0,1000);
    hist1d[66] = new TH1D("seedAmplitude_Glo1",addstr(fHTitle,"seedAmplitude_Glo1").c_str(),1000,0,1000);

    hist1d[67] = new TH1D("rhRHTimeUnCali",addstr(fHTitle,"rhRHTimeUnCali").c_str(),500,-25,25);
    hist1d[70] = new TH1D("rhRtTimeCali",addstr(fHTitle,"rhRtTimeCali").c_str(),500,-25,25);
    hist1d[72] = new TH1D("rhRtCali",addstr(fHTitle,"rhRtCali").c_str(),500,-25,25);

    hist1d[75] = new TH1D("phoCov2IEtaIEta_Glo0",addstr(fHTitle,"phoCov2IEtaIEta_Glo0").c_str(),1000,0,0.005);
    hist1d[76] = new TH1D("phoCov2IEtaIEta_Glo1",addstr(fHTitle,"phoCov2IEtaIEta_Glo1").c_str(),1000,0,0.005);

    hist1d[78] = new TH1D("phoCov2IEtaIPhi_Glo0",addstr(fHTitle,"phoCov2IEtaIPhi_Glo0").c_str(),1000,0,0.005);
    hist1d[79] = new TH1D("phoCov2IEtaIPhi_Glo1",addstr(fHTitle,"phoCov2IEtaIPhi_Glo1").c_str(),1000,0,0.005);

    hist1d[81] = new TH1D("phoCov2IPhiIPhi_Glo0",addstr(fHTitle,"phoCov2IPhiIPhi_Glo0").c_str(),1000,0,0.005);
    hist1d[82] = new TH1D("phoCov2IPhiIPhi_Glo1",addstr(fHTitle,"phoCov2IPhiIPhi_Glo1").c_str(),1000,0,0.005);

    hist1d[83] = new TH1D("rhEnergy",addstr(fHTitle,"rhEnergy").c_str(),1500,0,1500);
    hist1d[84] = new TH1D("rhRHisoot",addstr(fHTitle,"rhRtisOOT").c_str(),3,0,2);
    hist1d[85] = new TH1D("rhisweird",addstr(fHTitle,"rhisWeird").c_str(),3,0,2);
    hist1d[86] = new TH1D("rhisdiweird",addstr(fHTitle,"rhisDiWeird").c_str(),3,0,2);
    hist1d[87] = new TH1D("rhSwissCross",addstr(fHTitle,"rhSwissCross").c_str(),15200,-150,2);
    hist1d[88] = new TH1D("rhisgs6",addstr(fHTitle,"rhisGS6").c_str(),3,0,2);
    hist1d[89] = new TH1D("rhisgs1",addstr(fHTitle,"rhisGS1").c_str(),3,0,2);
    hist1d[90] = new TH1D("rhadctogev",addstr(fHTitle,"rhadcToGev").c_str(),500,0,5);
	hist1d[91] = new TH1D("rhpedrms12",addstr(fHTitle,"rhpedRMS12").c_str(),1000,0,10);

    hist1d[92] = new TH1D("phoDiEta",addstr(fHTitle,"phoDiEta").c_str(),180,-0.1,3.5);
    hist1d[93] = new TH1D("phoDiPhi",addstr(fHTitle,"phoDiPhi").c_str(),350,-3.5,3.5);

    hist1d[95] = new TH1D("phonrh_Glo0",addstr(fHTitle,"pho nrh Glo0").c_str(),100,0,100);
    hist1d[96] = new TH1D("phonrh_Glo1",addstr(fHTitle,"pho nrh Glo1").c_str(),100,0,100);

    hist1d[97] = new TH1D("seedRHCali_druLoc0",addstr(fHTitle,"seedRHCali DRU Loc0").c_str(),1000,-5,5);
    hist1d[98] = new TH1D("seedRHCali_druLoc1",addstr(fHTitle,"seedRHCali DRU Loc1").c_str(),1000,-5,5);

	hist1d[101] = new TH1D("difampoeffamp_druloc",addstr(fHTitle,"DRU Local Diff E / Ave E; DiffE/AveE").c_str(),300,-0.3,0.3);

    hist1d[117] = new TH1D("seedRHTime_druLoc0",addstr(fHTitle,"seedRHTime DRU Loc0").c_str(),500,-25,25);
    hist1d[118] = new TH1D("seedRHTime_druLoc1",addstr(fHTitle,"seedRHTime DRU Loc1").c_str(),500,-25,25);

    hist1d[121] = new TH1D("seedTOF_druLoc0",addstr(fHTitle,"seedTOF DRU Loc0").c_str(),250,-1.25,1.25);
    hist1d[122] = new TH1D("seedTOF_druLoc1",addstr(fHTitle,"seedTOF DRU Loc1").c_str(),250,-1.25,1.25);

    hist1d[123] = new TH1D("seedEnergy_druLoc0",addstr(fHTitle,"seedEnergy DRU Loc0").c_str(),500,0,500);
    hist1d[124] = new TH1D("seedEnergy_druLoc1",addstr(fHTitle,"seedEnergy DRU Loc1").c_str(),500,0,500);
    hist1d[125] = new TH1D("seedAmplitude_druLoc0",addstr(fHTitle,"seedAmplitude DRU Loc0").c_str(),1000,0,1000);
    hist1d[126] = new TH1D("seedAmplitude_druLoc1",addstr(fHTitle,"seedAmplitude DRU Loc1").c_str(),1000,0,1000);

    hist1d[127] = new TH1D("rhRHTimeOOT1UnCali",addstr(fHTitle,"kOOT True rhRHTimeUnCali").c_str(),500,-25,25);
    hist1d[129] = new TH1D("rhRHTimeOOT0UnCali",addstr(fHTitle,"kOOT False rhRHTimeUnCali").c_str(),500,-25,25);

	hist1d[131] = new TH1D("spikes",addstr(fHTitle,"spikes").c_str(),20,0,1000);
    hist1d[132] = new TH1D("spikesRHSel",addstr(fHTitle,"spikesRHSelect").c_str(),20,0,1000);
    hist1d[133] = new TH1D("spikesRHEff",addstr(fHTitle,"spikesRHEff").c_str(),20,0,1000);

	for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

    //------ 2D Hists --------------------------------------------------------------------------
    hist2d[0] = new TH2D("srudifampveffampl",addstr(fHTitle,"SRU Local Diff E V Ave E; DiffE; AveE").c_str(),101,-1,100,100,0,100);
    hist2d[1] = new TH2D("difampveffampg",addstr(fHTitle,"Global Diff E V Ave E; DiffE; AveE").c_str(),450,-450,450,500,0,500);	
    hist2d[2] = new TH2D("sruamp0vamp1l",addstr(fHTitle,"SRU Local E 0 V E 1;E0;E1").c_str(),400,0,400,400,0,400);
    hist2d[3] = new TH2D("amp0vamp1g",addstr(fHTitle,"Global E 0 V E 1;E0;E1").c_str(),400,0,400,400,0,400);

    hist2d[4] = new TH2D("phoEnergyVHadOverEM_Glo0",addstr(fHTitle,"phoEnergy V HadOverEM Glo0;E;HOEM").c_str(),500,0,250,100,0,0.5);
    hist2d[5] = new TH2D("phoEnergyVHadOverEM_Glo1",addstr(fHTitle,"phoEnergy V HadOverEM Glo1;E;HOEM").c_str(),500,0,250,100,0,0.5);

    hist2d[7] = new TH2D("phoEnergyVSigIEtaIEta_Glo0",addstr(fHTitle,"phoEnergy V SigmaIEtaIEta Glo0;E;sIEIE").c_str(),500,0,250,200,0,0.05);
    hist2d[8] = new TH2D("phoEnergyVSigIEtaIEta_Glo1",addstr(fHTitle,"phoEnergy V SigmaIEtaIEta Glo1;E;sIEIE").c_str(),500,0,250,200,0,0.05);

    hist2d[10] = new TH2D("phoEnergyVR9_Glo0",addstr(fHTitle,"phoEnergy V R9 Glo0;E;R9").c_str(),500,0,250,100,0,1);
    hist2d[11] = new TH2D("phoEnergyVR9_Glo1",addstr(fHTitle,"phoEnergy V R9 Glo1;E;R9").c_str(),500,0,250,100,0,1);
    hist2d[12] = new TH2D("phoEnergyVR9_sruLoc",addstr(fHTitle,"phoEnergy V R9 SRU Loc;E;R9").c_str(),500,0,250,100,0,1);

    hist2d[13] = new TH2D("scEnergyVAmp_sruLoc0",addstr(fHTitle,"scEnergy V Amp SRU Loc0;E;Amp").c_str(),500,0,250,1000,0,1000);
    hist2d[14] = new TH2D("scEnergyVRtTime_sruLoc0",addstr(fHTitle,"scEnergy V RtTime SRU Loc0;E;RtTime").c_str(),500,0,250,600,-15,15);

    hist2d[16] = new TH2D("scEnergyVAmp_sruLoc1",addstr(fHTitle,"scEnergy V Amp SRU Loc1;E;Amp").c_str(),500,0,250,1000,0,1000);
    hist2d[17] = new TH2D("scEnergyVRtTime_sruLoc1",addstr(fHTitle,"scEnergy V RtTime SRU Loc1;E;RtTime").c_str(),500,0,250,600,-15,15);

    hist2d[19] = new TH2D("scEnergyVAmp_Glo0",addstr(fHTitle,"scEnergy V Amp Glo0;E;Amp").c_str(),500,0,250,1000,0,1000);
    hist2d[20] = new TH2D("scEnergyVRtTime_Glo0",addstr(fHTitle,"scEnergy V RtTime Glo0;E;RtTime").c_str(),500,0,250,600,-15,15);

    hist2d[22] = new TH2D("scEnergyVAmp_Glo1",addstr(fHTitle,"scEnergy V Amp Glo1;E;Amp").c_str(),500,0,250,1000,0,1000);
    hist2d[23] = new TH2D("scEnergyVRtTime_Glo1",addstr(fHTitle,"scEnergy V RtTime Glo1;E;RtTime").c_str(),500,0,250,600,-15,15);

    hist2d[25] = new TH2D("rhEnergyVrhRtTime",addstr(fHTitle,"rhEnergy V rhRtTime;E;RtTime").c_str(),500,0,250,600,-15,15);

    hist2d[28] = new TH2D("phoDiMassvPt0",addstr(fHTitle,"phoDiMass V Pt Pho0;Mass;Pt").c_str(),140,55,125,500,0,500);
    hist2d[29] = new TH2D("phoDiMassvPt1",addstr(fHTitle,"phoDiMass V Pt Pho1;Mass;Pt").c_str(),140,55,125,500,0,500);
    hist2d[30] = new TH2D("phoDiMassvE0",addstr(fHTitle,"phoDiMass V E Pho0;Mass;E").c_str(),140,55,125,250,0,250);
    hist2d[31] = new TH2D("phoDiMassvE1",addstr(fHTitle,"phoDiMass V E Pho1;Mass;E").c_str(),140,55,125,250,0,250);
    hist2d[32] = new TH2D("phoPt0vPt1",addstr(fHTitle,"Pt Pho0 V Pt Pho1;Pt pho0;Pt pho1").c_str(),200,0,200,200,0,200);
    hist2d[33] = new TH2D("phoE0vE1",addstr(fHTitle,"E Pho0 V E Pho1;E pho0;E pho1").c_str(),500,0,250,250,0,250);
    hist2d[34] = new TH2D("phoDiMassvphoDiDr",addstr(fHTitle,"phoDiMass V phoDiDr;Mass;Dr").c_str(),140,55,125,500,0,5);

    hist2d[35] = new TH2D("drudifampveffampl",addstr(fHTitle,"DRU Local Diff E V Ave E; DiffE; AveE").c_str(),420,-1,20,100,0,100);
    hist2d[36] = new TH2D("druamp0vamp1l",addstr(fHTitle,"DRU Local E 0 V E 1;E0;E1").c_str(),500,0,250,250,0,250);

    hist2d[40] = new TH2D("scEnergyVAmp_druLoc0",addstr(fHTitle,"scEnergy V Amp DRU Loc0;E;Amp").c_str(),500,0,250,1000,0,1000);
    hist2d[41] = new TH2D("scEnergyVRtTime_druLoc0",addstr(fHTitle,"scEnergy V RtTime DRU Loc0;E;RtTime").c_str(),500,0,250,600,-15,15);

    hist2d[43] = new TH2D("scEnergyVAmp_druLoc1",addstr(fHTitle,"scEnergy V Amp DRU Loc1;E;Amp").c_str(),500,0,250,1000,0,1000);
    hist2d[44] = new TH2D("scEnergyVRtTime_druLoc1",addstr(fHTitle,"scEnergy V RtTime DRU Loc1;E;RtTime").c_str(),500,0,250,600,-15,15);

    hist2d[46] = new TH2D("evntvsieie",addstr(fHTitle,"Evnt V phoSigmaIEtaIEta;Event;Sieie").c_str(),24000,0,2400,400,0,0.04);
    hist2d[47] = new TH2D("sieievr9sruloc",addstr(fHTitle,"SigmaIEtaIeta V R9 SRU Loc;Sieie;r9").c_str(),400,0,0.04,100,0,1);
    hist2d[48] = new TH2D("sieievr9druloc",addstr(fHTitle,"SigmaIEtaIeta V R9 DRU Loc;Sieie;r9").c_str(),400,0,0.04,100,0,1);
    hist2d[49] = new TH2D("sieievr9glo0",addstr(fHTitle,"SigmaIEtaIeta V R9 Glo0;Sieie;r9").c_str(),400,0,0.04,100,0,1);
    hist2d[50] = new TH2D("sieievr9glo1",addstr(fHTitle,"SigmaIEtaIeta V R9 Glo1;Sieie;r9").c_str(),400,0,0.04,100,0,1);
    hist2d[51] = new TH2D("sieieglo1v2",addstr(fHTitle,"SigmaIEtaIeta Glo0 V Glo1;Glo0;Glo1").c_str(),400,0,0.04,400,0,0.04);

    hist2d[55] = new TH2D("scEfAmpVrtdt_Glo",addstr(fHTitle,"sc Eff Amp V Rt dT Glo;EffAmp;dt").c_str(),150,0,150,800,-4,4);
    hist2d[56] = new TH2D("scEfAmpVrtdt_sruLoc",addstr(fHTitle,"sc Eff Amp V Rt dT SRU Loc;EffAmp;dt").c_str(),150,0,150,800,-4,4);
    hist2d[57] = new TH2D("scEfAmpVtrdt_druLoc",addstr(fHTitle,"sc Eff Amp V Rt dT DRU Loc;EffAmp;dt").c_str(),150,0,150,800,-4,4);

    hist2d[90] = new TH2D("phoEVHOEM_gloLoc",addstr(fHTitle,"Energy V HadOverEM GLO Loc").c_str(),1500,0,750,250,0,0.25);
    hist2d[91] = new TH2D("phoSieieVHOEM_gloLoc",addstr(fHTitle,"SIeIe V HadOverEM GLO Loc").c_str(),400,0,0.04,250,0,0.25);
    hist2d[92] = new TH2D("phoEVCieie_gloLoc",addstr(fHTitle,"Energy V CovIEtaIEta GLO Loc").c_str(),1500,0,750,250,0,0.0025);
    hist2d[93] = new TH2D("phoSieieVCieie_gloLoc",addstr(fHTitle,"SIeIe V CovIEtaIEta GLO Loc").c_str(),400,0,0.04,250,0,0.0025);
    hist2d[94] = new TH2D("phoEVCieip_gloLoc",addstr(fHTitle,"Energy V CovIEtaIPhi GLO Loc").c_str(),1500,0,750,250,0,0.0025);
    hist2d[95] = new TH2D("phoSieieVCieip_gloLoc",addstr(fHTitle,"SIeIe V CovIEtaIPhi GLO Loc").c_str(),400,0,0.04,250,0,0.0025);
    hist2d[96] = new TH2D("phoEVCipip_gloLoc",addstr(fHTitle,"Energy V CovIPhiIPhi GLO Loc").c_str(),1500,0,750,250,0,0.0025);
    hist2d[97] = new TH2D("phoSieieVCipip_gloLoc",addstr(fHTitle,"SIeIe V CovIPhiIPhi GLO Loc").c_str(),400,0,0.04,250,0,0.0025);
    hist2d[98] = new TH2D("phoEVERHSEtCDR04_gloLoc",addstr(fHTitle,"Energy V EcalRHSumEtConeDR04 GLO Loc").c_str(),1500,0,750,100,0,50);
    hist2d[99] = new TH2D("phoSieieVERHSEtCDR04_gloLoc",addstr(fHTitle,"SIeIe V EcalRHSumEtConeDR04 GLO Loc").c_str(),400,0,0.04,100,0,50);
    hist2d[100] = new TH2D("phoEVHTSEtCDR04_gloLoc",addstr(fHTitle,"Energy V HcalTwrSumEtConeDR04 GLO Loc").c_str(),1500,0,750,100,0,50);
    hist2d[101] = new TH2D("phoSieieVHTSEtCDR04_gloLoc",addstr(fHTitle,"SIeIe V HcalTwrSumEtConeDR04 GLO Loc").c_str(),400,0,0.04,100,0,50);
    hist2d[102] = new TH2D("phoEVTSPtSCDR04_gloLoc",addstr(fHTitle,"Energy V TrkSumPtSolidConeDR04 GLO Loc").c_str(),1500,0,750,100,0,50);
    hist2d[103] = new TH2D("phoSieieVTSPtSCDR04_gloLoc",addstr(fHTitle,"SIeIe V TrkSumPtSolidConeDR04 GLO Loc").c_str(),400,0,0.04,100,0,50);
    hist2d[104] = new TH2D("phoEVTSPtHCDR04_gloLoc",addstr(fHTitle,"Energy V TrkSumPtHallowConeDR04 GLO Loc").c_str(),1500,0,750,100,0,50);
    hist2d[105] = new TH2D("phoSieieVTSPHCDR04_gloLoc",addstr(fHTitle,"SIeIe V TrkSumPtHallowConeDR04 GLO Loc").c_str(),400,0,0.04,100,0,50);

	hist2d[108] = new TH2D("SwcrVRtTime",addstr(fHTitle,"SwissCross V RH Time;SwCrs;Rt [ns]").c_str(),400,0.1,1.1,300,-15,15);
    hist2d[110] = new TH2D("SwcrVRtTimeTopo",addstr(fHTitle,"SwissCross V RH Time +Topo Cuts;SwCrs;Rt [ns]").c_str(),400,0.1,1.1,300,-15,15);
    hist2d[112] = new TH2D("SwcrVRtTimeOOT",addstr(fHTitle,"SwissCross V RH Time +kOOT Cut;SwCrs;Rt [ns]").c_str(),400,0.1,1.1,300,-15,15);


    hist2d[115] = new TH2D("SwcrVRhEnergyRtOOT",addstr(fHTitle,"SwissCross V RH Energy +RtkOOT Cut;SwCrs;Energy [GeV]").c_str(),400,0.1,1.1,1000,0,1000);
    hist2d[116] = new TH2D("SwcrVRhEnergyTopo",addstr(fHTitle,"SwissCross V RH Energy +kOOT Cut;SwCrs;Energy [GeV]").c_str(),400,0.1,1.1,1000,0,1000);
    hist2d[117] = new TH2D("SwcrVRhEnergy",addstr(fHTitle,"SwissCross V RH Energy +kOOT Cut;SwCrs;Energy [GeV]").c_str(),400,0.1,1.1,1000,0,1000);

    hist2d[118] = new TH2D("EnergyVRtTimeKfG16",addstr(fHTitle,"Energy V RtTime kOOT False hasGS1&6;Energy [GeV];Time [ns]").c_str(),250,0,1000,600,-30,30);
    hist2d[119] = new TH2D("EnergyVRtTimeG16",addstr(fHTitle,"Energy V RtTime;Energy [GeV] hasGS1&6;Time [ns]").c_str(),250,0,1000,600,-30,30);
    hist2d[120] = new TH2D("EnergyVRtTimeKfG6",addstr(fHTitle,"Energy V RtTime kOOT False hasGS6;Energy [GeV];Time [ns]").c_str(),250,0,1000,600,-30,30);
    hist2d[121] = new TH2D("EnergyVRtTimeG6",addstr(fHTitle,"Energy V RtTime;Energy [GeV] hasGS6;Time [ns]").c_str(),250,0,1000,600,-30,30);
    hist2d[122] = new TH2D("EnergyVRtTimeKfG1",addstr(fHTitle,"Energy V RtTime kOOT False hasGS1;Energy [GeV];Time [ns]").c_str(),250,0,1000,600,-30,30);
    hist2d[123] = new TH2D("EnergyVRtTimeG1",addstr(fHTitle,"Energy V RtTime;Energy [GeV] hasGS1;Time [ns]").c_str(),250,0,1000,600,-30,30);
    hist2d[124] = new TH2D("EnergyVRtTimeKfG0",addstr(fHTitle,"Energy V RtTime kOOT False noGS;Energy [GeV];Time [ns]").c_str(),250,0,1000,600,-30,30);
    hist2d[125] = new TH2D("EnergyVRtTimeG0",addstr(fHTitle,"Energy V RtTime;Energy [GeV] noGS;Time [ns]").c_str(),250,0,1000,600,-30,30);

}//<<>>void makehists::initHists()

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {

    	//auto indir = "ecalTiming/gammares_ttcc_140_v11_diag_mod1_exp3/EGamma1/";
        //auto indir = "ecalTiming/gammares_mc/DYto2L-4Jets_MLL-50_1J_TuneCP5_13p6TeV_madgraphMLM-pythia8/";
        //auto indir = "ecalTiming/gammares_mc/ZprimeToEE_M-6000_TuneCP5_13p6TeV_pythia8/";
        //auto indir = "ecalTiming/gammares_llpana/";
        auto indir = "/ecalTiming/gammares_llpana_pd/MET/";

        //auto infilename = "list_files/egammares_gammares_mc_DYto2L-4Jets_MLL-50_1J_MINIAODSIM_Run3Summer23MiniAODv4_v2.txt";
        //auto infilename = "list_files/egammares_gammares_mc_DYto2L-4Jets_MLL-50_1J_MINIAODSIM_Run3Winter24MiniAOD_v2.txt";
        //auto infilename = "list_files/egammares_gammares_mc_ZprimeToEE_M-6000_MINIAODSIM_Run3Winter24MiniAOD_v2.txt";
        //auto infilename = "list_files/egres_QCD_HT_1000to1500_R17.txt";
        auto infilename = "list_files/egammares_Met_PD_AOD_Run2017E-17Nov2017v2.txt";

		int brun = 0;
        int erun = 999999;
        auto califilename = "none";
        
        //auto outfilename = "egammares_diag_23D_370496_370580_tt_exp3_v12_EB.root";
        //auto outfilename = "egammares_diag_winter24_DYto2L-4Jets_MLL-50_1J_v12_EB.root";
        //auto outfilename = "egammares_diag_winter24_ZprimeToEE_M-6000_v12_EB_test.root";
        //auto outfilename = "egammares_diag_llpana_qcd_mc_diag.root";
        auto outfilename = "egammares_llpana_met_pd_diag.root";

        //auto fhtitle = "Run2024E 14_0_4 EB ";
		//auto fhtitle = "Winter24 DY EB";
        //auto fhtitle = "QCD MC 1k-1.5K";
        auto fhtitle = "MET AOD 2017E";

		makehists base;				
        base.llpgana_hist_maker( indir, infilename, outfilename, califilename, brun, erun, fhtitle, true );

    //}//<<>>if( argc != 4 ) 
    return 1;

}//<<>>int main ( int argc, char *argv[] )

