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
#include <TRandom.h>

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

#include "wc_ku_timefitter_wErr_func.cpp"
//#include "wc_ku_timefitter_wErr_func_new.cpp"
#include "KUCMSRootHelperFunctions.hh"
#include "KUCMSHelperFunctions.hh"
#include "ecal_config/KUCMSEcalDetIDFunctions.hh"

using namespace std;

/*
typedef unsigned int uInt;

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
        auto idinfo = DetIDMap[cmsswId];
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
        auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

void NormX(TH2F *& hist){ 	//, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY){
 
	std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;
    
    for (auto ibinX = 1; ibinX <= hist->GetXaxis()->GetNbins(); ibinX++){

    	//const auto binwidthX = hist->GetXaxis()->GetBinWidth(ibinX);
        const auto norm = hist->Integral(ibinX,ibinX,1,hist->GetYaxis()->GetNbins());
		if( norm == 0.0 ) continue;
        for (auto ibinY = 1; ibinY <= hist->GetYaxis()->GetNbins(); ibinY++){
        	//const auto binwidthY = hist->GetYaxis()->GetBinWidth(ibinY);
               
            //get multiplier/divisor
            //auto multiplier = 1.f;      
            //if (varBinsX) multiplier *= binwidthX;
            //if (varBinsY) multiplier *= binwidthY;
               
            // get content/error
            auto content = hist->GetBinContent(ibinX,ibinY);
            auto error   = hist->GetBinError  (ibinX,ibinY);
               
            // scale it
            //if (isUp){ 
            //   content *= multiplier;
            //   error   *= multiplier;
            //} else {
            content /= norm;
            error   /= norm;
            //}
               
            // set new contents
            hist->SetBinContent(ibinX,ibinY,content);
            hist->SetBinError  (ibinX,ibinY,error);
        }
	}
}

void scaleHist(TH2F *& hist, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY){

	std::cout << "Scaling " << (isUp?"up":"down") << " hist: " << hist->GetName() << std::endl;

    for (auto ibinX = 1; ibinX <= hist->GetXaxis()->GetNbins(); ibinX++){

    	const auto binwidthX = hist->GetXaxis()->GetBinWidth(ibinX);
        for (auto ibinY = 1; ibinY <= hist->GetYaxis()->GetNbins(); ibinY++){

        	const auto binwidthY = hist->GetYaxis()->GetBinWidth(ibinY);
            // get multiplier/divisor
            auto multiplier = 1.f;
            if( varBinsX ) multiplier *= binwidthX;
            if( varBinsY ) multiplier *= binwidthY;
			auto scale = ( not isUp ) ? multiplier : 1/multiplier; 
            // get content/error
            auto content = hist->GetBinContent(ibinX,ibinY)*scale;
            auto error = hist->GetBinError(ibinX,ibinY)*scale;
			// set new contents
            hist->SetBinContent(ibinX,ibinY,content);
            hist->SetBinError  (ibinX,ibinY,error);

        }//<<>>for (auto ibinY = 1; ibinY <= hist->GetYaxis()->GetNbins(); ibinY++)
	}//<<>>for (auto ibinX = 1; ibinX <= hist->GetXaxis()->GetNbins(); ibinX++)

//-------------------------------------------------------------------------------------------------------------------------------

}//<<>>void scaleHist(TH2D *& hist, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY)
*/

//void plot2dResolution( string indir, string infilelistname, string outfilename, string tvarname, string calimapname, string isd_type, int brun, int erun, int leta, int heta  ){
void plot2dResolution( std::string indir, std::string infilelistname, std::string outfilename, std::string tvarname, std::string calimapname, std::string isd_type, bool useAmp, std::string xbinstr, double noise, double constant ){

	bool debug = false;
    //bool debug = true;

	TRandom* getRandom = new TRandom();
	getRandom->SetSeed(0);
	unsigned int nRand(100);

    std::cout << "opening output file" << std::endl;
    string histoutfilename(outfilename+".root");
    TFile* fOutFile = new TFile( histoutfilename.c_str(), "RECREATE" );
    std::cout << "fOutFile : " << fOutFile << std::endl;

    double phoseedtimeCaliIc_0(0.0);
    double phoseedtimeCaliIcErr_0(0.0);
    double phoseedtimeCaliIc_1(0.0);
    double phoseedtimeCaliIcErr_1(0.0);

    // Declaration of leaf types
    UInt_t          run;
    vector<unsigned int> *resRhID;
    vector<float>   *resAmp;
    vector<float>   *resE;
    vector<float>   *resRtTime;
    vector<float>   *resTOF;

    // List of branches
    TBranch        *b_run;   //!
    TBranch        *b_resRhID;   //!
    TBranch        *b_resAmp;   //!
    TBranch        *b_resE;   //!
    TBranch        *b_resRtTime;   //!
    TBranch        *b_resTOF;   //!


    // >> calcs  <<
    std::cout << "Setting up 2D plot" << std::endl;

	string locsname("SRO_");
    string locdname("DRO_");
    string globname("ZEE_");

    string histname("Data_Hist");
    string histnameg("Data_GHist");
    string histnames("Data_SHist");
    string histname0("Data_OthHist");
    string histname1("td_Hist");
    string histname2("tdc_Hist");
    string histname3("teffa_Hist");
    string histname4("tdfc_Hist");
    string histname5("td_eta_Hist");
    string histname6("td_phi_Hist");
    string fTitle("#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)");
    string fTitle0("#Delta(Photon Seed Time) [ns] vs. A_{0}/#sigma_{n} (EBEB)");
    string fTitle1("Photon Seed Time [ns]");
    string fTitle2("Photon Seed Time Calibrated [ns]");
    string fTitle3("Amplitude 0 vs Amplitude 1");
    string fTitle4("Photon Seed Time Filtered Calibrated [ns]");
    string fTitle5("#Delta(Photon Seed Time) [ns] vs iEta");
    string fTitle6("#Delta(Photon Seed Time) [ns] vs iPhi");
    string fXTitle("A_{eff}/#sigma_{n} (EBEB)");
    std::vector<Double_t> fXBins;
    Bool_t fXVarBins = false;//dummy not used
    string fYTitle("#Delta(Photon Seed Time) [ns] (EBEB)");
    std::vector<Double_t> fYBins;
    Bool_t fYVarBins = false;//dummy not used;
    //string ybinstr("CONSTANT 1920 -3 3");
    //string ybinstr("CONSTANT 600 -3 3");
    //string ybinstr("CONSTANT 120 -3 3");
    string ybinstr("CONSTANT 720 -9 9");
    string fZTitle("");

    setBins(xbinstr,fXBins,fXVarBins);
    setBins(ybinstr,fYBins,fYVarBins);

    const auto xbins = &fXBins[0];
    const auto ybins = &fYBins[0];
    int nMyBins = fXBins.size()-1;

    int tdiv = 720;
    float tstart = -9.0;
    float tend = 9.0;

    int ediv = 720;
    float eend = 3600;

    std::string tehistname = "timeErrHist";
	auto timeErrBinHist = new TH1F(tehistname.c_str(), tehistname.c_str(),1000,0,1);
    TH1F * caliErrBinHist[nMyBins+1];
    std::string eahistname = "effAHist";
    caliErrBinHist[0] = new TH1F(eahistname.c_str(),eahistname.c_str(),225,0,2250);
    caliErrBinHist[0]->Sumw2();
    for( auto ibin = 1; ibin < nMyBins+1; ibin++ ){
		auto binhistname = "bin" + to_string(ibin) + "ErrHist";
		caliErrBinHist[ibin] = new TH1F(binhistname.c_str(),binhistname.c_str(),500,0,0.5);
        caliErrBinHist[ibin]->Sumw2();
    }//<<>>for( auto ibin = 1; ibin < nMyBins+1; ibin++ )

	// local same
	auto lochist = locsname+histname;
    auto theHistLS = new TH2F(lochist.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto lochistg = locsname+histnameg;
    auto theGHistLS = new TH2F(lochistg.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto lochists = locsname+histnames;
    auto theSHistLS = new TH2F(lochists.c_str(),fTitle.c_str(), ediv, 0, eend, tdiv, tstart, tend);
    auto lochist0 = locsname+histname0;
    auto theOthHistLS = new TH2F(lochist0.c_str(),fTitle0.c_str(), ediv, 0, eend, tdiv, tstart, tend);
    auto lochist1 = locsname+histname1;
    auto thetdHistLS = new TH1F(lochist1.c_str(),fTitle1.c_str(),tdiv,tstart,tend);
    auto lochist2 = locsname+histname2;
    auto thetdcHistLS = new TH1F(lochist2.c_str(),fTitle2.c_str(),tdiv,tstart,tend);

    //auto thetdfcHist = new TH1F(histname4.c_str(),fTitle4.c_str(),tdiv,tstart,tend);
    auto lochist3 = locsname+histname3;
    auto theEffaHistLS = new TH2F(lochist3.c_str(),fTitle3.c_str(),ediv,0,eend,ediv,0,eend);
    auto lochist5 = locsname+histname5;
    auto theetaHistLS = new TH2F(lochist5.c_str(),fTitle5.c_str(),182,-90.5,90.5,tdiv,tstart,tend);
    auto lochist6 = locsname+histname6;
    auto thephiHistLS = new TH2F(lochist6.c_str(),fTitle6.c_str(),363,-1.5,361.5,tdiv,tstart,tend);

    theHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistLS->GetZaxis()->SetTitle(fZTitle.c_str());
    theGHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
    theGHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
    theGHistLS->GetZaxis()->SetTitle(fZTitle.c_str());

	// local diffrnent
    auto locdhist = locdname+histname;
    auto theHistLD = new TH2F(locdhist.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    theHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistLD->GetZaxis()->SetTitle(fZTitle.c_str());

	// global
    auto globhist = globname+histname;
    auto theHistGB = new TH2F(globhist.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    theHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistGB->GetZaxis()->SetTitle(fZTitle.c_str());

    std::cout << "Setting up DetIDs." << std::endl;
    std::map<UInt_t,DetIDStruct> DetIDMap;
    SetupDetIDsEB( DetIDMap );
    SetupDetIDsEE( DetIDMap );

    double goodlev(0);
    double goodlin(0);
	double goodgev(0);
    double goodgin(0);
    double gevents(0);

    std::cout << "open input files list : " << infilelistname << std::endl;
    string disphotreename("tree/llpgtree");

    std::ifstream infilelist(infilelistname);
    std::string infiles;
    while (std::getline(infilelist,infiles)){

        std::stringstream ss(infiles);
        std::string infilename;
        std::string califilename;
        std::string srunstr;
        std::string erunstr;
        std::string subdir;
        ss >> infilename >> califilename >> srunstr >> erunstr >> subdir;
		char hash('#');
		char space(' ');
		char checkchar = ( not infilename.empty() ) ? infilename[0] : space;
		if( checkchar == hash || checkchar == space ) continue;
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "open input cali : " << califilename << std::endl;
		auto srun = std::stoi(srunstr);
        auto erun = std::stoi(erunstr);
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Subdir: " << subdir << std::endl;

        std::ifstream infile( infilename );
        std::string str;

        auto fInTree = new TChain(disphotreename.c_str());
        std::cout << "Adding files to TChain." << std::endl;
        while (std::getline(infile,str)){
            const std::string eosdir("root://cmseos.fnal.gov//store/user/");
            //const std::string eosdir("/uscms/home/jaking/nobackup/el8/ecaltiming/CMSSW_14_0_4/src/GammaResTool/GammaResTool/");
            auto tfilename = eosdir + indir + subdir + str;
            //auto tfilename = indir + "/" + str;
            //std::cout << "--  adding file: " << tfilename << std::endl;
			std::cout << "-";
        	fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
		std::cout << std::endl;

		//auto calidir = "/home/t3-ku/jaking/ecaltiming/skimmed_trees/local_chain/";
		string calidir = "cali_root_files/";
		TFile* fCaliFile(NULL);
		if( calimapname != "none" ){
			//string calidir = "";
			string califilein = calidir+califilename;
			std::cout << "Using caliFile : " << califilein << std::endl;
        	if( califilename != "none" ) fCaliFile = TFile::Open( califilein.c_str(), "read" );
		}//<<>>if( calimapname != "none" ){
        std::cout << "fInFile : " << infilename  << " fInTree : " << fInTree << " fCaliFile : " << fCaliFile << std::endl;

        std::cout << "set branches to get from fInFile : fInTree" << std::endl;

		run = 0;
		resRhID = 0;
		resAmp = 0;
		resE = 0;
		resRtTime = 0;
		resTOF = 0;

        fInTree->SetBranchAddress( "run", &run, &b_run );   //!
        fInTree->SetBranchAddress( "resRhID", &resRhID, &b_resRhID );   //!
        fInTree->SetBranchAddress( "resAmp", &resAmp, &b_resAmp);   //!
        fInTree->SetBranchAddress( "resE", &resE, &b_resE);   //!
        fInTree->SetBranchAddress( "resRtTime", &resRtTime, &b_resRtTime);   //!
        fInTree->SetBranchAddress( "resTOF" , &resTOF, &b_resTOF);   //!

		vector<TH2F*> calimaps{0,0,0,0,0,0,0,0};
        std::cout << "get maps from fCaliFile" << std::endl;
        if( fCaliFile ){
			fCaliFile->cd();
        	string cmbsrt("AveXtalRecTime");
        	string itcnt("");
        	string ebmaprt(cmbsrt+"EBMap"+itcnt);
        	string eberrmaprt(cmbsrt+"ErrEBMap"+itcnt);
        	string epmaprt(cmbsrt+"EPMap"+itcnt);
        	string emmaprt(cmbsrt+"EMMap"+itcnt);
        	calimaps[0] = (TH2F*)fCaliFile->Get(ebmaprt.c_str());
        	calimaps[1] = (TH2F*)fCaliFile->Get(eberrmaprt.c_str());
        	calimaps[2] = (TH2F*)fCaliFile->Get(epmaprt.c_str());
        	calimaps[3] = (TH2F*)fCaliFile->Get(emmaprt.c_str());
		}//<<>>if( calimapname != "none" )

        std::cout << "Getting calibration values and plotting" << std::endl;

        auto nEntries = fInTree->GetEntries();
		if( debug ) nEntries = ( nEntries < 1000 ) ? nEntries : 1000;
		if( debug ) std::cout << "Mf2d Proccessing " << nEntries << " entries : " << std::endl;
        for (auto centry = 0U; centry < nEntries; centry++){

			if( centry%10000000 == 0 or centry == 0){ 
				std::cout << "Proccessed " << centry << " of " << nEntries;
				std::cout << " " << (1000.0*static_cast<float>(centry)/static_cast<float>(nEntries))/10.0 << "%" << std::endl;
			}//<<>>if( centry%10000000 == 0 or centry == 0)

			auto entry = fInTree->LoadTree(centry);

			if(debug) std::cout << " - Start loop " << std::endl;

			gevents++;

            b_run->GetEntry(entry);   //!
            b_resRhID->GetEntry(entry);   //!
            b_resAmp->GetEntry(entry);   //!
            b_resE->GetEntry(entry);   //!
            b_resRtTime->GetEntry(entry);   //!
            b_resTOF->GetEntry(entry);   //!

			//int didx = 0;
			for( int didx = debug?0:4; didx < 4; didx++ ){
				if(debug) std::cout << "Run " << run << " id " << (*resRhID)[didx] << " Amp " << (*resAmp)[didx] << " E " << (*resE)[didx];
            	if(debug) std::cout << " Rt " << (*resRtTime)[didx]  << " TOF " << (*resTOF)[didx] << std::endl;
			}//<<>>for( int didx = 0; didx < 4; didx++ ){
			if( srun != 0 && ( run < srun || run > erun ) ) continue;

			if(debug) std::cout << " - Finshed Get Entry " << std::endl;

			auto idinfoL0 = DetIDMap[(*resRhID)[0]];
            auto idinfoL1 = DetIDMap[(*resRhID)[1]];
            auto idinfoG0 = DetIDMap[(*resRhID)[2]];
            auto idinfoG1 = DetIDMap[(*resRhID)[3]];

            auto i1L0 = idinfoL0.i1;
            auto i1L1 = idinfoL1.i1;
            auto i2L0 = idinfoL0.i2;     
            auto i2L1 = idinfoL1.i2;

            auto i1G0 = idinfoG0.i1;
            auto i1G1 = idinfoG1.i1;
            auto i2G0 = idinfoG0.i2;
            auto i2G1 = idinfoG1.i2;

			auto L0EB = idinfoL0.ecal;
            auto L1EB = idinfoL1.ecal;
            auto G0EB = idinfoG0.ecal;
            auto G1EB = idinfoG1.ecal;

            if(debug) std::cout << " - Finshed Setting IDInfo " << std::endl;

	        int bin_offset = 86;
	        int adjust = 0.0;

			vector<float> seedTimeIC0{0,0};
            vector<float> seedTimeICE0{0,0};
            vector<float> seedTimeIC1{0,0};
            vector<float> seedTimeICE1{0,0};

	        if( fCaliFile ){ //and califilename != "none" ){

				if ( (*resRhID)[0] ){
            		if ( L0EB == ECAL::EB ){
            			seedTimeIC0[0] = calimaps[0]->GetBinContent( i2L0 + bin_offset, i1L0 ) - adjust;
						seedTimeICE0[0] = calimaps[1]->GetBinContent( i2L0 + bin_offset, i1L0 );
            		}else if ( L0EB == ECAL::EP ){
            	    	seedTimeIC0[0] = calimaps[2]->GetBinContent( i2L0, i1L0 ) - adjust;
            		}else if (L0EB == ECAL::EM ){
            	    	seedTimeIC0[0] = calimaps[3]->GetBinContent( i2L0, i1L0 ) - adjust;
            		}//<<>>if ( L0EB == ECAL::EB )
				}//<<>>if ( (*resRhID)[0] )
            
                if ( (*resRhID)[1] ){
        			if ( L1EB == ECAL::EB ){
            			seedTimeIC1[0] = calimaps[0]->GetBinContent( i2L1 + bin_offset, i1L1 ) - adjust;
                		seedTimeICE1[0] = calimaps[1]->GetBinContent( i2L1 + bin_offset, i1L1 );
            		}else if ( L1EB == ECAL::EP ){
                		seedTimeIC1[0] = calimaps[2]->GetBinContent( i2L1, i1L1 ) - adjust;
            		}else if (L1EB == ECAL::EM ){
                		seedTimeIC1[0] = calimaps[3]->GetBinContent( i2L1, i1L1 ) - adjust;
					}//<<>>if ( L1EB == ECAL::EB )
                }//<<>>if ( (*resRhID)[1] )

                if ( (*resRhID)[2] ){
                	if ( G0EB == ECAL::EB ){
                    	seedTimeIC0[1] = calimaps[0]->GetBinContent( i2G0 + bin_offset, i1G0 ) - adjust;
                    	seedTimeICE0[1] = calimaps[1]->GetBinContent( i2G0 + bin_offset, i1G0 );
                	}else if ( G0EB == ECAL::EP ){
                    	seedTimeIC0[1] = calimaps[2]->GetBinContent( i2G0, i1G0 ) - adjust;
                	}else if (G0EB == ECAL::EM ){
                    	seedTimeIC0[1] = calimaps[3]->GetBinContent( i2G0, i1G0 ) - adjust;
                	}//<<>>if ( G0EB == ECAL::EB )
                }//<<>>if ( (*resRhID)[2] )

                if ( (*resRhID)[3] ){
                	if ( G1EB == ECAL::EB ){
                    	seedTimeIC1[1] = calimaps[0]->GetBinContent( i2G1 + bin_offset, i1G1 ) - adjust;
                    	seedTimeICE1[1] = calimaps[1]->GetBinContent( i2G1 + bin_offset, i1G1 );
                	}else if ( G1EB == ECAL::EP ){
                    	seedTimeIC1[1] = calimaps[2]->GetBinContent( i2G1, i1G1 ) - adjust;
                	}else if (G1EB == ECAL::EM ){
                    	seedTimeIC1[1] = calimaps[3]->GetBinContent( i2G1, i1G1 ) - adjust;
                	}//<<>>if ( G1EB == ECAL::EB )
                }//<<>>if ( (*resRhID)[3] )

	        }//<<>>if( calimapname != "none" )


//-------------------set for local, repo calcs for global --------------------------------

            if(debug) std::cout << " - Calc 2D Hist" << std::endl;

            double leffa0 = useAmp ? (*resAmp)[0] : (*resE)[0]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
            double leffa1 = useAmp ? (*resAmp)[1] : (*resE)[1]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
            double geffa0 = useAmp ? (*resAmp)[2] : (*resE)[2]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
            double geffa1 = useAmp ? (*resAmp)[3] : (*resE)[3]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
            double lxfill = (leffa0*leffa1)/sqrt(pow(leffa0,2)+pow(leffa1,2));
            double gxfill = (geffa0*geffa1)/sqrt(pow(geffa0,2)+pow(geffa1,2));

            double ldTOF = (*resTOF)[0]-(*resTOF)[1]; //phoseedTOF_0-phoseedTOF_1;
            double gdTOF = (*resTOF)[2]-(*resTOF)[3]; //phoseedTOF_0-phoseedTOF_1;
           	double lyf0 = (*resRtTime)[0]-seedTimeIC0[0];
           	double lyf1 = (*resRtTime)[1]-seedTimeIC1[0];
           	double gyf0 = (*resRtTime)[2]-seedTimeIC0[1];
           	double gyf1 = (*resRtTime)[3]-seedTimeIC1[1];
			//Double_t  Gaus(Double_t mean = 0, Double_t sigma = 1)
			if( constant != 0 ){
				double cterm = 2*pow(constant,2);
                lyf0 = getRandom->Gaus(lyf0,std::sqrt((pow(noise/leffa0,2)+cterm)/2));
                lyf1 = getRandom->Gaus(lyf1,std::sqrt((pow(noise/leffa1,2)+cterm)/2));
                gyf0 = getRandom->Gaus(gyf0,std::sqrt((pow(noise/geffa0,2)+cterm)/2));
                gyf1 = getRandom->Gaus(gyf1,std::sqrt((pow(noise/geffa1,2)+cterm)/2));
			}//<<>>if( cterm != 0 )
            double lyfill = lyf0-lyf1+ldTOF;
            double gyfill = gyf0-gyf1+gdTOF;

            double timeErr0 = 0; //phoseedtimeErr_0/25.0;
            double timeErr1 = 0; //phoseedtimeErr_1/25.0;
			double ldtserr0 = timeErr0*timeErr0+seedTimeICE0[0]*seedTimeICE0[0];
            double ldtserr1 = timeErr1*timeErr1+seedTimeICE1[0]*seedTimeICE1[0];
			double ldterr = sqrt(ldtserr0+ldtserr1);
            double gdtserr0 = timeErr0*timeErr0+seedTimeICE0[1]*seedTimeICE0[1];
            double gdtserr1 = timeErr1*timeErr1+seedTimeICE1[1]*seedTimeICE1[1];
            double gdterr = sqrt(gdtserr0+gdtserr1);

		    auto le_cut = ((*resE)[0]>=2)&&((*resE)[0]<=120)&&((*resE)[1]>=2)&&((*resE)[1]<=120);
            auto ge_cut = ((*resE)[2]>=2)&&((*resE)[2]<=120)&&((*resE)[3]>=2)&&((*resE)[3]<=120);
	        auto leta_cut = (L0EB == ECAL::EB)&&(L1EB == ECAL::EB);
            auto geta_cut = (G0EB == ECAL::EB)&&(G1EB == ECAL::EB);
            auto goodLocTime = (*resRtTime)[0] != 0 && (*resRtTime)[1] != 0; 
            auto goodGloTime = (*resRtTime)[2] != 0 && (*resRtTime)[3] != 0; 
			auto goodLocRHs = (*resRhID)[0] != 0 && (*resRhID)[1] != 0;
            auto goodGloRHs = (*resRhID)[2] != 0 && (*resRhID)[3] != 0;

			auto isd_cut = idinfoL0.TT == idinfoL1.TT; // true = same, fasle = different
	        auto levent_good = le_cut && leta_cut && goodLocRHs && goodLocTime;
            auto gevent_good = ge_cut && geta_cut && goodGloRHs && goodGloTime;

            if( levent_good ) goodlev++;
            if( goodLocRHs ) goodlin++;
            if( gevent_good ) goodgev++;
            if( goodGloRHs ) goodgin++;

			if(debug) std::cout << " - lxfill : " << lxfill << " lyfill : " << lyfill;
			if(debug) std::cout  << " flag : " << le_cut << " " << leta_cut << " " << goodLocRHs << std::endl; 
            if(debug) std::cout << " - gxfill : " << gxfill << " gyfill : " << gyfill;
            if(debug) std::cout  << " flag : " << ge_cut << " " << geta_cut << " " << goodGloRHs << std::endl;

//-----------------   set up fills for local same diffrent and global -----------------------

			if(debug) std::cout << " - Fill 2D Hist" << std::endl;
	        if( levent_good && isd_cut ){
				theHistLS->Fill(lxfill,lyfill);
				theSHistLS->Fill(lxfill,lyfill);
				theOthHistLS->Fill(leffa0,lyfill);
				//for( auto i = 0; i < nRand; i++ ){ 
				//	//virtual Double_t	Gaus(Double_t mean = 0, Double_t sigma = 1)
				//	auto gsfill = getRandom->Gaus(lyfill,ldterr);
                //  	theGHistLS->Fill(lxfill,gsfill);
				//}//<<>>for( int i = 0; i < 100; i++ )	
				if(debug) std::cout << " - Fill effA dist Hist" << std::endl;
				//string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 950 2250");
				caliErrBinHist[0]->Fill(lxfill);
				timeErrBinHist->Fill(timeErr0); timeErrBinHist->Fill(timeErr1);
               	if(debug) std::cout << " - Fill bin Hists" << std::endl;
				if( lxfill < 75 ) caliErrBinHist[1]->Fill(ldterr);
				else if( lxfill < 100 ) caliErrBinHist[2]->Fill(ldterr);
               	else if( lxfill < 125 ) caliErrBinHist[3]->Fill(ldterr);
               	else if( lxfill < 150 ) caliErrBinHist[4]->Fill(ldterr);
               	else if( lxfill < 175 ) caliErrBinHist[5]->Fill(ldterr);
               	else if( lxfill < 225 ) caliErrBinHist[6]->Fill(ldterr);
               	else if( lxfill < 275 ) caliErrBinHist[7]->Fill(ldterr);
               	else if( lxfill < 325 ) caliErrBinHist[8]->Fill(ldterr);
               	else if( lxfill < 375 ) caliErrBinHist[9]->Fill(ldterr);
               	else if( lxfill < 475 ) caliErrBinHist[10]->Fill(ldterr);
               	else if( lxfill < 600 ) caliErrBinHist[11]->Fill(ldterr);
               	else if( lxfill < 950 ) caliErrBinHist[12]->Fill(ldterr);
               	else if( lxfill < 2250 ) caliErrBinHist[13]->Fill(ldterr);
				//else std::cout << "Over 2250" << std::endl;
            }//<<>>if( levent_good )
            if(debug) std::cout << " - Fill 1D Hist" << std::endl;
            if( levent_good ){
				theEffaHistLS->Fill(leffa0,leffa1); 
				thetdHistLS->Fill((*resRtTime)[0]); thetdHistLS->Fill((*resRtTime)[1]);
                thetdcHistLS->Fill((*resRtTime)[0]-seedTimeIC0[0]); thetdcHistLS->Fill((*resRtTime)[1]-seedTimeIC1[0]);
				theetaHistLS->Fill(i2L0,lyfill); theetaHistLS->Fill(i2L1,lyfill);
                thephiHistLS->Fill(i1L0,lyfill); thephiHistLS->Fill(i1L1,lyfill);
			}//<<>>if( levent_good )
			if( levent_good && not isd_cut ){ theHistLD->Fill(lxfill,lyfill); }
            if( gevent_good ){ theHistGB->Fill(gxfill,gyfill); }

			if(debug) std::cout << " - Fill hists done" << std::endl;

        } // for (auto entry = 0U; entry < nEntries; entry++)
	 	//delete fInFile;

        //delete fCaliFile;  <<<<<<<<<<<<<<   delete califiles ????????????????
		if(debug) std::cout << " -------- Next Input file " << std::endl;

    } // while (std::getline(infilelist,infiles))

// --------------  process new histos ----------------------------------------------

	if(debug) std::cout << " - Scale 2D Hist" << std::endl;
    scaleHist(theHistLS,false,fXVarBins,fYVarBins);
    scaleHist(theGHistLS,false,fXVarBins,fYVarBins);
    scaleHist(theHistLD,false,fXVarBins,fYVarBins);
    scaleHist(theHistGB,false,fXVarBins,fYVarBins);
    normTH2F(theetaHistLS);
    normTH2F(thephiHistLS);

    fOutFile->cd();
    theHistLS->Write();
    theGHistLS->Write();
    theSHistLS->Write();
    theOthHistLS->Write();
    theEffaHistLS->Write();
    thetdHistLS->Write();
    thetdcHistLS->Write();
    theetaHistLS->Write();
    thephiHistLS->Write();
    timeErrBinHist->Write();
    for( auto ibin = 0; ibin < nMyBins; ibin++ ){ caliErrBinHist[ibin]->Write(); }
	theHistLD->Write();
    theHistGB->Write();

    delete theHistLS;
    delete theGHistLS;
    delete theSHistLS;
    delete theOthHistLS;
    delete theEffaHistLS;
    delete thetdHistLS;
    delete thetdcHistLS;
    delete theetaHistLS;
    delete thephiHistLS;
    delete timeErrBinHist;
    for( auto ibin = 0; ibin < nMyBins; ibin++ ){ delete caliErrBinHist[ibin]; }
    delete theHistLD;
    delete theHistGB;

    delete fOutFile;
	delete getRandom;

	auto passlev = 100.0*goodlev/gevents;
    auto passlin = 100.0*goodlin/gevents;
    auto passgev = 100.0*goodgev/gevents;
    auto passgin = 100.0*goodgin/gevents;

    if(gevents) std::cout << "Processed " << gevents << " with %" << passlev << " from %" << passlin << " of Local"; 
	if(gevents) std::cout << " and %" << passgev << " from %" << passgin << " of Global" << std::endl;
    std::cout << "Thats all Folks!" << std::endl;

}//<<>> void plot2dResolution( std::string indir, std::string infilelistname, 

int main ( int argc, char *argv[] ){

    	//std::string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 950 2250 9000"); 
   		//std::string xbinstr("VARIABLE 0 75 100 125 150 175 200 225 250 275 325 375 450 550 725 925 1325 1700 2250"); 
    	//std::string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 750 950 1275 1700 2250");
        //std::string xbinstr("VARIABLE 60 80 100 120 150 180 240 320 400 480 560 640 800 1600 3200");
        //std::string xbinstr("VARIABLE 60 120 180 240 300 360 420 480 540 600 690 780 960 1260 1740 2220 3200");
        //std::string xbinstr("VARIABLE 60 120 180 240 300 420 510 600 780 960 1260 1740 2220 3200");
        //std::string xbinstr("VARIABLE 25 50 75 100 125 150 175 200 225 250 300 400 600 1000 1800"); // ZEE
        std::string xbinstr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200 1800" ); // : xa
		bool useAmp(true);
        //std::string xbinstr("VARIABLE 0 10 12 14 16 18 20 22 24 26 30 34 40 48 62 78 120");
        //bool useAmp(false);

        //if( argc != 7 ) { std::cout << "Insufficent arguments." << std::endl; }
        //else {

            //std::string indir = "jaking/ecalTiming/";
            //std::string indir = "jaking/ecalTiming/gammares_ttcc_140_v11_diag_ebsf_ccgt/";
            //std::string indir = "jaking/ecalTiming/gammares_mc/";
            //std::string indir = "jaking/ecalTiming/gammares_llpana/MET/";
            //std::string indir = "jaking/ecalTiming/gammares_llpana_pd/MET/";
            //std::string indir = "jaking/ecalTiming/gammares_llpana_qcd/";
			//std::string indir = "jaking/ecalTiming/gammares_llpana_v2/MET/";
			//std::string indir = "jaking/ecalTiming/gammares_llpana_pd/";
            //std::string indir = "jaking/ecalTiming/gammares_llpana/";
            //std::string indir = "jaking/ecalTiming/gammares_llpana_mc/";
    		std::string indir("jaking/ecalTiming/gammares_ECAL_CC_HCAL_DI-v3/");
    		//std::string indir("jaking/ecalTiming/gammares_r24f_prompt/");

            //std::string infilelistname = "list_files/egammares_Met_PD_AOD_Run2017E-17Nov2017v2_reso_plotfilelist.txt";
            //std::string infilelistname = "list_files/egammares_MetPD_nocali_AOD_Run2017E-17Nov2017v2_reso_plotfilelist.txt";
            //std::string infilelistname = "list_files/kuntuple_QCDHT_Met75_R17_v20_infileslist.txt";
            //std::string infilelistname = "list_files/egammares_MetPD_MINIAOD_Run2017E_reso_plotfilelist.txt";
            //std::string infilelistname = "list_files/egammares_MetPD_MINIAOD_Run2017E_reso_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_Met_PD_AOD_Run2017E-17Nov2017v2_reso_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_Met_PD_AOD_Run2017E-17Nov2017v2_reso_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_MetPD_MINIAOD_Run2017E_304475_none_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_MetPD_MINIAOD_Run2017E_304475_Cali_plotfilelist.txt";
			//std::string infilelistname = "master_list_files/egammares_MetPD_MINIAOD_Run2017E_304475_CaliFltrd_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_MetPD_MINIAOD_Run2017E_v2_304475_Cali_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_DEGPD_AOD_Run2017F_v2_304475_noCali_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_DEGPD_AOD_Run2017F_v2_304475_Cali_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_DEGPD_AOD_Run2017EF_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/egammares_EGMPD_MINIAOD_Run2018D_v2_327238_Cali_plotfilelist.txt";
            //std::string infilelistname = "master_list_files/kuntuple_QCDHT_Met75_R17_v20_infileslist.txt";
            //std::string infilelistname = "master_list_files/kuntuple_QCDHT100200_Met75_R17_v20_plotfileslist.txt";
            //std::string infilelistname = "master_list_files/egammares_DYJetsToLL_Met75_R17_v20_plotfilelist.txt";
            std::string infilelistname = "egammares_DEG_CCvRT_R24_v20_plotfilelist.txt";

            //std::string outfilename = "egres_MetPD_AOD_Run2017E_Cali_Fltr_304475_14011_v12_resplots";
            //std::string outfilename = "egres_MetPD_AOD_Run2017E_noCali_Full_14011_v12_resplots";
            //std::string outfilename = "egres_QCDMC_AOD_R17_14011_v10_resplots";
            //std::string outfilename = "egres_MetPD_MINIAOD_Run2017E_14011_v12_resplots";
            //std::string outfilename = "egres_MetPD_MINIAOD_Run2017E_Cali_noFltr_304475_14011_v12_resplots";
            //std::string outfilename = "egres_MetPD_MINIAOD_Run2017E_Cali_Fltr_304475_14011_v12_resplots";
            //std::string outfilename = "egres_MetPD_AOD_Run2017E_v2_Cali_304475_14011_v12_resplots";
            //std::string outfilename = "egres_DEGPD_AOD_Run2017F_v2_Cali_FilteredFit_Full_14011_v12_resplots";
			//std::string outfilename = "egres_MetPD_AOD_Run2017E_Cali_FltrFit_304475_14011_v12_resplots";
            //std::string outfilename = "egres_EgmPD_MINIAOD_Run2018D_324305_327239_Filter_14011_v12_resplots";
            //std::string outfilename = "egres_QCD_AODSIM_Met75_R17_FilterFit_14011_v12_resplots";
            //std::string outfilename = "egres_QCDHT100200_AODSIM_Met75_R17_FilterFit_14011_v12_resplots";
            //std::string outfilename = "egres_DY1JetsToLL_AODSIM_R17_FilterFit_Smear_v3_ZBins_14011_v12_resplots";
            //std::string outfilename = "egres_DEGPD_AOD_Run2017F_v2_Cali_FilteredFit_ZBins_Full_14011_v12_resplots";
            //std::string outfilename = "egres_DEGPD_AOD_Run2017EF_v2_Cali_FilteredFit_ZBins_Full_14011_v12_resplots";
            std::string outfilename = "egres_DEGPD_AODSIM_R24_ECAL_v3_14011_v12_resplots";

            std::string tvarname = ""; //argv[4];
			std::string calimapname = "infile"; //argv[5];
            std::string isd_type = "yes"; //argv[6];
			double noise = 0.0;// no smearing term
			double constant = 0.0;// no smearing term
            //double noise = 14.0;// smearing term
            //double constant = 0.1083;// smearing term			
            //auto brun = std::stoi(argv[7]);
            //auto erun = std::stoi(argv[8]);
            //auto leta = std::stoi(argv[9]);
            //auto heta = std::stoi(argv[10]);
      		plot2dResolution( indir, infilelistname, outfilename, tvarname, calimapname, isd_type, useAmp, xbinstr, noise, constant );
            //float rnoise = 0.0;// no smearing term
            //float rconstant = 0.0;// no smearing term
            //float rnoise = 14.0;// smearing term
            //float rconstant = 0.1083;// smearing term  
			//std::string fitInFile = outfilename + ".root";
			//runTimeFitter( fitInFile, outfilename, "SRO_Data_Hist", xbinstr, rnoise, rconstant );
            //runTimeFitter( fitInFile, outfilename, "DRO_Data_Hist", xbinstr, rnoise, rconstant );
            //runTimeFitter( fitInFile, outfilename, "ZEE_Data_Hist", xbinstr, rnoise, rconstant );
//void runTimeFitter(const std::string fInFileName, const std::string fOutFileText, const std::string f2DHistName, const std::string & inplotname, std::string xbinstr, float noise, float constant ){

        //}//<<>>//if( argc != 7 )
        return 1;

}//<<>>int main ( int argc, char *argv[] )
