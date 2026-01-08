// ROOT inludes
#include "TStyle.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
//#include "TString.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TPaveText.h"

// STL includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <sys/stat.h>

#include "KUCMSRootHelperFunctions.hh"
#include "KUCMSHelperFunctions.hh"
#include "wc_ku_gausTimeFit.cpp"

//************** Primary function *************************************************************************************
void runTimeFitter(const std::string fInFileName, const std::string fOutFileText, const std::string f2DHistName, std::string xbinstr, float noise, float constant ){ 

  std::cout << "Initializing TimeFitter..." << std::endl;

  TFile * fInFile;
  TFile * fOutFile;
  TPaveText * fConfigPave;

  fInFile = TFile::Open(fInFileName.c_str());
  std::string fOutFileName = fOutFileText + ".root";
  fOutFile = TFile::Open(fOutFileName.c_str(),"UPDATE");

  std::string fTitle = "#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)";
  std::string fXTitle = "A_{eff}/#sigma_{n} (EBEB)";
  std::vector<Double_t> fXBins;
  bool fXVarBins = false;//dummy not used
  setBins(xbinstr,fXBins,fXVarBins);
  int fNBinsX = fXBins.size()-1;

  TimeFitType fTimeFitType = TimeFitType::Gaus1core;
  float fRangeLow = 2.0f;
  float fRangeUp = 2.0f;
  std::string fTimeText = "t_{1}-t_{2}";
  std::string fSigmaVarText = "A_{eff}/#sigma_{n}";
  std::string fSigmaVarUnit = "";

  // set fitter

  std::map<std::string,TH1F*> ResultsMap;
  const auto xbins = &fXBins[0];

  auto chi2ndfName = f2DHistName+"_chi2ndf";
  auto chi2ndfTitle = fTitle+" "+"#chi^{2}/NDF."+";"+fXTitle+";"+"#chi^{2}/NDF.";
  ResultsMap["chi2ndf"]  = new TH1F(chi2ndfName.c_str(),chi2ndfName.c_str(),fNBinsX,xbins);
  ResultsMap["chi2ndf"]->Sumw2();
  auto chi2probName = f2DHistName+"_chi2prob";
  auto chi2probTitle = fTitle+" "+"#chi^{2} Prob."+";"+fXTitle+";"+"#chi^{2} Prob.";
  ResultsMap["chi2prob"] = new TH1F(chi2probName.c_str(),chi2probTitle.c_str(),fNBinsX,xbins);
  ResultsMap["chi2prob"]->Sumw2();
  auto mutext = "#mu("+fTimeText+") [ns]";
  auto muName = f2DHistName+"_mu";
  auto muTitle = fTitle+" "+mutext+";"+fXTitle+";"+mutext;
  ResultsMap["mu"] = new TH1F(muName.c_str(),muTitle.c_str(),fNBinsX,xbins);
  ResultsMap["mu"]->Sumw2();
  auto sigtext = "#sigma("+fTimeText+") [ns]";
  auto sigName = f2DHistName+"_sigma";
  auto sigTitle = fTitle+" "+sigtext+";"+fXTitle+";"+sigtext;
  ResultsMap["sigma"] = new TH1F(sigName.c_str(),sigTitle.c_str(),fNBinsX,xbins);
  ResultsMap["sigma"]->Sumw2();
  auto occtext = "#occ("+fTimeText+")";
  auto occName = f2DHistName+"_occ";
  auto occTitle = fTitle+" "+occtext+";"+fXTitle+";"+occtext;
  ResultsMap["occ"] = new TH1F(occName.c_str(),occTitle.c_str(),fNBinsX,xbins);
  ResultsMap["occ"]->Sumw2();
  auto rmstext = "#rms("+fTimeText+")";
  auto rmsName = f2DHistName+"_rms";
  auto rmsTitle = fTitle+" "+rmstext+";"+fXTitle+";"+rmstext;
  ResultsMap["rms"] = new TH1F(rmsName.c_str(),rmsTitle.c_str(),fNBinsX,xbins);
  ResultsMap["rms"]->Sumw2();

  std::cout << "Getting input hist: Data" << std::endl;

  std::string inHistName = "Data";
  TH2D * Hist2D = (TH2D*)fInFile->Get(f2DHistName.c_str());

  for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++){ 

 	TimeFitStruct* TimeFit = new TimeFitStruct(fTimeFitType,fRangeLow,fRangeUp);

  	const std::string histname1 = Hist2D->GetName();
	std::string fullHistName = histname1 + "_ibin" + std::to_string(ibinX);
	float effamp = Hist2D->GetXaxis()->GetBinCenter(ibinX);
    TimeFit->hist = (TH1F*)Hist2D->ProjectionY(fullHistName.c_str(),ibinX,ibinX);
	TimeFit->hist->SetTitle(fullHistName.c_str());

    TimeFit->PrepFit();
    TimeFit->DoFit();
    TimeFit->hist->Write(TimeFit->hist->GetName(),TObject::kWriteDelete);
    TimeFit->form->Write(TimeFit->form->GetName(),TObject::kWriteDelete);
    TimeFit->fit->Write(TimeFit->fit->GetName(),TObject::kWriteDelete);

    TimeFit->GetFitResult();
    const auto & result = TimeFit->result;
	auto errHistName = "bin" + std::to_string(ibinX) + "ErrHist";
    auto errHist = (TH2F*)fInFile->Get(errHistName.c_str());
	auto merr = errHist->GetMean();
	auto nerr = errHist->Integral();
    auto sgerr = errHist->GetStdDev();
    auto err = sgerr/sqrt(nerr);
    std::cout << " Bin " << ibinX  <<  " : " << fXBins[ibinX-1] << "-" << fXBins[ibinX]; 
	std::cout << " : no guass: " << result.sigma << " err: " << result.esigma << std::endl; 
    ResultsMap["chi2ndf"]->SetBinContent(ibinX,result.chi2ndf);
    ResultsMap["chi2prob"]->SetBinContent(ibinX,result.chi2prob);
    ResultsMap["mu"]->SetBinContent(ibinX,result.mu);
    ResultsMap["mu"]->SetBinError(ibinX,result.emu);
	float smear = effamp ? sq2(noise/effamp)+2*sq2(constant) : 0;
	float smearer = result.esigma;//todo
	float smsig  = smear ? std::sqrt(sq2(result.sigma)+smear) : result.sigma;
	float smsigerr = smear ? smearer : result.esigma;
    ResultsMap["sigma"]->SetBinContent(ibinX,smsig);
    ResultsMap["sigma"]->SetBinError(ibinX,smsigerr);
    ResultsMap["occ"]->SetBinContent(ibinX,result.occ);
    ResultsMap["occ"]->SetBinError(ibinX,sqrt(result.occ));
    ResultsMap["rms"]->SetBinContent(ibinX,result.rms);
    ResultsMap["rms"]->SetBinError(ibinX,err);

	TimeFit->DeleteInternal();
	delete TimeFit;

  }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)

  for (const auto & ResultsPair : ResultsMap) ResultsPair.second->Write(ResultsPair.second->GetName(),TObject::kWriteDelete);

	//-------------------------------------------------------------------------------------------
	//  Do sigma fit
	//-------------------------------------------------------------------------------------------

  std::cout << "Prepping sigma fit for: Data" << std::endl;

  const auto & hist = ResultsMap["sigma"];
  const auto x_low = hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetFirst());
  const auto x_up  = hist->GetXaxis()->GetBinUpEdge (hist->GetXaxis()->GetLast());
  const std::string histname = hist->GetName();

  const auto formname = histname+"_form";
  const auto fitname  = histname+"_fit";
  auto fitformstr = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))";
  TFormula* form = new TFormula(formname.c_str(),fitformstr);

  std::cout << "Prepping sigma fit params for Data" << std::endl;

  TF1* fit = new TF1(fitname.c_str(),form->GetName(),x_low,x_up);
  fit->SetParName(0,"N"); fit->SetParameter(0,50); fit->SetParLimits(0,0,100);
  fit->SetParName(1,"C"); fit->SetParameter(1,0.5); fit->SetParLimits(1,0,1);
  fit->SetLineColor(hist->GetLineColor());

  std::cout << "Saving sigma form for Data" << std::endl;

  fOutFile->cd();
  form->Write(form->GetName(),TObject::kWriteDelete);

  std::cout << "Fitting sigma hist for Data" << std::endl;
  hist->Fit(fit->GetName(),"RBQM");
  fOutFile->cd();
  fit->Write(fit->GetName(),TObject::kWriteDelete);

	//---------------------------------------------------------------------------------------
	// Clean Up
	//---------------------------------------------------------------------------------------

  std::cout << "Deleting info for Data" << std::endl;

  for (auto & Pair : ResultsMap) delete Pair.second;
  //ResultsMap.clear();

  delete Hist2D;
  delete fOutFile;
  delete fInFile;

}//<<>>void runTimeFitter(......

