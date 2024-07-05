//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  6 05:46:49 2024 by ROOT version 6.30/03
// from TTree llpgtree/llpgtree
// found on file: ku_24E_diag_140_gammares_v11.root
//////////////////////////////////////////////////////////

#ifndef thistree_h
#define thistree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class thistree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   vector<unsigned int> *rhCaliID;
   vector<float>   *rhCaliEnergy;
   vector<float>   *rhCaliRtTime;
   vector<float>   *rhCaliCCTime;
   vector<unsigned int> *resRhID;
   vector<float>   *resAmp;
   vector<float>   *resE;
   vector<float>   *resRtTime;
   vector<float>   *resCCTime;
   vector<float>   *resTOF;
   vector<unsigned int> *res1RhID;
   vector<float>   *res1Amp;
   vector<float>   *res1E;
   vector<float>   *res1RtTime;
   vector<float>   *res1CCTime;
   vector<float>   *res1TOF;
   vector<unsigned int> *res2RhID;
   vector<float>   *res2Amp;
   vector<float>   *res2E;
   vector<float>   *res2RtTime;
   vector<float>   *res2CCTime;
   vector<float>   *res2TOF;
   vector<unsigned int> *resZ1RhID;
   vector<float>   *resZ1Amp;
   vector<float>   *resZ1E;
   vector<float>   *resZ1RtTime;
   vector<float>   *resZ1CCTime;
   vector<float>   *resZ1TOF;
   vector<unsigned int> *resZ2RhID;
   vector<float>   *resZ2Amp;
   vector<float>   *resZ2E;
   vector<float>   *resZ2RtTime;
   vector<float>   *resZ2CCTime;
   vector<float>   *resZ2TOF;
   vector<unsigned int> *rhID;
   vector<float>   *rhRtTime;
   vector<float>   *rhTOF;
   vector<float>   *rhEnergy;
   vector<bool>    *rhRtisOOT;
   vector<bool>    *rhisWeird;
   vector<bool>    *rhisDiWeird;
   vector<float>   *rhSwCross;
   vector<float>   *rhAmp;
   vector<bool>    *rhisGS6;
   vector<bool>    *rhisGS1;
   vector<float>   *rhadcToGeV;
   vector<float>   *rhpedrms12;
   vector<float>   *phoEnergy;
   vector<vector<unsigned int> > *phoRhIds;
   vector<float>   *phoPt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoHadOverEM;
   vector<float>   *phoSigmaIEtaIEta;
   vector<float>   *phoCov2IEtaIEta;
   vector<float>   *phoCov2IEtaIPhi;
   vector<float>   *phoCov2IPhiIPhi;
   vector<float>   *phoEcalRHSumEtConeDR04;
   vector<float>   *phoHcalTwrSumEtConeDR04;
   vector<float>   *phoTrkSumPtSolidConeDR04;
   vector<float>   *phoTrkSumPtHollowConeDR04;
   vector<float>   *phoR9;
   vector<int>     *phoSelType;
   Float_t         phoDiMass;
   Float_t         phoDiAngle;
   Float_t         phoDiDr;
   Float_t         phoDiPhi;
   Float_t         phoDiEta;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_rhCaliID;   //!
   TBranch        *b_rhCaliEnergy;   //!
   TBranch        *b_rhCaliRtTime;   //!
   TBranch        *b_rhCaliCCTime;   //!
   TBranch        *b_resRhID;   //!
   TBranch        *b_resAmp;   //!
   TBranch        *b_resE;   //!
   TBranch        *b_resRtTime;   //!
   TBranch        *b_resCCTime;   //!
   TBranch        *b_resTOF;   //!
   TBranch        *b_res1RhID;   //!
   TBranch        *b_res1Amp;   //!
   TBranch        *b_res1E;   //!
   TBranch        *b_res1RtTime;   //!
   TBranch        *b_res1CCTime;   //!
   TBranch        *b_res1TOF;   //!
   TBranch        *b_res2RhID;   //!
   TBranch        *b_res2Amp;   //!
   TBranch        *b_res2E;   //!
   TBranch        *b_res2RtTime;   //!
   TBranch        *b_res2CCTime;   //!
   TBranch        *b_res2TOF;   //!
   TBranch        *b_resZ1RhID;   //!
   TBranch        *b_resZ1Amp;   //!
   TBranch        *b_resZ1E;   //!
   TBranch        *b_resZ1RtTime;   //!
   TBranch        *b_resZ1CCTime;   //!
   TBranch        *b_resZ1TOF;   //!
   TBranch        *b_resZ2RhID;   //!
   TBranch        *b_resZ2Amp;   //!
   TBranch        *b_resZ2E;   //!
   TBranch        *b_resZ2RtTime;   //!
   TBranch        *b_resZ2CCTime;   //!
   TBranch        *b_resZ2TOF;   //!
   TBranch        *b_rhID;   //!
   TBranch        *b_rhRtTime;   //!
   TBranch        *b_rhTOF;   //!
   TBranch        *b_rhEnergy;   //!
   TBranch        *b_rhRtisOOT;   //!
   TBranch        *b_rhisWeird;   //!
   TBranch        *b_rhisDiWeird;   //!
   TBranch        *b_rhSwCross;   //!
   TBranch        *b_rhAmp;   //!
   TBranch        *b_rhisGS6;   //!
   TBranch        *b_rhisGS1;   //!
   TBranch        *b_rhadcToGeV;   //!
   TBranch        *b_rhpedrms12;   //!
   TBranch        *b_phoEnergy;   //!
   TBranch        *b_phoRhIds;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoHadOverEM;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoCov2IEtaIEta;   //!
   TBranch        *b_phoCov2IEtaIPhi;   //!
   TBranch        *b_phoCov2IPhiIPhi;   //!
   TBranch        *b_phoEcalRHSumEtConeDR04;   //!
   TBranch        *b_phoHcalTwrSumEtConeDR04;   //!
   TBranch        *b_phoTrkSumPtSolidConeDR04;   //!
   TBranch        *b_phoTrkSumPtHollowConeDR04;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoSelType;   //!
   TBranch        *b_phoDiMass;   //!
   TBranch        *b_phoDiAngle;   //!
   TBranch        *b_phoDiDr;   //!
   TBranch        *b_phoDiPhi;   //!
   TBranch        *b_phoDiEta;   //!

   thistree(TTree *tree=0);
   virtual ~thistree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef thistree_cxx
thistree::thistree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ku_24E_diag_140_gammares_v11.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ku_24E_diag_140_gammares_v11.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ku_24E_diag_140_gammares_v11.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

thistree::~thistree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t thistree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t thistree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void thistree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   rhCaliID = 0;
   rhCaliEnergy = 0;
   rhCaliRtTime = 0;
   rhCaliCCTime = 0;
   resRhID = 0;
   resAmp = 0;
   resE = 0;
   resRtTime = 0;
   resCCTime = 0;
   resTOF = 0;
   res1RhID = 0;
   res1Amp = 0;
   res1E = 0;
   res1RtTime = 0;
   res1CCTime = 0;
   res1TOF = 0;
   res2RhID = 0;
   res2Amp = 0;
   res2E = 0;
   res2RtTime = 0;
   res2CCTime = 0;
   res2TOF = 0;
   resZ1RhID = 0;
   resZ1Amp = 0;
   resZ1E = 0;
   resZ1RtTime = 0;
   resZ1CCTime = 0;
   resZ1TOF = 0;
   resZ2RhID = 0;
   resZ2Amp = 0;
   resZ2E = 0;
   resZ2RtTime = 0;
   resZ2CCTime = 0;
   resZ2TOF = 0;
   rhID = 0;
   rhRtTime = 0;
   rhTOF = 0;
   rhEnergy = 0;
   rhRtisOOT = 0;
   rhisWeird = 0;
   rhisDiWeird = 0;
   rhSwCross = 0;
   rhAmp = 0;
   rhisGS6 = 0;
   rhisGS1 = 0;
   rhadcToGeV = 0;
   rhpedrms12 = 0;
   phoEnergy = 0;
   phoRhIds = 0;
   phoPt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoHadOverEM = 0;
   phoSigmaIEtaIEta = 0;
   phoCov2IEtaIEta = 0;
   phoCov2IEtaIPhi = 0;
   phoCov2IPhiIPhi = 0;
   phoEcalRHSumEtConeDR04 = 0;
   phoHcalTwrSumEtConeDR04 = 0;
   phoTrkSumPtSolidConeDR04 = 0;
   phoTrkSumPtHollowConeDR04 = 0;
   phoR9 = 0;
   phoSelType = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("rhCaliID", &rhCaliID, &b_rhCaliID);
   fChain->SetBranchAddress("rhCaliEnergy", &rhCaliEnergy, &b_rhCaliEnergy);
   fChain->SetBranchAddress("rhCaliRtTime", &rhCaliRtTime, &b_rhCaliRtTime);
   fChain->SetBranchAddress("rhCaliCCTime", &rhCaliCCTime, &b_rhCaliCCTime);
   fChain->SetBranchAddress("resRhID", &resRhID, &b_resRhID);
   fChain->SetBranchAddress("resAmp", &resAmp, &b_resAmp);
   fChain->SetBranchAddress("resE", &resE, &b_resE);
   fChain->SetBranchAddress("resRtTime", &resRtTime, &b_resRtTime);
   fChain->SetBranchAddress("resCCTime", &resCCTime, &b_resCCTime);
   fChain->SetBranchAddress("resTOF", &resTOF, &b_resTOF);
   fChain->SetBranchAddress("res1RhID", &res1RhID, &b_res1RhID);
   fChain->SetBranchAddress("res1Amp", &res1Amp, &b_res1Amp);
   fChain->SetBranchAddress("res1E", &res1E, &b_res1E);
   fChain->SetBranchAddress("res1RtTime", &res1RtTime, &b_res1RtTime);
   fChain->SetBranchAddress("res1CCTime", &res1CCTime, &b_res1CCTime);
   fChain->SetBranchAddress("res1TOF", &res1TOF, &b_res1TOF);
   fChain->SetBranchAddress("res2RhID", &res2RhID, &b_res2RhID);
   fChain->SetBranchAddress("res2Amp", &res2Amp, &b_res2Amp);
   fChain->SetBranchAddress("res2E", &res2E, &b_res2E);
   fChain->SetBranchAddress("res2RtTime", &res2RtTime, &b_res2RtTime);
   fChain->SetBranchAddress("res2CCTime", &res2CCTime, &b_res2CCTime);
   fChain->SetBranchAddress("res2TOF", &res2TOF, &b_res2TOF);
   fChain->SetBranchAddress("resZ1RhID", &resZ1RhID, &b_resZ1RhID);
   fChain->SetBranchAddress("resZ1Amp", &resZ1Amp, &b_resZ1Amp);
   fChain->SetBranchAddress("resZ1E", &resZ1E, &b_resZ1E);
   fChain->SetBranchAddress("resZ1RtTime", &resZ1RtTime, &b_resZ1RtTime);
   fChain->SetBranchAddress("resZ1CCTime", &resZ1CCTime, &b_resZ1CCTime);
   fChain->SetBranchAddress("resZ1TOF", &resZ1TOF, &b_resZ1TOF);
   fChain->SetBranchAddress("resZ2RhID", &resZ2RhID, &b_resZ2RhID);
   fChain->SetBranchAddress("resZ2Amp", &resZ2Amp, &b_resZ2Amp);
   fChain->SetBranchAddress("resZ2E", &resZ2E, &b_resZ2E);
   fChain->SetBranchAddress("resZ2RtTime", &resZ2RtTime, &b_resZ2RtTime);
   fChain->SetBranchAddress("resZ2CCTime", &resZ2CCTime, &b_resZ2CCTime);
   fChain->SetBranchAddress("resZ2TOF", &resZ2TOF, &b_resZ2TOF);
   fChain->SetBranchAddress("rhID", &rhID, &b_rhID);
   fChain->SetBranchAddress("rhRtTime", &rhRtTime, &b_rhRtTime);
   fChain->SetBranchAddress("rhTOF", &rhTOF, &b_rhTOF);
   fChain->SetBranchAddress("rhEnergy", &rhEnergy, &b_rhEnergy);
   fChain->SetBranchAddress("rhRtisOOT", &rhRtisOOT, &b_rhRtisOOT);
   fChain->SetBranchAddress("rhisWeird", &rhisWeird, &b_rhisWeird);
   fChain->SetBranchAddress("rhisDiWeird", &rhisDiWeird, &b_rhisDiWeird);
   fChain->SetBranchAddress("rhSwCross", &rhSwCross, &b_rhSwCross);
   fChain->SetBranchAddress("rhAmp", &rhAmp, &b_rhAmp);
   fChain->SetBranchAddress("rhisGS6", &rhisGS6, &b_rhisGS6);
   fChain->SetBranchAddress("rhisGS1", &rhisGS1, &b_rhisGS1);
   fChain->SetBranchAddress("rhadcToGeV", &rhadcToGeV, &b_rhadcToGeV);
   fChain->SetBranchAddress("rhpedrms12", &rhpedrms12, &b_rhpedrms12);
   fChain->SetBranchAddress("phoEnergy", &phoEnergy, &b_phoEnergy);
   fChain->SetBranchAddress("phoRhIds", &phoRhIds, &b_phoRhIds);
   fChain->SetBranchAddress("phoPt", &phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoHadOverEM", &phoHadOverEM, &b_phoHadOverEM);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoCov2IEtaIEta", &phoCov2IEtaIEta, &b_phoCov2IEtaIEta);
   fChain->SetBranchAddress("phoCov2IEtaIPhi", &phoCov2IEtaIPhi, &b_phoCov2IEtaIPhi);
   fChain->SetBranchAddress("phoCov2IPhiIPhi", &phoCov2IPhiIPhi, &b_phoCov2IPhiIPhi);
   fChain->SetBranchAddress("phoEcalRHSumEtConeDR04", &phoEcalRHSumEtConeDR04, &b_phoEcalRHSumEtConeDR04);
   fChain->SetBranchAddress("phoHcalTwrSumEtConeDR04", &phoHcalTwrSumEtConeDR04, &b_phoHcalTwrSumEtConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04, &b_phoTrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04, &b_phoTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoSelType", &phoSelType, &b_phoSelType);
   fChain->SetBranchAddress("phoDiMass", &phoDiMass, &b_phoDiMass);
   fChain->SetBranchAddress("phoDiAngle", &phoDiAngle, &b_phoDiAngle);
   fChain->SetBranchAddress("phoDiDr", &phoDiDr, &b_phoDiDr);
   fChain->SetBranchAddress("phoDiPhi", &phoDiPhi, &b_phoDiPhi);
   fChain->SetBranchAddress("phoDiEta", &phoDiEta, &b_phoDiEta);
   Notify();
}

Bool_t thistree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void thistree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t thistree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef thistree_cxx
