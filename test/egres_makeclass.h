//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 12 06:23:14 2024 by ROOT version 6.30/03
// from TTree llpgtree/llpgtree
// found on file: ku_twotier_140_gammares_test.root
//////////////////////////////////////////////////////////

#ifndef egres_makeclass_h
#define egres_makeclass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class egres_makeclass {
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
   vector<float>   *unrhJitter;
   vector<float>   *unrhNonJitter;
   vector<float>   *unrhEncNonJitter;
   vector<float>   *unrhEnergy;

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
   TBranch        *b_unrhJitter;   //!
   TBranch        *b_unrhNonJitter;   //!
   TBranch        *b_unrhEncNonJitter;   //!
   TBranch        *b_unrhEnergy;   //!

   egres_makeclass(TTree *tree=0);
   virtual ~egres_makeclass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef egres_makeclass_cxx
egres_makeclass::egres_makeclass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ku_twotier_140_gammares_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ku_twotier_140_gammares_test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ku_twotier_140_gammares_test.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

egres_makeclass::~egres_makeclass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t egres_makeclass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t egres_makeclass::LoadTree(Long64_t entry)
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

void egres_makeclass::Init(TTree *tree)
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
   unrhJitter = 0;
   unrhNonJitter = 0;
   unrhEncNonJitter = 0;
   unrhEnergy = 0;
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
   fChain->SetBranchAddress("unrhJitter", &unrhJitter, &b_unrhJitter);
   fChain->SetBranchAddress("unrhNonJitter", &unrhNonJitter, &b_unrhNonJitter);
   fChain->SetBranchAddress("unrhEncNonJitter", &unrhEncNonJitter, &b_unrhEncNonJitter);
   fChain->SetBranchAddress("unrhEnergy", &unrhEnergy, &b_unrhEnergy);
   Notify();
}

Bool_t egres_makeclass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void egres_makeclass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t egres_makeclass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef egres_makeclass_cxx
