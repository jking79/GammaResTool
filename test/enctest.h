//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr  3 11:55:13 2023 by ROOT version 6.24/07
// from TTree llpgtree/llpgtree
// found on file: ku_KUCC_tt_enctest_126_gammares_v11.root
//////////////////////////////////////////////////////////

#ifndef enctest_h
#define enctest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class enctest {
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
   TBranch        *b_unrhJitter;   //!
   TBranch        *b_unrhNonJitter;   //!
   TBranch        *b_unrhEncNonJitter;   //!
   TBranch        *b_unrhEnergy;   //!

   enctest(TTree *tree=0);
   virtual ~enctest();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef enctest_cxx
enctest::enctest(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ku_KUCC_tt_enctest_126_gammares_v11.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ku_KUCC_tt_enctest_126_gammares_v11.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ku_KUCC_tt_enctest_126_gammares_v11.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

enctest::~enctest()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t enctest::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t enctest::LoadTree(Long64_t entry)
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

void enctest::Init(TTree *tree)
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
   fChain->SetBranchAddress("unrhJitter", &unrhJitter, &b_unrhJitter);
   fChain->SetBranchAddress("unrhNonJitter", &unrhNonJitter, &b_unrhNonJitter);
   fChain->SetBranchAddress("unrhEncNonJitter", &unrhEncNonJitter, &b_unrhEncNonJitter);
   fChain->SetBranchAddress("unrhEnergy", &unrhEnergy, &b_unrhEnergy);
   Notify();
}

Bool_t enctest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void enctest::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t enctest::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef enctest_cxx
