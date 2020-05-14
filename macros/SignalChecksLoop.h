//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr  6 16:05:27 2020 by ROOT version 6.18/04
// from TTree newtree/
// found on file: SkimFileMassSignalMC2016.root
//////////////////////////////////////////////////////////

#ifndef SignalChecksLoop_h
#define SignalChecksLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class SignalChecksLoop {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           WMatchedJet1;
   Int_t           WMatchedJet2;
   Double_t        JetEta1;
   Double_t        JetEta2;
   Double_t        JetPhi1;
   Double_t        JetPhi2;
   Double_t        JetPt1;
   Double_t        JetPt2;
   Double_t        PrunedMass1;
   Double_t        PrunedMass2;
   Double_t        Jet1_tau2overtau1;
   Double_t        Jet2_tau2overtau2;
   Double_t        Evtweight;
   Int_t           GenHadTau;
   Int_t           nAK8;
    Int_t          BTagsdeep;
   Double_t        HT;
   Double_t        DeltaRtoClosestB;
      Double_t        DeltaRtoClosestBLead;
   Double_t        MET;
   Double_t        ZpT1;
   Double_t        ZpT2;
    Int_t           ZHad;

   Int_t           NJets;
   Int_t           NVtx;
               
   // List of branches
   TBranch        *b_WMatchedJet1;   //!
   TBranch        *b_WMatchedJet2;   //!
   TBranch        *b_JetEta1;   //!
   TBranch        *b_JetEta2;   //!
   TBranch        *b_JetPhi1;   //!
   TBranch        *b_JetPhi2;   //!
   TBranch        *b_JetPt1;   //!
   TBranch        *b_JetPt2;   //!
   TBranch        *b_PrunedMass1;   //!
   TBranch        *b_PrunedMass2;   //!
   TBranch        *b_Jet1_tau2overtau1;   //!
   TBranch        *b_Jet2_tau2overtau1;   //!
   TBranch        *b_Evtweight;   //!
   TBranch        *b_GenHadTau;   //!
   TBranch        *b_nAK8;   //!
    TBranch        *b_BTagsdeep;   //!

   TBranch        *b_HT;   //!
   TBranch        *b_DeltaRtoClosestB;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_ZpT1;   //!
    TBranch        *b_ZpT2;   //!
    TBranch *b_ZHad;
   TBranch        *b_NJets;   //!
   TBranch        *b_NVtx;   //!
    TBranch        *b_DeltaRtoClosestBLead;   //!

    TString treename_;

   SignalChecksLoop(TString treename,TChain *tree=0);
   virtual ~SignalChecksLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SignalChecksLoop_cxx
SignalChecksLoop::SignalChecksLoop(TString treename,TChain *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
    treename_=treename;

    tree=new TChain(treename_);
    tree->Add("BVeto/SkimFileMassSignalMC2016.root");
    tree->Add("BVeto/SkimFileMassSignalMC2017.root");
    tree->Add("BVeto/SkimFileMassSignalMC2018.root");
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SkimFileMassSignalMC2016.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SkimFileMassSignalMC2016.root");
      }
      f->GetObject("newtree",tree);

   }
 */
   Init(tree);
}

SignalChecksLoop::~SignalChecksLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SignalChecksLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SignalChecksLoop::LoadTree(Long64_t entry)
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

void SignalChecksLoop::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
    fChain->SetBranchAddress("DeltaRtoClosestBLead", &DeltaRtoClosestBLead, &b_DeltaRtoClosestBLead);
    fChain->SetBranchAddress("DeltaRtoClosestB", &DeltaRtoClosestB, &b_DeltaRtoClosestB);
   fChain->SetBranchAddress("WMatchedJet1", &WMatchedJet1, &b_WMatchedJet1);
   fChain->SetBranchAddress("WMatchedJet2", &WMatchedJet2, &b_WMatchedJet2);
   fChain->SetBranchAddress("JetEta1", &JetEta1, &b_JetEta1);
   fChain->SetBranchAddress("JetEta2", &JetEta2, &b_JetEta2);
   fChain->SetBranchAddress("JetPhi1", &JetPhi1, &b_JetPhi1);
   fChain->SetBranchAddress("JetPhi2", &JetPhi2, &b_JetPhi2);
   fChain->SetBranchAddress("JetPt1", &JetPt1, &b_JetPt1);
   fChain->SetBranchAddress("JetPt2", &JetPt2, &b_JetPt2);
   fChain->SetBranchAddress("PrunedMass1", &PrunedMass1, &b_PrunedMass1);
   fChain->SetBranchAddress("PrunedMass2", &PrunedMass2, &b_PrunedMass2);
   fChain->SetBranchAddress("Jet1_tau2overtau1", &Jet1_tau2overtau1, &b_Jet1_tau2overtau1);
   fChain->SetBranchAddress("Jet2_tau2overtau2", &Jet2_tau2overtau2, &b_Jet2_tau2overtau1);
   fChain->SetBranchAddress("Evtweight", &Evtweight, &b_Evtweight);
   fChain->SetBranchAddress("GenHadTau", &GenHadTau, &b_GenHadTau);
   fChain->SetBranchAddress("nAK8", &nAK8, &b_nAK8);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("DeltaRtoClosestB", &DeltaRtoClosestB, &b_DeltaRtoClosestB);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("ZpT1", &ZpT1, &b_ZpT1);
    fChain->SetBranchAddress("ZpT2", &ZpT2, &b_ZpT2);
    fChain->SetBranchAddress("ZHad",&ZHad,&b_ZHad);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
    fChain->SetBranchAddress("BTagsdeep", &BTagsdeep, &b_BTagsdeep);

   Notify();
}

Bool_t SignalChecksLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SignalChecksLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SignalChecksLoop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SignalChecksLoop_cxx
