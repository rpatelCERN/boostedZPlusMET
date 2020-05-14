#define SignalChecksLoop_cxx
#include "SignalChecksLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void SignalChecksLoop::Loop()
{
//   In a ROOT session, you can do:
//      root> .L SignalChecksLoop.C
//      root> SignalChecksLoop t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
    TH1D*hSoftDrop=new TH1D("hSoftDrop",";Z-boson p_{T};Efficiency", 100,40,140);
    float METBins[7]={300,450,600,800,1000, 1200, 2000    };

    TH1D*hSoftDropLowPU=new TH1D("hSoftDropLowPU",";Z-boson p_{T};Efficiency", 100,40,140);
    TH1D*hSoftDropHighPU=new TH1D("hSoftDropHighPU",";Z-boson p_{T};Efficiency", 100,40,140);

    TH1D*METPlotSR_ANBinsHighPU=new TH1D("METPlotSR_ANBinsHighPU","",6, METBins);
    TH1D*METPlotSR_ANBinsLowPU=new TH1D("METPlotSR_ANBinsLowPU","",6, METBins);

    
    TH1D*hGenZ=new TH1D("hGenZ",";Z-boson p_{T};Efficiency", 28,200, 1600);
    TH1D*hGenZJetFound=new TH1D("hGenZJetFound",";Z-boson p_{T};Efficiency", 28,200, 1600);
    TH1D*hGenZJetMassWindow=new TH1D("hGenZJetMassWindow",";Z-boson p_{T};Efficiency", 28,200, 1600);
    TH1D*DeltaRtoBSub=new TH1D("DeltaRtoBSub","",31, 0,3.1);
    TH1D*DeltaRtoBLead=new TH1D("DeltaRtoBLead","",31, 0,3.1);
    TH1D*MinDeltaR=new TH1D("MinDeltaR","",31, 0,3.1);
    TH1D*NBtags=new TH1D("NBtags","", 5, 0, 5);
    TH2D*DeltaRPlane=new TH2D("DeltaRPlane", "", 31,0,3.1, 31, 0,3.1);
   Long64_t nentries = fChain->GetEntriesFast();
    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       if(nAK8<2 || MET<300)continue;
       
       
       if(DeltaRtoClosestB<3)DeltaRtoBSub->Fill(DeltaRtoClosestB);
        else DeltaRtoBSub->Fill(3);
        if(DeltaRtoClosestBLead<3)DeltaRtoBLead->Fill(DeltaRtoClosestBLead);
        else DeltaRtoBLead->Fill(3);
       
       if(DeltaRtoClosestB>3)DeltaRtoClosestB=3;
       if(DeltaRtoClosestBLead>3)DeltaRtoClosestBLead=3;
       
       if(DeltaRtoClosestBLead>DeltaRtoClosestB)MinDeltaR->Fill(DeltaRtoClosestB);
       else MinDeltaR->Fill(DeltaRtoClosestBLead);
       DeltaRPlane->Fill(DeltaRtoClosestB,DeltaRtoClosestBLead);
       NBtags->Fill(BTagsdeep);
       // if (Cut(ientry) < 0) continue;
       if(DeltaRtoClosestB>0.8 && PrunedMass2>70 && PrunedMass2<100 )hSoftDrop->Fill(PrunedMass1,Evtweight);
      // if(DeltaRtoClosestB>0.8)hSoftDrop->Fill(PrunedMass1);
       if(NVtx>20){
            hSoftDropHighPU->Fill(PrunedMass1);
           if(PrunedMass1>70 && PrunedMass1<100 && PrunedMass2>70 && PrunedMass2<100  )METPlotSR_ANBinsHighPU->Fill(MET);
       }
       else {
          hSoftDropLowPU->Fill(PrunedMass1);
             if(PrunedMass1>70 && PrunedMass1<100 && PrunedMass2>70 && PrunedMass2<100  )METPlotSR_ANBinsLowPU->Fill(MET);

       }
       if(ZHad==4 && ZpT1>0 && ZpT2>0){
           hGenZ->Fill(ZpT1);
           hGenZ->Fill(ZpT2);
           if(WMatchedJet1>0)hGenZJetFound->Fill(ZpT1);
           if(WMatchedJet2>0)hGenZJetFound->Fill(ZpT2);
           if(WMatchedJet1>0 && PrunedMass1>70&& PrunedMass1<100)hGenZJetMassWindow->Fill(ZpT1);
           if(WMatchedJet2>0 && PrunedMass2>70&& PrunedMass2<100)hGenZJetMassWindow->Fill(ZpT2);
           
       }
   }
    TFile*fout=new TFile("SignalEfficiencyChecks"+treename_+".root","RECREATE");
    hSoftDropHighPU->Write("SoftDropMassShapeHighPU");
    
    hSoftDropLowPU->Write("SoftDropMassShapeLowPU");
    hSoftDrop->Write("LeadJetSoftMass");
    METPlotSR_ANBinsHighPU->Write("METShapeHighPU");
    METPlotSR_ANBinsLowPU->Write("METShapeLowPU");

    hGenZJetFound->Write("ZMatchEffNum");
    hGenZ->Write("ZDenom");
    hGenZJetMassWindow->Write("ZMatchEffMassWindowNum");
    DeltaRtoBSub->Write(treename_+"DeltaRtoBSub");
    DeltaRtoBLead->Write(treename_+"DeltaRtoBLead");
    MinDeltaR->Write("MinDR");
    DeltaRPlane->Write("DeltaRPlane");
    NBtags->Write("Btags");
}
