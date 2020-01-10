#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPad.h"

#include <vector>
#include <map>
#include <iostream>
#include <assert.h>

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.cc"
#include "RA2bTree.cc"
#include "TriggerEfficiencySextet.cc"
#include "defaultArgs.h"
#include "TriggerCorrector.h"
using namespace std;

int main(int argc, char** argv){
    int reg_(0);
    skimSamples::region reg;
if( argc >= 2 ){
        reg_ = atoi(argv[1]);
std::cout<<"reg "<<reg_<<std::endl;
}
    bool looseCuts(false);

    defaultOptions options(argv[0],"");
    options.opts->add_options()("l,loose_cuts","apply loose jet pt cuts",cxxopts::value<bool>(looseCuts))("r,region","region to analyze",cxxopts::value<int>(reg_));
    options.opts->parse(argc, argv);

    reg = static_cast<skimSamples::region>(reg_);
    //int reg_(reg);
    TString Era(argv[2]);
    //std::cout<<"era "<<int(argv[1])<<std::endl;
    gROOT->ProcessLine(".L tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    
    skimSamples* skims_ = new skimSamples(reg,Era);

    typedef bool(*cuts)(RA2bTree*);
    vector<cuts> baselineCuts;

    if( looseCuts ){
        baselineCuts.push_back(*FiltersCut<RA2bTree>);
        if( reg == skimSamples::kLowDphi ){ 
            baselineCuts.push_back(*lowDPhiCuts<RA2bTree>);
        }else{
            baselineCuts.push_back(*DeltaPhiCuts<RA2bTree>);
        }
        if( reg == skimSamples::kSLm ){
            baselineCuts.push_back(*singleMuCut<RA2bTree>);
        }
        if( reg == skimSamples::kSLe ){
            baselineCuts.push_back(*singleEleCut<RA2bTree>);
        }
        baselineCuts.push_back(*METHTlooseCut<RA2bTree>);
        baselineCuts.push_back(*AK8MultCut<RA2bTree>);
    }else{
        if( reg == skimSamples::kSignal ){
            baselineCuts.push_back(*baselineCut<RA2bTree>);
        }else if( reg == skimSamples::kSLm ){
            baselineCuts.push_back(*singleMuBaselineCut<RA2bTree>);
        }else if( reg == skimSamples::kSLe ){
            baselineCuts.push_back(*singleEleBaselineCut<RA2bTree>);
        }else if( reg == skimSamples::kLowDphi ){ 
            baselineCuts.push_back(*lowDphiBaselineCut<RA2bTree>);
        }
	else if(reg==skimSamples::kPhoton){
            baselineCuts.push_back(*photonBaselineCut<RA2bTree>);
		//std::cout<<"Single Photon "<<std::endl;
	}
	else if(reg==skimSamples::kDYm){
            baselineCuts.push_back(*doubleMuBaselineCut<RA2bTree>);
		//std::cout<<"Single Photon "<<std::endl;
	}
	else if(reg==skimSamples::kDYe){
            baselineCuts.push_back(*doubleEleBaselineCut<RA2bTree>);
		//std::cout<<"Single Photon "<<std::endl;
	}
	else
            assert(1);
    }

    skimSamples skims = *skims_;
    TFile* outputFile;
    TString regionName;
    TString cutName="";
	TriggerCorrector trigcorror;
	TriggerCorrector trigcorrorHT;
	TriggerCorrector trigcorrorPhoBarrel;
	TriggerCorrector trigcorrorPhoEndcap;
        
	TriggerCorrector trigcorrorFakeMHT;
	if(Era=="MC2016"){
	trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2016");
	trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2016");
	trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2016");
        trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2016_JetHT");
        trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2016_JetHT");
 	}
	if(Era=="MC2017"){
	trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2017");
	trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2017");
	trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2017");
        trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2017_JetHT");
        trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2017_JetHT");
        }
	if(Era=="MC2018"){
	trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2018");
	trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2018");
	trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2018");
        trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2018_JetHT");
        trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2018_JetHT");
        }
    if( looseCuts )
        cutName="_looseCuts";
    if( reg == skimSamples::kSignal )
        regionName="Signal";
    if( reg == skimSamples::kSLm )
        regionName="_singleMu";
    if( reg == skimSamples::kSLe )
        regionName="_singleEle";
    if( reg == skimSamples::kLowDphi )
        regionName="_lowDphi";
    if(reg==skimSamples::kPhoton) regionName="_photon";
    if(reg==skimSamples::kDYe) regionName="_DYe";
    if(reg==skimSamples::kDYm) regionName="_DYm";
    outputFile = new TFile("SkimFileMass"+cutName+regionName+Era+".root","RECREATE");
    // background MC samples - 0 lepton regions
   for( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++){
   ///for( int iSample = 0 ; iSample < 0 ; iSample++){
// if(skims.sampleName[iSample]!="data" && skims.sampleName[iSample]!="data2017" && skims.sampleName[iSample]!="data2018" )continue; 
//	continue;
	//if(skims.sampleName[iSample]!="ZJets")continue;
	//if(!skims.sampleName[iSample].Contains("data2018"))continue;
       RA2bTree* ntuple = skims.ntuples[iSample];
 	//TTree*newtree=(TTree*)ntuple->fChain->CloneTree(0);
 	TTree*newtree=new TTree("newtree","");//(TTree*)ntuple->fChain->CloneTree(0);
	int BTags;	
        double MET,HT,Weight,JetPt1, JetPt2, JetPhi1, JetPhi2,PrunedMass1, PrunedMass2, Jet1_tau2overtau1, Jet2_tau2overtau1;
	double JetEta1,JetEta2;
	double DeltaRtoClosestB;
        int NJets;
        //TBranch*b_BTags, *b_Weight,*b_MET,*b_HT,*b_JetPt1, *b_JetPt2,*b_PrunedMass1, *b_PrunedMass2, *b_Jet1_tau2overtau1, *b_Jet2_tau2overtau1, *b_GenHadTau;
    	int WMatchedJet1, WMatchedJet2;
	int GenHadTau=0;
	int nAK8=0;
        newtree->Branch("WMatchedJet1", &WMatchedJet1, "WMatchedJet1/I");	
        newtree->Branch("WMatchedJet2", &WMatchedJet2, "WMatchedJet2/I");	
        newtree->Branch("JetPhi1", &JetPhi1, "JetPhi1/D");	
        newtree->Branch("JetPhi2", &JetPhi2, "JetPhi2/D");	
        newtree->Branch("JetPt1", &JetPt1, "JetPt1/D");	
        newtree->Branch("JetPt2", &JetPt2, "JetPt2/D");	
        newtree->Branch("JetEta1", &JetEta1, "JetEta1/D");	
        newtree->Branch("JetEta2", &JetEta2, "JetEta2/D");	
        newtree->Branch("PrunedMass1", &PrunedMass1, "PrunedMass1/D");	
        newtree->Branch("PrunedMass2", &PrunedMass2, "PrunedMass2/D");	
        newtree->Branch("Jet1_tau2overtau1", &Jet1_tau2overtau1, "Jet1_tau2overtau1/D");	
        newtree->Branch("Jet2_tau2overtau2", &Jet2_tau2overtau1, "Jet2_tau2overtau1/D");	
        newtree->Branch("Evtweight",&Weight, "Evtweight/D");  
        newtree->Branch("GenHadTau", &GenHadTau, "GenHadTau/I");
        newtree->Branch("nAK8", &nAK8, "nAK8/I");
	newtree->Branch("HT", &HT, "HT/D");
	newtree->Branch("NJets", &NJets, "NJets/I");
        //newtree->SetBranchAddress("BTags",&BTags,&b_BTags);
        newtree->Branch("MET",&MET, "MET/D");
        newtree->Branch("DeltaRtoClosestB", &DeltaRtoClosestB,"DeltaRtoClosestB/D");
        
        int numEvents = ntuple->fChain->GetEntries();
        ntupleBranchStatus<RA2bTree>(ntuple);
        int bin = -1;
        double weight=0.;
        float trigWeight=1.0;
	float trigunc=0.0;
        bool passBaseline;
        double jetMass1,jetMass2;
        TString filename;
        cout << skims.sampleName[iSample]<<numEvents <<endl;
    for( int iEvt = 0 ; iEvt < min(options.MAX_EVENTS,numEvents) ; iEvt++ ){
    //for( int iEvt = 0 ; iEvt < min(10,numEvents) ; iEvt++ ){
            ntuple->GetEntry(iEvt);
            if( iEvt % 100000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << min(options.MAX_EVENTS,numEvents) << endl;
             passBaseline=true;
             for( auto baselineCut : baselineCuts ){
             passBaseline&=baselineCut(ntuple);
             }
            if( ! passBaseline ) continue;
	    //std::cout<<"Pass Baseline "<<std::endl; 	    
            filename = ntuple->fChain->GetFile()->GetName();
            if( ( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. )continue;
            bin = -1;
            if(reg == skimSamples::kSignal ){
                //std::vector<double> EfficiencyCenterUpDown = Eff_MetMhtSextetReal_CenterUpDown(ntuple->HT, ntuple->MHT, ntuple->NJets);
                //trigWeight=EfficiencyCenterUpDown[0];
	        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
            }else if( reg == skimSamples::kSLm ){
                //trigWeight=singleMuonTrigWeights(ntuple);
	        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
            }else if( reg == skimSamples::kSLe ){
                //trigWeight=singleElectronTrigWeights(ntuple);
	        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
            }else if( reg == skimSamples::kLowDphi ){
                //trigWeight=lowDphiTrigWeights(ntuple);
	        trigWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
	    }
	    double prefireweight=1.0;
	    if( filename.Contains("2017") && !( skims.sampleName[iSample].Contains("data")))prefireweight=ntuple->NonPrefiringProb;
	    if(filename.Contains("2016"))lumi=35922.;
	    if(filename.Contains("2017"))lumi=41529.;
	    if(filename.Contains("2018"))lumi=59740.;
	   if(skims.sampleName[iSample]=="data" || skims.sampleName[iSample]=="data2017" || skims.sampleName[iSample]=="data2018")	trigWeight=1.0;
	    weight = ntuple->Weight*lumi*prefireweight*trigWeight;//*customPUweights(ntuple)*trigWeight;
	    //weight = ntuple->Weight *lumi*trigWeight*customPUweights(ntuple);    
	    //std::cout<<"Weight "<<ntuple->Weight<<std::endl;
            Weight=weight;
	    MET=ntuple->MET;
	    HT=ntuple->HT;
	    BTags=ntuple->BTags;
	    NJets=ntuple->NJets;
	    nAK8=ntuple->JetsAK8->size();
	    DeltaRtoClosestB=dRtoClosestB(ntuple);
	    if(nAK8>0){
            JetPhi1=ntuple->JetsAK8->at(0).Phi();  
            JetPt1=ntuple->JetsAK8->at(0).Pt();  
            JetEta1=ntuple->JetsAK8->at(0).Eta();  
	    PrunedMass1=ntuple->JetsAK8_softDropMass->at(0);
	   Jet1_tau2overtau1=ntuple->JetsAK8_NsubjettinessTau2->at(0)/ntuple->JetsAK8_NsubjettinessTau1->at(0);
	   }
	    if(nAK8>1){
            JetPhi2=ntuple->JetsAK8->at(1).Phi();
            JetPt2=ntuple->JetsAK8->at(1).Pt();
            JetEta2=ntuple->JetsAK8->at(1).Eta();
	    PrunedMass2=ntuple->JetsAK8_softDropMass->at(1);
	    Jet2_tau2overtau1=ntuple->JetsAK8_NsubjettinessTau2->at(1)/ntuple->JetsAK8_NsubjettinessTau1->at(1);
	    
	    }
	    //std::cout<<"MET"<<MET<<std::endl;
	WMatchedJet1=0;
	WMatchedJet2=0;
	if(skims.sampleName[iSample]!="data" && skims.sampleName[iSample]!="data2017" && skims.sampleName[iSample]!="data2018"){
	for( int i=0 ; i < ntuple->GenParticles->size() ; i++ ){
                        if(abs(ntuple->GenParticles_PdgId->at(i))>0 && abs(ntuple->GenParticles_PdgId->at(i))<5 && abs(ntuple->GenParticles_ParentId->at(i))==24){ 
      	//if( abs(ntuple->GenParticles_PdgId->at(i)) == 24 && ntuple->GenParticles->at(i).Pt()>200){ //&& ntuple->JetsAK8->at(0).DeltaR(ntuple->GenParticles->at(i))<0.4){
			if(nAK8<1)continue;
			float deta=ntuple->JetsAK8->at(0).Eta()-ntuple->GenParticles->at(i).Eta();
			float dphi=ntuple->JetsAK8->at(0).Phi()-ntuple->GenParticles->at(i).Phi();
			if(sqrt((deta*deta)+(dphi*dphi))<0.4)WMatchedJet1=1;
	    		if(nAK8<2)continue;
			deta=ntuple->JetsAK8->at(1).Eta()-ntuple->GenParticles->at(i).Eta();
			dphi=ntuple->JetsAK8->at(1).Phi()-ntuple->GenParticles->at(i).Phi();
			if(sqrt((deta*deta)+(dphi*dphi))<0.4)WMatchedJet2=1;
		}
    	     }
	}
	   if(skims.sampleName[iSample]!="data" && skims.sampleName[iSample]!="data2017" && skims.sampleName[iSample]!="data2018"){
		bool HadTau=false;
		GenHadTau=0;
	   	for(unsigned int t=0; t<ntuple->GenTaus_had->size(); ++t) if(ntuple->GenTaus_had->at(t))++GenHadTau ;
	   } 
		newtree->Fill();
        }// end event loop
	//newtree->SetName(filename);
	outputFile->cd();
	newtree->Write(skims.sampleName[iSample]);
	  
  }// end sample loop
if(reg == skimSamples::kSignal ){
    for( int iSample = 0 ; iSample < skims.signalNtuples.size() ; iSample++){

  
        RA2bTree* ntuple = skims.signalNtuples[iSample];
 	TTree*newtree=new TTree("newtree", "");//ntuple->fChain->CloneTree(0);
	int BTags;	
        double MET,HT,Weight,JetPt1, JetPt2,JetPhi1, JetPhi2,JetEta1, JetEta2,PrunedMass1, PrunedMass2, Jet1_tau2overtau1, Jet2_tau2overtau1;
        //TBranch*b_BTags, *b_Weight,*b_MET,*b_JetPt1, *b_JetPt2,*b_PrunedMass1, *b_PrunedMass2, *b_Jet1_tau2overtau1, *b_Jet2_tau2overtau1, *b_GenHadTau;
    	int WMatchedJet1, WMatchedJet2;
	int NJets;
	int GenHadTau=0;
	int nAK8=0;
	double DeltaRtoClosestB;
        double ZpT;
        newtree->Branch("WMatchedJet1", &WMatchedJet1, "WMatchedJet1/I");	
       newtree->Branch("WMatchedJet2", &WMatchedJet2, "WMatchedJet2/I");	
        newtree->Branch("JetEta1", &JetEta1, "JetEta1/D");	
        newtree->Branch("JetEta2", &JetEta2, "JetEta2/D");	
        newtree->Branch("JetPhi1", &JetPhi1, "JetPhi1/D");	
        newtree->Branch("JetPhi2", &JetPhi2, "JetPhi2/D");	
        newtree->Branch("JetPt1", &JetPt1, "JetPt1/D");	
        newtree->Branch("JetPt2", &JetPt2, "JetPt2/D");	
        newtree->Branch("PrunedMass1", &PrunedMass1, "PrunedMass1/D");	
        newtree->Branch("PrunedMass2", &PrunedMass2, "PrunedMass2/D");	
        newtree->Branch("Jet1_tau2overtau1", &Jet1_tau2overtau1, "Jet1_tau2overtau1/D");	
        newtree->Branch("Jet2_tau2overtau2", &Jet2_tau2overtau1, "Jet2_tau2overtau1/D");	
        newtree->Branch("Evtweight",&Weight, "Evtweight/D");  
        newtree->Branch("GenHadTau", &GenHadTau, "GenHadTau/I");
        newtree->Branch("nAK8", &nAK8, "nAK8/I");
	newtree->Branch("HT", &HT, "HT/D");
        newtree->Branch("DeltaRtoClosestB", &DeltaRtoClosestB,"DeltaRtoClosestB/D");
        //newtree->SetBranchAddress("BTags",&BTags,&b_BTags);
        newtree->Branch("MET",&MET, "MET/D");
	newtree->Branch("ZpT", &ZpT, "ZpT/D");
	newtree->Branch("NJets", &NJets, "NJets/I");
         
        //newtree->SetBranchAddress("BTags",&BTags,&b_BTags);
        //newtree->SetBranchAddress("MET",&MET, &b_MET);
        int numEvents = ntuple->fChain->GetEntries();
        ntupleBranchStatus<RA2bTree>(ntuple);
        int bin = -1;
        double weight=0.;
        float trigWeight=1.0;
        float trigunc;
        bool passBaseline;
        double jetMass1,jetMass2;
        TString filename;
    for( int iEvt = 0 ; iEvt < min(options.MAX_EVENTS,numEvents) ; iEvt++ ){
    //  for( int iEvt = 0 ; iEvt <1000; iEvt++ ){
            ntuple->GetEntry(iEvt);
            if( iEvt % 100000 == 0 ) cout << skims.signalSampleName[iSample] << ": " << iEvt << "/" << min(options.MAX_EVENTS,numEvents) << endl;
	                 passBaseline=true;
	    for( auto baselineCut : baselineCuts ){
             passBaseline&=baselineCut(ntuple);
            }
            if( ! passBaseline ) continue;
	    if(getNumGenHiggses(ntuple)>0)continue;
		//std::cout<<"Here ZZ "<<std::endl;	
               // std::vector<double> EfficiencyCenterUpDown = Eff_MetMhtSextetReal_CenterUpDown(ntuple->HT, ntuple->MHT, ntuple->NJets);
               // trigWeight=EfficiencyCenterUpDown[0];
	    if(Era=="MC2016")lumi=35922.;
	    if(Era=="MC2017")lumi=41529.;
	    if(Era=="MC2018")lumi=59740.;
	    double prefireweight=1.0;
	    if( Era=="MC2017")prefireweight=ntuple->NonPrefiringProb;
	     trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
	    double isrweight=SignalISRCorrection(ntuple);	    
	    weight=isrweight*ntuple->Weight*lumi*prefireweight*trigWeight/0.25;
            Weight=weight;
	    HT=ntuple->HT;
	    MET=ntuple->MET;
	    BTags=ntuple->BTags;
	    DeltaRtoClosestB=dRtoClosestB(ntuple);
	    nAK8=ntuple->JetsAK8->size();
	    NJets=ntuple->NJets;
	    if(nAK8>0){
            JetEta1=ntuple->JetsAK8->at(0).Eta();  
            JetPhi1=ntuple->JetsAK8->at(0).Phi();  
            JetPt1=ntuple->JetsAK8->at(0).Pt();  
	    PrunedMass1=ntuple->JetsAK8_softDropMass->at(0);
	   Jet1_tau2overtau1=ntuple->JetsAK8_NsubjettinessTau2->at(0)/ntuple->JetsAK8_NsubjettinessTau1->at(0);
	   }
	    if(nAK8>1){
            JetEta2=ntuple->JetsAK8->at(1).Eta();  
            JetPhi2=ntuple->JetsAK8->at(1).Phi();
            JetPt2=ntuple->JetsAK8->at(1).Pt();
	    PrunedMass2=ntuple->JetsAK8_softDropMass->at(1);
	    Jet2_tau2overtau1=ntuple->JetsAK8_NsubjettinessTau2->at(1)/ntuple->JetsAK8_NsubjettinessTau1->at(1);
	    }
	std::vector<unsigned int >ZBosonIndex;
    	for( int i=0 ; i < ntuple->GenParticles->size() ; i++ ){
        	if( ntuple->GenParticles_PdgId->at(i) == 23 &&
        	    ntuple->GenParticles_ParentId->at(i) == 1000023 &&
        	    ntuple->GenParticles_Status->at(i) == 22 )ZBosonIndex.push_back(i);
		
        	//    numZs++;
    		} 
	ZpT=ntuple->GenParticles->at(ZBosonIndex.at(0)).Pt();
        if(ZpT>ntuple->GenParticles->at(ZBosonIndex.at(1)).Pt())ZpT=ntuple->GenParticles->at(ZBosonIndex.at(1)).Pt();
        //std::cout<<"Z boson int "<<ntuple->GenParticles->at(ZBosonIndex.at(0)).Pt()<< " Mother " <<ntuple->GenParticles->at(ntuple->GenParticles_ParentIdx->at(ZBosonIndex.at(0))).Pt()<<std::endl;      
        newtree->Fill();
	}
        outputFile->cd();
        newtree->Write(skims.signalSampleName[iSample]);
   }
}
    outputFile->Close();

}
