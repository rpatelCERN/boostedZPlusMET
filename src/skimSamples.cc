// ROOT/custom libraries
#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"
#include "TH1F.h"
// STL libraries
#include <iostream>
#include <vector>
static const TString BASE_DIRMC="root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV17/";
static const TString BASE_DIR="root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV17/";

class skimSamples{

public : 

    TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*GJets,*GJets0p4,*Other,*DY,*TTinc;
    TChain *T5HH750, *T5HH1000, *T5HH1100,*T5HH1200,*T5HH1300,*T5HH1400,*T5HH1500,*T5HH1600,*T5HH1700,*T5HH1800,*T5HH1900,*T5HH2000,*T5HH2100;
    TChain*T5HH2200, *T5HH2300,*T5HH2400,*T5HH2500;
    TChain *data;
    TChain *data2017;
    TChain *data2018;
    std::vector<RA2bTree*> ntuples,signalNtuples;
    RA2bTree* dataNtuple; RA2bTree* dataNtuple2017; RA2bTree* dataNtuple2018;
    std::vector<TString> sampleName, signalSampleName;
    std::vector<TString> dataSampleName; 
    std::vector<int> fillColor, lineColor, sigLineColor;
    std::vector<int> NSignalEvents;
    enum region {kSignal,kSLm,kSLe,kLowDphi,kPhoton,kDYe, kDYm,kNumRegions};
//    enum eras {2016,2017,2018,3};
    TString regionNames[kNumRegions]={"signal","SLm","SLe","kLowDphi","photon","DYe", "DYm" };

    TString skimType;
    TString skimTypeMC;

    skimSamples(region r=kSignal, TString Era="MC2016"){

        skimType="";

        if( r == kSignal ){
            skimType=BASE_DIR+"/tree_signalUnblind/";
            skimTypeMC=BASE_DIRMC+"tree_signal/";
        }
        if( r == kPhoton ){
            skimType=BASE_DIR+"tree_GJet_CleanVars/";
            skimTypeMC=BASE_DIRMC+"tree_GJet_CleanVars/";
        }
        if( r == kSLm ){
            skimType=BASE_DIR+"tree_SLm/";
           skimTypeMC=BASE_DIRMC+"tree_SLm/";
        }
        if( r == kSLe ){
            skimType=BASE_DIR+"tree_SLe/";
            skimTypeMC=BASE_DIRMC+"tree_SLe/";
        }
        if(r==kLowDphi){
            skimType=BASE_DIR+"tree_LDP/";
            skimTypeMC=BASE_DIRMC+"tree_LDP/";
        }
        if(r==kDYe){
            skimType=BASE_DIR+"tree_DYe_CleanVars/";
            skimTypeMC=BASE_DIRMC+"tree_DYe_CleanVars/";
        }
        if(r==kDYm){
            skimType=BASE_DIR+"tree_DYm_CleanVars/";
            skimTypeMC=BASE_DIRMC+"tree_DYm_CleanVars/";
        }
        ///////////////////////////////////////////////////////////////////////
        // - - - - - - - - - - BACKGROUND INPUTS - - - - - - - - - - - - - - //
        ///////////////////////////////////////////////////////////////////////
//ONLY DATA IS AVAILABLE FOR NOW!!
        std::vector<TString> OtherFileNames;
        OtherFileNames.push_back("tree_WWTo1L1Nu2Q_");
        OtherFileNames.push_back("tree_WWTo2L2Nu_");
        OtherFileNames.push_back("tree_WWZ_");
        OtherFileNames.push_back("tree_WZTo1L1Nu2Q_");
        OtherFileNames.push_back("tree_WZTo1L3Nu_");
        OtherFileNames.push_back("tree_WZZ_");
        OtherFileNames.push_back("tree_ZZTo2L2Q_");
        OtherFileNames.push_back("tree_ZZTo2Q2Nu_");
        OtherFileNames.push_back("tree_ZZZ_");
        OtherFileNames.push_back("tree_TTTT_");
        OtherFileNames.push_back("tree_TTWJetsToLNu_");
        OtherFileNames.push_back("tree_TTWJetsToQQ_");
        OtherFileNames.push_back("tree_TTGJets_");
        OtherFileNames.push_back("tree_TTZToLLNuNu_");
        OtherFileNames.push_back("tree_TTZToQQ_");
        Other = new TChain("tree");
        for( int i = 0 ; i < OtherFileNames.size() ; i++ ){
            Other->Add(skimTypeMC+"/"+OtherFileNames[i]+Era+".root");
        }
        if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi || r == kPhoton ){
            //ntuples.push_back(new RA2bTree(Other));
            //sampleName.push_back("Other");
            fillColor.push_back(kRed+1);
            lineColor.push_back(1);
        }

        std::vector<TString> ZJetsFileNames;
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200_");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400_");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600_");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800_");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200_");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500_");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf_");
        ZJets = new TChain("tree");
        for( int i = 0 ; i < ZJetsFileNames.size() ; i++ ){
            ZJets->Add(skimTypeMC+"/"+ZJetsFileNames[i]+Era+".root");
        }
        if( r == kSignal || r == kLowDphi ){ 
            //ntuples.push_back(new RA2bTree(ZJets));
            //sampleName.push_back("ZJets");
            fillColor.push_back(kGreen+1);
            lineColor.push_back(1);
        }

        std::vector<TString> WJetsFileNames;
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200_");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500_");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400_");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf_");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600_");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800_");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200_");
        WJets = new TChain("tree");
        for( int i = 0 ; i < WJetsFileNames.size() ; i++ ){
            WJets->Add(skimTypeMC+"/"+WJetsFileNames[i]+Era+".root");
        }
        if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
            //ntuples.push_back(new RA2bTree(WJets));
           // sampleName.push_back("WJets");
            fillColor.push_back(kBlue);
            lineColor.push_back(1);
        }

        std::vector<TString> SnglTFileNames;
        SnglTFileNames.push_back("tree_ST_s-channel_");
        SnglTFileNames.push_back("tree_ST_t-channel_antitop_");
        SnglTFileNames.push_back("tree_ST_t-channel_top_");
        SnglTFileNames.push_back("tree_ST_tW_antitop_");
        SnglTFileNames.push_back("tree_ST_tW_top_");
        SnglT = new TChain("tree");
        for( int i = 0 ; i < SnglTFileNames.size() ; i++ ) {
            SnglT->Add(skimTypeMC+"/"+SnglTFileNames[i]+Era+".root");
        }
        if( r == kSignal || r == kSLm || r == kSLe ){
            //ntuples.push_back(new RA2bTree(SnglT));
            //sampleName.push_back("SnglT");
            fillColor.push_back(kOrange);
            lineColor.push_back(1);
        }

        std::vector<TString> TTincFileNames;
        TTincFileNames.push_back("tree_TTJets_");
        TTinc = new TChain("tree");
        for( int i = 0 ; i < TTincFileNames.size() ; i++ ){
            TTinc->Add(skimTypeMC+"/"+TTincFileNames[i]+Era+".root");
        }

        std::vector<TString> TTFileNames;
        TTFileNames.push_back("tree_TTJets_HT-600to800_");
        TTFileNames.push_back("tree_TTJets_HT-800to1200_");
        TTFileNames.push_back("tree_TTJets_HT-1200to2500_");
        TTFileNames.push_back("tree_TTJets_HT-2500toInf_");
        TTFileNames.push_back("tree_TTJets_SingleLeptFromT_");
        TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar_");
        TTFileNames.push_back("tree_TTJets_DiLept_");        
        TT = new TChain("tree");
        for( int i = 0 ; i < TTFileNames.size() ; i++ ){
            TT->Add(skimTypeMC+"/"+TTFileNames[i]+Era+".root");
        }
        if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
            //ntuples.push_back(new RA2bTree(TT));
            //sampleName.push_back("TT");
            fillColor.push_back(kCyan);
            lineColor.push_back(kCyan);
        }

        std::vector<TString> DYFileNames;
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-100to200_");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-200to400_");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-400to600_");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-600to800_");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-800to1200_");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-1200to2500_");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-2500toInf_");
        DY = new TChain("tree");
        for( int i = 0 ; i < DYFileNames.size() ; i++ ){
            DY->Add(skimTypeMC+"/"+DYFileNames[i]+Era+".root");
           // DY->Add(skimTypeLDP+"/"+DYFileNames[i]);
	    //std::cout<<DYFileNames[i]<<std::endl;
        }
 	if(r==kDYe || r==kDYm){
	
        ntuples.push_back(new RA2bTree(DY));
        sampleName.push_back("DY");
        fillColor.push_back(kGreen);
        lineColor.push_back(1);
	}
        std::vector<TString> GJets0p4FileNames;
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-100to200_");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-200to400_");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-400to600_");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-600toInf_");
        GJets0p4 = new TChain("tree");
        for( int i = 0 ; i < GJets0p4FileNames.size() ; i++ ){
            GJets0p4->Add(skimTypeMC+"/"+GJets0p4FileNames[i]+Era+".root");
        }
        if( r == kPhoton ){
            ntuples.push_back(new RA2bTree(GJets0p4));
            sampleName.push_back("GJets");
            fillColor.push_back(kGreen);
            lineColor.push_back(1);
        }

        std::vector<TString> QCDFileNames;
        QCDFileNames.push_back("tree_QCD_HT-200to300_");
        QCDFileNames.push_back("tree_QCD_HT-300to500_");
        QCDFileNames.push_back("tree_QCD_HT-500to700_");
        QCDFileNames.push_back("tree_QCD_HT-700to1000_");
        QCDFileNames.push_back("tree_QCD_HT-1000to1500_");
        QCDFileNames.push_back("tree_QCD_HT-1500to2000_");
        QCDFileNames.push_back("tree_QCD_HT-2000toInf_");
        QCD = new TChain("tree");
        for( int i = 0 ; i < QCDFileNames.size() ; i++ ){
            QCD->Add(skimTypeMC+"/"+QCDFileNames[i]+Era+".root");
        }
        if( r == kSignal || r == kPhoton || r == kLowDphi ){
           // ntuples.push_back(new RA2bTree(QCD));
           // sampleName.push_back("QCD");
            fillColor.push_back(kGray);
            lineColor.push_back(1);
        }

        ////////////////////////////////////////////////////////////
        // - - - - - - - - - - - DATA INPUTS - - - - - - - - - -  //
        ////////////////////////////////////////////////////////////

        std::vector<TString> METFileNames;
        METFileNames.push_back("tree_MET_2016B.root");
        METFileNames.push_back("tree_MET_2016C.root");
        METFileNames.push_back("tree_MET_2016D.root");
        METFileNames.push_back("tree_MET_2016E.root");
        METFileNames.push_back("tree_MET_2016F.root");
        METFileNames.push_back("tree_MET_2016G.root");
        METFileNames.push_back("tree_MET_2016H.root");
       // METFileNames.push_back("tree_HTMHT_re2016H3.root");
        if( r == kSignal || r == kLowDphi ){
            data = new TChain("tree");
            for( int i = 0 ; i < METFileNames.size() ; i++ ){
                data->Add(skimType+"/"+METFileNames[i]);
            }    
            dataNtuple = new RA2bTree(data);
	    //ntuples.push_back(dataNtuple);
	    //sampleName.push_back("data"); 
	    fillColor.push_back(kWhite);
	    lineColor.push_back(1);
        }
        
	METFileNames.resize(0);
        METFileNames.push_back("tree_MET_2017B.root");
        METFileNames.push_back("tree_MET_2017C.root");
        METFileNames.push_back("tree_MET_2017D.root");
        METFileNames.push_back("tree_MET_2017E.root");
        METFileNames.push_back("tree_MET_2017F.root");
       // METFileNames.push_back("tree_HTMHT_re2016H3.root");
        if( r == kSignal || r == kLowDphi ){
            data2017 = new TChain("tree");
            for( int i = 0 ; i < METFileNames.size() ; i++ ){
                data2017->Add(skimType+"/"+METFileNames[i]);
            }    
            dataNtuple = new RA2bTree(data2017);
	    //ntuples.push_back(dataNtuple);
	    //sampleName.push_back("data2017"); 
	    fillColor.push_back(kWhite);
	    lineColor.push_back(1);
        }

	METFileNames.resize(0);
        METFileNames.push_back("tree_MET_2018A.root");
        METFileNames.push_back("tree_MET_2018B.root");
        METFileNames.push_back("tree_MET_2018C.root");
        METFileNames.push_back("tree_MET_2018D.root");
       // METFileNames.push_back("tree_HTMHT_re2016H3.root");
        if( r == kSignal || r == kLowDphi ){
            data2018 = new TChain("tree");
            for( int i = 0 ; i < METFileNames.size() ; i++ ){
                data2018->Add(skimType+"/"+METFileNames[i]);
            }    
            dataNtuple = new RA2bTree(data2018);
	    //ntuples.push_back(dataNtuple);	    
	    //sampleName.push_back("data2018"); 
	    fillColor.push_back(kWhite);
	    lineColor.push_back(1);
        }
        std::vector<TString> SingleElectronNames;
        SingleElectronNames.push_back("tree_SingleElectron_2016B.root");
        SingleElectronNames.push_back("tree_SingleElectron_2016C.root");
        SingleElectronNames.push_back("tree_SingleElectron_2016D.root");
        SingleElectronNames.push_back("tree_SingleElectron_2016E.root");
        SingleElectronNames.push_back("tree_SingleElectron_2016F.root");
        SingleElectronNames.push_back("tree_SingleElectron_2016G.root");
        SingleElectronNames.push_back("tree_SingleElectron_2016H.root");
        if( r == kSLe || r==kDYe ){
            data = new TChain("tree");
            for( int i = 0 ; i < SingleElectronNames.size() ; i++ ){
                data->Add(skimType+"/"+SingleElectronNames[i]);
	
            }
            dataNtuple = new RA2bTree(data);
	    ntuples.push_back(dataNtuple);
	    sampleName.push_back("data"); 
	    fillColor.push_back(kBlack);
	    lineColor.push_back(1);

        }

        std::vector<TString> SingleMuonNames;
        SingleMuonNames.push_back("tree_SingleMuon_2016B.root");
        SingleMuonNames.push_back("tree_SingleMuon_2016C.root");
        SingleMuonNames.push_back("tree_SingleMuon_2016D.root");
        SingleMuonNames.push_back("tree_SingleMuon_2016E.root");
        SingleMuonNames.push_back("tree_SingleMuon_2016F.root");
        SingleMuonNames.push_back("tree_SingleMuon_2016G.root");
        SingleMuonNames.push_back("tree_SingleMuon_2016H.root");
        if( r == kSLm || r==kDYm){
            data = new TChain("tree");
            for( int i = 0 ; i < SingleMuonNames.size() ; i++ ){
                data->Add(skimType+"/"+SingleMuonNames[i]);
            }
            dataNtuple = new RA2bTree(data);
	    ntuples.push_back(dataNtuple);
	    sampleName.push_back("data"); 
	    fillColor.push_back(kBlack);
	    lineColor.push_back(1);
        }
        SingleElectronNames.resize(0);
        SingleElectronNames.push_back("tree_SingleElectron_2017B.root");
        SingleElectronNames.push_back("tree_SingleElectron_2017C.root");
        SingleElectronNames.push_back("tree_SingleElectron_2017D.root");
        SingleElectronNames.push_back("tree_SingleElectron_2017E.root");
        SingleElectronNames.push_back("tree_SingleElectron_2017F.root");
        if( r == kSLe || r==kDYe ){
            data2017 = new TChain("tree");
            for( int i = 0 ; i < SingleElectronNames.size() ; i++ ){
                data2017->Add(skimType+"/"+SingleElectronNames[i]);
	
            }
            dataNtuple2017 = new RA2bTree(data2017);
	    ntuples.push_back(dataNtuple2017);
	    sampleName.push_back("data2017"); 
	    fillColor.push_back(kBlack);
	    lineColor.push_back(1);

        }
        SingleElectronNames.resize(0);
        SingleElectronNames.push_back("tree_EGamma_2018A.root");
        SingleElectronNames.push_back("tree_EGamma_2018B.root");
        SingleElectronNames.push_back("tree_EGamma_2018C.root");
        SingleElectronNames.push_back("tree_EGamma_2018D.root");
        if( r == kSLe || r==kDYe ){
            data2018 = new TChain("tree");
            for( int i = 0 ; i < SingleElectronNames.size() ; i++ ){
                data2018->Add(skimType+"/"+SingleElectronNames[i]);

            }
            dataNtuple2018 = new RA2bTree(data2018);
            ntuples.push_back(dataNtuple2018);
            sampleName.push_back("data2018");
            fillColor.push_back(kBlack);
            lineColor.push_back(1);

        }
        SingleMuonNames.resize(0);
        SingleMuonNames.push_back("tree_SingleMuon_2017B.root");
        SingleMuonNames.push_back("tree_SingleMuon_2017C.root");
        SingleMuonNames.push_back("tree_SingleMuon_2017D.root");
        SingleMuonNames.push_back("tree_SingleMuon_2017E.root");
        SingleMuonNames.push_back("tree_SingleMuon_2017F.root");
        if( r == kSLm || r==kDYm){
            data2017 = new TChain("tree");
            for( int i = 0 ; i < SingleMuonNames.size() ; i++ ){
                data2017->Add(skimType+"/"+SingleMuonNames[i]);
            }
            dataNtuple2017 = new RA2bTree(data2017);
	    ntuples.push_back(dataNtuple2017);
	    sampleName.push_back("data2017"); 
	    fillColor.push_back(kBlack);
	    lineColor.push_back(1);
        }
        SingleMuonNames.resize(0);
        SingleMuonNames.push_back("tree_SingleMuon_2018A.root");
        SingleMuonNames.push_back("tree_SingleMuon_2018B.root");
        SingleMuonNames.push_back("tree_SingleMuon_2018C.root");
        SingleMuonNames.push_back("tree_SingleMuon_2018D.root");
        if( r == kSLm || r==kDYm){
            data2018 = new TChain("tree");
            for( int i = 0 ; i < SingleMuonNames.size() ; i++ ){
                data2018->Add(skimType+"/"+SingleMuonNames[i]);
            }
            dataNtuple2018 = new RA2bTree(data2018);
	    ntuples.push_back(dataNtuple2018);
	    sampleName.push_back("data2018"); 
	    fillColor.push_back(kBlack);
	    lineColor.push_back(1);
        }
        std::vector<TString> SinglePhotonFileNames;
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2016B.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2016C.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2016D.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2016E.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2016F.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2016G.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2016H.root");
        if( r == kPhoton ){
            data = new TChain("tree");
            for( int i = 0 ; i < SinglePhotonFileNames.size() ; i++ ){
                data->Add(skimType+"/"+SinglePhotonFileNames[i]);
            }
            dataNtuple = new RA2bTree(data);
	    ntuples.push_back(dataNtuple);
	    sampleName.push_back("data"); 
	    fillColor.push_back(kBlack);
	    lineColor.push_back(1);
        }
        SinglePhotonFileNames.resize(0);
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2017B.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2017C.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2017D.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2017E.root");
        SinglePhotonFileNames.push_back("tree_SinglePhoton_2017F.root");
        if( r == kPhoton ){
            data2017 = new TChain("tree");
            for( int i = 0 ; i < SinglePhotonFileNames.size() ; i++ ){
                data2017->Add(skimType+"/"+SinglePhotonFileNames[i]);
            }
            dataNtuple = new RA2bTree(data2017);
	    ntuples.push_back(dataNtuple);
	     sampleName.push_back("data2017"); 
	    fillColor.push_back(kBlack);
	    lineColor.push_back(1);
        }
        SinglePhotonFileNames.resize(0);
        SinglePhotonFileNames.push_back("tree_EGamma_2018A.root");
        SinglePhotonFileNames.push_back("tree_EGamma_2018B.root");
        SinglePhotonFileNames.push_back("tree_EGamma_2018C.root");
        SinglePhotonFileNames.push_back("tree_EGamma_2018D.root");
        if( r == kPhoton ){
            data2018 = new TChain("tree");
            for( int i = 0 ; i < SinglePhotonFileNames.size() ; i++ ){
                data2018->Add(skimType+"/"+SinglePhotonFileNames[i]);
            }
            dataNtuple = new RA2bTree(data2018);
            ntuples.push_back(dataNtuple);
            sampleName.push_back("data2018");
            fillColor.push_back(kBlack);
            lineColor.push_back(1);
        }
        T5HH1000 = new TChain("tree");
        T5HH1100 = new TChain("tree");
        T5HH1200 = new TChain("tree");
        T5HH1300 = new TChain("tree");
        T5HH1400 = new TChain("tree");
        T5HH1500 = new TChain("tree");
        T5HH1600 = new TChain("tree");
        T5HH1700 = new TChain("tree");
        T5HH1800 = new TChain("tree");
        T5HH1900 = new TChain("tree");
        T5HH2000 = new TChain("tree");
        T5HH2100 = new TChain("tree");
        T5HH2200 = new TChain("tree");
        T5HH2300 = new TChain("tree");
        T5HH2400 = new TChain("tree");
        T5HH2500 = new TChain("tree");
	TFile*fsig=TFile::Open("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1000_1_"+Era+".root","READ");
        TH1F*temp=(TH1F*)fsig->Get("nEventProc");
	fsig->Close();

	T5HH1000->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1000_1_"+Era+".root");
	T5HH1100->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1100_1_"+Era+".root");
	T5HH1200->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1200_1_"+Era+".root");
	T5HH1300->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1300_1_"+Era+".root");
	T5HH1400->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1400_1_"+Era+".root");
	T5HH1500->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1500_1_"+Era+".root");
	T5HH1600->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1600_1_"+Era+".root");
	T5HH1700->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1700_1_"+Era+".root");
	T5HH1800->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1800_1_"+Era+".root");
	T5HH1900->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_1900_1_"+Era+".root");
	T5HH2000->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_2000_1_"+Era+".root");	
	T5HH2100->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_2100_1_"+Era+".root");
	T5HH2200->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_2200_1_"+Era+".root");
	T5HH2300->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_2300_1_"+Era+".root");
	T5HH2400->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_2400_1_"+Era+".root");
	T5HH2500->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal/tree_T5qqqqZH_2500_1_"+Era+".root");
      if( r == kSignal ){

            signalNtuples.push_back(new RA2bTree(T5HH1000));
            signalNtuples.push_back(new RA2bTree(T5HH1100));
            signalNtuples.push_back(new RA2bTree(T5HH1200));
            signalNtuples.push_back(new RA2bTree(T5HH1300));
            signalNtuples.push_back(new RA2bTree(T5HH1400));
            signalNtuples.push_back(new RA2bTree(T5HH1500));
            signalNtuples.push_back(new RA2bTree(T5HH1600));
            signalNtuples.push_back(new RA2bTree(T5HH1700));
            signalNtuples.push_back(new RA2bTree(T5HH1800));
            signalNtuples.push_back(new RA2bTree(T5HH1900));
            signalNtuples.push_back(new RA2bTree(T5HH2000));
            signalNtuples.push_back(new RA2bTree(T5HH2100));
            signalNtuples.push_back(new RA2bTree(T5HH2200));
            signalNtuples.push_back(new RA2bTree(T5HH2300));
            signalNtuples.push_back(new RA2bTree(T5HH2400));
            signalNtuples.push_back(new RA2bTree(T5HH2500));
    for( int i = 0 ; i < signalNtuples.size() ; i++ ){
    //       for( int i = 0 ; i < 1 ; i++ ){
		int mass=1000+100*i;
            	signalSampleName.push_back(TString::Format("T5ZZ%d" ,mass));
            	sigLineColor.push_back(kRed);
       	        fsig=TFile::Open(TString::Format("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoZEvents/tree_signal_JECup/tree_T5qqqqZH_%d_1_",mass)+Era+".root","READ");	
            	temp=(TH1F*)fsig->Get("nEventProc");
		NSignalEvents.push_back(temp->GetBinContent(1));
		fsig->Close();
	    }

        }
    };

    RA2bTree* findNtuple(TString name){
        for( int iSam = 0 ; iSam < sampleName.size() ; iSam++ ){
            if( sampleName[iSam] == name )
                return ntuples[iSam] ;
        }
        for( int iSam = 0 ; iSam < signalSampleName.size() ; iSam++ ){
            if( signalSampleName[iSam] == name )
                return signalNtuples[iSam] ;
        }
        return NULL;
    };

};
