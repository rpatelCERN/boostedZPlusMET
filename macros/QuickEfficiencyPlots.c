//
//  QuickEfficiencyPlots.c
//  
//
//  Created by Rishi Patel on 4/7/20.
//

#include <stdio.h>
#include"tdrstyle.C"
#include "CMS_lumi.C"
void QuickEfficiencyPlots(){
    TFile*fin=new TFile("SignalEfficiencyChecksT5ZZ1700.root","READ");
   TH1D*ZMatchEffNum=(TH1D*) fin->Get("ZMatchEffNum");
    TH1D*ZMatchEffNumMass=(TH1D*) fin->Get("ZMatchEffMassWindowNum");
    TH1D*ZDenom=(TH1D*) fin->Get("ZDenom");
    TGraphAsymmErrors *eff=new TGraphAsymmErrors (ZMatchEffNum,ZDenom);
    TGraphAsymmErrors *eff2=new TGraphAsymmErrors (ZMatchEffNumMass,ZDenom);
    eff->SetTitle("; Z-boson p_{T} [GeV]; Efficiency (%)");
    TCanvas*c1=new TCanvas("c1","", 800,800);
    writeExtraText = true;       // if extra text
    extraText  = " Simulation";
    lumi_sqrtS = "137 fb^{-1}(13 TeV)";


    
    for(unsigned int i=0; i<eff->GetN(); ++i){
        double x,y;
        eff->GetPoint(i, x,y);
        eff->SetPoint(i, x,y*100);
        eff->SetPointEYhigh(i,eff->GetErrorYhigh(i)*100);
        eff->SetPointEYlow(i,eff->GetErrorYlow(i)*100);

        eff2->GetPoint(i, x,y);
        eff2->SetPoint(i, x,y*100);
        eff2->SetPointEYhigh(i,eff2->GetErrorYhigh(i)*100);
        eff2->SetPointEYlow(i,eff2->GetErrorYlow(i)*100);
    }
    eff->SetMarkerStyle(kOpenCircle);
    eff->SetMarkerColor(kBlack);
    eff->SetLineColor(kBlack);
    eff->GetXaxis()->SetRangeUser(200, 1500);
    eff->GetYaxis()->SetRangeUser(0, 105);
    eff->Draw("APE");
    eff2->SetMarkerStyle(kFullCircle);
     eff2->SetMarkerColor(kBlack);
    eff2->SetLineColor(kBlack);

    eff2->Draw("PESame");
    CMS_lumi( c1,0,0 );
 TLegend *leg=new TLegend(0.2669173,0.1576873,0.679198,0.3942181,NULL,"brNDC NB");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);
    leg->SetHeader("T5qqqqZZ M_{#tilde{g}} = 1.7 TeV");
    leg->AddEntry(eff, "AK8 Jets Match Gen Z", "PLE");
    leg->AddEntry(eff2, "AK8 Jets Match Gen Z and m_{J} [70,100]", "PLE");

    leg->Draw();
    c1->Update();
    return;
    /*
    TCanvas*c2=new TCanvas("c2","", 800,800);
    TFile*fin1300=new TFile("SignalEfficiencyChecksT5ZZ1300.root","READ");
    TH1D*LeadJetMass1300=(TH1D*)fin1300->Get("LeadJetSoftMass");
    TFile*fin1700=new TFile("SignalEfficiencyChecksT5ZZ1700.root","READ");
    TH1D*LeadJetMass1700=(TH1D*)fin1700->Get("LeadJetSoftMass");

    TFile*fin2100=new TFile("SignalEfficiencyChecksT5ZZ2100.root","READ");
    TH1D*LeadJetMass2100=(TH1D*)fin2100->Get("LeadJetSoftMass");
    LeadJetMass2100->Rebin(2);
    LeadJetMass1700->Rebin(2);
    LeadJetMass1300->Rebin(2);
    LeadJetMass2100->DrawNormalized();
    LeadJetMass1700->DrawNormalized("same");
    LeadJetMass1300->DrawNormalized("same");
*/
     }
