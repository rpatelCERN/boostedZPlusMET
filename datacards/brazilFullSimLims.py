#! /usr/bin/env python

import os
import glob
import math
import array
import sys
import time
import ROOT
from array import array


import tdrstyle
tdrstyle.setTDRStyle()
# ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
# ROOT.setTDRStyle();
# ROOT.gStyle.SetPadLeftMargin(0.16);
# ROOT.gStyle.SetPadRightMargin(0.10);
# ROOT.gStyle.SetPadTopMargin(0.10);
# ROOT.gStyle.SetPalette(1);

## ===========================================================================================
## ===========================================================================================
## ===========================================================================================

def columnToList(fn,col):
	f = open(fn,'r');

	olist = [];
	for line in f: 
		linelist = line.strip().split()
		olist.append( linelist[col] );
	return olist

def ExtractFile(iname, tag):
	f = ROOT.TFile(iname);
	t = f.Get("limit");
	lims = [];
	lims.append(tag);
	for i in range(6):
		t.GetEntry(i);
		lims.append( t.limit )
	return lims;

if __name__ == '__main__':

	idir = "./";
	results = [];			  
	#results.append( ExtractFile(idir+'/higgsCombineT5ZZ750.Asymptotic.mH750.root','750') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1000.txt.Asymptotic.mH1000.root','1000') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1100.txt.Asymptotic.mH1100.root','1100') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1200.txt.Asymptotic.mH1200.root','1200') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1300.txt.Asymptotic.mH1300.root','1300') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1400.txt.Asymptotic.mH1400.root','1400') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1500.txt.Asymptotic.mH1500.root','1500') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1600.txt.Asymptotic.mH1600.root','1600') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1700.txt.Asymptotic.mH1700.root','1700') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1800.txt.Asymptotic.mH1800.root','1800') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ1900.txt.Asymptotic.mH1900.root','1900') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ2000.txt.Asymptotic.mH2000.root','2000') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ2100.txt.Asymptotic.mH2100.root','2100') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ2200.txt.Asymptotic.mH2200.root','2200') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ2300.txt.Asymptotic.mH2300.root','2300') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ2400.txt.Asymptotic.mH2400.root','2400') );
	results.append( ExtractFile(idir+'/higgsCombineT5ZZ2500.txt.Asymptotic.mH2500.root','2500') );
	#xsecs=[2.26585, 0.325388 , 0.163491, 0.0856418, 0.0460525, 0.0252977, 0.0141903, 0.00810078, 0.00470323, 0.00276133, 0.00163547, 0.000981077,0.000591918]
	xsecs=[0.385,0.191,0.0985,0.0522,0.0284,0.0157,0.00887,0.00507,0.00293,0.00171,0.00101,0.000598,0.000356,0.000213,0.000128,0.0000768]
	#xsecs=[0.325388 , 0.163491, 0.0856418, 0.0460525, 0.0252977, 0.0141903, 0.00810078, 0.00470323, 0.00276133, 0.00163547, 0.000981077,0.000591918]
	names   = [];
	l_obs   = [];
	l_m2sig = [];
	l_m1sig = [];
	l_exp   = [];
	l_p1sig = [];
	l_p2sig = [];
	count=0;
	for r in results:
		names.append(r[0]);
		l_m2sig.append(r[1]*xsecs[count]);
		l_m1sig.append(r[2]*xsecs[count]);
		l_exp.append(r[3]*xsecs[count]);
		l_p1sig.append(r[4]*xsecs[count]);
		l_p2sig.append(r[5]*xsecs[count]);
		l_obs.append(r[6]*xsecs[count]);
		
		count=count+1
	print "l_exp = ", l_exp
	print "l_obs = ", l_obs

	a_xax = array('d', []);
	a2_xax = array('d', []);
	a_exp = array('d', []);
	a_obs = array('d', []);
	a_1sig = array('d', []);
	a_2sig = array('d', []);
	#Need to do this a bit more clever
	for i in range(len(names)): a_xax.append( float(names[i]) );
	for i in range(len(names)): a2_xax.append( float(names[i]) );
	for i in range(len(names)-1,-1,-1): a2_xax.append( float(names[i]));
	for i in range(len(l_obs)): a_obs.append( float(l_obs[i]) );
	for i in range(len(l_exp)): a_exp.append( float(l_exp[i]) );
	
	for i in range(len(l_m2sig)): a_2sig.append( float(l_m2sig[i]) );
	for i in range(len(l_p2sig)-1,-1,-1): a_2sig.append( float(l_p2sig[i]) );
	
	for i in range(len(l_m1sig)): a_1sig.append( float(l_m1sig[i]) );
	for i in range(len(l_p1sig)-1,-1,-1): a_1sig.append( float(l_p1sig[i]) );
	#print a_2sig, len(a_2sig)
	#print a2_xax, len(a2_xax)
	a_2sig.append(results[0][6])
	a2_xax.append(0.5)

	g_exp = ROOT.TGraph(len(a_xax), a_xax, a_exp)
	g_obs = ROOT.TGraph(len(a_xax), a_xax, a_obs)
	g_1sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_1sig)
	g_2sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_2sig)

	#print g_2sig;

	can = ROOT.TCanvas("can","can",1200,1000);
	hrl = ROOT.TH1F("hrl","hrl",30,1000,2500);
	# hrl = can.DrawFrame(0,0,6,15);
	hrl.GetYaxis().SetTitle("Cross section [pb] ");
	hrl.GetXaxis().SetTitle("M_{#tilde{g}} [GeV]");
	#for i in range(0,15):
	#	if i%3==0:
	#		hrl.GetXaxis().SetBinLabel(i+1,names[i])
	#	if i==14:hrl.GetXaxis().SetBinLabel(i+1, names[i])
	#hrl.GetXaxis().SetBinLabel(2,names[1])
	#hrl.GetXaxis().SetBinLabel(3,names[2])
	#hrl.GetXaxis().SetBinLabel(4,names[3])
	#hrl.GetXaxis().SetBinLabel(5,names[4])
	#hrl.GetXaxis().SetBinLabel(6,names[5])
	hrl.SetMaximum(0.1);
	hrl.SetMinimum(0.0002);
	hrl.Draw();

	#can.SetGrid(); 
	can.SetLogy();
	hrl.GetXaxis().SetRangeUser(1300,2500);
        hrl.SetMinimum(0.0001);
	txta = ROOT.TLatex(0.17,0.96,"CMS");
	txta.SetNDC();
	txtb = ROOT.TLatex(0.75,0.905,"Preliminary");
	txtb.SetNDC(); txtb.SetTextFont(52);
	
	txtc = ROOT.TLatex(0.75,0.96,"137 fb^{-1} (13 TeV)");
	txtc.SetNDC(); txtc.SetTextFont(42); txtc.SetTextSize(0.04);
        txtd = ROOT.TLatex(0.45,0.86,"pp#rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow q #bar{q} #tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2} #rightarrow Z #tilde{#chi}^{0}_{1} ");
	txte = ROOT.TLatex(0.45,0.8,"m_{#tilde{g}} #font[122]{-} m_{#tilde{#chi}^{0}_{2}} = 50 GeV, m_{#tilde{#chi}^{0}_{1}} = 1 GeV  ");
        txte.SetNDC(); txte.SetTextFont(42); txte.SetTextSize(0.05);	
	#txtd = ROOT.TLatex(0.65,0.85,"#tilde{#chi}^{2}_{0}#rightarrow Z #tilde{#chi}^{1}_{0}, #tilde{#chi}^{2}_{0}#rightarrow Z #tilde{#chi}^{1}_{0}");
	txtd.SetNDC(); txtd.SetTextFont(42); txtd.SetTextSize(0.05);
	f=open("LatestGluGluNNLO.txt", 'r')
	a_stop = array('d', []);
	a_xsec = array('d', []);
	a_xsecUP = array('d', []);
	a_xsecDown = array('d', []);
	for line in f:
		parse=line.split("|")
		a_stop.append(float(parse[1]));
		a_xsec.append(float(parse[2]))
                a_xsecUP.append(float(parse[2])+(float(parse[2])*float(parse[3])/100.))
                a_xsecDown.append(float(parse[2])-(float(parse[2])*float(parse[3])/100.))
	g_xsec=ROOT.TGraph(len(a_stop), a_stop, a_xsec)
        g_xsecSup=ROOT.TGraph(len(a_stop), a_stop, a_xsecUP)
        g_xsecSdn=ROOT.TGraph(len(a_stop), a_stop, a_xsecDown)
	#leg = ROOT.TLegend(0.20,0.15,0.35,0.30);
	leg = ROOT.TLegend(0.6477462,0.5292308,0.9474124,0.7333333);
	leg.SetFillStyle(1001);
	leg.SetFillColor(0);    
	leg.SetBorderSize(1);  
	# leg.SetNColumns(2);
	leg.AddEntry(g_exp,"expected","l")
	leg.AddEntry(g_obs,"observed","l")
	leg.AddEntry(g_2sig,"expected 2#sigma_{exp}","f")
	leg.AddEntry(g_1sig,"expected 1#sigma_{exp}","f")
  	leg.AddEntry(g_xsec, "Theory #pm #sigma_{theory}", "l")

        leg2 = ROOT.TLegend(0.6477462,0.4564103,0.9482471,0.6605128); 
        leg2.SetBorderSize(0);
        leg2.SetLineColor(1);
        leg2.SetLineStyle(1);
        leg2.SetLineWidth(1);
        leg2.SetFillColor(0);
        leg2.SetFillStyle(0);
        entry=leg2.AddEntry("NULL"," ","l");
        entry.SetLineColor(4);
        entry.SetLineStyle(2);
        entry.SetLineWidth(1);
        entry.SetMarkerColor(1);
   	entry.SetMarkerStyle(21);
   	entry.SetMarkerSize(1);
   	entry.SetTextFont(42);
        leg3 = ROOT.TLegend(0.648581,0.4489744,0.9490818,0.6330769);
        leg3.SetBorderSize(0);
        leg3.SetLineColor(1);
        leg3.SetLineStyle(1);
        leg3.SetLineWidth(1);
        leg3.SetFillColor(0);
        leg3.SetFillStyle(0);
        entry=leg3.AddEntry("NULL"," ","l");
        entry.SetLineColor(4);
        entry.SetLineStyle(2);
        entry.SetLineWidth(1);
        entry.SetMarkerColor(1);
   	entry.SetMarkerStyle(21);
   	entry.SetMarkerSize(1);
   	entry.SetTextFont(42);

	#oneLine = ROOT.TF1("oneLine","1",175,550);
	#oneLine.SetLineColor(ROOT.kRed+2);
	#oneLine.SetLineWidth(2);
	#oneLine.SetLineStyle(1);
	
	g_1sig.SetFillColor(ROOT.kGreen);
	g_1sig.SetFillStyle(1000);
	g_2sig.SetFillColor(ROOT.kYellow);
	g_2sig.SetFillStyle(1000);
	g_exp.SetLineStyle(2);
	g_exp.SetLineWidth(2);
	g_obs.SetLineWidth(2);
	g_2sig.Draw('f');
	g_1sig.Draw('fsames');
	#g_1sig.Draw('f');
	g_obs.Draw('lsames');
	g_exp.Draw('lsames');
	for i in range(0,100):
		#print "Mass %d  Exp Excl %g " %(2000+i,g_exp.Eval(2000+i))
		#print "Theory Xsec %g " %g_xsec.Eval(2000+i)
		if(g_obs.Eval(1900+i)<g_xsec.Eval(1900+i)):print "Mass %d" %(1900+i)
		#if(g_exp.Eval(2000+i)<g_xsec.Eval(2000+i)):print "Mass %d" %(2000+i)
	#oneLine.Draw("LSAMES");
	txta.Draw();
	#txtb.Draw();
	txtc.Draw();
	txtd.Draw();	
	txte.Draw();	
	leg.Draw();
   	leg2.Draw();
   	leg3.Draw();
	g_xsecSdn.SetLineStyle(2);
	g_xsecSup.SetLineStyle(2);
	g_xsec.SetLineStyle(1);
	g_xsecSup.SetLineColor(ROOT.kBlue);
	g_xsecSdn.SetLineColor(ROOT.kBlue);
	g_xsec.SetLineColor(ROOT.kBlue);
	g_xsec.Draw("lsame")
	g_xsecSup.Draw("lsame")
	g_xsecSdn.Draw("lsame")
	ROOT.gPad.Update();
   	ROOT.gPad.RedrawAxis();
	#hrl.Draw("same")
	can.SaveAs('T5ZZCombResults.pdf');
	can.SaveAs('T5ZZCombResults.C');
	fout=ROOT.TFile("BoostedZPlots.root","RECREATE")
        g_xsec.Write("GluinoTheoryXsec");
        g_xsecSup.Write("GluinoTheoryXsecSigmaUp");
        g_xsecSdn.Write("GluinoTheoryXsecSigmaDn");
        g_exp.Write("ExpectedLimit");
        g_obs.Write("ObservedLimit");
        g_2sig.Write("ExpectedLimitTwoSigmaBand"); 
        g_1sig.Write("ExpectedLimitSigmaBand"); 
        hrl.Write("CombResult");

