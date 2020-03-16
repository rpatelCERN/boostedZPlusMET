import os
from array import array
from ROOT import *
f=open("realistic-counting-experimentoptimization.txt", "r")
#Signal yields in MET bins for each of the mass points
#Bkg Yields in the same MET bins for the mass points
#finputs=TFile("BkgEstimates.root", "READ");
#fsig=TFile("BoostedZZSignalCorrectedNominalwHEM.root", "READ");
fsig=TFile("SignalNominalResync.root", "READ");
mGo=[1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200, 2300,2400,2500]
ContaminScales=[0.0, 0.0, 0.0, 0.092,0.090,0.085,0.083,0.082,0.083, 0.081,0.076,0.079,0.079,0.078,0.079, 0.078]
#ContaminScales=[208.629,111.446,60.249,32.349,18.318,10.254,5.825,3.432,2.069,1.211,0.705,0.432,0.263,0.159,0.098,0.059]
#ContaminScales=[0.552,0.295,0.159,0.086,0.048,0.027,0.015,0.009,0.005,0.003,0.002,0.001,0.001,0,0,0]
#mGo=[1900,2000,2100]
#for i in range(0,16):
#	signal=fsig.Get("T5ZZ%dMETShapeSRMETBins" %mGo[i])
	#print ContaminScales[i]/signal.Integral()
#	ContaminScales[i]=ContaminScales[i]/signal.Integral()
#	print ContaminScales[i]
fpre=TFile("Prefire.root", "READ");
fisr=TFile("ISRSys.root", "READ");
fstat=TFile("MCStat.root","READ");
fJEC=TFile("JECSys.root", "READ");
fJER=TFile("JERSys.root", "READ");
fScale=TFile("ScaleSys.root", "READ");
fSmear=TFile("SmearSys.root", "READ");
count=0
for m in mGo:		
        signal=fsig.Get("T5ZZ%dMETShapeSRMETBins" %m)
	signal.Scale(1.0-ContaminScales[count])
	Pre=fpre.Get("T5ZZ%dPrefire" %m)
        ISRUp=fisr.Get("T5ZZ%dISRDown" %m)
        ISRDown=fisr.Get("T5ZZ%dISRUp" %m)
        MCStat=fstat.Get("T5ZZ%dMCStat" %m)
	JECup=fJEC.Get("T5ZZ%dJECUp" %m)
	JECdown=fJEC.Get("T5ZZ%dJECDown" %m)
	JERup=fJER.Get("T5ZZ%dJERUp" %m)
	JERdown=fJER.Get("T5ZZ%dJERDown" %m)
	Smear=fSmear.Get("T5ZZ%dSmear" %m);
	ScaleUp=fScale.Get("T5ZZ%dScaleUp" %m)
	ScaleDown=fScale.Get("T5ZZ%dScaleDown" %m)
        fout=open("T5ZZ%dBtags.txt" %m,"w")
	f.seek(0);
	count=count+1
	for line in f:
                if "rate" in line:line=line.replace("sig1", "%g" %signal.GetBinContent(1));
                if "rate" in line:line=line.replace("sig2", "%g" %signal.GetBinContent(2));
                if "rate" in line:line=line.replace("sig3", "%g" %signal.GetBinContent(3));
                if "rate" in line:line=line.replace("sig4", "%g" %signal.GetBinContent(4));
                if "rate" in line:line=line.replace("sig5", "%g" %signal.GetBinContent(5));
                if "rate" in line:line=line.replace("sig6", "%g" %signal.GetBinContent(6));
		if "ISRSys" in line:line="ISRSys lnN %g/%g - %g/%g - %g/%g - %g/%g - %g/%g - %g/%g -\n" %(ISRUp.GetBinContent(1),ISRDown.GetBinContent(1), ISRUp.GetBinContent(2),ISRDown.GetBinContent(2),ISRUp.GetBinContent(3),ISRDown.GetBinContent(3),ISRUp.GetBinContent(4),ISRDown.GetBinContent(4),ISRUp.GetBinContent(5),ISRDown.GetBinContent(5),ISRUp.GetBinContent(6),ISRDown.GetBinContent(6))
		#if "ISRSys" in line:line="ISRSys lnN %g - %g - %g - %g - %g - %g -\n" %(ISRUp.GetBinContent(1), ISRUp.GetBinContent(2),ISRUp.GetBinContent(3),ISRUp.GetBinContent(4),ISRUp.GetBinContent(5),ISRUp.GetBinContent(6))
		if "SignalStat" in line:
			line="SignalStat1 lnN %g - - - - - - - - - - -\n" %(MCStat.GetBinContent(1))#, MCStat.GetBinContent(2),MCStat.GetBinContent(3),MCStat.GetBinContent(4),MCStat.GetBinContent(5),MCStat.GetBinContent(6))
                	fout.write(line)
			line="SignalStat2 lnN  - - %g - - - - - - - - -\n" %(MCStat.GetBinContent(2))#, MCStat.GetBinContent(2),MCStat.GetBinContent(3),MCStat.GetBinContent(4),MCStat.GetBinContent(5),MCStat.GetBinContent(6))
                	fout.write(line)
			line="SignalStat3 lnN - - - - %g - - - - - - -\n" %(MCStat.GetBinContent(3))#, MCStat.GetBinContent(2),MCStat.GetBinContent(3),MCStat.GetBinContent(4),MCStat.GetBinContent(5),MCStat.GetBinContent(6))
                	fout.write(line)
			line="SignalStat4 lnN - - - - - - %g - - - - -\n" %(MCStat.GetBinContent(4))#, MCStat.GetBinContent(2),MCStat.GetBinContent(3),MCStat.GetBinContent(4),MCStat.GetBinContent(5),MCStat.GetBinContent(6))
                	fout.write(line)
			line="SignalStat5 lnN - - - - - - - - %g - - -\n" %(MCStat.GetBinContent(5))#, MCStat.GetBinContent(2),MCStat.GetBinContent(3),MCStat.GetBinContent(4),MCStat.GetBinContent(5),MCStat.GetBinContent(6))
                	fout.write(line)
			line="SignalStat6 lnN - - - - - - - - - - %g -\n" %(MCStat.GetBinContent(6))#, MCStat.GetBinContent(2),MCStat.GetBinContent(3),MCStat.GetBinContent(4),MCStat.GetBinContent(5),MCStat.GetBinContent(6))
                	fout.write(line)
			continue	
		if "SignalJEC" in line:line="SignalJEC lnN %g/%g - %g/%g - %g/%g - %g/%g - %g/%g - %g/%g -\n" %(JECup.GetBinContent(1),JECdown.GetBinContent(1), JECup.GetBinContent(2),JECdown.GetBinContent(2),JECup.GetBinContent(3),JECdown.GetBinContent(3), JECup.GetBinContent(4),JECdown.GetBinContent(4),JECup.GetBinContent(5),JECdown.GetBinContent(5), JECup.GetBinContent(6),JECdown.GetBinContent(6))
		if "SignalJER" in line and m!=2000:line="SignalJER lnN %g/%g - %g/%g - %g/%g - %g/%g - %g/%g - %g/%g -\n" %(JERup.GetBinContent(1),JERdown.GetBinContent(1), JERup.GetBinContent(2),JERdown.GetBinContent(2),JERup.GetBinContent(3),JERdown.GetBinContent(3), JERup.GetBinContent(4),JERdown.GetBinContent(4),JERup.GetBinContent(5),JERdown.GetBinContent(5), JERup.GetBinContent(6),JERdown.GetBinContent(6))
		if "SignalJER" in line and m==2000:line="SignalJER lnN %g - %g - %g - %g - %g - %g -\n" %(2-JERup.GetBinContent(1),2- JERup.GetBinContent(2),2-JERup.GetBinContent(3),2-JERup.GetBinContent(4),2-JERup.GetBinContent(5),2-JERup.GetBinContent(6))
		if "SignalRes" in line:line="SignalRes lnN %g - %g - %g - %g - %g - %g -\n" %(Smear.GetBinContent(1),Smear.GetBinContent(2),Smear.GetBinContent(3),Smear.GetBinContent(4),Smear.GetBinContent(5),Smear.GetBinContent(6))
		if "SignalPrefire" in line: line="SignalPrefire lnN %g - %g - %g - %g - %g - %g -\n" %(Pre.GetBinContent(1), Pre.GetBinContent(2),Pre.GetBinContent(3),Pre.GetBinContent(4),Pre.GetBinContent(5),Pre.GetBinContent(6))
		if "SignalScale" in line: line="SignalScale lnN %g/%g - %g/%g - %g/%g - %g/%g - %g/%g - %g/%g -\n" %(ScaleDown.GetBinContent(1),ScaleUp.GetBinContent(1), ScaleDown.GetBinContent(2),ScaleUp.GetBinContent(2),ScaleDown.GetBinContent(3),ScaleUp.GetBinContent(3),ScaleDown.GetBinContent(4),ScaleUp.GetBinContent(4),ScaleDown.GetBinContent(5),ScaleUp.GetBinContent(5),ScaleDown.GetBinContent(6),ScaleUp.GetBinContent(6))
                fout.write(line)
        fout.close()
	#os.system("combineCards.py T5ZZ%d.txt T5ZZ%dSingleJet.txt T5ZZ%dBtags.txt  > T5ZZCombine%d.txt" %(m,m,m,m) ) 
	#os.system("combineCards.py T5ZZ%d.txt T5ZZ%dSingleJet.txt   > T5ZZCombine%d.txt" %(m,m,m) ) 
	print m
	os.system("combine -M Asymptotic T5ZZ%dBtags.txt  -m %d -n T5ZZ%d.txt " %(m,m,m))
