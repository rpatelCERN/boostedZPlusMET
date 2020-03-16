import os
mgo=[1900]

for m in mgo:
	os.system("text2workspace.py --X-allow-no-signal --X-allow-no-background T5ZZ%dBtags.txt -o T5ZZ%dBtags.root" %(m,m))
  	os.system("combineTool.py -M Impacts -d T5ZZ%dBtags.root -m 125 --doInitialFit --expectSignal 1 -t -1 --rMin=-10 " %m)
    	os.system("combineTool.py -M Impacts --doFits --expectSignal 1 -t -1 -m 125 -d T5ZZ%dBtags.root --rMin=-10 --robustFit 1 " %m)
   	os.system("rm -rf higgsCombine_paramFit_Test_PoisMETUncStat6.MultiDimFit.mH125.root")
  	os.system("rm -rf rm -rf higgsCombine_paramFit_Test_MassShape*.MultiDimFit.mH125.root")
  	os.system("rm -rf rm -rf higgsCombine_paramFit_Test_MET*.MultiDimFit.mH125.root")
  	os.system("combineTool.py -M Impacts -d T5ZZ%dBtags.root -m 125 -o impactsT5ZZ%d.json" %(m,m))
   	os.system("plotImpacts.py -i impactsT5ZZ%d.json -o impactsT5ZZ%d" %(m,m))
