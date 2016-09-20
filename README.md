Repo for studying boosted Higgs tools with the RA2b trees.  This package relies on another package for setting up some basic 
classes to facilitate some things.  

#General Stuff:

## Setting up the code:
git clone https://github.com/awhitbeck/AnalysisTools.git
git clone https://github.com/awhitbeck/boostedHiggsPlusMET.git

## Compiling:

There is currently no makefile for this repository.  However, to compile any of the code, one can use <pre>build.sh</pre>
by passing the name of the target .cc minus the file extension as an arugment.  The executable will will have the same
root as the target file name, but will be append with '.exe'
 
## Batch:

Before submitting any batch jobs, be sure to retar your working area and copy it into the <pre>src</pre> directory.

<pre>
cd $CMSSW_BASE/../
tar -cf workingArea.tar CMSSW_7_4_2 --exclude='*.dag.*' --exclude='*tar' --exclude='*root' --exclude='*png' --exclude='*pdf' --exclude='*stdout' --exclude='*stderr' --exclude='*condor'
mv workingArea.tar $CMSSW_BASE/src/boostedHiggsHeppy/src/
cd $CMSSW_BASE/src/boostedHiggsHeppy/src/
</pre>

# Baseline Skims

baselineSkim.cc is used to reduce the heppy trees down to a managable size by
selecting one the events that pass a set of baseline cuts.  These cuts are
defined in <pre>src/selecBaseline.cc</pre>. This script take two arguments, 
the input file name and the output eos directory.  The input file must exist 
in the input list found in data/RA2bInputs.txt.  

Example:
<pre>
./baselineSkim.exe /store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV10/Spring16.TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root ./
</pre>

### Batch submit:

To submit the code to condor @ the LPC, the simplest method is to first create a list of DAG files.  A script exists
to generate these dag files automatically.  The input arguments to the bash script (baslineSkim.sh) which is executed 
on the worker nodes are passed to the JDL files via the DAG files.  An example job configuration is:
<pre>
JOB baselineSkim_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0 baselineSkim.jdl
VARS baselineSkim_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0 arguments="/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV10/Spring16.QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_9_RA2AnalysisTree.root root://cmseos.fnal.gov//store/user/awhitbe1/RA2bSkims_V10_v0/ "
</pre>
The jdl file for run the skims is <code>baselineSkim.jdl</code>.  One should specify the three command line arguments to <code>baselineSkim.sh</code> 
when submitting the job.  e.g.:

<pre>
condor_submit baselineSkim.jdl arguments='root://cmseos.fnal.gov//store/user/lpchbb/HeppyNtuples/V14/ root://cmseos.fnal.gov//store/user/awhitbe1/heppySkims/ TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root'
</pre>