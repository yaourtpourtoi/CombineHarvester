###Before combine, prefarably in lx01:
#Producing datacards (in working area, not here):
#e.g.:
$./scripts/makeDatacards_cpdecay_2018_tauID_SF.py --cfg=scripts/plot_cpdecays_2018.cfg -c 'mt' scripts/params_2018.json -s 'cpdecay' --batch --embedding --output_folder="et_DCembedding2018/"



###Now coming to combine(should be in lx02):
###Use Protocol_Decays.md  for more commands
###Use src/HttSystematics_SMRun2.cc  to find systematics



#combine help:
#https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/

#Modify the morphing code, i.e. bin/TauIDSF.cpp, according to the needs then:
$scramv1 b clean;
$scramv1 b;

#e.g. of command to run to produce .txt datacards from .root datacards: (be careful of input_folder and output_folder in the bin/TauIDSF.cpp as sometimes they are set ad-hoc)
$TauIDSF --output_folder="tauIDSF_output/embed/2018" --embed=true --era=2018

#make workspace from .txt datacards
$combineTool.py  -M T2W -i output/tauIDSF_output/embed/2018/*/* -o ws.root --parallel 8

#do the fits using workspace:
$for i in {1..18}; do combineTool.py -m 125 -M MultiDimFit -d output/tauIDSF_output/embed/2018/htt_mt_${i}_13TeV/125/ws.root --there -n .r --saveFitResult  --cminDefaultMinimizerTolerance=0.1 --cminDefaultMinimizerStrategy=0; done;
#--cminDefaultMinimizerTolerance=0.1 is the default. --cminDefaultMinimizerStrategy=1 is the default.

#to see fit result:
$root -l output/tauIDSF_output/embed/2018/htt_mt_9_13TeV/125/multidimfit.r.root
root[0] fit_mdf->Print()

#produce postfit shapes:
$for i in {1..18}; do PostFitShapesFromWorkspace -m 125 -d output/et_tauIDSF_output/MC/2018/htt_mt_${i}_13TeV/125/combined.txt.cmb -w output/et_tauIDSF_output/MC/2018/htt_mt_${i}_13TeV/125/ws.root -o output_shapes/TightVsEle/MC/2018/shapes_${i}.root --print --sampling --postfit -f output/et_tauIDSF_output/MC/2018/htt_mt_${i}_13TeV/125/multidimfit.r.root:fit_mdf; done;

