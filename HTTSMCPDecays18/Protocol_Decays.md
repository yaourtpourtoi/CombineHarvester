# Always do the following commend before running anything (it seems to prevent random seg faults )

ulimit -s unlimited

# Creating datacard

    MorphingSMCPDecays18 --output_folder="pas_1206_v2" --mergeXbbb=true 

the option --no_shape_systs=true can be used as well to remove all shape uncertainties except for bbb's

# Building the workspaces:
(combined sub folder only)

    combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/pas_1206_v2/cmb/* -o ws.root --parallel 8

(or to build all subdirectories)

    combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/pas_1206_v2/*/* -o ws.root --parallel 8

# Run maximum likelihood scan

    combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90 --points 21 --redefineSignalPOIs alpha  -d output/test_cp/cmb/125/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1 --cminDefaultMinimizerStrategy=0 

    If want to scale to some lumi X, include the rate parameter lumi_scale=X in the --setParameters option (scaling 2017+2018 to full run2 = 1.35)

    To run on IC batch use (1 point per job):
    `--job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" --split-points 1`
    To run on lx batch use:
    `--job-mode lxbatch --sub-opts '-q 1nh --split-points 1'

    Useful option to save all nuisance parameter values when performing MultiDimFit (and doesn't seem to cost extra runtime):
    --saveSpecifiedNuis all

for this fit and others when running on data it helps to define fall back algorithms with higher tolerance, e.g:

--cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10

# Plot scan

1D scans can be plotted using scripts/plot1DScan.py script.
To plot alpha:

   python scripts/plot1DScan.py --main=output/newnew_mt_bins/cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root --POI=alpha --output=alpha_cmb --no-numbers --no-box --x-min=-90 --x-max=90 --y-max=8

## New scan plotting script

Plot 1D scans using `scripts/draw_nll_scans.py`, see instructions below.

# do 2D scans of kappas
build workspace with 

  combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/merge_sig/cmb/* --PO do_kappas -o ws_kappas.root --parallel 8

run files (note i would use batch jobs):

  combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,muggH=1,kappaH=1,kappaA=0  --redefineSignalPOIs kappaH,kappaA --points 2000  -d output/merge_sig/cmb/125/ws_kappas.root --algo grid -t -1 --there -n .kappas --alignEdges 1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=1


# GOF tests for unblinding

KS tests

(note when unblinding the data in stages change cmb accordingly!)
Also change job-mode options if you are not running on IC batch!

Run toys for all cats seperatly:

   combineTool.py -M GoodnessOfFit --algorithm KS  --there -d output/pas_1206_v2/cmb/125/ws.root -n ".KS.toys" --fixedSignalStrength=1 --there --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 1 -s 0:500:1 

(this runs 500 jobs with 1 toys each)

Run observed
   combineTool.py -M GoodnessOfFit --algorithm KS  --there -d output/pas_1206_v2/cmb/125/ws.root -n ".KS" --fixedSignalStrength=1 --setParameters muV=1,alpha=0,muggH=1,mutautau=1

Collect results and make plots

    combineTool.py -M CollectGoodnessOfFit --input output/pas_1206_v2/cmb/125/higgsCombine.KS.GoodnessOfFit.mH125.root output/pas_1206_v2/cmb/125/higgsCombine.KS.toys.GoodnessOfFit.mH125.*.root --there -o cmb_KS.json

   python ../CombineTools/scripts/plotGof.py --statistic KS --mass 125.0 cmb_KS.json --title-right="137 fb^{-1} (13 TeV)" --output='-KS'

saturated model tests:

For saturated model we always run seperatly for each category

Run toys for all cats seperatly:

   combineTool.py -M GoodnessOfFit --algorithm saturated  --there -d output/pas_1206_v2/htt_*_13TeV/125/ws.root -n ".saturated.toys" --fixedSignalStrength=1 --there --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 50 -s 0:10:1

(this runs 10 jobs per channel/year with 50 toys each)
 
and for combined
   combineTool.py -M GoodnessOfFit --algorithm saturated  --there -d output/pas_1206_v2/cmb/125/ws.root -n ".saturated.toys" --fixedSignalStrength=1 --there --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 1 -s 0:500:1

(this runs 500 jobs with 1 toys each) 

Run observed
   combineTool.py -M GoodnessOfFit --algorithm saturated  --there -d output/pas_1206_v2/htt_*_13TeV/125/ws.root -n ".saturated" --fixedSignalStrength=1 --setParameters muV=1,alpha=0,muggH=1,mutautau=1
and for combined 
  combineTool.py -M GoodnessOfFit --algorithm saturated  --there -d output/pas_1206_v2/cmb/125/ws.root -n ".saturated" --fixedSignalStrength=1 --setParameters muV=1,alpha=0,muggH=1,mutautau=1


Collect output and make plots:

    bins=( "mt_1" "mt_2" "mt_3" "mt_4" "mt_5" "mt_6" "tt_1" "tt_2" "tt_3" "tt_4" "tt_5" "tt_6" "tt_7" "tt_8" "tt_9" "tt_10" "tt_11" )
    years=( "2016" "2017" "2018" )

    for i in "${bins[@]}"
    do
      for j in "${years[@]}"
      do
        combineTool.py -M CollectGoodnessOfFit --input output/pas_1206_v2/htt_"$i"_"$j"_13TeV/125/higgsCombine.saturated.GoodnessOfFit.mH125.root output/gof/htt_"$i"_"$j"_13TeV/125/higgsCombine.saturated.toys.GoodnessOfFit.mH125.*.root --there -o "$i"_"$j"_saturated.json
      python ../CombineTools/scripts/plotGof.py --statistic saturated --mass 125.0 "$i"_"$j"_saturated.json --title-right="137 fb^{-1} (13 TeV)" --output="$i"_"$j"'-saturated'
      done
    done

then do same for cmb if you ran this as well

# perform ZTT validation

Morphing step

    ZTTValidation

T2W

    combineTool.py -M T2W  -i output/ztt_validation/*/* -o ws.root --parallel 8

run fits

    combineTool.py -m 125 -M MultiDimFit  -d output/ztt_validation/cmb/125/ws.root  --there -n .r_ztt --saveFitResult --setParameterRanges r=0.999,1.001 --expectSignal 1 --cminDefaultMinimizerStrategy=0


make plots

do GOF

We will concider KS and saturated model tests.

KS:
Initially we will perform seperate fits for each category - when we have a full systematic template taking into account correlations ect between categories we can just use one fit then there will be no need to do the next steps for every category individually but rather just the cmb 

Run toys for all cats seperatly:

   combineTool.py -M GoodnessOfFit --algorithm KS  --there -d output/ztt_validation/htt_*_13TeV/125/ws.root -n ".KS.toys" --fixedSignalStrength=1 --there --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 100 -s 0:5:1

Run observed
   combineTool.py -M GoodnessOfFit --algorithm KS  --there -d output/ztt_validation/htt_*_13TeV/125/ws.root -n ".KS" --fixedSignalStrength=1 

Collect output and make plots:

    bins=( "mt_3" "mt_4" "mt_5" "tt_3" "tt_5" "tt_6" "tt_7" "tt_8" "tt_9" )

    for i in "${bins[@]}"
    do
      combineTool.py -M CollectGoodnessOfFit --input output/ztt_validation/htt_"$i"_13TeV/125/higgsCombine.KS.GoodnessOfFit.mH125.root output/ztt_validation/htt_"$i"_13TeV/125/higgsCombine.KS.toys.GoodnessOfFit.mH125.*.root --there -o "$i"_KS.json
      python ../CombineTools/scripts/plotGof.py --statistic KS --mass 125.0 "$i"_KS.json --title-right="60 fb^{-1} (13 TeV)" --output='-KS'
    done

When just using the combined category:

    combineTool.py -M CollectGoodnessOfFit --input output/ztt_validation/cmb/125/higgsCombine.KS.GoodnessOfFit.mH125.root output/ztt_validation/cmb/125/higgsCombine.KS.toys.GoodnessOfFit.mH125.*.root --there -o cmb_KS.json

   python ../CombineTools/scripts/plotGof.py --statistic KS --mass 125.0 cmb_KS.json --title-right="60 fb^{-1} (13 TeV)" --output='-KS'


Saturated model:

For saturated model we always run seperatly for each category

Run toys for all cats seperatly:

   combineTool.py -M GoodnessOfFit --algorithm saturated  --there -d output/ztt_validation/htt_*_13TeV/125/ws.root -n ".saturated.toys" --fixedSignalStrength=1 --there --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 100 -s 0:50:1

Run observed
   combineTool.py -M GoodnessOfFit --algorithm saturated  --there -d output/ztt_validation/htt_*_13TeV/125/ws.root -n ".saturated" --fixedSignalStrength=1 

Collect output and make plots:

    bins=( "mt_3" "mt_4" "mt_5" "tt_3" "tt_5" "tt_6" "tt_7" "tt_8" "tt_9" )

    for i in "${bins[@]}"
    do
      combineTool.py -M CollectGoodnessOfFit --input output/gof/htt_"$i"_13TeV/125/higgsCombine.saturated.GoodnessOfFit.mH125.root output/gof/htt_"$i"_13TeV/125/higgsCombine.saturated.toys.GoodnessOfFit.mH125.*.root --there -o "$i"_saturated.json
      python ../CombineTools/scripts/plotGof.py --statistic saturated --mass 125.0 "$i"_saturated.json --title-right="60 fb^{-1} (13 TeV)" --output="$i"'-saturated'
    done


# Run impacts

First create workspace using top instructions.

Then perform initial fit:

    combineTool.py -M Impacts -d cmb/125/ws.root -m 125 --robustFit 1 -t -1  --doInitialFit --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0

To run impacts for each systematic on crab (RECOMMENDED):
First open `custom_crab.py` and edit the workarea name. 
Make sure you have a valid grid proxy.
Then run:

    combineTool.py -M Impacts -d ws.root -m 125 --robustFit 1 -t  -1  --doFits --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0 --merge 1 --job-mode crab3 --task-name grid-test-impacts --custom-crab custom_crab.py

Otherwise for SGE batch use:

    combineTool.py -M Impacts -d ws.root -m 125 --robustFit 1 -t  -1  --doFits --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0 --merge 1 --job-mode 'SGE'  --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0"

Less recommended:
For lxplus batch use `--job-mode condor --sub-opts='+JobFlavour = "longlunch"` but this is not fully tested (eg. might run out of time).

Collect results:

    combineTool.py -M Impacts -d cmb/125/ws.root -m 125 -o impacts.json

Make impact plot:

    plotImpacts.py -i impacts.json -o impacts

# this seems to help convergence:

combineTool.py -M Impacts -d output/pas_2206/htt_stage2/125/ws.root -m 125 --robustFit 1 --doInitialFit --setParameters alpha=0,muV=1,muggH=1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --setParameterRanges muV=-10,10:muggH,-10,10:alpha=-180,180 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP

combineTool.py -M Impacts -d ../output/pas_2206/cmb/125/ws.root -m 125 --robustFit 1 --doFits --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0,muV=1,muggH=1  --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --setParameterRanges muV=-10,10:muggH,-10,10:alpha=-180,180 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --merge 1 --job-mode 'SGE'  --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0"

# Perform fits plots/fits/GOF of background only categories unrolled in phiCP bins

This is useful if you want to compare data/MC agreement in these completly unblinded categories

first run morphing (use --backgroundOnly=1 for ZTT categories or =2 for jetFakes category) 
'MorphingSMCPDecays18 --output_folder="ztt_checks" --mergeXbbb=true --backgroundOnly=1'

run T2W:

    combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/ztt_checks/cmb/* -o ws.root --parallel 8

Then make plots in the usual way

We could also perform fits / GOF tests for these categories - but this may count as unblinding so might be better to not do this before unblinding. 

To perform KS test:

   combineTool.py -M GoodnessOfFit --algorithm KS  --there -d output/ztt_checks/cmb/125/ws.root -n ".KS.toys" --freezeParameters muggH,alpha,muV --there --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 100 -s 0:5:1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=1

Run observed

   combineTool.py -M GoodnessOfFit --algorithm KS  --there -d output/ztt_checks/cmb/125/ws.root -n ".KS" --freezeParameters muggH,alpha,muV --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=1

Collect outputs and make plots

    combineTool.py -M CollectGoodnessOfFit --input output/ztt_checks/cmb/125/higgsCombine.KS.GoodnessOfFit.mH125.root output/ztt_checks/cmb/125/higgsCombine.KS.toys.GoodnessOfFit.mH125.*.root --there -o cmb_KS.json

    python ../CombineTools/scripts/plotGof.py --statistic KS --mass 125.0 cmb_KS.json --title-right="60 fb^{-1} (13 TeV)" --output='-KS'

# New (Apr 2020) plotting scripts using PostFitShapesFromWorkspace output and MultiDimFit result (1D scan)

Run following set-up commands (after sourcing cmsenv as usual):

    pip3 install --user --upgrade numpy
    pip3 install --user dftools pysge oyaml uproot
    pip3 install --user --upgrade dftools

Alternatively, use your private python3 conda environment.

For CMS plotting style (required) do:

    export MPLCONFIGDIR=./scripts/mpl_configdir/

Also create the output directory called `plots`:

    mkdir plots

## Scans of alpha
To see the full list of options with examples can use `python3 scripts/draw_nll_scans.py --help` (like usual).

To plot 1D scan of alpha using MultiDimFit output (ie. run fit first using above commands):

    python3 scripts/draw_nll_scans.py --input-folder output/01042020/ --channel tt --mode single --plot-name alpha_cmb

There are some options available to plot multiple by category and channel 
(will add by year as well soon).

## 2D scans of kappa
To plot 2D scans of kappa (related to Yukawa couplings) use option `--mode 2d_kappa`, eg:

    python3 scripts/draw_nll_scans.py --input-folder output/01042020 --mode 2d_kappa

## Prefit/postfit distributions
First create workspace for category/merged category of interest using above command.
Then run PostFitShapesFromWorkspace using the workspace.

    combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/test_cp/cmb/* -o ws_ps.root --parallel 8


### Producing prefit shapes
Just use one workspace and produce shapes using `--freeze` option:

For `alpha=0` prefit:

    PostFitShapesFromWorkspace -m 125 -d output/merge/cmb/125/combined.txt.cmb -w output/merge/cmb/125/ws.root --print --total-shapes-bin=true  -o shapes_eff_sm.root

For `alpha=90` prefit:

    PostFitShapesFromWorkspace -m 125 -d output/merge/cmb/125/combined.txt.cmb -w output/merge/cmb/125/ws.root --print --total-shapes-bin=true --freeze alpha=90 -o shapes_eff_ps.root

For multiple channels, can accelerate prefit shapes by looping over folders and
create workspace + PostFitShapes for separate bins.

### Producing postfit shapes

Run MultiDimFit and save fit result:

(the -t -1 has been removed from the command below so be careful if you are not unblinding!)

combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d output/pas_1206_v2/cmb/125/ws.root --algo none  --there -n .bestfit --cminDefaultMinimizerStrategy=0 --saveFitResult

This will create multidimfit.bestfit.root

Note in case of partial unblinding replace cmb with the subdirectory corresponding to the unblinding step

Add `--postfit --sampling -f <fit_result>` to PostFitShapesFromWorkspace command:

this will produce all plots at once in one root file but it tends to be slow:

    PostFitShapesFromWorkspace -m 125 -d output/merge/cmb/125/combined.txt.cmb -w output/merge/cmb/125/ws.root -o shapes_eff.root --print --postfit --sampling --total-shapes-bin=true -f output/30032020_ps/cmb/125/multidimfit.bestfit.root:fit_mdf 

it is better to produce shapes seperate for each year to speed this up (e.g for tt in 2016):

    PostFitShapesFromWorkspace -m 125 -d output/pas_1202/tt_2016/125/combined.txt.cmb -w output/pas_1202/tt_2016/125/ws.root -o shapes_tt_2016.root --print --postfit --sampling --total-shapes-bin=true -f output/pas_1206_v2/cmb/125/multidimfit.bestfit.root:fit_mdf

If this is still too slow you can produce the shapes for each channel seperatly instead

### Drawing distributions

Use `scripts/run_draw1d_cpdecays.py` (with option `--mode prefit` or `--mode posfit`, eg:
    
    python3 scripts/run_draw1d_cpdecays.py --channel tt --year 2016 --draw-signals --signal-scale 50 --mode prefit --datacard shapes_eff_sm_prefitonly_tt_2016.root --alt-datacard shapes_eff_ps_prefitonly_tt_2016.root

If not wanting to draw signals, use `--draw-signals False` and no need to specify `alt-datacard`.
