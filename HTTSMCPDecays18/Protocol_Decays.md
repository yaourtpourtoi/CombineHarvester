# Always do the following commend before running anything (it seems to prevent random seg faults )

ulimit -s unlimited

# Creating datacard

    MorphingSMCPDecays18 --output_folder="cpdecay2017" --postfix="-2D" 

the option --no_shape_systs=true can be used as well to remove all shape uncertainties except for bbb's

# Building the workspaces:

    combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/test_cp/cmb/* -o ws.root --parallel 8

# Run maximum likelihood scan

    combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90 --points 21 --redefineSignalPOIs alpha  -d output/test_cp/cmb/125/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1 --cminDefaultMinimizerStrategy=0 

    If want to scale to some lumi X, include the rate parameter lumi_scale=X in the --setParameters option (scaling 2017+2018 to full run2 = 1.35)

    To run on IC batch use (1 point per job):
    `--job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" --split-points 1`
    To run on lx batch use:
    `--job-mode lxbatch --sub-opts '-q 1nh --split-points 1'

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


# perform ZTT validation

Morphing step

    ZTTValidation

T2W

    combineTool.py -M T2W  -i output/ztt_validation/*/* -o ws.root --parallel 8

run fits

    combineTool.py -m 125 -M MultiDimFit  -d output/ztt_validation/htt_tt_3_13TeV/125/ws.root  --there -n .r_ztt --saveFitResult --setParameterRanges r=0.999,1.001

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

first perform initial fit:

  'combineTool.py -M Impacts -d cmb/125/ws.root -m 125 --robustFit 1 -t -1  --doInitialFit --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0'

then run impact with:

  'combineTool.py -M Impacts -d cmb/125/ws.root -m 125 --robustFit 1 -t -1  --doFits --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0 --job-mode 'SGE'  --prefix-file ic --sub-opts "-q hep.q -l h_rt=0:180:0" --merge=1'

Collect results:

  `combineTool.py -M Impacts -d cmb/125/ws.root -m 125 -o impacts.json`

Make impact plot:

  `plotImpacts.py -i impacts.json -o impacts`

Perform fits plots/fits/GOF of background only categories unrolled in phiCP bins
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

Alternatively, use your private python3 conda environment.

For CMS plotting style (required) do:

    export MPLCONFIGDIR=./scripts/mpl_configdir/

Also create the output directory called `plots`:

    mkdir plots

## Scans of alpha
To plot 1D scan of alpha using MultiDimFit output (ie. run fit first using above commands):

    python3 scripts/draw_nll_scans.py --input-folder output/01042020/ --channel tt --mode single --plot-name alpha_cmb

There are some options available to plot multiple by category and channel 
(will add by year as well soon).

## Prefit/postfit distributions
First create workspace for category/merged category of interest using above command.
Then run PostFitShapesFromWorkspace using the workspace.

For prefit both of these steps need to be done with `alpha=0` and `alpha=90` as initial values
such that we can plot both SM and PS distributions on prefit plots.
To change these need to do (using vim or any other text editor):

    vim ../CombinePdfs/python/CPMixtureDecays.py

Change `alpha[0,-90,90] --> alpha[90,-90,90]` for PS and rerun command for workspace creation with different workspace name (`ws.root --> ws_ps.root`), eg:

    combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/test_cp/cmb/* -o ws_ps.root --parallel 8

Otherwise, only using Asmiov of SM (`alpha=0`).

### Producing prefit shapes

For `alpha=0` prefit:

    PostFitShapesFromWorkspace -m 125 -d output/merge/cmb/125/combined.txt.cmb -w output/merge/cmb/125/ws.root --print -o shapes_eff_sm.root

For `alpha=90` prefit:

    PostFitShapesFromWorkspace -m 125 -d output/merge/cmb/125/combined.txt.cmb -w output/merge/cmb/125/ws_ps.root --print -o shapes_eff_ps.root

(You can just use the same workspace for the above and use --freeze alpha=90)

For multiple channels, can accelerate prefit shapes by looping over folders and
create workspace + PostFitShapes for separate bins.

### Producing postfit shapes

Run MultiDimFit and save fit result:

    combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90 --points 21 --redefineSignalPOIs alpha  -d output/test_cp/cmb/125/ws.root --algo none -t -1 --there -n .alpha --saveFitResult

This will create multidimfit.alpha.root

Add `--postfit --sampling -f <fit_result>` to PostFitShapesFromWorkspace command:

    PostFitShapesFromWorkspace -m 125 -d output/merge/cmb/125/combined.txt.cmb -w output/merge/cmb/125/ws.root -o shapes_eff.root --print --postfit --sampling -f output/30032020_ps/cmb/125/multidimfit.alpha.root:fit_mdf 


### Drawing distributions

Use `scripts/run_draw1d_cpdecays.py` (with option `--mode prefit` or `--mode posfit`, eg:
    
    python3 scripts/run_draw1d_cpdecays.py --channel tt --year 2016 --draw-signals --signal-scale 50 --mode prefit --datacard shapes_eff_sm_prefitonly_tt_2016.root --alt-datacard shapes_eff_ps_prefitonly_tt_2016.root

If not wanting to draw signals, use `--draw-signals False` and no need to specify `alt-datacard`.
