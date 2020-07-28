# Always do the following commend before running anything (it seems to prevent random seg faults )

ulimit -s unlimited


# creating datacards
    
To create datacards first run morphing:    

    `MorphingSMCP2016 --output_folder="cp261217" --postfix="-2D"`

By default this only uses 2016 data at the moment. To run on 2017 data use option --era=2017. To combine 2016, 2017 and 2018 use option --era="2016,2017,2018"

To not include shape systematics use the option:
    `--no_shape_systs=true`

# Building the workspaces:

Build the to fit alpha with ggH rate floating workspace using

    `combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixture:CPMixture -i output/test/cmb/* -o ws.root --parallel 8`

Build the to fit mu vs alpha build workspace using

    `combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixture:CPMixture -i output/test/cmb/* -o ws.root --parallel 8 --PO fit_2D`

Build the workspace using SM template only in the 0-jet and boosted categories:

    `combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixture:CPMixture -i output/test/cmb/* -o ws.root --parallel 8 --PO sm_fix`

# Run maximum liklihood scan

Worspaces are created with 4 free parameters kappag, kappaW, kappaZ, mutautau and alpha. The first 4 parameters are used to float the ggH, VBF and VH rates.
By using the option --freezeParameters these parameters can be prevented from floating in the fit this can be used to change how the VBF/VH background is treated in the fit e.g. VBF cross-section taken as SM or not.
To take VBF/VH cross-section as SM in both rate and shape use:
    `--freezeParameters kappaW,kappaZ`
To take the H->tautau BR as SM use:
    `--freezeParameters mutautau`

To run with no systematics: '--freezeParameters allConstrainedNuisances'

To scan alpha:

combineTool.py -m 125 -M MultiDimFit --setParameters alpha=0 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d output/cp210720/cmb/125/ws.root --algo grid  --there -n .alpha --floatOtherPOIs 1 --points=21 --alignEdges 1 -t -1 --cminDefaultMinimizerStrategy 0 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10  --cminDefaultMinimizerTolerance=1

(Note we may want to change the --cminDefaultMinimizerStrategy 0 at some point but the fit does not always converget for --cminDefaultMinimizerStrategy 1 - could adjust the number of function calls and/or tolerance to get this to work)

To free all systematics add '--freezeParameters allConstrainedNuisances'

To scan mu:
We need to freeze tautau BR to 1 as the fit has no sensitivity to muggH it only has sensitivity to muggH*mutautau so by setting mutautau to 1 we effecitivly just get one rate parameter that scales the XS*BR
    `combineTool.py -m 125 -M MultiDimFit --setParameters muggH=1 --setParameterRanges muggH=-0.1,3 --points 20 --redefineSignalPOIs muggH --freezeParameters mutautau -d ws.root --algo grid -t -1 --there -n .mu --floatOtherPOIs 1 --cminDefaultMinimizerStrategy 0`

To run on IC batch use (1 point per job):
 `--job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=0:180:0" --split-points 1`
To run on lx batch use:
  `--job-mode lxbatch --sub-opts '-q 1nh --split-points 1'

Run 2--cminDefaultMinimizerStrategy=0D liklihood scan of mu vs alpha using:
    `combineTool.py -m 125 -M MultiDimFit --setParameters muggH=1,alpha=0 --freezeParameters mutautau --setParameterRanges alpha=-90,90:muggH=0,2.5 --redefineSignalPOIs alpha,muggH -d output/cp261118_nobbb/cmb/125/ws.root --there -n ".2DScan" --points 2000 --algo grid -t -1 --parallel=8 --alignEdges`

# Plot scan

1D scans can be plotted using scripts/plot1DScan.py script.
To plot alpha:

  python scripts/plot1DScan.py --main output/cp160720/cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root --POI alpha --output alpha --no-numbers --no-box --x_title "#alpha_{gg} (#circ)" --y-max 5 --x-max=90 --x-min=-90 --rezero --improve

plot fa3 instead (using alpha points)

python scripts/plot1DScanfa3.py --main output/cp210720/cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root --POI fa3 --output fa3_cutbased --no-numbers --no-box --x_title "f_{a3}^{ggH}cos(\phi_{a3}^{ggH})" --y-max 3.5 --x-max=1 --x-min=-1 --rezero --improve

or by channel use:

scripts/plot1DScan.py --main output/cp210720/cmb/125/higgsCombine.alpha.v3.MultiDimFit.mH125.root --POI alpha --output alpha_bychan --no-numbers --no-box --x_title "#alpha_{gg} (#circ)" --x-max=90 --x-min=-90 --others output/cp210720/em_cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root:e#mu:1 output/cp210720/et_cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root:"e#tau":2 output/cp210720/tt_cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root:"#tau#tau":3 output/cp210720/mt_cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root:"#mu#tau":6 --rezero --improve --y-max=2.6 --rezero --improve

Plot 1D scan of mu:

  python scripts/plot1DScan.py --main output/test/cmb/125/higgsCombine.mu.MultiDimFit.mH125.root --POI muggH --output muggH --no-numbers --no-box --x_title "#mu_{ggH}^{#tau#tau}" --y-max 12 --rezero --improve

2D scans can be plotted using scripts/plotMultiDimFit.py script:

    `python scripts/plotMultiDimFit.py --title-right "77.8 fb^{-1} (13 TeV)" --cms-sub "" --mass 125 -o mu_vs_alpha output/cp281118_nobbb/cmb/125/higgsCombine.2DScan.MultiDimFit.mH125.root --x-min -90 --x-max 90`

nclusivly

# Run impacts

First create workspace using top instructions.

Then perform initial fit:

    combineTool.py -M Impacts -d cmb/125/ws.root -m 125 --robustFit 1 -t -1  --doInitialFit --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10

this one seems to work but is slow:

#combine -M MultiDimFit -n blah_initialFit_Test --algo singles --redefineSignalPOIs alpha -t -1 --setParameterRanges alpha=-90,90:muV=-2,4:muggH=0.5,1.5 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd MINIMIZER_analytic --robustFit 1 -v 2 -m 125 -d output/cp210720/cmb/125/ws.root --setParameters alpha=0
    
combineTool.py -M Impacts -d cmb/125/ws.root -m 125 --robustFit 1 -t -1  --doInitialFit --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd MINIMIZER_analytic  --setParameters alpha=0 --setParameterRanges alpha=-90,90:muV=-2,4:muggH=0.5,1.5  --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10

To run impacts for each systematic on crab (RECOMMENDED):
First open `custom_crab.py` and edit the workarea name.
Make sure you have a valid grid proxy.
Then run:

    combineTool.py -M Impacts -d ws.root -m 125 --robustFit 1 -t  -1  --doFits --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10 --merge 1 --job-mode crab3 --task-name grid-test-impacts --custom-crab custom_crab.py

Otherwise for SGE batch use:

    combineTool.py -M Impacts -d ws.root -m 125 --robustFit 1 -t  -1  --doFits --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP  --setParameters alpha=0 --setParameterRanges alpha=-90,90  --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10 --merge 1 --job-mode 'SGE'  --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" 

Less recommended:
For lxplus batch use `--job-mode condor --sub-opts='+JobFlavour = "longlunch"` but this is not fully tested (eg. might run out of time).

Collect results:

    combineTool.py -M Impacts -d cmb/125/ws.root -m 125 -o impacts.json

Make impact plot:

    plotImpacts.py -i impacts.json -o impacts


## Make post-fit plots

Run fit:
   combineTool.py -m 125 -M MultiDimFit --redefineSignalPOIs muggH  -d output/cp110219/htt_01jet/125/ws.root --algo none --there -n .muggH.plots --floatOtherPOIs 1  --saveFitResult

Run postfitshapes:
   PostFitShapesFromWorkspace -m 125 -d output/cp110219/htt_01jet/125/combined.txt.cmb -w output/cp110219/htt_01jet/125/ws.root -o shapes_unblinding_01jets.root --print --sampling --postfit -f output/cp110219/htt_01jet/125/multidimfit.muggH.plots.freeze.root:fit_mdf

you can also freeze a parameter for plotting purposes e.g if you want to show CP-odd signal use --freeze alpha=90

Plot can be produced using the postFitPlotJetFakes.py script, e.g:

`python scripts/postFitPlotJetFakes.py --mode prefit --file_dir htt_tt_2018_6_ -f shapes_prefit_tt_6.root --file_alt=shapes_prefit_tt_6_ps.root --ratio  --proper_errors_uniform --manual_blind --log_y`

Run with (some) systematics frozen
Do fit and store workspace
  'combineTool.py -m 125 -M MultiDimFit  --redefineSignalPOIs alpha -d output/cp260219/cmb/125/ws.root --algo none  --there -n .saveWS  --saveWorkspace'

Now run scans freezing all systematics with "--freezeParameters all" - can also freeze individual systematics or groups
  'combineTool.py -m 125 -M MultiDimFit --setParameters alpha=0 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d output/cp260219/cmb/125/higgsCombine.saveWS.MultiDimFit.mH125.root --algo grid --there -n .alpha.nosyst2 --floatOtherPOIs 1 --points=37 --alignEdges 1 --freezeParameters all --snapshot MultiDimFit'

can also do the following to freeze theory uncertanties only:

combineTool.py -m 125 -M MultiDimFit --setParameters alpha=0 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d output/cp260219/cmb/125/higgsCombine.saveWS.MultiDimFit.mH125.root --algo grid --there -n .alpha.notheory --floatOtherPOIs 1 --points=37 --alignEdges 1  --snapshot MultiDimFit --freezeParameters CMS_scale_gg_13TeV,CMS_FiniteQuarkMass_13TeV,CMS_PS_ggH_13TeV,CMS_UE_ggH_13TeV,BR_htt_THU,BR_htt_PU_mq,BR_htt_PU_alphas,QCDScale_ggH,QCDScale_qqH,QCDScale_WH,QCDScale_ZH,pdf_Higgs_WH,pdf_Higgs_ZH,pdf_Higgs_gg,pdf_Higgs_qq,CMS_ggH_mig01,CMS_ggH_mig12 --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=0:180:0" --split-points 1 --task-name alpha.notheory
