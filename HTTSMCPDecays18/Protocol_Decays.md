# Creating datacard


    MorphingSMCPDecays18 --output_folder="cpdecay2017" --postfix="-2D" 

the option --no_shape_systs=true can be used as well to remove all shape uncertainties except for bbb's

# Building the workspaces:

    combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/test_cp/cmb/* -o ws.root --parallel 8

# Run maximum liklihood scan

    combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1,lumi_scale=1 --setParameterRanges alpha=-90,90 --points 20 --redefineSignalPOIs alpha  -d output/test_cp/cmb/125/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1 

    If want to scale to some lumi X, include the rate parameter lumi_scale=X in the --setParameters option (scaling 2017+2018 to full run2 = 1.35)

    To run on IC batch use (1 point per job):
    `--job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" --split-points 1`
    To run on lx batch use:
    `--job-mode lxbatch --sub-opts '-q 1nh --split-points 1'

# Plot scan

1D scans can be plotted using scripts/plot1DScan.py script.
To plot alpha:

    python scripts/plot1DScan.py --main=output/test_cp/cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root --POI=alpha --output=alpha --no-numbers --no-box --x_title="#alpha (#circ)" --y-max=0.7

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
      combineTool.py -M CollectGoodnessOfFit --input output/ztt_validation/htt_"$i"_13TeV/125/higgsCombine.saturated.GoodnessOfFit.mH125.root output/ztt_validation/htt_"$i"_13TeV/125/higgsCombine.saturated.toys.GoodnessOfFit.mH125.*.root --there -o "$i"_saturated.json
      python ../CombineTools/scripts/plotGof.py --statistic saturated --mass 125.0 "$i"_saturated.json --title-right="60 fb^{-1} (13 TeV)" --output="$i"'-saturated'
    done

