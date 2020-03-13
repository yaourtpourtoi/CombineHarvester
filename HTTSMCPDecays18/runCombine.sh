#!/bin/sh

ERA=all
DIRECTORY=test_desy_Run2
CHANNEL=Combined

cp /nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/Jan20/CMSSW_10_2_16/src/HiggsCP/output_cards_newDNN_2016.root ./shapes/DESY/2016/htt_mt.inputs-sm-13TeV.root
cp /nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/Jan20/CMSSW_10_2_16/src/HiggsCP/output_cards_newDNN_2017.root ./shapes/DESY/2017/htt_mt.inputs-sm-13TeV.root
cp /nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/Jan20/CMSSW_10_2_16/src/HiggsCP/output_cards_newDNN_2018.root ./shapes/DESY/2018/htt_mt.inputs-sm-13TeV.root

MorphingSMCPDecays18 --output_folder=$DIRECTORY --input_folder_mt="DESY" --era=$ERA --do_jetfakes=true

combineTool.py -M T2W -P CombineHarvester.CombinePdfs.CPMixtureDecays:CPMixtureDecays -i output/$DIRECTORY/htt_mt_${CHANNEL}_13TeV/* -o ws.root --parallel 8 

combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1,lumi_scale=1 --setParameterRanges alpha=-90,90 --points 20 --redefineSignalPOIs alpha -d output/$DIRECTORY/htt_mt_${CHANNEL}_13TeV/125/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1

python scripts/plot1DScan.py --main=output/$DIRECTORY/htt_mt_${CHANNEL}_13TeV/125/higgsCombine.alpha.MultiDimFit.mH125.root --POI=alpha --output=alpha --no-numbers --no-box --x_title="#alpha (#circ)" --y-max=0.7
