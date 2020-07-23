# run combined years plots using:
#python3 scripts/run_draw1d_cpdecays_ztt.py --channel mt --year 2018  --mode postfit --datacard shapes_ztt_cmb.root --cmb_years

import oyaml as yaml
import pandas as pd
import numpy as np
import argparse

from plotting_ztt import (
    draw_1d, 
    create_df,
    var_kw,
    process_kw,
    nbins_kw,
)

def parse_arguments():
    epilog = (
        "Example:\n \n"
        "python3 scripts/run_draw1d_cpdecays.py --channel mt --year 2016 "
        "--draw-signals --signal-scale 50 --mode prefit " 
        "--datacard shapes_sm_eff.root --alt-datacard shapes_ps_eff.root "
    )
    parser = argparse.ArgumentParser(epilog=epilog)


    parser.add_argument(
        "--channel", default="tt", choices=["tt", "mt"],
        help="Which channel to use",
    )
    parser.add_argument(
        "--year", default=2018,
        help="Which year to use",
    )
    parser.add_argument(
        "--draw-signals", action='store_true', default=False,
        help="Draw signals?",
    )
    parser.add_argument(
        "--signal-scale", default=1.,
        help="Scale the signal by this value (not for 'unrolled' plots)",
    )
    parser.add_argument(
        "--mode", default="prefit", choices=["prefit", "postfit"],
        help="Which histograms to use",
    )
    parser.add_argument(
        "--datacard", default="shapes_eff.root",
        help="Path to PostFitShapes datacard",
    )
    parser.add_argument(
        "--alt-datacard", default=None,
        help="If set, use alternative template (for signal) found at this path",
    )
    parser.add_argument(
        "--no-ff", action="store_true", default=False,
        help="Don't use fake factors",
    )
    parser.add_argument(
        "--no-embedding", action="store_true", default=False,
        help="Don't use embedded samples",
    )
    parser.add_argument(
        "--cmb_years", action="store_true", default=False,
        help="Do combined years",
    )

    arguments = parser.parse_args()

    # Use ff and embedding by default, but making it less confusing
    # for users
    arguments.ff = not arguments.no_ff
    del arguments.no_ff

    arguments.embedding = not arguments.no_embedding
    del arguments.no_embedding

    return arguments

def draw1d_cpdecays(
    channel, year, draw_signals, signal_scale, ff, embedding, mode,
    datacard, alt_datacard, cmb_years
):

    # Plotting SM and PS template
    signals = []
    if draw_signals:
        # signals = ["H_sm", "H_ps"]
        signals = ["Bestfit", "H_ps"]

    leg_kw = {"offaxis": True, "fontsize": 9, "labelspacing":0.12,}

    # The following channel kwargs config file define the 
    # merging for the plotting, eg. merge all SM higgs signals into H_sm.
    # Should always use scripts/plot_kw_postfit.yaml when using ff and embedding
    ch_kw = {}
    with open("scripts/plot_kw.yaml", "r") as f:
        ch_kw = yaml.safe_load(f)
    if ff: # always use FF
        ch_kw = {}
        with open("scripts/plot_kw_postfit_ztt.yaml", "r") as f:
            ch_kw = yaml.safe_load(f)
    #if embedding: # always use embedding
    #    for ch, proc in ch_kw.items():
    #        if ch in ["tt", "mt", "et", "em",]:
    #            proc["EmbedZTT"] = ["EmbedZTT"]
    #            del proc["ZTT"]

    # Histogram processes to load in
    # By default we use fake factors and embedding
    if embedding and ff:
        if channel == "tt":
            processes = ['data_obs', 'EmbedZTT', 'ZL', 'TTT', 'VVT', 'jetFakes','Wfakes']
        elif channel == "mt":
            processes = ['data_obs', 'ZTT', 'ZL', 'TTT', 'VVT', 'jetFakes']
    elif ff:
        processes = ['data_obs', 'ZTT', 'ZL', 'TTT', 'VVT', 'jetFakes', 'EWKZ',]
    elif embedding:
        processes = [
            'data_obs', 'EmbedZTT', 'ZL', 'TTT', 'VVT', 'VVJ', 'W', 'QCD', 'ZJ'
        ]
    else:
        processes = [
            'data_obs', 'ZTT', 'ZL', 'ZJ', 'TTT', 'TTJ', 'VVT', 'VVJ', 
            'W', 'QCD', 'EWKZ',
        ]

    if len(signals) > 0:
        processes.extend([
            "TotalSig",
            # "ggH_sm_htt", "qqH_sm_htt", "ZH_sm_htt", "WH_sm_htt",
            # "ggH_ps_htt", "qqH_ps_htt", "ZH_ps_htt", "WH_ps_htt",
        ])
        
    # Draw categories (defined in nbins_kw in plotting.py):
    # 1-2: backgrounds, 3+: signal (higgs) categories
    # correspond to CH bins defined in Morphing scripts
    if channel == "tt":
        bins_to_plot = list(range(1,12))
        bins_to_plot = [1,2,3,7]
    elif channel == "mt":
        bins_to_plot = list(range(1,7))
        bins_to_plot = [3,4,5,6,30,40,50,60]
    for bin_number in bins_to_plot:
        #category = nbins_kw[channel][bin_number][3]
        category = 'ztt'
        # Initialise empty and change depending on category bellow
        plot_var = ""

        # Making use of python3 f-strings
        directory = f"htt_{channel}_{year}_{bin_number}_13TeV_{mode}"

        print(f"Doing category {category}")
        if category in ["embed", "fakes"]:
            # MVA score plots for background categories
            if channel == "tt":
                plot_var = "BDT_score"
            elif channel == "mt":
                plot_var = "NN_score"
            partial_blind = False
            blind = False
            unrolled = False
            norm_bins = True
        elif "signal" in category:
            # signal inclusive category, added blind option
            if channel == "tt":
                plot_var = "BDT_score"
            elif channel == "mt":
                plot_var = "NN_score"
            partial_blind = False
            blind = True # blind all of data for signal category
            unrolled = False
            norm_bins = True
        else:
            # 'unrolled' category plots
            plot_var = "Bin_number"
            partial_blind = False # unblind only first window of 'unrolled'
            blind = False
            unrolled = False
            norm_bins = False


        signal_scale = 1. # no need to scale on log plot

        # Create dataframe to plot
        df_plot = create_df(datacard, directory, channel, processes, ch_kw)

        # In order to have alternative (PS) hypothesis create dataframe for this
        # and then replace PS hypothesis in df_plot by PS entries here.
        # This is because, with PostFitShapesFromWorkspace, we don't have any
        # entries for PS signals when SM (alpha=0) 
        if draw_signals:
            if cmb_years: df_plot_alt = create_df(alt_datacard, directory, channel, ['TotalSigPS'], ch_kw)
            else: df_plot_alt = create_df(alt_datacard, directory, channel, processes, ch_kw)
            if df_plot_alt is not None:
                df_plot_alt.reset_index(inplace=True)
                df_plot_alt.loc[df_plot_alt["parent"] == "Bestfit", "parent"] = "H_ps"
                df_plot_alt.set_index(["parent","binvar0","binvar1"], inplace=True)
                df_plot = pd.concat([
                    df_plot,
                    df_plot_alt.loc[
                        df_plot_alt.index.get_level_values("parent") == "H_ps"
                    ]
                ], axis='index', sort=False)
        if partial_blind:
            # Unblind first window of unrolled bins only (for now)
            data_mask = df_plot.index.get_level_values("parent") == "data_obs"
            blind_mask = df_plot.index.get_level_values("binvar0") >= \
                nbins_kw[channel][bin_number][1]
            df_plot.loc[data_mask & blind_mask, "sum_w"] = np.nan

        # Always use mcstat=True and mcsyst=True when plotting systematic unc.
        category_=category
        if 'fakes' not in category and 'ztt' not in category: category_='higgs'
        year_ = year
        if cmb_years: year_='cmb'
        draw_1d(
            df_plot, plot_var, channel, category_, year_, blind=False, sigs=signals, 
            signal_scale=signal_scale, ch_kw=ch_kw, process_kw=process_kw, 
            var_kw=var_kw, leg_kw=leg_kw, unrolled=unrolled, norm_bins=norm_bins,
            nbins=nbins_kw[channel][bin_number], mcstat=True, mcsyst=True,
            logy=False, sm_bkg_ratio=True, postfix=mode,
        )

if __name__ == "__main__":
    # Unpack arguments from argparse so that the options can be passed directly 
    # into the function.
    draw1d_cpdecays(**vars(parse_arguments()))
