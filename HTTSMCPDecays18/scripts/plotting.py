import glob
import os
import oyaml as yaml
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import dftools
import scipy
import uproot
from tqdm.auto import tqdm
mpl.use('pdf')
plt.style.use('cms')

def create_df(
    datacard, directory, channel, processes, ch_kw={},
):
    df = pd.DataFrame()
    _file = uproot.open(datacard)
    for process in processes:
        if process not in _file[directory]:
            # in CH we remove processes with 0 yield
            # can add them back here in this case (using near zero event weights)
            # TO DO: can be implemented directly in CH (Morphing script)
            df = pd.concat([df, pd.DataFrame({
                "varname0": ["var"] * nbins,
                "binvar0": bin_edges[:-1],
                "binvar1": bin_edges[1:],
                "sum_w": [1e-10]*nbins,
                "sum_ww": [1e-10]*nbins, 
                "parent": process,
            })], axis='index', sort=False)
            continue
        hist = _file["{}/{}".format(directory, process)]
        # print(dir(hist))
        bins = hist.bins
        bin_edges = hist.edges
        nbins = hist.numbins
        weights = hist.values
        variance = hist.variances
        
        df = pd.concat([df, pd.DataFrame({
            "varname0": ["var"] * nbins,
            "binvar0": bin_edges[:-1],
            "binvar1": bin_edges[1:],
            "sum_w": weights,
            "sum_ww": variance,
            "parent": hist.name.decode(),
        })], axis='index', sort=False)
    df.set_index(["parent","binvar0","binvar1"], inplace=True)
    
    df = dftools.transform.merge(df, ch_kw[channel])
    
    return df

def draw_1d(
    df_, plot_var, channel, year, blind, sigs=[], 
    signal_scale=1., ch_kw={}, process_kw={}, var_kw={}, 
    leg_kw={}, sig_kw={}, norm_mc=False, fractions=False,
    unrolled=False, nbins=[[4], 14, "", ""], mcstat=True,
    sig_ratio=False, norm_bins=False, mcsyst=True,
):
    # to keep the original one the same
    df = df_.copy(deep=True)

    # low edges
    binning = df.index.get_level_values("binvar0").unique()
    if norm_bins:
        df["binwidth"] = df.reset_index().eval("binvar1-binvar0").values
        df["sum_w"] = df.eval("sum_w/binwidth")
        df["sum_ww"] = df.eval("sum_ww/(binwidth**2)")
    
    with mpl.backends.backend_pdf.PdfPages(
        f"plots/{plot_var}_{nbins[3]}_{channel}_{year}.pdf",
        keep_empty=False,
    ) as pdf:
        fig, ax = plt.subplots(
            figsize=(2.9, 3.2), dpi=200,
            nrows=2, ncols=1,
            sharex=True, sharey=False,
            gridspec_kw={"height_ratios": (3, 1), "hspace": 0.1, "wspace": 0.1},
        )
        if unrolled:
            fig.set_size_inches(5.3, 3.2)
    
        if year == "2016":
            lumi = 35.9
        elif year == "2017": 
            lumi = 41.5
        elif year == "2018": 
            lumi = 59.7
        
        if unrolled:
            dftools.draw.cms_label(ax[0], "Preliminary", lumi=lumi, extra_label=nbins[2])
        else:
            dftools.draw.cms_label(ax[0], "Preliminary", lumi=lumi)
        
        # To fix when y axis is too large, scientific notation starts showing 
        # up in top left and overlaps with CMS logo
        scale_by_tenthou = False
        if not unrolled and (df["sum_w"] > 1e5).any():
            scale_by_tenthou = True
            df["sum_w"] = df.eval("sum_w/1e5")
            df["sum_ww"] = df.eval("sum_ww/(1e5**2)")


        data_mask = df.index.get_level_values("parent") != "data_obs"
        
        df_data = df.loc[~data_mask,:]
        df_mc = df.loc[data_mask,:]
        if norm_mc:
            df_mc = df_data.sum() * df_mc / df_mc.sum() 

        # scale signals if set
        sig_mask = ~df_mc.index.get_level_values("parent").isin(sigs)
        df_mc.loc[~sig_mask, "sum_w"] = df_mc.loc[~sig_mask, "sum_w"].copy(deep=True) * signal_scale
        df_mc.loc[~sig_mask, "sum_ww"] = df_mc.loc[~sig_mask, "sum_ww"].copy(deep=True) * signal_scale**2

        ymax = max([
            df_data["sum_w"].max(), 
            df_mc["sum_w"].max() 
        ])
        ymc_max = max([
            df_mc["sum_w"].max(), 
        ])
            
        if len(sigs) == 0:
            ymax *= 1.6
            leg_kw = {
                "offaxis": False, "fontsize": 9, "labelspacing":0.12,
                "ncol": 2, "loc": 9, 
            }
        elif len(sigs) > 0:
            if channel == "tt":
                ymax *= 1.6
            elif channel == "mt":
                ymax *= 1.8
            leg_kw = {
                "offaxis": False, "fontsize": 7, "labelspacing":0.12,
                "ncol": 2, "loc": 9, "framealpha": 0.6,
            }
        if unrolled:
            leg_kw = {
                "offaxis": True, "fontsize": 9, "labelspacing":0.14,
            }

        process_kw["labels"]["H_sm"] = r'$\mathrm{SM\ H} \rightarrow\tau\tau$'
        process_kw["labels"]["H_ps"] = r'$\mathrm{PS\ H} \rightarrow\tau\tau$'
        if signal_scale != 1.:
            process_kw["labels"]["H_sm"] = f'${signal_scale}\\times \\mathrm{{SM\ H}} \\rightarrow\\tau\\tau$'
            process_kw["labels"]["H_ps"] = f'${signal_scale}\\times \\mathrm{{PS\ H}} \\rightarrow\\tau\\tau$'
         
        # For MC use poisson errors by default
        # If using Gaussian errors (like the one from CH PostFitShapes) use symmetric errors
        # can use KW and set these only when using mcsyst
        if mcsyst:
            interval_func = lambda x, variance: (x-np.sqrt(variance), x+np.sqrt(variance))

        dftools.draw.data_mc(
            ax, df_data, df_mc, "binvar0", binning, 
            log=unrolled, legend=True,
            proc_kw=process_kw, legend_kw=leg_kw,
            sigs=sigs,
            blind=blind,
            add_ratios=fractions, 
            mcstat=mcstat, 
            mcstat_top=mcstat,
            interval_func=interval_func,
        )
        
        if not unrolled:
            ax[0].set_ylim(0., ymax)
            if blind:
                ax[0].set_ylim(0., ymc_max*2.)
        ax[0].set_ylabel(r'Events')
        if norm_bins and scale_by_tenthou:
            ax[0].set_ylabel(r'$10^{5}\ \mathrm{events}\ /\ \mathrm{bin}$')
        elif norm_bins:
            ax[0].set_ylabel(r'Events / bin')
        ax[1].set_ylabel(r'Ratio')
        
        if unrolled:
            #ax[0].set_ylim(1e-1, 1e3)
            ax[0].set_ylim(1e-1, ymc_max*10)
            ax[1].set_ylim(0, 2)
            ax[1].set_yticks([0.5, 1., 1.5])
            
            binvar0_bins = len(binning)
            vert_lines = list(range(nbins[1], binvar0_bins, int(nbins[1])))
            ax[0].vlines(vert_lines, *ax[0].get_ylim(), linestyles='--', colors='black', zorder=20)
            ax[1].vlines(vert_lines, *ax[0].get_ylim(), linestyles='--', colors='black', zorder=20)
            
            ax[1].set_xticks(list(range(0, binvar0_bins, int(nbins[1]))))

            # annotate windows (BDT score)
            widebin_cents = list(range(int(nbins[1]/2), binvar0_bins, int(nbins[1])))
            _, ypos = ax[0].transData.inverted().transform(ax[0].transAxes.transform((0, 0.95)))
            for idx, xpos in enumerate(widebin_cents):
                if channel == "tt":
                    ftsize = 9
                elif channel == "mt":
                    ftsize = 6
                ax[0].text(
                    xpos, ypos, f"({nbins[0][idx]}, {nbins[0][idx+1]})", 
                    ha='center', va='top', fontsize=ftsize,
                )
                
        try:
            ax[1].set_xlabel(var_kw[plot_var])
        except KeyError:
            print(f"{plot_var} not defined in var_kw")
            ax[1].set_xlabel(plot_var.replace("_"," "))
        
        if sig_ratio:
            denom_mask = df_mc.index.get_level_values("parent") != "H_sm"
            bin_edges, bin_cents = dftools.draw.bin_lows_to_edges_cents(binning)
            for sig_name in ["H_sm", "H_ps"]:
                numer_mask = df_mc.index.get_level_values("parent") != sig_name
                sig_ratio = df_mc.loc[~numer_mask, "sum_w"].values / df_mc.loc[~denom_mask, "sum_w"].values
                sig_ratio_ww = df_mc.loc[~numer_mask, "sum_ww"].values / df_mc.loc[~denom_mask, "sum_w"].values**2
                ax[1].hist(
                    binning, bins=list(binning)+[binning[-1] + df["binwidth"].values[0]], 
                    weights=sig_ratio, histtype='step', lw=1,
                    color=process_kw["colours"][sig_name],
                    zorder=-1,
                )
                up = sig_ratio + np.sqrt(sig_ratio_ww)
                down = sig_ratio - np.sqrt(sig_ratio_ww)
                ax[1].fill_between(
                    bin_edges, list(up)+[up[-1]], list(down)+[down[-1]],
                    step='post', color=process_kw["colours"][sig_name],
                    alpha=0.2, zorder=-1,
                )

        ax[1].set_yticks([0.6, 0.8, 1., 1.2, 1.4])
        #ax[1].axhline(1.2, ls='--', color='gray')
        #ax[1].axhline(0.8, ls='--', color='gray')
        ax[1].set_ylim(0.6, 1.4)
        if unrolled:
            ax[1].set_ylim(0., 2.)
            ax[1].set_yticks([0., 0.5, 1., 1.5, 2.])
        fig.align_labels(ax)

        print(f"Saving plots/{plot_var}_{nbins[3]}_{channel}_{year}")
        pdf.savefig(fig, bbox_inches='tight')

##################### KWs

var_kw = {
    "m_vis": r'$m_{\tau_{h}\tau_{h}} (\mathrm{GeV})$',
    "svfit_mass": r'$m_{\tau\tau} (\mathrm{GeV})$',
    "svfit_mass_err": r'$m_{\tau\tau}^{\mathrm{error}} (\mathrm{GeV})$',
    "n_jets": r'$n_{\mathrm{jets}}$',
    "n_btag": r'$n_{\mathrm{b-tag}}^{\mathrm{medium}}$',
    "n_loose_btag": r'$n_{\mathrm{b-tag}}^{\mathrm{loose}}$',
    "mjj": r'$m_{jj} (\mathrm{GeV})$',
    "jdeta": r'$\Delta\eta(\vec{p}_{j_{1}}, \vec{p}_{j_{2}})$',
    "jpt_1": r'$p_{\mathrm{T}}^{j_{1}} (\mathrm{GeV})$',
    "jpt_2": r'$p_{\mathrm{T}}^{j_{2}} (\mathrm{GeV})$',
    "jeta_1": r'$\eta_{j_{1}}$',
    "jeta_2": r'$\eta_{j_{2}}$',
    "shifted_dphi_jtt": r'$\Delta\phi(\vec{p}_j,\vec{p}_{\mathrm{Z}})$',
    "shifted_dphi_jtt_smear": r'$\Delta\phi(\vec{p}_j,\vec{p}_{\mathrm{Z}})$',
    "shifted_dphi_j20tt": r'$\Delta\phi(\vec{p}_j,\vec{p}_{\mathrm{Z}})$',
    "shifted_dphi_j20tt_smear": r'$\Delta\phi(\vec{p}_j,\vec{p}_{\mathrm{Z}})$',
    "sjdphi": r'$\Delta\phi(\vec{p}_{j_1}, \vec{p}_{j_2})$',
    "sjdphi_smear": r'$\Delta\phi(\vec{p}_{j_1}, \vec{p}_{j_2})$',
    "pt_tt": r'$p_{\mathrm{T}}^{\tau\tau} (\mathrm{GeV})$',
    "pt_vis": r'$p_{\mathrm{T}}^{\tau_h\tau_h} (\mathrm{GeV})$',
    "pt_1": r'$p_{\mathrm{T}}^{\tau_{h_{1}}} (\mathrm{GeV})$',
    "pt_2": r'$p_{\mathrm{T}}^{\tau_{h_{2}}} (\mathrm{GeV})$',
    "eta_1": r'$\eta_{\tau_{h_{1}}}$',
    "eta_2": r'$\eta_{\tau_{h_{2}}}$',
    "tau_decay_mode_1": r'$\tau_{h_{1}}\ \mathrm{decay\ mode}$',
    "tau_decay_mode_2": r'$\tau_{h_{2}}\ \mathrm{decay\ mode}$',
    "mva_dm_1": r'$\tau_{h_{1}}\ \mathrm{MVA\ decay\ mode}$',
    "mva_dm_2": r'$\tau_{h_{2}}\ \mathrm{MVA\ decay\ mode}$',
    "met": r'$p_{\mathrm{T}}^{\mathrm{miss}} (\mathrm{GeV})$',
    "residual_pt": r'$(\vec{p}_{\mathrm{MET}} + \vec{p}_\mathrm{jet} + \vec{p}_\mathrm{Z})_{\mathrm{T}}\ (\mathrm{GeV})$',
    "IC_15Mar2020_max_score": r'BDT score',
    "NN_score": r'NN score',
    "Bin_number": r'Bin number',
}

process_kw={
    "labels": {
        "SMTotal": "SM Total", 
        #"Backgrounds": "Bkgs", 
        "Backgrounds": "Minors", 
        "Minors": "Minors", 
        "ZL": r'$\mathrm{Z}\rightarrow \mu\mu$',
        "QCD": "QCD",
        "TT": r'$t\bar{t}$',
        "Electroweak": "Electroweak",
        "ZTT": r'$\mathrm{Z}\rightarrow\tau\tau$',
        "ZMM": r'$\mathrm{Z}\rightarrow \mu\mu$',
        "ZEE": r'$\mathrm{Z}\rightarrow ee$',
        "jetFakes": r'$\mathrm{jet}\rightarrow \tau_{h}$',
        "EmbedZTT": r'$\mu\rightarrow\tau \ \mathrm{Embed.}$',
        "ggH": r'$gg\mathrm{H} \rightarrow\tau\tau$',
        "qqH": r'$qq\mathrm{H} \rightarrow\tau\tau$',
        "VH": r'$\mathrm{VH} \rightarrow\tau\tau$',
        "H_sm": r'$\mathrm{SM\ H} \rightarrow\tau\tau$',
        "H_ps": r'$\mathrm{PS\ H} \rightarrow\tau\tau$',
    },
    "colours": {
        "SMTotal": 'black', 
        "Backgrounds": "#d9d9d9", 
        "Minors": "#d9d9d9",
        "ZL": "#64C0E8",
        "QCD": "#ffb8c9",
        "TT": "#9B98CC",
        "Electroweak": "#DE5A6A",
        "ZTT": "#E8AD46",
        "ZMM": "#64C0E8",
        "ZEE": "#64C0E8",
        "jetFakes": "#addd8e",
        "EmbedZTT": "#E8AD46",
        "ggH": "#ef3b2c",
        "qqH": "#2171b5",
        "VH": "#c994c7",
        "H_sm": "#4292c6",
        "H_ps": "#2ca25f", # dark green
    },
}

nbins_kw = {
    "tt": {
        1: [[None], 1, "embed", "embed"], # embed
        2: [[None], 1, "fakes", "fakes"], # fakes
        3: [[0., 0.7, 0.8, 0.9, 1.], 16, r'$\rho\rho$', "rho-rho"], # rho-rho
        4: [[0., 0.7, 0.8, 0.9, 1.], 4, r'$a_{1}^{1\mathrm{pr}}\rho + a_{1}^{1\mathrm{pr}}a_{1}^{1\mathrm{pr}}$', "0a1-rho_0a1-0a1"], # 0a1-rho + 0a1-0a1
        5: [[0., 0.7, 0.8, 0.9, 1.], 8, r'$a_{1}^{3\mathrm{pr}}\rho$', "a1-rho"], # a1-rho
        6: [[0., 0.7, 0.8, 0.9, 1.], 4, r'$a_{1}^{3\mathrm{pr}}a_{1}^{3\mathrm{pr}}$', "a1-a1"], # a1-a1
        7: [[0., 0.7, 0.8, 0.9, 1.], 16, r'$\pi\rho$', "pi-rho"], # pi-rho
        8: [[0., 0.7, 0.8, 1.], 6, r'$\pi\pi$', "pi-pi"], # pi-pi
        9: [[0., 0.7, 0.8, 0.9, 1.], 4, r'$\pi a_{1}^{3\mathrm{pr}}$', "pi-a1"], # pi-a1
        10: [[0., 0.7, 0.8, 0.9, 1.], 4, r'$\pi a_{1}^{1\mathrm{pr}}$', "pi-0a1"], # pi-0a1
        11: [[0., 0.7, 0.8, 0.9, 1.], 4, r'$a_{1}^{3\mathrm{pr}} a_{1}^{1\mathrm{pr}}$', "a1-0a1"], # a1-0a1
    },
    "mt": {
        1: [[None], 1, "embed", "embed"], # embed
        2: [[None], 1, "fakes", "fakes"], # fakes
        3: [[0.0, 0.45, 0.6, 0.7, 0.8, 0.9, 1.0], 16, r'$\mu\rho$', "mu-rho"], # mu-rho
        4: [[0.0, 0.45, 0.6, 0.7, 0.8, 0.9, 1.0], 12, r'$\mu\pi$', "mu-pi"], # mu-pi
        5: [[0.0, 0.45, 0.6, 0.7, 0.8, 0.9, 1.0], 8, r'$\mu a_{1}^{3\mathrm{pr}}$', "mu-a1"], # mu-a1
        6: [[0.0, 0.45, 0.6, 0.7, 0.8, 0.9, 1.0], 4, r'$\mu a_{1}^{1\mathrm{pr}}$', "mu-0a1"], # mu-0a1
    },
}

nllscan_kw = {
    "tt": {
        0: [r"$\tau_h\tau_h\ \mathrm{combined}$", "combined", "#DE5A6A"],
        1: ["embed", "embed", ""], # embed
        2: ["fakes", "fakes", ""], # fakes
        3: [r'$\rho\rho$', "rho-rho", "#9B98CC"], # rho-rho
        4: [r'$a_{1}^{1\mathrm{pr}}\rho + a_{1}^{1\mathrm{pr}}a_{1}^{1\mathrm{pr}}$', "0a1-rho_0a1-0a1", "#d9d9d9"], # 0a1-rho + 0a1-0a1
        5: [r'$a_{1}^{3\mathrm{pr}}\rho$', "a1-rho", "#64C0E8"], # a1-rho
        6: [r'$a_{1}^{3\mathrm{pr}}a_{1}^{3\mathrm{pr}}$', "a1-a1", "#ffb8c9"], # a1-a1
        7: [ r'$\pi\rho$', "pi-rho", "#E8AD46"], # pi-rho
        8: [ r'$\pi\pi$', "pi-pi", "#addd8e"], # pi-pi
        9: [r'$\pi a_{1}^{3\mathrm{pr}}$', "pi-a1", "#c994c7"], # pi-a1
        10: [r'$\pi a_{1}^{1\mathrm{pr}}$', "pi-0a1", "#2ca25f"], # pi-0a1
        11: [r'$a_{1}^{3\mathrm{pr}} a_{1}^{1\mathrm{pr}}$', "a1-0a1", "#2171b5"], # a1-0a1
    },
    "mt": {
        0: [r"$\tau_{\mu}\tau_h\ \mathrm{combined}$", "combined", "#DE5A6A"],
        1: ["embed", "embed", ""], # embed
        2: ["fakes", "fakes", ""], # fakes
        3: [r'$\mu\rho$', "mu-rho", "#9B98CC"], # mu-rho
        4: [r'$\mu\pi$', "mu-pi", "#E8AD46"], # mu-pi
        5: [r'$\mu a_{1}^{3\mathrm{pr}}$', "mu-a1", "#addd8e"], # mu-a1
        6: [r'$\mu a_{1}^{1\mathrm{pr}}$', "mu-0a1", "#c994c7"], # mu-0a1
    },
}
