import oyaml as yaml
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import dftools
import scipy
import argparse

from plotting import (
    draw_1d, 
    create_df,
    var_kw,
    process_kw,
    nllscan_kw,
)
mpl.use('pdf')
plt.style.use('cms')

def parse_arguments():

    epilog = (
        "Examples:\n "
        "python3 scripts/draw_nll_scans.py --input-folder output/01042020 "
        "--channel tt --mode single --plot-name alpha_cmb \n "

        "python3 scripts/draw_nll_scans.py --input-folder output/01042020 "
        "--channel tt --mode split_by_category --plot-name alpha_tt_split "
        "--y-scale linear \n "
        
        "python3 scripts/draw_nll_scans.py --input-folder output/01042020 "
        "--mode 2d_kappa \n "
    )
    parser = argparse.ArgumentParser(epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--input-folder",
        default="output/01042020/cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root",
        help="Path to top directory of output folder",
    )
    parser.add_argument(
        "--channel",
        default="tt", choices=["tt", "mt"],
        help="Which channel to use, only required for combined scans",
    )
    parser.add_argument(
        "--cat",
        default="cmb",
        help="Which sub-category to use",
    )
    parser.add_argument(
        "--nsigmas",
        default=[],
        help="How many sigmas to annotate on plot",
    )
    parser.add_argument(
        "--mode",
        default="single",
        choices=["single", "split_by_category", "2d_kappa", "muV", "muggH", "alpha", "alpha_obs", "mutautau"],
        help="What scans to plot",
    )
    parser.add_argument(
        "--y-scale",
        default="linear", choices=["linear", "sqrt"],
        help="What type of scale to use on y-axis",
    )
    parser.add_argument(
        "--plot-name",
        default="nllscan",
        help="Name of figure to save",
    )
    parser.add_argument(
        "--add-significance", action="store_true",
        default=False,
        help="If set, add significance label between bestfit (SM for Asimov) and maximum (PS for Asimov)",
    )
    parser.add_argument(
        "--observed", action="store_true",
        default=False,
        help="If set, add observed parameter scan",
    )

    return parser.parse_args()

def custom_cms_label(ax, label, lumi=35.9, energy=13):
    ax.text(
        0, 1, r'$\mathbf{CMS}\ \mathit{'+label+'}$',
        ha='left', va='bottom', transform=ax.transAxes,
    )
    if lumi != 137:
        ax.text(
            1, 1, r'${:.1f}\ \mathrm{{fb}}^{{-1}}$ ({:.0f} TeV)'.format(lumi, energy),
            ha='right', va='bottom', transform=ax.transAxes,
        )
    else:
        ax.text(
            1, 1, r'${:.0f}\ \mathrm{{fb}}^{{-1}}$ ({:.0f} TeV)'.format(lumi, energy),
            ha='right', va='bottom', transform=ax.transAxes,
        )

def prepare_scan(scan, parameter="alpha",):
    """
    Helper function to read in MultiDimFit ROOT file and convert to
    x- and y-values to use for plotting the NLL scan


    Paramters
    =========
    scan: str
        Path to MultiDimFit ROOT output file

    parameter: str
        Which parameter to scan. Default is 'alpha'


    Return
    ======
    xs: pd.Series
        parameter values from scan ROOT file

    ys: pd.Series
        -2deltaLL values constructed scan ROOT file
    """

    # Load file with uproot and get 'limit' branch
    # as is done for CH MultiDimFit outputs
    f = uproot.open(scan)["limit"]

    # Get the parameter, 'deltaNLL', 'quantileExpected' TBranches
    df = f.pandas.df([parameter,"deltaNLL","quantileExpected"],
        namedecode="utf-8")

    # Some preselection as is usually done
    df = df.query("quantileExpected > -0.5")
    df = df.loc[~df.duplicated(),:]
    df = df.sort_values(by=parameter)

    xs = df[parameter]
    ys = 2*df["deltaNLL"]

    return xs, ys


def single_scan(input_folder, cat, nsigmas, plot_name, add_significance=False):
    """
    Function to plot NLL scan using ROOT output file from MultiDimFit


    Paramters
    =========
    input_folder: str
        Path to top of output directory within which ROOT file output 
        from MultiDimFit is stored

    cat: str
        Sub-category to plot. Usually 'cmb' for full

    nsigmas: str
        How many sigmas to annotate

    plot_name: str
        Name of plot to be saved as pdf
    """
    nsigs = []
    if nsigmas == 1:
        nsigs = [1]
    elif nsigmas == 2:
        nsigs = [1, 2]

    # Plot single scan (for combined scan for instance or any other)
    with mpl.backends.backend_pdf.PdfPages(f"plots/{plot_name}.pdf", keep_empty=False,) as pdf:
        fig, ax = plt.subplots(
            figsize=(4, 3), dpi=200,
        )

        path = f"{input_folder}/{cat}/125/higgsCombine.alpha.MultiDimFit.mH125.root"
        xvalues, yvalues = prepare_scan(path)

        # Helper function from dftools to return DataFrame with spline fit
        results = dftools.draw.nllscan(
            xvalues, yvalues, ax=ax, nsigs=nsigs,
            left_bracket=(-90,0), right_bracket=(0,90),
        )
        
        # Aesthetics and lines 1sigma, 2sigma lines
        custom_cms_label(ax, "Preliminary", lumi=59.7)
        ax.set_xticks([-90, -45, 0, 45, 90])
        ax.set_xlim(-90., 90)
        ax.set_ylim(0., None)

        # Add label of sigma here
        if len(nsigs) >= 1:
            ax.text(
                -85, 1.01, r'$1\sigma$', 
                ha='left', va='bottom', color='gray'
            )
            if len(nsigs) >= 2:
                ax.text(
                    -85, 4.01, r'$2\sigma$', 
                    ha='left', va='bottom', color='gray'
                )

        bestfit = results.query("nsig == 0.")["xval"].values
        result = np.abs(results.query("abs(nsig) == 1.")["xval"].values)
        hi_string = f"{result[1]:.0f}"
        lo_string = f"{result[0]:.0f}"
        if hi_string == lo_string:
            full_string = f"{bestfit[0]:.0f} \\pm {hi_string}\\ {{}}^{{\circ}}$"
        else:
            full_string = f"{bestfit[0]:.0f}_{{-{lo_string}}}^{{+{hi_string}}} {{}}^{{\circ}}$"

        result_label = (
            r"$\hat{\phi}_{\tau} = " + full_string
        )
        ax.text(
            0.5, 0.8, 
            result_label,
            ha='center', va='bottom',
            transform=ax.transAxes, 
        )

        if add_significance:
            # Add significance of SM vs PS discrimination
            significance = np.sqrt(
                scipy.stats.chi2.ppf(
                    scipy.stats.chi2.cdf(
                        max(yvalues)-bestfit,
                    1),
                1)
            )[0]
            sig_label = r"$0^+ \mathrm{vs}\ 0^- =$\ "+f"{significance:.2f}$\sigma$"
            print(f"SM vs PS significance is {sig_label}")
            ax.text(
                0.5, 0.7, 
                sig_label,
                ha='center', va='bottom',
                transform=ax.transAxes, 
            )

        ax.set_xlabel(r"$\phi_{\tau} (\mathrm{degrees})$")
        ax.set_ylabel(r"$-2\Delta\log\mathcal{L}$")

        print(f"Saving figure as {plot_name}.pdf")
        pdf.savefig(fig, bbox_inches='tight')


def prepare_results(ax, results, pc_level, parameter, pos, observed):
    ''' Helper function '''

    nsigs = []
    for pc in pc_level:
        nsigs.append(scipy.stats.norm.ppf((1+pc)/2))
    bestfit = results.query("nsig == 0.")["xval"].values
    print(bestfit)
    result_values = []
    label_strings = []
    for nsig in nsigs:
        result_values.append(np.abs(results.query(f"abs(nsig) == {nsig}")["xval"].values))
    for result in result_values:
        print(result)
        unc_hi = result[1] - bestfit[0]
        unc_lo = bestfit[0] + result[0]
        hi_string = f"{unc_hi:.0f}"
        lo_string = f"{unc_lo:.0f}"
        if parameter not in ["alpha"]:
            hi_string = f"{unc_hi:.2f}"
            lo_string = f"{unc_lo:.2f}"
        string = (lo_string, hi_string)
        label_strings.append(string)
    #if hi_string == lo_string:
    #    
    #    full_string = f"{bestfit[0]:.0f} \\pm {hi_string}\\ {{}}^{{\circ}}$"
    #else:
    full_strings = []
    for idx, string in enumerate(label_strings):
        if parameter not in ["alpha"]:
            if abs(float(string[0])) == abs(float(string[1])):
               full_string = f"{bestfit[0]:.2f} \\pm {abs(float(string[0]))}\\ "
            else:
               full_string = (f"{bestfit[0]:.2f}_{{-{string[0]}}}^{{+{string[1]}}}\\ ")
        else:
            if abs(float(string[0])) == abs(float(string[1])):
                full_string = f"{bestfit[0]:.0f} \\pm {abs(float(string[0])):.0f}\\ "
            else:
               full_string = (f"{bestfit[0]:.0f}_{{{string[0]}}}^{{+{string[1]}}}\\ ")
        # full_string = (f"{bestfit[0]:.2f}_{{-{string[0]}}}^{{+{string[1]}}}\\ ")
        # if parameter not in ["alpha"]:
        #     full_string = (f"{bestfit[0]:.2f}_{{-{string[0]}}}^{{+{string[1]}}}\\ ")
        full_strings.append(full_string)
        
        if parameter == "alpha":
            result_label = (
                r"$\hat{\phi}^{\mathrm{prefit\ exp.}}_{\tau\tau} = " + full_string + "{{}}^{{\circ}}"
            )
            if observed:
                result_label = (
                    r"$\hat{\phi}^{\mathrm{obs.}}_{\tau\tau} = " + full_string + "{{}}^{{\circ}}"
                )
        elif parameter == "muggH":
            result_label = (
                r"$\hat{\mu}_{gg\mathrm{H}}^{\tau\tau} = " + full_string
            )
        elif parameter == "muV":
            result_label = (
                r"$\hat{\mu}_{\mathrm{V}}^{\tau\tau} = " + full_string
            )
        elif parameter == "mutautau":
            result_label = (
                r"$\hat{\mu}^{\tau\tau} = " + full_string
            )
            
        xpos = pos[0]
        ypos = pos[1]
        if idx == 0:
            ax.text(
                # 0.5, 0.85, 
                xpos, ypos,
                result_label + r"(68\%\ \mathrm{CL})$",
                ha='center', va='bottom',
                transform=ax.transAxes, 
            )
        elif idx == 1:
            ax.text(
                # 0.5, 0.75, 
                xpos, ypos-0.1,
                result_label + r"(95\%\ \mathrm{CL})$",
                ha='center', va='bottom',
                transform=ax.transAxes, 
            )
        
        print(f"{bestfit[0]:.0f}_{{{string[0]}}}^{{+{string[1]}}} {{}}^{{\circ}}$")

def single_parameter_scan(input_folder, parameter, cat, plot_name, observed, add_significance=False):
    ''' New function to do scans of all parameters including observed scans '''

    boundaries = []
    if parameter == "alpha":
        boundaries = [(-90,0), (0,90)]
    elif parameter in ["muggH", "muV", "mutautau"]:
        boundaries = [(0,1), (1,2)]
    else:
        boundaries = [(0,1), (1,2)]
    with mpl.backends.backend_pdf.PdfPages(f"plots/{plot_name}.pdf", keep_empty=False,) as pdf:
        fig, ax = plt.subplots(
            figsize=(4,3), dpi=200,
        )

        path = f"{input_folder}/{cat}/125/higgsCombine.{parameter}.MultiDimFit.mH125.root"
        xvalues, yvalues = prepare_scan(path, parameter)
        if observed:
            path = f"{input_folder}/{cat}/125/higgsCombine.{parameter}.observed.MultiDimFit.mH125.root"
            xvalues_obs, yvalues_obs = prepare_scan(path, parameter)

        scan_kw = dict(color=nllscan_kw["tt"][0][2])
        # either specify sigma or do using percent level
        pc_level = [0.68, 0.95]
        if parameter not in ["alpha"]:
            pc_level = [0.68]
        nsigs = []
        for pc in pc_level:
            nsigs.append(scipy.stats.norm.ppf((1+pc)/2))
        if observed:
            scan_kw = dict(color=nllscan_kw["mt"][3][2])
            scan_kw.update(ls='--')
            line_kw = dict(lw=0)
        results = dftools.draw.nllscan(
            xvalues, yvalues, ax=ax, nsigs=nsigs, 
            left_bracket=boundaries[0], right_bracket=boundaries[1],
            marker_kw=scan_kw, spline_kw=scan_kw, line_kw=line_kw,
        )
        if observed:
            scan_kw = dict(color=nllscan_kw["tt"][0][2])
            results_obs = dftools.draw.nllscan(
                xvalues_obs, yvalues_obs, ax=ax, nsigs=nsigs, 
                left_bracket=boundaries[0], right_bracket=boundaries[1],
                marker_kw=scan_kw, spline_kw=scan_kw,
            )
        custom_cms_label(ax, "Preliminary", lumi=137)
        if parameter == "alpha":
            ax.set_xticks([-90, -45, 0, 45, 90])
            ax.set_xlim(-90., 90)
            ax.set_ylim(0., None)
            
            ax.text(-85, 1.01, r'$68\%$', ha='left', va='bottom', color='gray')
            ax.text(-85, 3.92, r'$95\%$', ha='left', va='bottom', color='gray')
        elif parameter == "mutautau":
            ax.set_xlim(0, 2)
            ax.set_ylim(0., None)
            ax.text(0.05, 1.01, r'$68\%$', ha='left', va='bottom', color='gray')
            ax.text(0.05, 3.92, r'$95\%$', ha='left', va='bottom', color='gray')
            #ax.text(0.55, 1.01, r'$68\%$', ha='left', va='bottom', color='gray')
            #ax.text(0.55, 3.92, r'$95\%$', ha='left', va='bottom', color='gray')
        else:
            ax.set_xlim(0, 2)
            ax.set_ylim(0., None)
            ax.text(0.05, 1.01, r'$68\%$', ha='left', va='bottom', color='gray')
            # ax.text(0.05, 3.92, r'$95\%$', ha='left', va='bottom', color='gray')

        prepare_results(ax, results, pc_level, parameter, pos=[0.5,0.85], observed=False)
        if observed:
            prepare_results(ax, results_obs, pc_level, parameter, pos=[0.5,0.65], observed=True)
        # Add significance of SM vs PS discrimination
        significance = np.sqrt(
            scipy.stats.chi2.ppf(
                scipy.stats.chi2.cdf(
                    [max(yvalues)],
                1),
            1)
        )[0]
        sig_label = r"$0^+ \mathrm{vs}\ 0^- =$\ "+f"{significance:.2f}$\sigma$"
        print(f"SM vs PS significance is {sig_label}")
        if observed:
            significance = np.sqrt(
                scipy.stats.chi2.ppf(
                    scipy.stats.chi2.cdf(
                        [max(yvalues_obs)],
                    1),
                1)
            )[0]
            sig_label = r"$0^+ \mathrm{vs}\ 0^- =$\ "+f"{significance:.2f}$\sigma$"
            print(f"SM vs PS significance is {sig_label}")

        if add_significance:
            ax.text(
                0.5, 0.65,
                sig_label,
                ha='center', va='bottom',
                transform=ax.transAxes,
            )
    
        if parameter == "alpha":
            ax.set_xlabel(r"$\phi_{\tau\tau} (\mathrm{degrees})$")
        elif parameter == "muggH":
            ax.set_xlabel(r"$\mu_{gg\mathrm{H}}^{\tau\tau}$")
        elif parameter == "muV":
            ax.set_xlabel(r"$\mu_{\mathrm{V}}^{\tau\tau}$")
            ax.set_ylim(0, 3.2)
        elif parameter == "mutautau":
            ax.set_xlabel(r"$\mu^{\tau\tau}$")
            #ax.set_ylim(0, 6.2)
            #ax.set_xlim(0.5, 1.5)
        #ax.set_ylim(0, 6)
        #ax.set_xlim(-2, 2)
        ax.set_ylabel(r"$-2\Delta\log\mathcal{L}$")
        pdf.savefig(fig, bbox_inches='tight')

def split_by_category_scan(input_folder, channel, plot_name, y_scale="linear"):
    """
    Function to plot NLL scan using multiple ROOT output file from MultiDimFit.
    The combined channel scan will be plotted, included all the sub-categories
    of that channel.


    Paramters
    =========
    input_folder: str
        Path to top of output directory within which ROOT file output 
        from MultiDimFit is stored

    channel: str
        Channel to plot for

    plot_name: str
        Name of plot to be saved as pdf

    y_scale: str
        The type of scale to use on y-axis. Can be useful for plotting multiple
        scans of minor significance.
    """

    with mpl.backends.backend_pdf.PdfPages(f"plots/{plot_name}.pdf", keep_empty=False,) as pdf:
        fig, ax = plt.subplots(
            figsize=(4, 3), dpi=200,
        )

        if channel == "tt":
            bins = [0] + [3,7,5,9,8,4,10,6,11] # this is ordered by hand
        elif channel == "mt":
            bins = [0] + list(np.arange(3, 7, 1))
        labels = []
        for (category, bin_) in zip([f"htt_{channel}"] + [f"htt_{channel}_{x}_13TeV" for x in list(np.arange(3, 12, 1))], bins):

            # Load scans in category loop
            path = f"{input_folder}/{category}/125/higgsCombine.alpha.MultiDimFit.mH125.root"
            xvalues, yvalues = prepare_scan(path)

            # Draw 1sigma line and label only for combined channel scan
            if bin_ == 0:
                nsigmas = [1,2]
            else:
                nsigmas = []

            scan_kw = dict(color=nllscan_kw[channel][bin_][2],)

            results = dftools.draw.nllscan(
                xvalues, yvalues, ax=ax, nsigs=nsigmas, 
                left_bracket=(-90,0), right_bracket=(0,90),
                marker_kw=scan_kw, spline_kw=scan_kw,
            )
            labels.append(nllscan_kw[channel][bin_][0])

            # Add label of sigma here
            if len(nsigmas) >= 1:
                ax.text(
                    -85, 1.01, r'$1\sigma$', 
                    ha='left', va='bottom', color='gray'
                )
                if len(nsigmas) >= 2:
                    ax.text(
                        -85, 4.01, r'$2\sigma$', 
                        ha='left', va='bottom', color='gray'
                    )
                bestfit = results.query("nsig == 0.")["xval"].values

        custom_cms_label(ax, "Preliminary", lumi=137)
        ax.set_xticks([-90, -45, 0, 45, 90])
        ax.set_xlim(-90., 90)
        ax.set_ylim(0., None)


        full_handles, full_labels = ax.get_legend_handles_labels()
        handles = [
            x for idx, x in enumerate(full_handles) 
            if "Spline" in full_labels[idx]
        ]
        ax.legend(
            handles, labels, bbox_to_anchor=(1, 1.03), 
            labelspacing=0.3, borderpad=0.2
        )
        ax.set_xlabel(r"$\phi_{\tau} (\mathrm{degrees})$")

        if y_scale == "sqrt":
            ax.set_yscale(
                'function', 
                functions=(lambda x: np.maximum(x, 0)**0.5, lambda x: x**2)
            )
        ax.set_ylabel(r"$-2\Delta\log\mathcal{L}$")

        print(f"Saving figure as {plot_name}.pdf")
        pdf.savefig(fig, bbox_inches='tight')

def split_by_year_scan(input_folder, plot_name, y_scale="linear"):
    """
    Function to plot NLL scan using multiple ROOT output file from MultiDimFit.
    The NLL scan for each year will be shown


    Paramters
    =========
    input_folder: str
        Path to top of output directory within which ROOT file output 
        from MultiDimFit is stored

    plot_name: str
        Name of plot to be saved as pdf

    y_scale: str
        The type of scale to use on y-axis. Can be useful for plotting multiple
        scans of minor significance.
    """
    channel = "years"
    with mpl.backends.backend_pdf.PdfPages(f"plots/{plot_name}.pdf", keep_empty=False,) as pdf:
        fig, ax = plt.subplots(
            figsize=(4, 3), dpi=200,

        )
        labels = []
        significances = []
        for idx, category in enumerate(["cmb","htt_2016", "htt_2017", "htt_2018"]):
            # Load scans in category loop
            path = f"{input_folder}/{category}/125/higgsCombine.alpha.MultiDimFit.mH125.root"
            xvalues, yvalues = prepare_scan(path)

            if idx == 0:
                nsigmas = [scipy.stats.norm.ppf((1+0.68)/2),scipy.stats.norm.ppf((1+0.95)/2)]
            else:
                nsigmas = []

            scan_kw = dict(color=nllscan_kw[channel][idx][2],)

            results = dftools.draw.nllscan(
                xvalues, yvalues, ax=ax, nsigs=nsigmas, 
                left_bracket=(-90,0), right_bracket=(0,90),
                marker_kw=scan_kw, spline_kw=scan_kw,
            )
            labels.append(nllscan_kw[channel][idx][0])
            
            # print significance of each category
            significance = np.sqrt(
                scipy.stats.chi2.ppf(
                    scipy.stats.chi2.cdf(
                        max(yvalues),
                    1),
                1)
            )
            sig_label = r"$0^+ \mathrm{vs}\ 0^- =$\ "+f"{significance:.2f}$\sigma$"
            print(f"SM vs PS significance for year {nllscan_kw[channel][idx][1]} is {sig_label}")
            significances.append(sig_label)
        custom_cms_label(ax, "Preliminary", lumi=137)
        ax.set_xticks([-90, -45, 0, 45, 90])
        ax.set_xlim(-90., 90)
        ax.set_ylim(0., None)

        full_handles, full_labels = ax.get_legend_handles_labels()
        handles = [x for idx, x in enumerate(full_handles) if "Spline" in full_labels[idx]]
        ax.legend(handles, labels, loc=9, labelspacing=0.3, borderpad=0.2, framealpha=0.7)
        ax.set_xlabel(r"$\phi_{\tau\tau} (\mathrm{degrees})$")

        ax.set_ylabel(r"$-2\Delta\log\mathcal{L}$")

        pdf.savefig(fig, bbox_inches='tight')
        pass

def scan_2d_kappa(input_folder, category="cmb", plot_name="scan_2d_kappa",):
    """
    Function to plot NLL scan using multiple ROOT output file from MultiDimFit.
    This is specifically for 2D scans of kappas (ie. reduced Yukawa couplings)


    Paramters
    =========
    input_folder: str
        Path to top of output directory within which ROOT file output 
        from MultiDimFit is stored

    plot_name: str
        Name of plot to be saved as pdf

    category: str
        Category name as in CH, usually 'cmb' for these kind of scans
    """

    with mpl.backends.backend_pdf.PdfPages(f"plots/{plot_name}.pdf", keep_empty=False,) as pdf:
        fig, ax = plt.subplots(
            figsize=(4, 3), dpi=200,
        )
        path = f"{input_folder}/{category}/125/higgsCombine.kappas.MultiDimFit.mH125.root"
        parameter0 = "kappaH"
        parameter1 = "kappaA"
        f = uproot.open(path)["limit"]
        df = f.pandas.df([parameter0, parameter1, "deltaNLL","quantileExpected"],
            namedecode="utf-8")
        df = df.query("quantileExpected > -0.5 and deltaNLL < 1000 ")
        df = df.loc[~df.duplicated(),:]
        df = df.sort_values(by=[parameter0, parameter1])
        custom_cms_label(ax, "Preliminary", lumi=137)
        
        xbins = df[parameter0].unique()
        ybins = df[parameter1].unique()
        df["deltaNLL"] = 2*df["deltaNLL"]
        # print(df)
        z = df.set_index([parameter0, parameter1])["deltaNLL"].unstack().values.T
        # some nans...remove by setting to high value (high NLL)
        # this is only a temp. fix, hopefully fix to Physics model will remove these
        z[np.isnan(z)] = 300
        # print(z)
        
        pos = ax.imshow(
            z, origin='lower', interpolation='bicubic',
            extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]],
            aspect='auto', cmap="Blues_r",
            vmin=0., vmax=25.,
        )
        cbar = fig.colorbar(pos, ax=ax)
        cbar.set_label(r'$-2\Delta\log \mathcal{L}$')
        
        X, Y = np.meshgrid(xbins, ybins)
        ax.contour(
            scipy.ndimage.zoom(X, 4),
            scipy.ndimage.zoom(Y, 4),
            scipy.ndimage.zoom(z, 4),
            #z,
            levels=[scipy.stats.chi2.ppf(0.68, df=2)],
            colors=['black'],
        )
        ax.contour(
            scipy.ndimage.zoom(X, 4),
            scipy.ndimage.zoom(Y, 4),
            scipy.ndimage.zoom(z, 4),
            levels=[scipy.stats.chi2.ppf(0.95, df=2)],
            colors=['black'], linestyles='dashed',
        )
        bf = (
            df.loc[df["deltaNLL"]==df["deltaNLL"].min(), parameter0],
            df.loc[df["deltaNLL"]==df["deltaNLL"].min(), parameter1],
        )
        ax.plot(
            *bf, 'P', color='black',
            ms=5, label="Best fit",
        )
        ax.plot(
            1, 0, '*', color='#e31a1c',
            ms=4, label="SM",
        )
        handles, labels = ax.get_legend_handles_labels()
        handles = handles[::-1] + [
            mpl.lines.Line2D([0], [0], color='black', lw=1),
            mpl.lines.Line2D([0], [0], color='black', lw=1, ls='--'),
        ]
        labels = labels[::-1] + [r'$68\%$ CI', r'$95\%$ CI']
        ax.legend(
            handles, labels,
            loc=3, labelspacing=0.1, borderpad=0.2,
            fancybox=True, edgecolor='#d9d9d9',
            framealpha=0., handlelength=1.,
        )

        ax.text(
            0.75, 0.05, r"$\mu_{gg\mathrm{H}}^{\tau\tau} = \mu_{V}^{\tau\tau} = 1$",
            ha='center', va='bottom', transform=ax.transAxes,
        )

        ax.set_xlabel(r'$\kappa_{\tau}$')
        ax.set_ylabel(r'$\tilde{\kappa}_{\tau}$')
        # ax.set_ylim(-1.5, 1.5)
        # ax.set_xlim(-1.5, 1.5)
        print(f"Saving figure as plots/{plot_name}.pdf")
        pdf.savefig(fig, bbox_inches='tight')

def main(input_folder, channel, cat, nsigmas, mode, plot_name, y_scale, observed, add_significance):
    if mode == "single":
        single_scan(input_folder, cat, nsigmas, plot_name, add_significance)
    elif mode == "split_by_category":
        split_by_category_scan(input_folder, channel, plot_name, y_scale)
    elif mode == "split_by_year":
        split_by_year_scan(input_folder, plot_name, y_scale)
    elif mode == "mutautau":
        single_parameter_scan(input_folder, "mutautau", cat, plot_name, observed, add_significance)
    elif mode == "muV":
        single_parameter_scan(input_folder, "muV", cat, plot_name, observed, add_significance)
    elif mode == "muggH":
        single_parameter_scan(input_folder, "muggH", cat, plot_name, observed, add_significance)
    elif mode == "alpha":
        single_parameter_scan(input_folder, "alpha", cat, plot_name, observed, add_significance)
    elif mode == "2d_kappa":
        scan_2d_kappa(input_folder, cat, plot_name)

if __name__ == "__main__":
    main(**vars(parse_arguments()))
