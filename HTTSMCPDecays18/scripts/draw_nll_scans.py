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
        "--y-scale linear "
    )
    parser = argparse.ArgumentParser(epilog=epilog)

    parser.add_argument(
        "--input-folder",
        default="output/01042020/cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root",
        help="Path to top directory of output folder",
    )
    parser.add_argument(
        "--channel",
        default="tt", choices=["tt", "mt"],
        help="Which channel to use",
    )
    parser.add_argument(
        "--mode",
        default="single",
        choices=["single", "split_by_category"],
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

    return parser.parse_args()

def custom_cms_label(ax, label, lumi=35.9, energy=13):
    ax.text(
        0, 1, r'$\mathbf{CMS}\ \mathit{'+label+'}$',
        ha='left', va='bottom', transform=ax.transAxes,
    )
    ax.text(
        1, 1, r'${:.0f}\ \mathrm{{fb}}^{{-1}}$ ({:.0f} TeV)'.format(lumi, energy),
        ha='right', va='bottom', transform=ax.transAxes,
    )

def prepare_scan(scan):
    """
    Helper function to read in MultiDimFit ROOT file and convert to
    x- and y-values to use for plotting the NLL scan


    Paramters
    =========
    scan: str
        Path to MultiDimFit ROOT output file


    Return
    ======
    xs: pd.Series
        'alpha' values from scan ROOT file

    ys: pd.Series
        -2deltaLL values constructed scan ROOT file
    """

    # Load file with uproot and get 'limit' branch
    # as is done for CH MultiDimFit outputs
    f = uproot.open(scan)["limit"]

    # Get the 'alpha', 'deltaNLL', 'quantileExpected' TBranches
    df = f.pandas.df(["alpha","deltaNLL","quantileExpected"],
        namedecode="utf-8")

    # Some preselection as is usually done
    df = df.query("quantileExpected > -0.5")
    df = df.loc[~df.duplicated(),:]
    df = df.sort_values(by="alpha")

    xs = df["alpha"]
    ys = 2*df["deltaNLL"]

    return xs, ys


def single_scan(input_folder, channel, plot_name, add_significance=False):
    """
    Function to plot NLL scan using ROOT output file from MultiDimFit


    Paramters
    =========
    input_folder: str
        Path to top of output directory within which ROOT file output 
        from MultiDimFit is stored

    channel: str
        Channel to plot for

    plot_name: str
        Name of plot to be saved as pdf
    """

    # Plot single scan (for combined scan for instance or any other)
    with mpl.backends.backend_pdf.PdfPages(f"plots/{plot_name}.pdf", keep_empty=False,) as pdf:
        fig, ax = plt.subplots(
            figsize=(3.2, 2.8), dpi=200,
        )

        path = f"{input_folder}/cmb/125/higgsCombine.alpha.MultiDimFit.mH125.root"
        xvalues, yvalues = prepare_scan(path)

        # Helper function from dftools to return DataFrame with spline fit
        results = dftools.draw.nllscan(
            xvalues, yvalues, ax=ax, nsigs=[1,2], 
            left_bracket=(-90,0), right_bracket=(0,90),
        )
        
        # Aesthetics and lines 1sigma, 2sigma lines
        custom_cms_label(ax, "Preliminary", lumi=137)
        ax.set_xticks([-90, -45, 0, 45, 90])
        ax.set_xlim(-90., 90)
        ax.set_ylim(0., None)
        ax.text(-85, 1.01, r'$1\sigma$', ha='left', va='bottom', color='gray')
        ax.text(-85, 4.01, r'$2\sigma$', ha='left', va='bottom', color='gray')

        bestfit = results.query("nsig == 0.")["xval"].values
        result = np.abs(results.query("abs(nsig) == 1.")["xval"].values)
        hi_string = f"{result[1]:.0f}"
        lo_string = f"{result[0]:.0f}"
        if hi_string == lo_string:
            full_string = f"{bestfit[0]:.0f} \\pm {hi_string}\\ {{}}^{{\circ}}$"
        else:
            full_string = f"{bestfit[0]:.0f}_{{{lo_string}}}^{{+{hi_string}}} {{}}^{{\circ}}$"

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
            figsize=(3.2, 2.8), dpi=200,
        )

        if channel == "tt":
            bins = [0] + list(np.arange(3, 12, 1))
        elif channel == "mt":
            bins = [0] + list(np.arange(3, 7, 1))
        labels = []
        for (category, bin_) in zip([f"htt_{channel}"] + [f"htt_{channel}_{x}_13TeV" for x in list(np.arange(3, 12, 1))], bins):

            # Load scans in category loop
            path = f"{input_folder}/{category}/125/higgsCombine.alpha.MultiDimFit.mH125.root"
            xvalues, yvalues = prepare_scan(path)

            # Draw 1sigma line and label only for combined channel scan
            if bin_ == 0:
                nsigmas = [1]
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

def main(input_folder, channel, mode, plot_name, y_scale, add_significance):
    if mode == "single":
        single_scan(input_folder, channel, plot_name, add_significance)
    elif mode == "split_by_category":
        split_by_category_scan(input_folder, channel, plot_name, y_scale)

if __name__ == "__main__":
    main(**vars(parse_arguments()))
