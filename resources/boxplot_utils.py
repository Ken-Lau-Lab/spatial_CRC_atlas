import os
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns; sns.set_style("white")

from math import ceil
from matplotlib.markers import TICKDOWN
from scipy import stats


def significance_bar(
    start,
    end,
    height,
    displaystring,
    linewidth=1.2,
    markersize=8,
    boxpad=0.3,
    fontsize=14,
    color="k",
):
    """
    Draw significance bracket on matplotlib figure
    """
    # draw a line with downticks at the ends
    plt.plot(
        [start, end],
        [height] * 2,
        "-",
        color=color,
        lw=linewidth,
        marker=TICKDOWN,
        markeredgewidth=linewidth,
        markersize=markersize,
    )
    # draw the text with a bounding box covering up the line
    plt.text(
        0.5 * (start + end),
        height,
        displaystring,
        ha="center",
        va="center",
        bbox=dict(
            facecolor="1.", edgecolor="none", boxstyle="Square,pad=" + str(boxpad)
        ),
        size=fontsize,
    )


def boxplots_group(
    a,
    outdir,
    obs=None,
    colors=None,
    figsize=(3, 3),
    ncols=6,
    sig=True,
    cmap_dict=None,
    titles=None,
):
    """
    Plot trends from  `adata.obs` metadata. Save all plots in
    grid of single `.png` file.

    Parameters
    ----------
    file : str
        The annotated data matrix of shape `n_obs` by `n_vars`. Rows correspond to
        samples and columns to genes.
    outdir : str
        Path to output directory for saving plots
    obs : list of str, optional (default=None)
        Columns from `adata.obs` to plot ImSig scores against
    colors : list of str, optional (default=None)
        Columns from `imsig_results` to plot
    figsize : tuple of float, optional (default=(5,5))
        Size of output figure in inches
    ncols : int, optional (default=5)
        Number of columns in gridspec
    sig : bool, optional (default=True)
        Perform significance testing (2-way t-test) between all groups and add
        significance bars to plot(s)
    cmap_dict : dictionary, optional (default=None)
        Dictionary of group, color pairs from `obs` to color boxes and points by

    Returns
    -------
    Saves plots as `.png` files to `outdir`
    """
    if isinstance(colors, str):
        # coerce single string to list for looping
        colors = [colors]
    n_panels = len(colors) if isinstance(colors, list) else 1
    if isinstance(obs, str):
        # coerce single string to list for looping
        obs = [obs]
    for x in obs:
        print("Saving ImSig boxplots for {}".format(x))
        if n_panels <= ncols:
            n_rows, n_cols = 1, n_panels
        else:
            n_rows, n_cols = ceil(n_panels / ncols), ncols
        fig = plt.figure(figsize=(n_cols * figsize[0], n_rows * figsize[1]))
        left, bottom = 0.1 / n_cols, 0.1 / n_rows
        gs = gridspec.GridSpec(
            nrows=n_rows,
            ncols=n_cols,
            wspace=0.1,
            left=left,
            bottom=bottom,
            right=1 - (n_cols - 1) * left - 0.01 / n_cols,
            top=1 - (n_rows - 1) * bottom - 0.1 / n_rows,
        )
        if sig:
            # initialize dictionary of p values
            sig_out = {
                "signature": [],
                "group1": [],
                "group2": [],
                "pvalue": [],
                "pvalue_adj": [],
            }
        for ic, c in enumerate(colors):
            plt.subplot(gs[ic])
            sns.boxplot(
                data=a.obs,
                x=x,
                y=c,
                hue=x,
                palette=cmap_dict if cmap_dict is not None else None,
                dodge=False,
                saturation=0.4,
                fliersize=0,
            )
            sns.stripplot(
                data=a.obs,
                x=x,
                y=c,
                hue=x,
                palette=cmap_dict if cmap_dict is not None else None,
                jitter=True,
                dodge=False,
                s=7,
                edgecolor="k",
                linewidth=0.5,
                alpha=0.7,
            )
            plt.legend([],[],frameon=False)
            if sig:
                sig_count = 0  # initiate significant count for bar height
                for i_sig in range(len(a.obs[x].unique())):
                    for i_sig_2 in [
                        x
                        for x in range(i_sig, len(a.obs[x].unique()))
                        if x != i_sig
                    ]:
                        # if label has more than two classes, automatically perform t-tests
                        # and add significance bars to plots
                        _, p_value = stats.ttest_ind(
                            a.obs.loc[a.obs[x] == a.obs[x].unique()[i_sig], c].dropna(),
                            a.obs.loc[a.obs[x] == a.obs[x].unique()[i_sig_2], c].dropna(),
                        )
                        # Bonferroni correction
                        pvalue_adj = p_value * len(a.obs[x].unique())
                        # dump results into dictionary
                        sig_out["signature"].append(c)
                        sig_out["group1"].append(a.obs[x].unique()[i_sig])
                        sig_out["group2"].append(a.obs[x].unique()[i_sig_2])
                        sig_out["pvalue"].append(p_value)
                        sig_out["pvalue_adj"].append(pvalue_adj)
                        if pvalue_adj <= 0.05:
                            sig_count += 1  # increment significant count
                            if pvalue_adj < 0.0001:
                                displaystring = r"***"
                            elif pvalue_adj < 0.001:
                                displaystring = r"**"
                            else:
                                displaystring = r"*"
                            # offset by 10 percent for each significant pair
                            height = (
                                a.obs[c].max() + 0.1 * a.obs[c].max() * sig_count
                            )
                            # set up significance bar
                            bar_centers = np.array([i_sig, i_sig_2])
                            significance_bar(
                                bar_centers[0],
                                bar_centers[1],
                                height,
                                displaystring,
                            )
            if titles is not None:
                plt.title(titles[ic])
            plt.xticks(rotation=90)
            plt.xlabel("")
            plt.ylabel(c)
        gs.tight_layout(fig)
        if sig:
            plt.savefig(
                os.path.abspath(os.path.join(outdir, "{}_boxplot_trends_sig.png".format(x))),
            )
            # save statistics to .csv file
            pd.DataFrame(sig_out).to_csv(
                os.path.abspath(
                    os.path.join(outdir, "{}_boxplot_trends_pvals.csv".format(x))
                ),
                index=False,
            )
        else:
            plt.savefig(
                os.path.abspath(os.path.join(outdir, "{}_boxplot_trends.png".format(x))),
            )
