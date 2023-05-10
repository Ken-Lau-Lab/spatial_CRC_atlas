# -*- coding: utf-8 -*-
"""
Custom plotting wrapper functions

@author: C Heiser
"""
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from scipy.cluster.hierarchy import linkage, leaves_list


def cody_heatmap(
    adata,
    groupby,
    features,
    layer=None,
    cluster_vars=False,
    vars_dict=None,
    groupby_order=None,
    groupby_colordict=None,
    cluster_obs=False,
    cmap="Greys",
    figsize=(5, 5),
    save=None,
    dpi=400,
    **kwargs,
):
    """
    Custom wrapper around `sc.pl.dotplot`

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object to plot from
    groupby : str
        Categorical column of `adata.obs` to group dotplot by
    features : list of str
        List of features from `adata.obs.columns` or `adata.var_names` to plot
    layer : str
        Key from `adata.layers` to use for plotting gene values
    cluster_vars : bool, optional (default=`False`)
        Hierarchically cluster `features` for a prettier dotplot. If `True`, return
        `features` in their new order (first return variable).
    vars_dict : dict, optional (default=`None`)
        Dictionary of groups of vars to highlight with brackets on dotplot. Keys are
        variable group names and values are list of variables found in `features`.
    groupby_order : list, optional (default=`None`)
        Explicit order for groups of observations from `adata.obs[groupby]`
    groupby_colordict : dict, optional (default=`None`)
        Dictionary mapping `groupby` categories (keys) to colors (values). Black
        outline will be added to provide contrast to light colors.
    cluster_obs : bool, optional (default=`False`)
        Hierarchically cluster `groupby` observations and show dendrogram
    cmap : str, optional (default="Greys")
        matplotlib colormap for dots
    figsize : tuple of float, optional (default=(5,5))
        Size of the figure in inches
    save : str or `None`, optional (default=`None`)
        Path to file to save image to. If `None`, return figure object (second return
        variable if `cluster_vars == True`).
    dpi : float, optional (default=400)
        Resolution in dots per inch for saving figure. Ignored if `save` is `None`.
    **kwargs
        Keyword arguments to pass to `sc.pl.dotplot`

    Returns
    -------
    features_ordered : list of str
        If `cluster_vars == True`, reordered `features` based on hierarchical
        clustering
    figure : sc._plotting.DotPlot
        If `save == None`, scanpy plotting object is returned
    """
    if np.all([x in adata.obs.columns for x in features]):
        same_origin = True
        print("Using {} features from adata.obs".format(len(features)))
        a_comb_sig = sc.AnnData(
            adata.obs[features].values,
            obs=adata.obs[[x for x in adata.obs.columns if x not in features]],
        )
        a_comb_sig.X = np.nan_to_num(a_comb_sig.X, 0)
        a_comb_sig.layers["raw_counts"] = a_comb_sig.X.copy()
        a_comb_sig.var_names = features
    elif np.all([x in adata.var_names for x in features]):
        same_origin = True
        print("Using {} features from adata.X".format(len(features)))
        a_comb_sig = adata[:, features].copy()
    else:
        same_origin = False
        print("Using {} features from adata.X and adata.obs".format(len(features)))
        a_comb_sig = adata[
            :, list(set(features).intersection(set(adata.var_names)))
        ].copy()
    if cluster_vars and vars_dict is None:
        assert (
            same_origin
        ), "In order to hierarchically cluster features, they must all reside in .X or .obs"
        print("Hierarchically clustering features")
        # first get hierchically-clustered indices of variables
        link = linkage(a_comb_sig.X.T)
        leaves = leaves_list(link)
        myplot = sc.pl.dotplot(
            a_comb_sig,
            a_comb_sig.var_names[leaves],  # use indices from clustermap
            groupby=groupby,
            dendrogram=cluster_obs if groupby_order is None else False,
            categories_order=groupby_order if groupby_order is not None else None,
            layer=layer,
            standard_scale="var",
            swap_axes=True,
            var_group_rotation=90,
            return_fig=True,
            figsize=figsize,
            **kwargs,
        )
    elif vars_dict is not None:
        print("Using vars_dict for ordering features")
        myplot = sc.pl.dotplot(
            a_comb_sig,
            vars_dict,
            groupby=groupby,
            dendrogram=cluster_obs if groupby_order is None else False,
            categories_order=groupby_order if groupby_order is not None else None,
            layer=layer,
            standard_scale="var",
            swap_axes=True,
            var_group_rotation=90,
            return_fig=True,
            figsize=figsize,
            **kwargs,
        )
    else:
        print("Features ordered as given")
        myplot = sc.pl.dotplot(
            a_comb_sig,
            features,
            groupby=groupby,
            dendrogram=cluster_obs if groupby_order is None else False,
            categories_order=groupby_order if groupby_order is not None else None,
            layer=layer,
            standard_scale="var",
            swap_axes=True,
            var_group_rotation=90,
            return_fig=True,
            figsize=figsize,
            **kwargs,
        )
    # style options to plot
    myplot.style(cmap=cmap, dot_edge_color="k", dot_edge_lw=1)
    myplot.ax_dict = myplot.get_axes()
    if groupby_colordict:
        myplot.ax_dict["mainplot_ax"].set_xticklabels(
            myplot.ax_dict["mainplot_ax"].get_xticklabels(),
            path_effects=[pe.withStroke(linewidth=0.2, foreground="k")],
        )
        [
            t.set_color(i)
            for (i, t) in zip(
                [
                    groupby_colordict[x.get_text()]
                    for x in myplot.ax_dict["mainplot_ax"].get_xticklabels()
                ],
                myplot.ax_dict["mainplot_ax"].xaxis.get_ticklabels(),
            )
        ]
    if save is None:
        if cluster_vars:
            return a_comb_sig.var_names[leaves], myplot
        else:
            return myplot
    else:
        print("Saving to {}".format(save))
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        if cluster_vars:
            return a_comb_sig.var_names[leaves]


def cody_violinmap(
    adata,
    groupby,
    features,
    layer=None,
    cluster_vars=False,
    vars_dict=None,
    groupby_order=None,
    groupby_colordict=None,
    cluster_obs=False,
    cmap="Greys",
    figsize=(5, 5),
    save=None,
    dpi=400,
    **kwargs,
):
    """
    Custom wrapper around `sc.pl.dotplot`

    Parameters
    ----------
    adata : sc.AnnData
        AnnData object to plot from
    groupby : str
        Categorical column of `adata.obs` to group dotplot by
    features : list of str
        List of features from `adata.obs.columns` or `adata.var_names` to plot
    layer : str
        Key from `adata.layers` to use for plotting gene values
    cluster_vars : bool, optional (default=`False`)
        Hierarchically cluster `features` for a prettier dotplot. If `True`, return
        `features` in their new order (first return variable).
    vars_dict : dict, optional (default=`None`)
        Dictionary of groups of vars to highlight with brackets on dotplot. Keys are
        variable group names and values are list of variables found in `features`.
    groupby_order : list, optional (default=`None`)
        Explicit order for groups of observations from `adata.obs[groupby]`
    groupby_colordict : dict, optional (default=`None`)
        Dictionary mapping `groupby` categories (keys) to colors (values). Black
        outline will be added to provide contrast to light colors.
    cluster_obs : bool, optional (default=`False`)
        Hierarchically cluster `groupby` observations and show dendrogram
    cmap : str, optional (default="Greys")
        matplotlib colormap for dots
    figsize : tuple of float, optional (default=(5,5))
        Size of the figure in inches
    save : str or `None`, optional (default=`None`)
        Path to file to save image to. If `None`, return figure object (second return
        variable if `cluster_vars == True`).
    dpi : float, optional (default=400)
        Resolution in dots per inch for saving figure. Ignored if `save` is `None`.
    **kwargs
        Keyword arguments to pass to `sc.pl.dotplot`

    Returns
    -------
    features_ordered : list of str
        If `cluster_vars == True`, reordered `features` based on hierarchical
        clustering
    figure : sc._plotting.DotPlot
        If `save == None`, scanpy plotting object is returned
    """
    if np.all([x in adata.obs.columns for x in features]):
        same_origin = True
        print("Using {} features from adata.obs".format(len(features)))
        a_comb_sig = sc.AnnData(
            adata.obs[features].values,
            obs=adata.obs[[x for x in adata.obs.columns if x not in features]],
        )
        a_comb_sig.X = np.nan_to_num(a_comb_sig.X, 0)
        a_comb_sig.layers["raw_counts"] = a_comb_sig.X.copy()
        a_comb_sig.var_names = features
    elif np.all([x in adata.var_names for x in features]):
        same_origin = True
        print("Using {} features from adata.X".format(len(features)))
        a_comb_sig = adata[:, features].copy()
    else:
        same_origin = False
        print("Using {} features from adata.X and adata.obs".format(len(features)))
        a_comb_sig = adata[
            :, list(set(features).intersection(set(adata.var_names)))
        ].copy()
    if cluster_vars and vars_dict is None:
        assert (
            same_origin
        ), "In order to hierarchically cluster features, they must all reside in .X or .obs"
        print("Hierarchically clustering features")
        # first get hierchically-clustered indices of variables
        link = linkage(a_comb_sig.X.T)
        leaves = leaves_list(link)
        myplot = sc.pl.stacked_violin(
            a_comb_sig,
            a_comb_sig.var_names[leaves],  # use indices from clustermap
            groupby=groupby,
            dendrogram=cluster_obs if groupby_order is None else False,
            categories_order=groupby_order if groupby_order is not None else None,
            layer=layer,
            standard_scale="var",
            cmap=cmap,
            linewidth=1,
            swap_axes=True,
            var_group_rotation=90,
            return_fig=True,
            figsize=figsize,
            **kwargs,
        )
    elif vars_dict is not None:
        print("Using vars_dict for ordering features")
        myplot = sc.pl.stacked_violin(
            a_comb_sig,
            vars_dict,
            groupby=groupby,
            dendrogram=cluster_obs if groupby_order is None else False,
            categories_order=groupby_order if groupby_order is not None else None,
            layer=layer,
            standard_scale="var",
            cmap=cmap,
            linewidth=1,
            swap_axes=True,
            var_group_rotation=90,
            return_fig=True,
            figsize=figsize,
            **kwargs,
        )
    else:
        print("Features ordered as given")
        myplot = sc.pl.stacked_violin(
            a_comb_sig,
            features,
            groupby=groupby,
            dendrogram=cluster_obs if groupby_order is None else False,
            categories_order=groupby_order if groupby_order is not None else None,
            layer=layer,
            standard_scale="var",
            cmap=cmap,
            linewidth=1,
            swap_axes=True,
            var_group_rotation=90,
            return_fig=True,
            figsize=figsize,
            **kwargs,
        )
    # style options to plot
    myplot.ax_dict = myplot.get_axes()
    if groupby_colordict:
        myplot.ax_dict["mainplot_ax"].set_xticklabels(
            myplot.ax_dict["mainplot_ax"].get_xticklabels(),
            path_effects=[pe.withStroke(linewidth=0.2, foreground="k")],
        )
        [
            t.set_color(i)
            for (i, t) in zip(
                [
                    groupby_colordict[x.get_text()]
                    for x in myplot.ax_dict["mainplot_ax"].get_xticklabels()
                ],
                myplot.ax_dict["mainplot_ax"].xaxis.get_ticklabels(),
            )
        ]
    if save is None:
        if cluster_vars:
            return a_comb_sig.var_names[leaves], myplot
        else:
            return myplot
    else:
        print("Saving to {}".format(save))
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        if cluster_vars:
            return a_comb_sig.var_names[leaves]
