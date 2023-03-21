# -*- coding: utf-8 -*-
"""
Functions and classes for analyzing 10X Visium spatial transcriptomic and imaging data

@author: C Heiser
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq

sc.set_figure_params(dpi=100, dpi_save=400)
sns.set_style("white")
plt.rcParams["font.family"] = "sans-serif"

from cNMF.cnmf import cnmf_markers
from colormath.color_objects import sRGBColor, XYZColor
from colormath.color_conversions import convert_color
from glob import glob
from matplotlib.markers import TICKDOWN
from mycolorpy import colorlist as mcp
from scipy import stats
from scipy.cluster import hierarchy
from shutil import rmtree
from sklearn.decomposition import non_negative_factorization
from MILWRM import show_pita


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


def read_10x_visium(
    directory,
    count_file="raw_feature_bc_matrix.h5",
    filter_labels=True,
    cleanup=True,
    save_to=None,
    verbose=True,
):
    """
    Read Visium data from 10x Genomics and compile AnnData object

    Parameters:
        directory (str): directory containing Visium files.
            Minimum requirements are 'spatial.tar.gz' and '*_raw_feature_bc_matrix.h5'.
        count_file (str): file within dir that has counts.
            Default '*_raw_feature_bc_matrix.h5'.
        filter_labels (bool): add 'in_tissue' label to adata.obs from
            '*_filtered_gene_bc_matrix.h5'
        cleanup (bool): remove files after reading into AnnData.
        save_to (str or None): name of .h5ad file to write to. If None, don't write.
        verbose (bool): print progress updates to console.

    Returns:
        a (AnnData.anndata): AnnData object with Visium data
    """
    # read visium data into scanpy AnnData object
    if verbose:
        print(
            "Reading Visium data from {}".format(
                glob("{}/*{}".format(directory, count_file))
            )
        )
    a = sc.read_visium(
        path=directory,
        count_file=os.path.basename(glob("{}/*{}".format(directory, count_file))[0]),
    )
    a.var_names_make_unique()
    # add in_tissue labels from filtered file if desired
    if filter_labels:
        if verbose:
            print(
                "Adding 'in_tissue' label from {}".format(
                    glob("{}/*filtered_feature_bc_matrix.h5".format(directory))[0]
                ),
                end="",
            )
        b = sc.read_10x_h5(
            glob("{}/*filtered_feature_bc_matrix.h5".format(directory))[0]
        )
        a.obs["in_tissue"] = "False"
        a.obs.loc[b.obs_names, "in_tissue"] = "True"
        if verbose:
            print(
                " - {} spots within tissue area".format(
                    len(
                        a.obs.loc[
                            a.obs.in_tissue == "True",
                        ]
                    )
                )
            )
    # remove unnecessary files
    if cleanup:
        if verbose:
            print("Cleaning up workspace")
        rmtree("{}/spatial".format(directory))
        for f in glob("{}/*".format(directory)):
            os.remove(f)
    # save AnnData as .h5ad file
    if save_to is not None:
        if verbose:
            print("Saving AnnData to {}/{}".format(directory, save_to))
        a.write("{}/{}".format(directory, save_to), compression="gzip")
    return a


def embed_3D_color(adata, use_rep="X_pca"):
    """
    embed first 3 dimensions of a low-dimensional representation (i.e. PCA)
    into 3D color space (RGB)

    Parameters:
        adata (AnnData.anndata): the data
        use_rep (str): representation from .obsm to embed in RGB space; default "X_pca"

    Returns:
        adata is edited in place, adding .obsm["{use_rep}_RGB"]
    """
    rep = adata.obsm[use_rep].copy()  # make copy of representation matrix
    # normalize to positive 3D cube on [0.0, 1.0]
    rep -= np.min(rep)
    rep /= np.max(rep, axis=0)
    # convert to RGB using colormath
    adata.obsm["{}_RGB".format(use_rep)] = np.array(
        [
            convert_color(
                XYZColor(rep[x, 0], rep[x, 1], rep[x, 2]), sRGBColor
            ).get_value_tuple()
            for x in range(len(rep))
        ]
    )


def upsample_pita(
    adata,
    regressor,
    features=None,
    use_rep=None,
    layer=None,
    plot_out=True,
    **kwargs,
):
    """
    use regression to upsample Visium features to pixel space

    Parameters:
        adata (AnnData.anndata): the data
        regressor (sklearn model): untrained regressor object from sklearn
        features (list of int or str): names or indices of features to cast onto spot
            image if None, cast all features. if plot_out, first feature in list will
            be plotted. if not specified and plot_out, first feature (index 0) will
            be plotted.
        use_rep (str): key from adata.obsm to use for plotting
            if None, use self.adata.X
        layer (str): key from adata.layers to use for plotting
            ignored if use_rep is not None
        plot_out (bool): show resulting image?
        **kwargs: arguments to pass to show_pita() function

    Returns:
        regressor (sklearn model): trained regressor object from sklearn
        X (np.ndarray):
        y (np.ndarray):
        upsampled (np.array): image of desired expression in upsampled pixel space
    """
    # coerce features to list if only single string
    if features and not isinstance(features, list):
        features = [features]

    if use_rep is None:
        # use HVGs if no gene features specified
        if not features:
            features = adata.var_names[adata.var.highly_variable == 1].tolist()
        if layer is None:
            print("Upsampling pita with {} features from adata.X".format(len(features)))
            mapper = pd.DataFrame(
                adata.X[:, [adata.var_names.get_loc(x) for x in features]],
                index=adata.obs_names,
            )
        else:
            print(
                "Upsampling pita with {} features from adata.layers['{}']".format(
                    len(features), layer
                )
            )
            mapper = pd.DataFrame(
                adata.layers[layer][:, [adata.var_names.get_loc(x) for x in features]],
                index=adata.obs_names,
            )
    else:
        if not features:
            print(
                "Upsampling pita with {} features from adata.obsm['{}']".format(
                    adata.obsm[use_rep].shape[1], use_rep
                )
            )
            mapper = pd.DataFrame(adata.obsm[use_rep], index=adata.obs_names)
        else:
            assert all(
                isinstance(x, int) for x in features
            ), "features must be integer indices if using rep from adata.obsm"
            print(
                "Upsampling pita with {} features from adata.obsm['{}']".format(
                    len(features), use_rep
                )
            )
            mapper = pd.DataFrame(
                adata.obsm[use_rep][:, features], index=adata.obs_names
            )

    feats = list(mapper.columns)  # get names of features from training data
    # merge mapper df with pixel_map_df to align observations with RGB channels
    mapper = adata.uns["pixel_map_df"][["barcode", "ch_0", "ch_1", "ch_2"]].merge(
        mapper, left_on="barcode", right_index=True, how="inner"
    )
    y = mapper[feats].values  # gene expression features are output variables
    X = mapper[["ch_0", "ch_1", "ch_2"]].values  # microscopy features are predictors
    regressor.fit(X, y)  # fit model
    y_pred = regressor.predict(X)  # predict at pixel resolution
    # put values into dataframe and merge with x,y coordinates for casting
    y_pred_df = pd.DataFrame(y_pred, index=mapper.index)
    upsampled_df = adata.uns["pixel_map_df"][["x", "y"]].merge(
        y_pred_df, how="outer", left_index=True, right_index=True
    )
    # convert to np.ndarray
    upsampled = upsampled_df.pivot(
        index="y",
        columns="x",
        values=[x for x in upsampled_df.columns if x not in ["x", "y"]],
    ).values
    upsampled = upsampled.reshape(
        (
            upsampled.shape[0],
            int(upsampled.shape[1] / (len(upsampled_df.columns) - 2)),
            (len(upsampled_df.columns) - 2),
        ),
        order="F",
    )

    if plot_out:
        show_pita(pita=upsampled, **kwargs)
    print("Done!")
    return regressor, X, y, upsampled


def deconvolve_cnmf(adata, cnmf_dir, k, raw_layer=None, update_H=False):
    """
    use matrix factorization to deconvolve Visium spot expression

    Parameters:
        adata (AnnData.anndata): the (ST) data
        cnmf_dir (str): directory containing cNMF results from scRNA-seq reference data
        k (int): value of k to pull consensus for
        raw_layer (str): key from adata.layers containing raw counts. if None, use adata.X
        update_H (bool): set to True, both W and H will be estimated from initial guesses.
            set to False, only W will be estimated.

    Returns:
    """
    # read in reference counts and cNMF results from cnmf_dir
    ref_counts_file = glob("{}/cnmf_tmp/*.norm_counts.h5ad".format(cnmf_dir))
    assert (
        len(ref_counts_file) == 1
    ), "Found more than one norm_counts file in {}/cnmf_tmp/".format(cnmf_dir)
    ref_counts_file = ref_counts_file[0]  # get name of norm counts file
    print("Reading cNMF reference from {}".format(ref_counts_file))
    ref_counts = sc.read(ref_counts_file)
    spectra_file = glob(
        "{}/*.spectra.k_{}.*.consensus.txt".format(cnmf_dir, k)
    )  # get name of gene spectra file
    assert (
        len(spectra_file) < 2
    ), "Found more than one consensus spectra file with k={} in {}".format(k, cnmf_dir)
    assert (
        len(spectra_file) > 0
    ), "Did not find consensus spectra file with k={} in {}".format(k, cnmf_dir)
    spectra_file = spectra_file[0]  # get name of consensus spectra file
    ref_counts.varm["cnmf_spectra"] = pd.read_csv(
        spectra_file, sep="\t", index_col=0
    ).T.values

    # add cnmf_markers to .uns of ST AnnData
    spectra_score_file = glob(
        "{}/*.gene_spectra_score.k_{}.*.txt".format(cnmf_dir, k)
    )  # get name of gene spectra file
    assert (
        len(spectra_score_file) < 2
    ), "Found more than one spectra score file with k={} in {}".format(k, cnmf_dir)
    assert (
        len(spectra_score_file) > 0
    ), "Did not find spectra score file with k={} in {}".format(k, cnmf_dir)
    cnmf_markers(adata, spectra_score_file[0])

    # subset to union of genes from scRNA and ST datasets
    adata = adata[
        :, list(set(adata.var_names).intersection(set(ref_counts.var_names)))
    ].copy()
    ref_counts = ref_counts[
        :, list(set(adata.var_names).intersection(set(ref_counts.var_names)))
    ].copy()

    # normalize and scale ST data as in cNMF preprocessing
    if raw_layer is not None:
        print("Moving layer '{}' to .X for normalization and scaling".format(raw_layer))
        adata.X = adata.layers[raw_layer].copy()
    sc.pp.normalize_total(
        adata, target_sum=1e6, inplace=True
    )  # transcript per million (TPM)
    sc.pp.scale(
        adata, zero_center=False, copy=False
    )  # do not zero-center to avoid negatives

    # perform NMF to factorize adata using ref_counts.varm["cnmf_spectra"]
    print("Factorizing counts using reference gene spectra")
    adata.X = adata.X.astype(ref_counts.varm["cnmf_spectra"].dtype)  # match dtypes
    (usages, spectra, n_iter) = non_negative_factorization(
        X=adata.X,
        H=ref_counts.varm["cnmf_spectra"].T,
        n_components=ref_counts.varm["cnmf_spectra"].shape[1],
        update_H=update_H,
        solver="cd",
        beta_loss="frobenius",
        max_iter=400,
        regularization=None,
        alpha=0.0,
        l1_ratio=0.0,
    )

    # add usages to ST anndata object
    print("Adding NMF usages to anndata object")
    usage = pd.DataFrame(
        usages,
        columns=[x for x in range(1, usages.shape[1] + 1)],
        index=adata.obs_names,
    )
    usage.columns = ["usage_" + str(col) for col in usage.columns]
    # normalize usages to total for each cell
    usage_norm = usage.div(usage.sum(axis=1), axis=0)
    usage_norm.index = usage_norm.index.astype(str)
    # add usages to .obs for visualization
    adata.obs = pd.merge(
        left=adata.obs, right=usage_norm, how="left", left_index=True, right_index=True
    )
    # replace missing values with zeros for all factors
    adata.obs.loc[:, usage_norm.columns].fillna(value=0, inplace=True)
    # add usages as array in .obsm for dimension reduction
    adata.obsm["cnmf_usages"] = adata.obs.loc[:, usage_norm.columns].values

    print("Done!")
    # return anndata, updated spectra, reference spectra, and number of iterations
    return adata, spectra, ref_counts.varm["cnmf_spectra"].T, n_iter


def blur_features_st(adata, features, use_rep="obs", spatial_graph_key=None, n_rings=1):
    """
    Blur values in an `AnnData` object using spatial nearest neighbors

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing Visium data
    use_rep : str
        Representation from `adata.obsm` to use as clustering data (e.g. "X_pca")
    features : list of int or None, optional (default=`None`)
        List of features to use from `adata.obsm[use_rep]` (e.g. [0,1,2,3,4] to
        use first 5 principal components when `use_rep`="X_pca"). If `None`, use
        all features from `adata.obsm[use_rep]`
    spatial_graph_key : str, optional (default=`None`)
        Key in `adata.obsp` containing spatial graph connectivities (i.e.
        `"spatial_connectivities"`). If `None`, compute new spatial graph using
        `n_rings` in `squidpy`.
    n_rings : int, optional (default=1)
        Number of hexagonal rings around each spatial transcriptomics spot to blur
        features by for capturing regional information. Assumes 10X Genomics Visium
        platform.

    Returns
    -------
    adata.obs is edited in place with new blurred columns
    """
    if spatial_graph_key is not None:
        # use existing spatial graph
        assert (
            spatial_graph_key in adata.obsp.keys()
        ), "Spatial connectivities key '{}' not found.".format(spatial_graph_key)
    else:
        # create spatial graph
        print("Computing spatial graph with {} hexagonal rings".format(n_rings))
        sq.gr.spatial_neighbors(adata, coord_type="grid", n_rings=n_rings)
        spatial_graph_key = "spatial_connectivities"  # set key to expected output
    # determine where to pull data values from
    if use_rep == "obs":
        tmp = adata.obs[features].copy()
    else:
        tmp = pd.DataFrame(adata.obsm[use_rep][:, features])
        tmp.columns = [use_rep + "_{}".format(x) for x in features]
    tmp2 = tmp.copy()  # copy of temporary dataframe for dropping blurred features into
    cols = tmp.columns  # get column names
    # perform blurring by nearest spot neighbors
    for x in range(len(tmp)):
        vals = tmp.iloc[
            list(
                np.argwhere(
                    adata.obsp[spatial_graph_key][
                        x,
                    ]
                )[:, 1]
            )
            + [x],
            :,
        ].mean()
        tmp2.iloc[x, :] = vals.values
    # add blurred features to anndata object
    adata.obs[["blur_" + x for x in cols]] = tmp2.loc[:, cols].values


def spa_corr(
    adata,
    features,
    use_rep="obs",
    spatial_graph_key=None,
    n_rings=1,
    mode="moran",
    **kwargs,
):
    """
    Perform spatial autocorrelation analysis (Moran's I or Geary's C) on non-gene
    features from an `AnnData` object (i.e. `.obs` or `.obsm` values)

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing Visium data
    use_rep : str
        Representation from `adata.obsm` to use as clustering data (e.g. "X_pca")
    features : list of int or None, optional (default=`None`)
        List of features to use from `adata.obsm[use_rep]` (e.g. [0,1,2,3,4] to
        use first 5 principal components when `use_rep`="X_pca"). If `None`, use
        all features from `adata.obsm[use_rep]`
    spatial_graph_key : str, optional (default=`None`)
        Key in `adata.obsp` containing spatial graph connectivities (i.e.
        `"spatial_connectivities"`). If `None`, compute new spatial graph using
        `n_rings` in `squidpy`.
    n_rings : int, optional (default=1)
        Number of hexagonal rings around each spatial transcriptomics spot to blur
        features by for capturing regional information. Assumes 10X Genomics Visium
        platform.
    mode : literal ["moran","geary"], optional (default="moran")
        Statistic to use. "moran" for Moran's I, "geary" for Geary's C.
    **kwargs : optional
        Keyword args to pass to `sq.gr.spatial_autocorr`, such as `jobs=4` to
        parallelize over 4 CPUs

    Returns
    -------
    stats : pd.DataFrame
        Results of spatial autocorrelation statistical testing
    """
    if spatial_graph_key is not None:
        # use existing spatial graph
        assert (
            spatial_graph_key in adata.obsp.keys()
        ), "Spatial connectivities key '{}' not found.".format(spatial_graph_key)
    else:
        # create spatial graph
        print("Computing spatial graph with {} hexagonal rings".format(n_rings))
        sq.gr.spatial_neighbors(adata, coord_type="grid", n_rings=n_rings)
        spatial_graph_key = "spatial_connectivities"  # set key to expected output
    # create AnnData object with specified non-gene features
    if use_rep == "obs":
        tmp = sc.AnnData(adata.obs[features].values, obs=adata.obs)
        tmp.var_names = features
    else:
        tmp = sc.AnnData(adata.obsm[use_rep][:, features], obs=adata.obs)
        tmp.var_names = [use_rep + "_{}".format(x) for x in features]
    # transfer spatial graph to new object
    tmp.obsp[spatial_graph_key] = adata.obsp[spatial_graph_key]
    # perform spatial autocorrelation analysis
    sq.gr.spatial_autocorr(tmp, genes=tmp.var_names, mode=mode, **kwargs)
    outname = "moranI" if mode == "moran" else "gearyC"
    return tmp.uns[outname]


def covar_matrix(
    X,
    save_to=None,
    correlate=True,
    cbar_label="Correlation",
    cmap="bwr",
    cluster_cmap="plasma",
    figsize=(10, 10),
):
    """
    Plot covariance matrix as heatmap and identify hierarchical clusters

    Parameters
    ----------
    X : pd.DataFrame
        The data in format `obs x var`. Each `var` will be correlated with every other
        `var` using all available `obs`.
    save_to : str, optional (default=`None`)
        Path to `.png` file to save heatmap to. If `None`, do not save plot.
    correlate : bool, optional (default=`True`)
        Whether or not to correlate X first

    Returns
    -------
    heatmap : sns.ClusterGrid
        Seaborn clustermap object containing plot
    clusters : pd.Series
        Index will match columns of `X`, values will be cluster IDs corresponding to
        heatmap.
    """
    corr = X.corr() if correlate else X  # pandas correlation function

    heatmap1 = sns.clustermap(corr)  # first clustermap to establish relationships
    plt.close()

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool))
    # apply the inverse permutation to the mask
    mask = mask[np.argsort(heatmap1.dendrogram_row.reordered_ind), :]
    mask = mask[:, np.argsort(heatmap1.dendrogram_col.reordered_ind)]

    # create dendrogram scipy object to define clusters for coloring
    _, dendroax = plt.subplots(1, 1)
    dendro = hierarchy.dendrogram(heatmap1.dendrogram_row.linkage, ax=dendroax)
    plt.close()

    lut = dict(
        zip(
            np.unique(dendro["leaves_color_list"]),
            mcp.gen_color(
                cmap=cluster_cmap, n=len(np.unique(dendro["leaves_color_list"]))
            ),
        )
    )
    row_colors = pd.Series(dendro["leaves_color_list"]).map(lut)

    # get clusters out of dendrogram object
    clusters = pd.Series(dendro["leaves_color_list"]).iloc[
        np.argsort(heatmap1.dendrogram_row.reordered_ind)
    ]
    clusters.index = corr.index  # get var names associated with each cluster

    # final clustermap generation, overwrite heatmap
    heatmap = sns.clustermap(
        corr,
        mask=mask,
        cmap=cmap,
        robust=True,
        xticklabels=True,
        yticklabels=True,
        figsize=figsize,
        row_colors=list(
            row_colors.iloc[np.argsort(heatmap1.dendrogram_row.reordered_ind)]
        ),
        col_colors=list(
            row_colors.iloc[np.argsort(heatmap1.dendrogram_row.reordered_ind)]
        ),
        dendrogram_ratio=0.05,
        colors_ratio=0.02,
        row_cluster=True,
        col_cluster=True,
        cbar_kws={"label": cbar_label, "orientation": "horizontal", "shrink": 0.5},
        cbar_pos=(0.6, 0.6, 0.15, 0.03),
        center=0.0,
    )
    # make it nasty
    heatmap.ax_heatmap.set_xticklabels(heatmap.ax_heatmap.get_xticklabels(), fontsize=7)
    heatmap.ax_heatmap.set_yticklabels(heatmap.ax_heatmap.get_yticklabels(), fontsize=7)
    heatmap.ax_cbar.set_xticklabels(heatmap.ax_cbar.get_xticklabels(), fontsize=7)
    heatmap.ax_cbar.set_xlabel("Correlation", fontsize=7)
    heatmap.ax_heatmap.yaxis.set_label_position("left")
    heatmap.ax_heatmap.yaxis.tick_left()
    heatmap.ax_heatmap.tick_params(axis="y", pad=18, left=False)
    heatmap.ax_heatmap.tick_params(axis="x", pad=18)
    heatmap.ax_row_dendrogram.set_visible(False)
    heatmap.ax_col_dendrogram.set_visible(False)
    col_pos = heatmap.ax_col_colors.get_position()
    row_pos = heatmap.ax_row_colors.get_position()
    heatmap.ax_col_colors.set_position(
        (
            col_pos.bounds[0],
            row_pos.bounds[1] - (col_pos.height + (col_pos.height * 0.2)),
            col_pos.width,
            col_pos.height,
        )
    )

    # save if you want
    if save_to is not None:
        heatmap.savefig(save_to)

    return heatmap, clusters


def covar_matrix_asym(
    corr,
    save_to=None,
    cbar_label="Correlation",
    cmap="bwr",
    cluster_cmap="plasma",
    figsize=(10, 10),
):
    """
    Plot asymmetrical covariance matrix as heatmap and identify hierarchical clusters

    Parameters
    ----------
    corr : pd.DataFrame
        The correlation matrix as a dataframe
    save_to : str, optional (default=`None`)
        Path to `.png` file to save heatmap to. If `None`, do not save plot.

    Returns
    -------
    heatmap : sns.ClusterGrid
        Seaborn clustermap object containing plot
    row_clusters : pd.Series
        Index will match rows of `corr`, values will be cluster IDs corresponding to
        heatmap.
    col_clusters : pd.Series
        Index will match columns of `corr`, values will be cluster IDs corresponding
        to heatmap.
    """
    heatmap1 = sns.clustermap(corr)  # first clustermap to establish relationships
    plt.close()

    # create dendrogram scipy object to define row clusters for coloring
    _, dendroax = plt.subplots(1, 1)
    dendro = hierarchy.dendrogram(heatmap1.dendrogram_row.linkage, ax=dendroax)
    plt.close()
    lut = dict(
        zip(
            np.unique(dendro["leaves_color_list"]),
            mcp.gen_color(
                cmap=cluster_cmap, n=len(np.unique(dendro["leaves_color_list"]))
            ),
        )
    )
    row_colors = pd.Series(dendro["leaves_color_list"]).map(lut)
    # get clusters out of dendrogram object
    row_clusters = pd.Series(dendro["leaves_color_list"]).iloc[
        np.argsort(heatmap1.dendrogram_row.reordered_ind)
    ]
    row_clusters.index = corr.index  # get var names associated with each cluster

    # create dendrogram scipy object to define column clusters for coloring
    _, dendroax = plt.subplots(1, 1)
    dendro = hierarchy.dendrogram(heatmap1.dendrogram_col.linkage, ax=dendroax)
    plt.close()
    lut = dict(
        zip(
            np.unique(dendro["leaves_color_list"]),
            mcp.gen_color(
                cmap=cluster_cmap, n=len(np.unique(dendro["leaves_color_list"]))
            ),
        )
    )
    col_colors = pd.Series(dendro["leaves_color_list"]).map(lut)
    # get clusters out of dendrogram object
    col_clusters = pd.Series(dendro["leaves_color_list"]).iloc[
        np.argsort(heatmap1.dendrogram_col.reordered_ind)
    ]
    col_clusters.index = corr.columns  # get var names associated with each cluster

    # final clustermap generation, overwrite heatmap
    heatmap = sns.clustermap(
        corr,
        cmap=cmap,
        robust=True,
        xticklabels=True,
        yticklabels=True,
        figsize=figsize,
        row_colors=list(
            row_colors.iloc[np.argsort(heatmap1.dendrogram_row.reordered_ind)]
        ),
        col_colors=list(
            col_colors.iloc[np.argsort(heatmap1.dendrogram_col.reordered_ind)]
        ),
        dendrogram_ratio=0.05,
        colors_ratio=0.02,
        row_cluster=True,
        col_cluster=True,
        cbar_kws={"label": cbar_label, "orientation": "horizontal", "shrink": 0.5},
        cbar_pos=(0.85, 0.85, 0.15, 0.03),
        center=0.0,
    )
    # make it nasty
    heatmap.ax_heatmap.set_xticklabels(heatmap.ax_heatmap.get_xticklabels(), fontsize=7)
    heatmap.ax_heatmap.set_yticklabels(heatmap.ax_heatmap.get_yticklabels(), fontsize=7)
    heatmap.ax_cbar.set_xticklabels(heatmap.ax_cbar.get_xticklabels(), fontsize=7)
    heatmap.ax_cbar.set_xlabel("Correlation", fontsize=7)
    heatmap.ax_heatmap.yaxis.set_label_position("left")
    heatmap.ax_heatmap.yaxis.tick_left()
    heatmap.ax_heatmap.tick_params(axis="y", pad=18, left=False)
    heatmap.ax_heatmap.tick_params(axis="x", pad=18)
    heatmap.ax_row_dendrogram.set_visible(False)
    heatmap.ax_col_dendrogram.set_visible(False)
    col_pos = heatmap.ax_col_colors.get_position()
    row_pos = heatmap.ax_row_colors.get_position()
    heatmap.ax_col_colors.set_position(
        (
            col_pos.bounds[0],
            row_pos.bounds[1] - (col_pos.height + (col_pos.height * 0.2)),
            col_pos.width,
            col_pos.height,
        )
    )

    # save if you want
    if save_to is not None:
        heatmap.savefig(save_to)

    return heatmap, row_clusters, col_clusters


def community_detection(X, min_abundance=0.05, mean_corr_percentile=0.5, **kwargs):
    """
    Plot covariance matrix as heatmap and identify hierarchically clustered communities

    Parameters
    ----------
    X : pd.DataFrame
        The data in format `obs x var`. Each `var` will be correlated with every other
        `var` using all available `obs`.
    min_abundance : float, optional (default=0.05)
        Minimum value to accept for correlation, below which values are set to 0.0
    mean_corr_percentile : float, optional (default=0.5)
        Minimum average correlation percentile (compared to global maximum cell state
        correlation needed) to keep a community
    **kwargs : optional
        Keyword args to pass to `covar_matrix`

    Returns
    -------
    heatmap : sns.ClusterGrid
        Seaborn clustermap object containing plot
    clusters : pd.Series
        Index will match columns of `X`, values will be cluster IDs corresponding to
        heatmap.
    """
    # set values below min_abundance to 0.0 and correlate, dropping all NaN cols/rows
    X = (
        X.mask(X < min_abundance, other=0.0)
        .corr()
        .dropna(axis=0, how="all")
        .dropna(axis=1, how="all")
    )
    # perform covariance analysis and clustering on correlated X
    heatmap, clusters = covar_matrix(
        X=X,
        correlate=False,
        **kwargs,
    )
    # get global maximum correlation for reference
    corr_global_max = X.values[np.triu_indices_from(X, k=1)].max()
    # loop through heirarchical clusters and determine percentile of average
    # correlation value in each cluster
    out = {}
    for cluster in clusters.unique():
        print(cluster, end=": ")
        tmp = X.loc[
            clusters.loc[
                clusters == cluster,
            ].index,
            clusters.loc[
                clusters == cluster,
            ].index,
        ]
        out[cluster] = (
            tmp.values[np.triu_indices_from(tmp, k=1)].mean() / corr_global_max
        )
        print(out[cluster])
    # create dataframe of cell state communities
    communities = pd.DataFrame(clusters.map(out), columns=["mean_corr"])
    communities["community"] = clusters
    # threshold communities for those with mean correlation percentile larger than
    # mean_corr_percentile
    communities = communities.loc[
        communities.mean_corr >= mean_corr_percentile
    ].sort_values("community")
    # return plot and communities dataframe
    return heatmap, communities


def community_detection_asym(
    X, col_obs, row_obs, min_abundance=0.05, mean_corr_percentile=0.5, **kwargs
):
    """
    Plot covariance matrix as heatmap and identify hierarchically clustered communities

    Parameters
    ----------
    X : pd.DataFrame
        The data with features to correlate in columns. All features from `col_obs`
        will be correlated with features from `row_obs`.
    col_obs : list of str
        Names of columns in `X` that you want to be columns in the correlation
    row_obs : list of str
        Names of columns in `X` that you want to be rows in the correlation
    min_abundance : float, optional (default=0.05)
        Minimum value to accept for correlation, below which values are set to 0.0
    mean_corr_percentile : float, optional (default=0.5)
        Minimum average correlation percentile (compared to global maximum cell state
        correlation needed) to keep a community
    **kwargs : optional
        Keyword args to pass to `covar_matrix_asym`

    Returns
    -------
    heatmap : sns.ClusterGrid
        Seaborn clustermap object containing plot
    clusters : pd.Series
        Index will match columns of `X`, values will be cluster IDs corresponding to
        heatmap.
    """
    # set values below min_abundance to 0.0 and correlate, dropping all NaN cols/rows
    X = X[row_obs + col_obs].copy()
    X = X.mask(X < min_abundance, other=0.0)
    corr = {}
    for c in col_obs:
        corr[c] = {}
        for r in row_obs:
            corr[c][r] = X[c].corr(X[r])
    corr = pd.DataFrame(corr)  # coerce to df
    corr.dropna(axis=0, how="all")
    corr.dropna(axis=1, how="all")
    # perform covariance analysis and clustering on correlated X
    heatmap, row_clusters, col_clusters = covar_matrix_asym(
        corr=corr,
        **kwargs,
    )
    # get global maximum correlation for reference
    corr_global_max = corr.values.max()
    # loop through heirarchical clusters and determine percentile of average
    # correlation value in each row cluster
    out = {}
    for cluster in row_clusters.unique():
        print(cluster, end=" (rows): ")
        tmp = corr.loc[
            row_clusters.loc[
                row_clusters == cluster,
            ].index,
        ]
        out[cluster] = tmp.values.mean() / corr_global_max
        print(out[cluster])
    # create dataframe of cell state communities
    row_communities = pd.DataFrame(row_clusters.map(out), columns=["mean_corr"])
    row_communities["community"] = row_clusters
    # threshold communities for those with mean correlation percentile larger than
    # mean_corr_percentile
    row_communities = row_communities.loc[
        row_communities.mean_corr >= mean_corr_percentile
    ].sort_values("community")
    # correlation value in each column cluster
    out = {}
    for cluster in col_clusters.unique():
        print(cluster, end=" (cols): ")
        tmp = corr.loc[
            :,
            col_clusters.loc[
                col_clusters == cluster,
            ].index,
        ]
        out[cluster] = tmp.values.mean() / corr_global_max
        print(out[cluster])
    # create dataframe of cell state communities
    col_communities = pd.DataFrame(col_clusters.map(out), columns=["mean_corr"])
    col_communities["community"] = col_clusters
    # threshold communities for those with mean correlation percentile larger than
    # mean_corr_percentile
    col_communities = col_communities.loc[
        col_communities.mean_corr >= mean_corr_percentile
    ].sort_values("community")
    # return plot and communities dataframe
    return heatmap, row_communities, col_communities
