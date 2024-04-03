"""Plotting for graph functions."""
from __future__ import annotations

from scipy.sparse import csr_matrix
from types import MappingProxyType
from typing import Union  # noqa: F401
from typing import Any, Mapping, Sequence, TYPE_CHECKING
from pathlib import Path
from typing_extensions import Literal

from anndata import AnnData

import numpy as np
import pandas as pd

import seaborn as sns
import squidpy as sq
import matplotlib.pyplot as plt

from squidpy.gr._utils import (
    _get_valid_values,
    _assert_categorical_obs,
    _assert_non_empty_sequence,
)
from squidpy.pl._utils import save_fig
from squidpy.pl._color_utils import Palette_t, _get_palette
from squidpy._constants._pkg_constants import Key


def threshold_expression(adata, thresh_dict, fractions=False, percentiles=False):
    """
    Bin expression into (+) and (-) labels for given cell states, archetypes, or genes

    Parameters
    ----------
    adata : anndata.AnnData
        ST data object
    thresh_dict : dictionary
        Dictionary with names of `adata.obs` columns or `adata.var_names` containing
        fractional abundances of cell states or archetypes as keys and threshold value
        above which to label positive pixels as values (e.g. {"MYE4_VUMCrefNMF30":0.1}
        to call MYE4+ above 10%)
    fractions : bool, optional (default=`False`)
        Values in `thresh_dict` are fractions on [0.0, 1.0] rather than counts values.
        Expression is min-max scaled before thresholding. ONLY FOR GENES
    percentiles : bool, optional (default=`False`)
    """
    for i in list(thresh_dict.keys()):
        if i in adata.obs.columns:
            print("Thresholding .obs column {} at {}".format(i, thresh_dict[i]))

            if fractions:
                adata.obs["{}_expr".format(i)] = adata.obs[i] - adata.obs[i].min()
                adata.obs["{}_expr".format(i)] = (
                    adata.obs["{}_expr".format(i)]
                    / adata.obs["{}_expr".format(i)].max()
                )
            else:
                adata.obs["{}_expr".format(i)] = adata.obs[i].values

        elif i in adata.var_names:
            print("Thresholding gene {} at {}".format(i, thresh_dict[i]))
            if isinstance(adata.X, csr_matrix):
                adata.obs["{}_expr".format(i)] = adata[:, i].X.todense()
            else:
                adata.obs["{}_expr".format(i)] = adata[:, i].X

            if fractions:
                adata.obs["{}_expr".format(i)] = (
                    adata.obs["{}_expr".format(i)]
                    - adata.obs["{}_expr".format(i)].min()
                )
                adata.obs["{}_expr".format(i)] = (
                    adata.obs["{}_expr".format(i)]
                    / adata.obs["{}_expr".format(i)].max()
                )

        else:
            print("{} not detected in adata.obs.columns or adata.var_names!".format(i))

        if percentiles:
            if np.percentile(adata.obs["{}_expr".format(i)], thresh_dict[i]) == adata.obs["{}_expr".format(i)].min():
                print("Not enough signal in {} to proceed; skipping.".format(i))
            else:
                adata.obs["{}_thresh".format(i)] = pd.cut(
                    x=adata.obs["{}_expr".format(i)],
                    bins=[
                        adata.obs["{}_expr".format(i)].min(),
                        np.percentile(adata.obs["{}_expr".format(i)], thresh_dict[i]),
                        adata.obs["{}_expr".format(i)].max()+0.1
                    ],
                    # label low values as '-', high values as "+"
                    labels=["{}-".format(i), "{}+".format(i)],
                )
        else:
            adata.obs["{}_thresh".format(i)] = pd.cut(
                x=adata.obs["{}_expr".format(i)],
                # assume values are in [0.0, 1.0]
                bins=[-1, thresh_dict[i], 1.1],
                # label low values as '-', high values as "+"
                labels=["{}-".format(i), "{}+".format(i)],
            )
        adata.obs.drop(columns="{}_expr".format(i), inplace=True)


def _get_data(adata: AnnData, cluster_key: str, func_name: str, **kwargs: Any) -> Any:
    key = getattr(Key.uns, func_name)(cluster_key, **kwargs)
    try:
        return adata.uns[key]
    except KeyError:
        raise KeyError(
            f"Unable to get the data from `adata.uns[{key!r}]`. "
            f"Please run `squidpy.gr.{func_name}(..., cluster_key={cluster_key!r})` first."
        ) from None


def co_occurrence(
    adata: AnnData,
    cluster_key: str,
    palette: Palette_t = None,
    categories: str | Sequence[str] | None = None,
    clusters: str | Sequence[str] | None = None,
    dist_max: float | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int | None = None,
    save: str | Path | None = None,
    legend_kwargs: Mapping[str, Any] = MappingProxyType({}),
    **kwargs: Any,
) -> None:
    """
    Plot co-occurrence probability ratio for each cluster.

    The co-occurrence is computed by :func:`squidpy.gr.co_occurrence`.

    Parameters
    ----------
    %(adata)s
    %(cluster_key)s
    categories
        Category instances for which to compare conditional probability.
    clusters
        Cluster instances for which to plot conditional probability.
    dist_max
        Maximum distance to consider for plotting co-occurrence results.
    %(cat_plotting)s
    legend_kwargs
        Keyword arguments for :func:`matplotlib.pyplot.legend`.
    kwargs
        Keyword arguments for :func:`seaborn.lineplot`.

    Returns
    -------
    %(plotting_returns)s
    """
    _assert_categorical_obs(adata, key=cluster_key)
    occurrence_data = _get_data(
        adata, cluster_key=cluster_key, func_name="co_occurrence"
    )

    legend_kwargs = dict(legend_kwargs)
    if "loc" not in legend_kwargs:
        legend_kwargs["loc"] = "center left"
        legend_kwargs.setdefault("bbox_to_anchor", (1, 0.5, 0.3, 0.3))

    # define categories
    categories = (
        adata.obs[cluster_key].cat.categories if categories is None else categories
    )
    # subset occurrence_data to categories given
    cat_locs = [adata.obs[cluster_key].cat.categories.get_loc(x) for x in categories]
    occurrence_data["occ"] = occurrence_data["occ"][cat_locs, :, :]
    occurrence_data["occ"] = occurrence_data["occ"][:, cat_locs, :]
    # subset occurrence_data once again based on dist_max
    if dist_max:
        i_max = np.where(occurrence_data["interval"] > dist_max)[0].min()
        occurrence_data["occ"] = occurrence_data["occ"][:, :, :i_max]
        occurrence_data["interval"] = occurrence_data["interval"][: i_max + 1]

    out = occurrence_data["occ"]
    interval = occurrence_data["interval"][1:]

    clusters = categories if clusters is None else clusters
    clusters = _assert_non_empty_sequence(clusters, name="clusters")
    clusters = sorted(_get_valid_values(clusters, categories))

    palette = _get_palette(
        adata,
        cluster_key=cluster_key,
        categories=adata.obs[cluster_key].cat.categories,
        palette=palette,
    )
    palette = dict((k, palette[k]) for k in categories)

    fig, axs = plt.subplots(
        1,
        len(clusters),
        figsize=(5 * len(clusters), 5) if figsize is None else figsize,
        dpi=dpi,
        constrained_layout=True,
    )
    axs = np.ravel(axs)  # make into iterable

    for g, ax in zip(clusters, axs):
        idx = np.where(pd.Index(categories) == g)[0][0]
        df = pd.DataFrame(out[idx, :, :].T, columns=categories).melt(
            var_name=cluster_key, value_name="probability"
        )
        df["distance"] = np.tile(interval, len(categories))
        plt.axhline(y=1, color="k", linestyle="--", linewidth=1.8)
        sns.lineplot(
            x="distance",
            y="probability",
            data=df,
            dashes=False,
            hue=cluster_key,
            hue_order=categories,
            palette=palette,
            ax=ax,
            **kwargs,
        )
        ax.legend(**legend_kwargs)
        ax.set_xlabel("Distance")
        ax.set_ylabel(
            rf"$\frac{{p(exp|{g})}}{{p(exp)}}$", rotation="horizontal", labelpad=28
        )

    if save is not None:
        save_fig(fig, path=save)


def co_occurrence_refNMF(
    adata,
    cluster_key="CNV Clone",
    ref_clusters="1",
    cluster_key_add="MYE4_thresh",
    cat_add="MYE4+",
    **kwargs,
):
    """
    Perform co-occurrence testing with multiple categorical spot labels

    Parameters
    ----------
    adata : anndata.AnnData
        ST data object
    cluster_key : str
        Key from `adata.obs` containing primary labels for co-occurrence comparison (reference)
    ref_clusters : str or list of str
        Cluster instances for which to plot conditional probability
    cluster_key_add : str
        Key from `adata.obs` containing secondary labels for co-occurrence testing (comparison)
    cat_add : str
        Category from `cluster_key_add` to add to `cluster_key` reference for co-occurrence testing
    """
    tmp = adata.copy()
    tmp.obs["{} vs. {}".format(cluster_key, cat_add)] = tmp.obs[cluster_key]
    tmp.obs["{} vs. {}".format(cluster_key, cat_add)] = tmp.obs[
        "{} vs. {}".format(cluster_key, cat_add)
    ].cat.add_categories(cat_add)
    tmp.obs.loc[
        tmp.obs[cluster_key_add] == cat_add, "{} vs. {}".format(cluster_key, cat_add)
    ] = cat_add
    tmp.uns["{} vs. {}_colors".format(cluster_key, cat_add)] = np.append(
        tmp.uns["{}_colors".format(cluster_key)], "#000000"
    )  # make new line black
    sq.gr.co_occurrence(tmp, cluster_key="{} vs. {}".format(cluster_key, cat_add))
    co_occurrence(
        tmp,
        cluster_key="{} vs. {}".format(cluster_key, cat_add),
        clusters=ref_clusters,
        **kwargs,
    )


def immune_excl_cooccurrence(a_ws, pat, block, markers_thresh_dict, TMA_dist=None, **kwargs):
    # threshold expression first
    threshold_expression(
        adata = a_ws,
        thresh_dict = markers_thresh_dict,
        **kwargs,
    )

    # perform analysis per clone region
    for clone in [x for x in a_ws.obs["CNV Clone"].cat.categories if x not in ["S","E"]]:
        for marker in list(markers_thresh_dict.keys()):
            if "{}_thresh".format(marker) in a_ws.obs.columns:
                try:
                    co_occurrence_refNMF(
                        adata=a_ws,
                        cluster_key="CNV Clone",
                        ref_clusters=clone,
                        cluster_key_add="{}_thresh".format(marker),
                        cat_add="{}+".format(marker),
                        save="cooccurrence_{}_{}_clone{}_{}{}.png".format(pat, block, clone, marker, "_TMA" if TMA_dist is not None else ""),
                        # which categories to include in plot as comparison to ref_cluster
                        categories=[clone,"S","{}+".format(marker)],
                        figsize=(5,3),
                        legend_kwargs={"loc":"best", "fontsize":10},
                    )
                except:
                    print("Failed on {}_thresh for clone {} of {}".format(marker, clone, pat))

    # combine major clones into 'Epi.' region and analyze again
    a_ws.obs["Compartment"] = "Epi."
    a_ws.obs.loc[a_ws.obs["CNV Clone"].isin(["S","E"]), "Compartment"] = "S"
    a_ws.obs.Compartment = a_ws.obs.Compartment.astype("category")
    a_ws.uns["Compartment_colors"] = pd.Index(['#1f77b4', '#e377c2'])
    for marker in list(markers_thresh_dict.keys()):
        if "{}_thresh".format(marker) in a_ws.obs.columns:
            try:
                co_occurrence_refNMF(
                    adata=a_ws,
                    cluster_key="Compartment",
                    ref_clusters="Epi.",
                    cluster_key_add="{}_thresh".format(marker),
                    cat_add="{}+".format(marker),
                    save="cooccurrence_{}_{}_EPI_{}{}.png".format(pat, block, marker, "_TMA" if TMA_dist is not None else ""),
                    # which categories to include in plot as comparison to ref_cluster
                    categories=["Epi.","S","{}+".format(marker)],
                    # restrict xlim to reasonable value for TMA (within one core distance)
                    dist_max=TMA_dist,
                    figsize=(5,3),
                    legend_kwargs={"loc":"best", "fontsize":10},
                )
            except:
                print("Failed on {}_thresh for EPI of {}".format(marker, pat))
