import infercnvpy as cnv
import numpy as np
import scanpy as sc
import seaborn as sns

# define color dictionary for stroma and tumor edge
base_colordict = {
    "S": "#e377c2",  # pink to match pathology_annotation:smooth_muscle
    "E": "#d62728",  # red to match pathology_annotation:carcinoma_border
}


def load_cnv(
        pat,
        sample_key,
        CNV_group="patient_name",
        dataset_dir="datasets/",
        infercnv_dir="infercnv/",
    ):
    """
    pat: str
        Patient ID
    sample_key: pd.DataFrame
        Visium sample key dataframe
    """
    # read in patient anndatas and concatenate
    key_tmp = sample_key.loc[sample_key[CNV_group] == pat, :].copy()
    outs = []
    for s in key_tmp.index:
        a = sc.read("{}/{}_master.h5ad".format(dataset_dir, s))
        print("Read adata from {}/{}_master.h5ad".format(dataset_dir, s))

        # read in CNV matrix and put in a.obsm slot
        tmp = np.load("{}/{}_cnv.npz".format(infercnv_dir, s), allow_pickle="TRUE")
        a.obsm["X_cnv"] = tmp.f.arr_0.item()

        # append to patient adata list
        outs.append(a)

    # concatenate anndata objects
    a_comb = outs[0].concatenate(
        outs[1:],
        join="outer",
        batch_categories=list(key_tmp.index),
        fill_value=0,
    )
    del a_comb.var

    # read in CNV genomic partitions
    a_comb.uns["cnv"] = np.load(
        "{}/uns_cnv_{}.npy".format(
            infercnv_dir,
            sample_key.loc[sample_key[CNV_group] == pat, "CNV_group"][0]
        ),
        allow_pickle="TRUE",
    ).item()

    return a_comb


def curate_cnv_ST(
        pat,
        d,
        sample_key,
        CNV_group="CNV_group",
        dataset_dir="datasets/",
        infercnv_dir="infercnv/",
    ):
    """
    pat: str
        Patient ID
    d: dict
        Dictionary of cnv_leiden cluster remapping
    sample_key: pd.DataFrame
        Visium sample key dataframe
    """
    # define clone color dictionary
    clones = [x for x in d.values() if x not in ["S", "E"]]
    clones_colordict = dict(
        zip(clones, sns.color_palette("tab10", len(clones)).as_hex())
    )
    colordict = {**base_colordict, **clones_colordict}
    # save CNV cluster colors
    np.save("{}/cnv_clone_colors_{}.npy".format(infercnv_dir, pat), colordict)

    # read in patient anndatas and concatenate
    key_tmp = sample_key.loc[sample_key[CNV_group] == pat, :].copy()
    outs = []
    for s in key_tmp.index:
        a = sc.read("{}/{}_master.h5ad".format(dataset_dir, s))
        print("Read adata from {}/{}_master.h5ad".format(dataset_dir, s))

        # map new CNV clone values into adata.obs
        print("Mapping new 'CNV clone' values into adata.obs")
        a.obs["CNV clone"] = a.obs.cnv_leiden.astype(str).replace(d)
        a.obs["CNV clone"] = a.obs["CNV clone"].astype("category")
        # create colormap as well
        a.uns["CNV clone_colors"] = [
            colordict[x] for x in a.obs["CNV clone"].cat.categories
        ]

        # save anndata back to master
        a.write("{}/{}_master.h5ad".format(dataset_dir, s), compression="gzip")
        print("Saved adata to {}/{}_master.h5ad".format(dataset_dir, s))

        # read in CNV matrix and put in a.obsm slot
        tmp = np.load("{}/{}_cnv.npz".format(infercnv_dir, s), allow_pickle="TRUE")
        a.obsm["X_cnv"] = tmp.f.arr_0.item()

        # plot spatial
        sc.pl.spatial(
            a,
            color=list(
                set(
                    ["cnv_score", "CNV score", "CNV clone", "pathology_annotation"]
                ).intersection(set(a.obs.columns))
            ),
            size=1.7,
            ncols=3,
            img_key=None,
            frameon=False,
            vmin=0.0,
            vmax=0.12,
            cmap="viridis",
            save="_{}_{}_CNV_curated.png".format(pat, s),
        )

        # append to patient adata list
        outs.append(a)

    # concatenate anndata objects
    a_comb = outs[0].concatenate(
        outs[1:],
        join="outer",
        batch_categories=list(key_tmp.index),
        fill_value=0,
    )
    del a_comb.var

    # read in CNV genomic partitions
    a_comb.uns["cnv"] = np.load(
        "{}/uns_cnv_{}.npy".format(
            infercnv_dir,
            sample_key.loc[sample_key[CNV_group] == pat, "CNV_group"][0]
        ),
        allow_pickle="TRUE",
    ).item()
    # create colormap as well
    a_comb.obs["CNV clone"] = a_comb.obs["CNV clone"].astype("category")
    a_comb.uns["CNV clone_colors"] = [
        colordict[x] for x in a_comb.obs["CNV clone"].cat.categories
    ]

    # plot heatmap with CNV Leiden clusters
    print("Plotting CNV heatmap with 'CNV clone' clusters for {}".format(pat))
    cnv.pl.chromosome_heatmap(
        a_comb,
        groupby="CNV clone",
        save="_{}.png".format(pat),
        dendrogram=True,
        figsize=(12, 8),
    )


def curate_cnv_scRNA(
        adata,
        pat,
        d,
        infercnv_dir="infercnv/",
    ):
    """
    adata: anndata.AnnData
        Anndata object containing scRNA-seq counts
    pat: str
        Patient ID
    d: dict
        Dictionary of cnv_leiden cluster remapping
    """
    # define clone color dictionary
    clones = [x for x in d.values()]
    clones_colordict = dict(
        zip(clones, sns.color_palette("tab10", len(clones)).as_hex())
    )
    colordict = {**base_colordict, **clones_colordict}
    # save CNV cluster colors
    np.save("{}/cnv_clone_colors_{}.npy".format(infercnv_dir, pat), colordict)
    # map new CNV clone values into adata.obs
    print("Mapping new 'CNV clone' values into adata.obs")
    tmp = adata[
        adata.obs.Patient == pat, :
    ].copy()  # make copy of adata from patient ID
    tmp.obs["CNV clone"] = "S"  # set background to stroma
    tmp.obs.loc[tmp.obs.cnv_leiden.isin(list(d.keys())), "CNV clone"] = (
        tmp.obs.loc[tmp.obs.cnv_leiden.isin(list(d.keys())), "cnv_leiden"]
        .astype(str)
        .replace(d)
    )
    tmp.obs["CNV clone"] = tmp.obs["CNV clone"].astype("category")
    # create colormap as well
    tmp.uns["CNV clone_colors"] = [
        colordict[x] for x in tmp.obs["CNV clone"].cat.categories
    ]
    # plot heatmap with CNV Leiden clusters
    print("Plotting CNV heatmap with 'CNV clone' clusters for {}".format(pat))
    cnv.pl.chromosome_heatmap(
        tmp,
        groupby="CNV clone",
        save="_VUMC_{}_curated.png".format(pat),
        dendrogram=True,
        figsize=(12, 8),
    )
    # transfer 'CNV clone' to master adata object
    if "CNV clone" not in adata.obs:
        # initialize 'CNV clone' column in adata.obs
        adata.obs["CNV clone"] = np.nan
    tmp.obs["CNV clone"] = (
        tmp.obs.Patient.astype(str) + " " + tmp.obs["CNV clone"].astype(str)
    )
    adata.obs.loc[adata.obs.Patient == pat, "CNV clone"] = tmp.obs["CNV clone"]


def remap_cnv(
    X_cnv,
    uns_cnv,
    target_size=[17, 10, 9, 6, 7, 8, 7, 5, 6, 6, 9, 9, 2, 5, 5, 7, 9, 2, 12, 4, 2, 3],
):
    X_cnv_out = None  # initialize output numpy array
    uns_cnv_out = {"chr_pos": {}}  # initialize output uns dict
    for i in range(len(uns_cnv["chr_pos"])):
        # doing this one chromosome at a time
        current_chr = "chr{}".format(i + 1)
        next_chr = "chr{}".format(i + 2)
        print(current_chr, end=": ")  # status update
        if i + 2 < 23:
            # calculate number of chr partitions
            start_pos = uns_cnv["chr_pos"][current_chr]
            end_pos = uns_cnv["chr_pos"][next_chr]
            size = end_pos - start_pos
        else:
            # if chr22, use size of X_cnv to calculate partition size
            start_pos = uns_cnv["chr_pos"][current_chr]
            end_pos = X_cnv.shape[1]
            size = end_pos - start_pos
        print("orig. size = {}".format(size), end=" - ")  # status update
        divisor = size // target_size[i]  # floor division to get reshape divisor
        remainder_len = size - (
            divisor * target_size[i]
        )  # calculate remainder to aggregate first
        print(
            "taking {} off the top and reshaping by {}".format(remainder_len, divisor)
        )  # status update
        # now manipulate X_cnv
        tmp = np.array(
            X_cnv[:, start_pos:end_pos].todense().copy()
        )  # copy chromosome columns into tmp
        # hold the remainder out first, then aggregate columns
        if remainder_len > 0:
            remainder_vec = tmp[:, 0:remainder_len].mean(axis=1)
            # aggregate by averaging consecutive cols
            agg = np.reshape(
                tmp[:, remainder_len:], (tmp.shape[0], target_size[i], divisor)
            ).mean(axis=2)
        else:
            remainder_vec = tmp[
                :, 0
            ].copy()  # if no remainder, just duplicate the first position
            # aggregate by averaging consecutive cols
            agg = np.reshape(tmp, (tmp.shape[0], target_size[i], divisor)).mean(axis=2)
        if X_cnv_out is None:
            X_cnv_out = np.append(remainder_vec[np.newaxis].T, agg, 1)
        else:
            X_cnv_out = np.append(X_cnv_out, remainder_vec[np.newaxis].T, 1)
            X_cnv_out = np.append(X_cnv_out, agg, 1)
        # update uns dict
        uns_cnv_out["chr_pos"][current_chr] = X_cnv_out.shape[1] - (target_size[i] + 1)
    # return structures
    return X_cnv_out, uns_cnv_out
