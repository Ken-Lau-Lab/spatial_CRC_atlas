import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt





def normalization( dat_ct):
    """this function normalize the data so that each cell has the same 
    number of total counts as the median value of the total counts among all cells.
    The data will also be log-like transformed
    Count values will also be transformed to z-scores for each gene"""
    sc.pp.normalize_total(dat_ct) 
    dat_ct.X = np.arcsinh(dat_ct.X).copy()
    sc.pp.scale(dat_ct)
    
    return dat_ct




def save_h5ad_compresseed(dat_ls, name_ls, compression_opt = None):
    """This function save a list of anndata objs with names/paths specified in name_ls
    Files will be saved in compressed form, and the compression option is 'gzip for all files unless specified 
    The function will return the number of files saved successfully'"""
    
    if (not compression_opt):
        compression_opt = ['gzip' for i in range(len(dat_ls))]
    saved_num = 0
    for j in range(len(dat_ls)):
        try:
            dat_ls[j].write(name_ls[j], compression = compression_opt[j])
        except BaseException as err:
            print(err)
            continue
        else:
            saved_num +=1
        
    return saved_num



def clustering(dat, n_pcs = 50, n_neighbors= None,  use_highly_variable = False, highly_variable_layer = None, flavor = 'seurat', neighbor_distance_metrics='euclidean', raw_layer_name = 'raw_counts', verbose = False, norm_raw_before_arcsinh = True):
    """This function calculate pcs, and get the umap of the dataset
    dat will be edited in place. NOTE: the dataset need to have a raw count layer if want to run highly variable gene
    
    param dat: the anndata object that has transformed data
    param n_pcs: number of PCs to pass in to sc.tl.pca() function
    param n_neighbors: number of neightbors to pass in to the sc.pp.neighbors() function
    param use_highly_variable: if run PCA with only highly variable genes
    param highly_variable_layer: the count matrix layer used to run highly variable gene, if None, arcsinh used 
    param flavor: flavor for highly variable gene function, see scanpy highly variable gene function for more detail, need to change layer accordingly
    param neighbor_distance_metrics: distance type when running the neighbor function, default euclidean distance (can use cosine which is seurat's default)
    param raw_layer_name the name of the layer that has the raw count
    param use_raw_before_arcsinh: whether normalize the raw count library size before arcsinh transformation or not, default set to True 
    
    return None
    """
    
    if(n_neighbors == None):
        n_neighbors = int( np.sqrt(dat.n_obs) )
        
    if(use_highly_variable):
        if ('highly_variable' in dat.var ):
            if(verbose): print("Running PCA with hvg")
            sc.tl.pca(dat, return_info=False , use_highly_variable = True)
        else: 
            # run highly variable gene from here 
            if( raw_layer_name not in dat.layers):
                #need a raw layer, else will give error
                print(" Error: need to have/specify a raw count layer")
                return 
            if(verbose): print("Transforming data")
            dat.layers["curr_layer" ] = dat.X.copy() # save the current layer
            dat.X = dat.layers[raw_layer_name].copy() # change to the raw layer 
            add_highly_variable_gene_col(dat, layer = highly_variable_layer, flavor = flavor, use_normed= norm_raw_before_arcsinh) # run highly variable gene function
            
            if(verbose): print("Running PCA with hvg")
            sc.tl.pca(dat, return_info=False , use_highly_variable = True)
            dat.X = dat.layers["curr_layer" ].copy() # change back to current layer 
            del dat.layers["curr_layer" ] #clean up
    else:
        if(verbose): print("Running PCA without hvg")
        sc.tl.pca(dat, return_info=False, use_highly_variable = False )
    
    if(verbose): print("Calculating Neighborhood")
    sc.pp.neighbors( dat, n_neighbors= n_neighbors , n_pcs=n_pcs, metric=neighbor_distance_metrics)
    if(verbose): print("Computing UMAP")
    sc.tl.umap(dat)
    
    return 




def pct_dropout(dat):
    dat.var["pct_dropout_by_counts"] = np.array(
            (1 - (dat.X.astype(bool).sum(axis=0) / dat.n_obs)) * 100
        ).squeeze()
    
def add_highly_variable_gene_col(dat, layer = None, flavor = 'seurat', use_normed = True):
    """ This function find highly variable genes for a dataset and add the variable in-place
    @param dat: anndata object with raw count
    @param layer: the layer of X want to use, use arcsinh is not specified
    @param flavor: flavor of method, see scanpy's documentation of the highly variable gene function
    
    @return None"""
    
    if (not layer):
        if(use_normed):
            X_norm = sc.pp.normalize_total(dat, inplace=False)['X']
            dat.layers["arcsinh"] = np.arcsinh(X_norm )
        else:
            dat.layers["arcsinh"] = np.arcsinh(dat.X.copy() )
        
    sc.pp.highly_variable_genes(dat, layer="arcsinh", flavor=flavor, inplace = True)
    
    return 






def parse_rank_gene_group_dict(rank_gene_dict, fields_to_keep = ['logfoldchanges', 'pvals_adj', 'pvals', 'names'], 
                               condition = 'fetal', sort_by = 'pvals_adj', ascending = True):
    """this function parse the dictionary output by the scanpy function scanpy.tl.rank_gene_group (differential gene expression) into a pd dataframe.
    this function only parse one condition (eg. one cell type's DGE compare with another condition )
    
    @param rank_gene_dict: the dictionary output by the rank_gene_group function (stored in the anndata's uns field, default 'rank_gene_groups'
    @param field to keep: fields in the dictionary keys to keep as columns in the resultant dataframe
    @param condition: if the rank gene group function is performed over >2 groups (eg. DGE between different cell types), which condition's DGE do you want to extract from the dictionary )
    
    return: the dataframe parsed from the dictionary that contains the fields_to_keep terms as columns"""
    
    result_df = pd.DataFrame({fields_to_keep[i] : rank_gene_dict[fields_to_keep[i]][condition] for i in range(len(fields_to_keep))
              })
    result_df.sort_values(by = sort_by, inplace=True, ascending=ascending)
    
    return result_df




#copied from ~/immune_exclusion/2_immune_geen... notebook
def compute_otsu_criteria(im, th):
    """ This function calculates a threshold that separates a bimodal distribution. Copied from wikipedia otsu method page: https://en.wikipedia.org/wiki/Otsu%27s_method
    @param im: image, an nxn matrix or an 1D array
    @param th: threshold
    
    return a otsu criterion value """
    # create the thresholded image
    thresholded_im = np.zeros(im.shape)
    thresholded_im[im >= th] = 1

    # compute weights
    nb_pixels = im.size
    nb_pixels1 = np.count_nonzero(thresholded_im)
    weight1 = nb_pixels1 / nb_pixels
    weight0 = 1 - weight1

    # if one of the classes is empty, eg all pixels are below or above the threshold, that threshold will not be considered
    # in the search for the best threshold
    if weight1 == 0 or weight0 == 0:
        return np.inf

    # find all pixels belonging to each class
    val_pixels1 = im[thresholded_im == 1]
    val_pixels0 = im[thresholded_im == 0]

    # compute variance of these classes
    var1 = np.var(val_pixels1) if len(val_pixels1) > 0 else 0
    var0 = np.var(val_pixels0) if len(val_pixels0) > 0 else 0

    return weight0 * var0 + weight1 * var1


#copied from ~/immune_exclusion/2_immune_geen... notebook
def get_otsu_threshold( data, threshold_range = None, threshold_step = 0.01):
    """This function computes a best_threshold that separates 2 peaks in a bimodal distribution
    @param data: vector (array or matrix) that has the bimodal distribution
    @param threshold_range: a list of values of thresholds where the best threshold will be selected from, or could be a range/arange object. If None, numpy.arange(data min, data max, step =threshold_step ) will be used
    @param threshold_step: step for the default threshold_range, will be implemented only if threshold_range is None
    
    return a threshold  minimizing the Otsu criteria"""

    # testing all thresholds from 0 to the maximum of the image
    if(not threshold_range):
        threshold_range = np.arange(int( np.min(data) ), int( np.max(data) ) , threshold_step)
    
    criterias = [compute_otsu_criteria(data, th) for th in threshold_range]

    # best threshold is the one minimizing the Otsu criteria
    best_threshold = threshold_range[np.argmin(criterias)]
    
    return best_threshold


    
    

    
# this is a random note, but how to convert dict with different length of values to a pd dataframe?
# pd.DataFrame(dict([(col_name,pd.Series(values)) for col_name,values in my_dict.items() ]))

def dictToDf(the_dict):
    df = pd.DataFrame( dict([(col_name,pd.Series(values)) for col_name,values in the_dict.items() ]) )
    return df
    


    
# from ~/plasticity/9_fetal_atlas... notebook    
def get_up_down_dge_from_rankGenDict(adata, rank_gene_group_key, 
                                     group_by, sub_group = None, num_genes = 20, sort_by = ['pvals_adj'], ascending = [True]):
    """this function get top up and down regulated genes from each group in output by scanpy's rank_gene_group (differential gene expression) function, return a dataframe as described below
    @param adata: the anndata object
    @param rank_gene_group_key: the name of rank_gene_group output stores adata's uns field (the names input as argument to 'key_added' parameter in the rank_gene_group function )
    @param group_by: the obs column name input as the 'group_by' argument in the sc.tl.rank_gene_group function
    @param sub_group: a subset of adata.obs[group_by] list (if interested in only some group's comparison
    @param num_genes: number of top significant genes output from this function
    
    @return: a dataframe with <group_name>_up and <group_name>_down as column names and significant differential genes in each column"""
    
    
    rank_gene_dict = adata.uns[rank_gene_group_key]
    if(not sub_group):
        sub_group = adata.obs[group_by].unique()
    
    all_df = pd.DataFrame()
    for g in sub_group:
        dge_df = parse_rank_gene_group_dict(adata.uns[rank_gene_group_key], fields_to_keep = ['logfoldchanges', 'pvals_adj', 'pvals', 'names', 'scores'], 
                               condition = g, sort_by = sort_by, ascending = ascending)
        g_up = dge_df[dge_df["logfoldchanges"] > 0].sort_values(by = sort_by, ascending = ascending)
        g_down = dge_df[dge_df["logfoldchanges"] < 0].sort_values(by = sort_by, ascending = ascending)
        
        all_df[f"{g}_up"] = g_up["names"][0:num_genes].values
        all_df[f"{g}_down"] = g_down["names"][0:num_genes].values
        
    return all_df


        