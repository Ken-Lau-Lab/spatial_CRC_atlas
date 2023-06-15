# Step 1

## scRNA

`scRNA_cytotrace.Rmd` - R notebook for running [`CytoTRACE`](https://cytotrace.stanford.edu/) on scRNA datasets in `.h5ad` format.

`scRNA_infercnv_all.ipynb` - Python notebook for running [`infercnvpy`](https://github.com/icbi-lab/infercnvpy) on scRNA datasets in `.h5ad` format. Performed on "all" scRNA data including Broad Institute and [Chen, _et al._](https://doi.org/10.1016/j.cell.2021.11.031) cohorts.

`scRNA_infercnv_patients.ipynb` - Python notebook for running [`infercnvpy`](https://github.com/icbi-lab/infercnvpy) on scRNA datasets in `.h5ad` format. Performed on scRNA data from patient samples with matched ST, LCM-WES, and MxIF data for comparison.

## ST

`ST_cytotrace.Rmd` - R notebook for running [`CytoTRACE`](https://cytotrace.stanford.edu/) on Visium ST datasets in `.h5ad` format. __NOTE:__ in order for this notebook to run properly, the `.h5ad` files must be _dense_, i.e. `.X` cannot be a `scipy.sparse` matrix. Conversion might be necessary prior to running this notebook in the future.

`ST_MILWRM_refNMF.Rmd` - R notebook for running reference NMF (refNMF) deconvolution followed by [`MILWRM`](https://github.com/Ken-Lau-Lab/MILWRM) tissue domain detection on Visium ST datasets in `.h5ad` format.

## WES

`WES_pseudobulk_polyps.r` - R script for visualizing "pseudobulk" and bulk mutational profiling on a patient-level basis.

`WES_phylogenetics_LCM.r` - R script for building phylogenetic trees from LCM-WES data on a patient-level basis.
