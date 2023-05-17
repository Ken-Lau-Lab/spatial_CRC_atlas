# Step 2

## scRNA

`scRNA_gene_signatures.ipynb` - Python notebook for scoring literature-derived gene signatures on scRNA datasets in `.h5ad` format.

`scRNA_infercnv_patients_curation.ipynb` - Python notebook for "curating" [`infercnvpy`](https://github.com/icbi-lab/infercnvpy) results on scRNA datasets in `.h5ad` format. Performed on scRNA data from patient samples with matched ST, LCM-WES, and MxIF data for comparison. "Curates" CNV Leiden clusters into interpretable "CNV Clones" for downstream analysis and PPT ordering.

## ST

`ST_compile_master_AnnData.ipynb` - Python notebook for compiling prior analyses (CytoTRACE, gene signature scores, refNMF cell state scores, MILWRM, LCM ROI masks, manual pathology annotations) on Visium ST datasets into  `_master.h5ad` files per specimen for downstream analyses.
