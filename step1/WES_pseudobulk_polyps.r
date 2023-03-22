rm(list = ls())

# SET CORRECT WORKING DIRECTORY
setwd("/home/cody/git/spatial_CRC_atlas/step1/")

# make directory for outputs
dir.create("./WES_out/", showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# packages and functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(stringr)
library(data.table)
library(maftools)
library(MesKit)
library(clusterProfiler)
library(cDriver)
library(ggpubr)
library(grid)
# load the genome reference
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
# source updated maftools functions
source("../resources/WES/maftools/R/oncomatrix.R")
source("../resources/WES/maftools/R/oncoplot.R")
source("../resources/WES/maftools/R/plotMutProfile.r")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# settings and metadata
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read clinical feature (cancer subtype)
clinical <- read.table(
  "../resources/WES/clinical_labels.txt",
  header=T,
  sep="\t"
)

# read genes with pathway annotations
pways <- read.table(
  "../resources/WES/pathway_genelists.txt",
  header=T,
  sep="\t",
  stringsAsFactors=F
)

# Driver genes of CRC collected from IntOGen v.2020.2
driverGene <- as.character(
  read.table(
    system.file("extdata/", "IntOGen-DriverGenes_COREAD.tsv", package="MesKit")
  )$V1
)
# add a couple genes to get PAT06439 on our oncoplot
driverGene <- append(driverGene, "SYMPK")
driverGene <- append(driverGene, "USH2A")
# add a couple genes to get PAT71397 ROI1 on our oncoplot
driverGene <- append(driverGene, "MUC16")
# add
driverGene <- append(driverGene, "POLE")

# read in genes to subset to from MSKCC cancer hotspots
# (https://www.cancerhotspots.org/#/home)
hotspots <- unique(
  read.table("../resources/WES/cancerhotspots.txt", sep = "\t", header = T
)$Gene)

# read in genes to subset to from COSMIC 
# (https://cancer.sanger.ac.uk/cosmic/)
COSMIC.CRC.genes <- read.table(
  "../resources/WES/COSMIC_genes_colon.csv",
  sep = ",",
  header = T
) %>%
  filter(Mutated.samples > 100)  # filter to COSMICs detected in > 100 samples
# split out unique gene names
COSMIC.CRC.genes <- unique(
  data.frame(do.call('rbind', strsplit(COSMIC.CRC.genes$Gene.name,"_",fixed=T))
)$X1)

# define minimum VAF
freq <- 0.1
# define minimum read depth
depth <- 15

# define list of patient IDs
patient_ids <- c(
  "HTA11_08622_A",
  "HTA11_07663",
  "HTA11_08622_B",
  "HTA11_10711",
  "HTA11_07862",
  "HTA11_01938",
  "WD86055",
  "HTA11_06134",
  "PAT00222",
  "PAT01586",
  "PAT01587",
  "PAT05785",
  "PAT06439",
  "PAT15211",
  "PAT30884",
  "PAT33430",
  "PAT40364",
  "PAT54273",
  "PAT59460",
  "PAT59600",
  "PAT59667",
  "PAT71397",
  "PAT71662",
  "PAT73458",
  "PAT73899",
  "PAT74143"
)

tumor_colors <- c(
  "#c4a4e1", # SSL
  "#fee799", # TA
  "#7a4fa3", # MSI
  "#ffc101" # MSS
)
names(tumor_colors) <- c(
  "SSL/HP",
  "TA/TVA",
  "MSI",
  "MSS"
)
tumor_colors = list(Subtype = tumor_colors)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read master WES data object
data.polyps <- read.maf(
  paste0("../data/WES/master_depth", depth, "_maftools.maf"),
  clinicalData = clinical
)
# subset to LCM samples only
data.polyps <- subsetMaf(
  data.polyps,
  clinQuery="Tumor_ID == 'B' & Subtype %in% c('SSL/HP','TA/TVA')"
)
dim(data.polyps@data)

# read master WES data object
data.lcm <- read.maf(
  paste0("../data/WES/master_depth", depth, "_maftools.maf"),
  clinicalData = clinical
)
# subset to LCM samples only
data.lcm <- subsetMaf(
  data.lcm,
  clinQuery="Tumor_ID != 'B'"
)
dim(data.lcm@data)

# remap tumor_sample_barcode to patient names in LCM samples
data.lcm@data$Tumor_Sample_Barcode <- unlist(lapply(
  data.lcm@data$Tumor_Sample_Barcode,
  FUN = function(x){
    clinical[clinical$Tumor_Sample_Barcode==x, "Patient_ID"]
  }
))
# update clinical slot of MAF object
data.lcm@clinical.data$Tumor_Sample_Barcode <- data.lcm@clinical.data$Patient_ID
data.lcm@clinical.data <- data.lcm@clinical.data[
  !duplicated(data.lcm@clinical.data$Tumor_Sample_Barcode),
]

# merge into master 'pseudobulk' MAF
data <- merge_mafs(c(data.polyps, data.lcm))

# manipulate some data attributes so maftools doesn't choke
data@data$Tumor_Sample_Barcode <- droplevels(data@data$Tumor_Sample_Barcode)
data@variant.type.summary <- data@variant.type.summary %>%
  filter(!grepl("_ROI", Tumor_Sample_Barcode))
data@variant.type.summary$Tumor_Sample_Barcode <- droplevels(
  data@variant.type.summary$Tumor_Sample_Barcode
)
data@variant.classification.summary <- data@variant.classification.summary %>%
  filter(!grepl("_ROI", Tumor_Sample_Barcode))
data@variant.classification.summary$Tumor_Sample_Barcode <- droplevels(
  data@variant.classification.summary$Tumor_Sample_Barcode
)

data@variants.per.sample

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save oncoplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot figure with custom altered pathways
png(
  filename=paste0(
    "WES_out/pseudobulk_oncoplot_pathway_depth",
    depth,
    ".png"
  ),
  width=3000,
  height=6000,
  res=500,
  type="cairo"
)
oncoplot(
  maf=data,
  gene_mar=12,
  barcode_mar=7,
  showTumorSampleBarcodes=T,
  drawRowBar = F,
  sampleOrder = c(
    "MAP.01938_polyp",
    "MAP.07862_polyp",
    "MAP.10711_polyp",
    "WD86055",
    "PAT05785",
    "PAT06439",
    "PAT15211",
    "PAT30884",
    "PAT59600",
    "PAT59667",
    "PAT71397",
    "PAT71662",
    "PAT73899",
    "PAT74143",
    "SG00003",
    "SG00004",
    "MAP.06134_polyp",
    "MAP.07663_polyp",
    "MAP.08622_polyp_1",
    "MAP.08622_polyp_2",
    "PAT00222",
    "PAT01586",
    "PAT01587",
    "PAT33430",
    "PAT40364",
    "PAT54273",
    "PAT59460",
    "PAT73458",
    "SG00001",
    "SG00002"
  ),
  TumorSampleLabels = c(
    "HTA11_01938",
    "HTA11_07862",
    "HTA11_10711",
    "PAT71397",
    "PAT05785",
    "PAT06439",
    "PAT15211",
    "PAT30884",
    "PAT59600",
    "PAT59667",
    "PAT71397",
    "PAT71662",
    "PAT73899",
    "PAT74143",
    "SG00003",
    "SG00004",
    "HTA11_06134",
    "HTA11_07663",
    "HTA11_08622_A",
    "HTA11_08622_B",
    "PAT00222",
    "PAT01586",
    "PAT01587",
    "PAT33430",
    "PAT40364",
    "PAT54273",
    "PAT59460",
    "PAT73458",
    "SG00001",
    "SG00002"
  ),
  barcodeSrt=90,
  clinicalFeatures=c("Subtype"),
  #sortByAnnotation=T,
  legend_height=1.5,
  pathways=pways,
  anno_height = 0.2,
  #logColBar=TRUE,
  sepwd_samples=0.5,
  annotationColor=tumor_colors,
  showTitle = FALSE,
  colbar_pathway = FALSE,
  removeNonMutated = FALSE,
  #fill = TRUE,
); dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MesKit prep
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in order to pass to MesKit, edits must be made to the above MAF file:
# 1. rename column 't_ref_count' to 'Ref_allele_depth'
# 2. rename column 't_alt_count' to 'Alt_allele_depth'
t <- read.table(
  paste0("../data/WES/master_depth", depth, "_maftools.maf"),
  sep="\t",
  header=T
) %>%
  dplyr::rename(
    Ref_allele_depth=t_ref_count,
    Alt_allele_depth=t_alt_count,
    VAF=AF
  )
write.table(
  t,
  file=paste0("../data/WES/pseudobulk_depth", depth, "_maftools.maf"),
  sep="\t",
  row.names=F,
  quote=F
)
# read in MAF using MesKit
maf2 <- readMaf(
  mafFile=paste0("../data/WES/pseudobulk_depth", depth, "_maftools.maf"),
  clinicalFile="../resources/WES/clinical_labels.txt",
  refBuild = "hg19"
)

# calculate MATH score of each sample
polyp.math <- mathScore(
  maf2,
  use.tumorSampleLabel = TRUE
)
write.csv(
  polyp.math$`Bulk`,
  file=paste0("WES_out/polyp_MATHscore_depth", depth, ".csv")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build phylogenetic tree(s)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NJ.all <- getPhyloTree(
  maf2,
  patient.id = NULL,
  method = "NJ",
  min.vaf = freq
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mutational signatures analysis (COSMIC v2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trimatrix.all <- triMatrix(NJ.all, level = 3)
# Visualize the 96 trinucleodide mutational profile
plotMutSigProfile(
  trimatrix.all,
  use.tumorSampleLabel = T
)
fit.all <- fitSignatures(trimatrix.all, signaturesRef = "cosmic_v2")

png(
  filename=paste0(
    "WES_out/polyp_COSMICsignaturesfit_depth",
    depth,
    ".png"
  ),
  width=8000,
  height=6000,
  res=400,
  type="cairo"
); cowplot::plot_grid(plotlist = plotMutSigProfile(fit.all)); dev.off()

# cosine similarities between mutational profile of each group and COSMIC v2
library(ComplexHeatmap)

tmp <- data.frame(fit.all$`Bulk`$cosine.similarity)
rownames(tmp) <- unlist(lapply(rownames(tmp), FUN = function(x) {
  clinical[clinical$Tumor_Sample_Barcode == x, "Tumor_Sample_Label"]
}))
colnames(tmp) <- str_replace_all(colnames(tmp), "Signature.", "")
png(
  filename=paste0(
    "WES_out/polyp_COSMICsignatures_depth",
    depth,
    ".png"
  ),
  width=3500,
  height=1000,
  res=500,
  type="cairo"
)
ComplexHeatmap::Heatmap(
  tmp,
  name = "Cosine\nsimilarity"
); dev.off()

# write master cosine similarity matrix to file
write.csv(
  fit.all$`Bulk`$cosine.similarity,
  file=paste0("WES_out/polyp_COSMICsignatures_depth", depth, ".csv")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mutational signatures analysis (COSMIC v3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.all.v3 <- fitSignatures(trimatrix.all, signaturesRef = "exome_cosmic_v3")

png(
  filename=paste0(
    "WES_out/polyp_COSMICsignaturesv3fit_depth",
    depth,
    ".png"
  ),
  width=8000,
  height=6000,
  res=400,
  type="cairo"
); cowplot::plot_grid(plotlist = plotMutSigProfile(fit.all.v3)); dev.off()

# cosine similarities between mutational profile of each group and COSMIC v3
tmp <- data.frame(fit.all.v3$`Bulk`$cosine.similarity)
rownames(tmp) <- unlist(lapply(rownames(tmp), FUN = function(x) {
  clinical[clinical$Tumor_Sample_Barcode == x, "Tumor_Sample_Label"]
}))
png(
  filename=paste0(
    "WES_out/polyp_COSMICsignaturesv3_depth",
    depth,
    ".png"
  ),
  width=6000,
  height=1000,
  res=500,
  type="cairo"
)
ComplexHeatmap::Heatmap(
  tmp,
  name = "Cosine\nsimilarity"
); dev.off()

# write master cosine similarity matrix to file
write.csv(
  fit.all.v3$`Bulk`$cosine.similarity,
  file=paste0("WES_out/polyp_COSMICsignaturesv3_depth", depth, ".csv")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save driver mutations and TMB to files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(
  x=data@variants.per.sample,
  file="../resources/WES/pseudobulk_variants_per_sample.csv"
)

muts <- data.frame(matrix(
  ncol = 4,
  nrow = length(
    levels(genesToBarcodes(data, genes="BRAF")$BRAF$Tumor_Sample_Barcode))
  ),
  row.names = levels(
    genesToBarcodes(data, genes="BRAF")$BRAF$Tumor_Sample_Barcode
  )
)
colnames(muts) <- c("APC", "KRAS", "TP53", "BRAF")
muts$BRAF <- "F"
muts$TP53 <- "F"
muts$APC <- "F"
muts$KRAS <- "F"
muts[
  as.character(genesToBarcodes(data, genes="BRAF")$BRAF$Tumor_Sample_Barcode),
  "BRAF"] <- "T"
muts[
  as.character(genesToBarcodes(data, genes="KRAS")$KRAS$Tumor_Sample_Barcode),
  "KRAS"] <- "T"
muts[
  as.character(genesToBarcodes(data, genes="TP53")$TP53$Tumor_Sample_Barcode),
  "TP53"] <- "T"
muts[
  as.character(genesToBarcodes(data, genes="APC")$APC$Tumor_Sample_Barcode),
  "APC"] <- "T"
write.csv(x=muts, file="../resources/WES/CRC_mutations_pseudobulk.csv")
