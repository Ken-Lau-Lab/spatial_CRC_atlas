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
library(BSgenome.Hsapiens.UCSC.hg19)  # TODO: hg38
library(org.Hs.eg.db)
# source updated maftools functions
source("../resources/WES/maftools/R/oncomatrix.R")
source("../resources/WES/maftools/R/oncoplot.R")
source("../resources/WES/maftools/R/plotMutProfile.r")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# settings and metadata
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read clinical feature (cancer subtype)
clinical <- read.table("../resources/WES/clinical_labels.txt", header=T, sep="\t")

# read genes with pathway annotations
pways <- read.table(
  "../resources/WES/pathway_genelists.txt",
  header=T,
  sep="\t",
  stringsAsFactors=F
)

# read genes with pathway annotations
pways_short <- read.table(
  "../resources/WES/pathway_genelists_curated_short_Feb23.txt",
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

# define list of tumor classes
tumor_colors <- c(
  "#c4a4e1", # SSL
  "#c4a4e1", # SSL
  "#c4a4e1", # SSL
  "#fee799", # TA
  "#fee799", # TA
  "#fee799", # TVA
  "#fee799", # TVA
  "#c4a4e1",  # HP
  "#7a4fa3", # MSI
  "#7a4fa3", # MSI
  "#7a4fa3", # MSI
  "#ffc101", # MSS
  "#ffc101", # MSS
  "#ffc101", # MSS
  "#ffc101", # MSS
  "#7a4fa3", # MSI
  "#7a4fa3", # MSI
  "#7a4fa3", # MSI
  "#7a4fa3", # MSI
  "#ffc101", # MSS
  "#ffc101", # MSS
  "#ffc101", # MSS
  "#ffc101", # MSS
  "#7a4fa3", # MSI
  "#ffc101",# MSS
  "#ffc101" # MSS
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
# Read all patients
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
data.polyps@variants.per.sample

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
data.lcm@variants.per.sample
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

data <- merge_mafs(c(data.polyps, data.lcm))
data@data$Tumor_Sample_Barcode <- droplevels(data@data$Tumor_Sample_Barcode)
data@variant.type.summary$Tumor_Sample_Barcode <- droplevels(data@variant.type.summary$Tumor_Sample_Barcode)
data@variant.type.summary <- data@variant.type.summary %>%
  filter(!grepl("_ROI", Tumor_Sample_Barcode))
data@variant.classification.summary$Tumor_Sample_Barcode <- droplevels(data@variant.classification.summary$Tumor_Sample_Barcode)
data@variant.classification.summary <- data@variant.classification.summary %>%
  filter(!grepl("_ROI", Tumor_Sample_Barcode))

data@variants.per.sample

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save oncoplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png(
  filename=paste0(
    "WES_out/polyp_oncoplot_nopathway_depth",
    depth,
    ".png"
  ),
  width=3000,
  height=3000,
  res=300,
  type="cairo"
)
oncoplot(
  maf=data,
  #genes=driverGene,
  top=50,
  drawRowBar=T,
  gene_mar=7,
  barcode_mar=7,
  showTumorSampleBarcodes=T,
  #TumorSampleLabels=c(
  #  "PAT15211",
  #  "PAT05785",
  #  "HTA11_08622_B",
  #  "HTA11_10711",
  #  "HTA11_07862",
  #  "HTA11_01938",
  #  "WD86055",
  #  "HTA11_06134"
  #),
  barcodeSrt=45,
  clinicalFeatures=c("Subtype"),
  sortByAnnotation=T,
  legend_height=0.7,
); dev.off()

# plot figure with custom altered pathways
png(
  filename=paste0(
    "WES_out/polyp_oncoplot_pathway_depth",
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

# plot shortened figure with custom altered pathways
png(
  filename=paste0(
    "WES_out/polyp_oncoplot_pathway_short_depth",
    depth,
    ".png"
  ),
  width=3000,
  height=3000,
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
  pathways=pways_short,
  anno_height = 0.2,
  #logColBar=TRUE,
  sepwd_samples=0.5,
  annotationColor=tumor_colors,
  showTitle = FALSE,
  colbar_pathway = F,
  removeNonMutated = FALSE,
  #fill = TRUE,
); dev.off()

# plot figure with top 3 automatically-determined altered pathways
png(
  filename=paste0(
    "WES_out/polyp_oncoplot_autopathway_depth",
    depth,
    ".png"
  ),
  width=3000,
  height=3000,
  res=300,
  type="cairo"
)
oncoplot(
  maf=data,
  gene_mar=12,
  barcode_mar=7,
  showTumorSampleBarcodes=T,
  barcodeSrt=45,
  clinicalFeatures=c("Subtype"),
  sortByAnnotation=T,
  legend_height=0.5,
  pathways="auto",
  colbar_pathway=T,
); dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# maftools analyses, cont'd.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create matrix of mutation counts for each sample and write to .csv file
mat <- mutCountMatrix(data)
fwrite(
  mat,
  file=paste0(
    "WES_out/polyp_mutCountMatrix_depth",
    depth,
    ".csv"
  ),
  row.names=TRUE
)

# plot oncogenic pathways affected across samples
pathways <- OncogenicPathways(data)
dev.copy(
  png,
  paste0(
    "WES_out/polyp_pathways_depth",
    depth,
    ".png"
  ),
  res=100
); dev.off()

#exclusive/co-occurance event analysis on top 10 mutated genes. 
interactions <- somaticInteractions(
  maf = data,
  top = 40,
  pvalue = c(0.05, 0.1)
)
dev.copy(
  png,
  paste0(
    "WES_out/polyp_interactions_depth",
    depth,
    ".png"
  ),
  res=400,
  width=2000,
  height=2000,
); dev.off()

# clinical enrichment
ce = clinicalEnrichment(maf = data, clinicalFeature = "Subtype")
ce$groupwise_comparision[p_value < 0.05]
plotEnrichmentResults(
  enrich_res = ce,
  pVal = 0.05,
  geneFontSize = 0.5,
  annoFontSize = 0.6
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MesKit prep
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write MAF file for MesKit analysis
maftools::write.mafSummary(
  data,
  basename=paste0("../data/WES/polyp_depth", depth)
)
# in order to pass to MesKit, edits must be made to the above MAF file:
# 1. rename column 't_ref_count' to 'Ref_allele_depth'
# 2. rename column 't_alt_count' to 'Alt_allele_depth'
t <- read.table(
  paste0("../data/WES/polyp_depth", depth, "_maftools.maf"),
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
  file=paste0("../data/WES/polyp_depth", depth, "_maftools.maf"),
  sep="\t",
  row.names=F,
  quote=F
)
# read in MAF using MesKit
maf2 <- readMaf(
  mafFile=paste0("../data/WES/polyp_depth", depth, "_maftools.maf"),
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
# cluster mutations on VAF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mut.cluster <- mutCluster(
  maf2, patient.id = "Bulk", use.ccf = F, use.tumorSampleLabel = TRUE
)
png(
  filename=paste0(
    "WES_out/polyp_MATHclusters_depth",
    depth,
    ".png"
  ),
  width=6000,
  height=4000,
  res=400,
  type="cairo"
); cowplot::plot_grid(plotlist = mut.cluster$cluster.plot); dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MesKit oncoplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hotspot_oncoplot <- grid.grabExpr(plotMutProfile_CH(
  maf2,
  class = "SP",
  geneList = hotspots,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 60,
  removeEmptyCols = F
))
COSMIC_oncoplot <- grid.grabExpr(plotMutProfile_CH(
  maf2,
  class = "SP",
  geneList = COSMIC.CRC.genes,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 60,
  removeEmptyCols = F
))
driverGene_oncoplot <- grid.grabExpr(plotMutProfile_CH(
  maf2,
  class = "SP",
  geneList = driverGene,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 32,
  removeEmptyCols = F
))
png(
  filename=paste0("WES_out/polyp_oncoplots_depth", depth, ".png"),
  width=9000,
  height=4500,
  res=500,
  type="cairo"
)
ggarrange(
  plotlist = list(driverGene_oncoplot, hotspot_oncoplot, COSMIC_oncoplot),
  labels = c("Drivers", "Hotspots (MSKCC)", "COSMIC (CRC)"),
  hjust = 0.5,
  vjust = 0.05,
  ncol = 3,
  widths = c(1, 1, 1)
) + theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 45, unit = "pt"))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build phylogenetic tree(s)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NJ.all <- getPhyloTree(
  maf2,
  patient.id = NULL,
  method = "NJ",
  min.vaf = freq
)

# combine plots
tree_polyps <- plotPhyloTree(
  NJ.all,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL
)

onco <- grid.grabExpr(plotMutProfile_CH(
  maf2,
  class = "SP",
  geneList = driverGene,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 32,
  removeEmptyCols = F
))
# combine into master plot
png(
  filename=paste0(
    "WES_out/polyp_trees_depth",
    depth,
    ".png"
  ),
  width=5000,
  height=4500,
  res=500,
  type="cairo"
)
ggarrange(
  plotlist = list(tree_polyps, onco),
  labels = c("A", "B"),
  ncol = 2,
  widths = c(1, 1)
)
dev.off()

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
# plotting oncoplot with driverGene
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# split out gene id
mut.class <- classifyMut(maf2, class="SP")
mut.class <- cbind(
  mut.class,
  data.frame(
    do.call(
      "rbind",
      strsplit(mut.class$Mut_ID,":",fixed=T)
    )
  )
)
mut.class$gene_id <- mut.class$X1

# define list of tumor classes
tumor_class_colors <- c(
  "#c4a4e1", # SSL
  "#c4a4e1", # HP
  "#fee799", # TA
  "#fee799" # TVA
)
names(tumor_class_colors) <- c(
  "SSL",
  "HP",
  "TA",
  "TVA"
)
tumor_class_colors = list(Subtype = tumor_class_colors)

# plot fancy oncoplot
png(
  filename=paste0(
    "WES_out/polyp_oncoplot_drivers_depth",
    depth,
    ".png"
  ),
  width=2000,
  height=2000,
  res=300,
  type="cairo"
)
oncoplot(
  maf=data,
  genes=intersect(mut.class$gene_id, driverGene),
  drawRowBar=T,
  gene_mar=7,
  barcode_mar=5,
  showTumorSampleBarcodes=T,
  TumorSampleLabels=c(
    "MAP08622_A",
    "MAP07663",
    "MAP08622_B",
    "MAP10711",
    "MAP07862",
    "MAP01938",
    "WD86055",
    "MAP06134"
  ),
  barcodeSrt=45,
  clinicalFeatures=c("Subtype"),
  annotationColor=tumor_class_colors,
  sortByAnnotation=T,
  legend_height=0.7,
); dev.off()

# plot fancy oncoplot
png(
  filename=paste0(
    "WES_out/polyp_oncoplot_hotspots_depth",
    depth,
    ".png"
  ),
  width=2000,
  height=4000,
  res=300,
  type="cairo"
)
oncoplot(
  maf=data,
  genes=intersect(mut.class$gene_id, hotspots),
  drawRowBar=T,
  gene_mar=7,
  barcode_mar=5,
  showTumorSampleBarcodes=T,
  TumorSampleLabels=c(
    "MAP08622_A",
    "MAP07663",
    "MAP08622_B",
    "MAP10711",
    "MAP07862",
    "MAP01938",
    "WD86055",
    "MAP06134"
  ),
  barcodeSrt=45,
  clinicalFeatures=c("Subtype"),
  annotationColor=tumor_class_colors,
  sortByAnnotation=T,
  legend_height=0.7,
); dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(
  x=data@variants.per.sample,
  file="../resources/WES/pseudobulk_variants_per_sample.csv"
)

muts <- data.frame(matrix(ncol = 4, nrow = length(levels(genesToBarcodes(data, genes="BRAF")$BRAF$Tumor_Sample_Barcode))), row.names = levels(genesToBarcodes(data, genes="BRAF")$BRAF$Tumor_Sample_Barcode))
colnames(muts) <- c("APC", "KRAS", "TP53", "BRAF")
muts$BRAF <- "F"
muts$TP53 <- "F"
muts$APC <- "F"
muts$KRAS <- "F"
muts[as.character(genesToBarcodes(data, genes="BRAF")$BRAF$Tumor_Sample_Barcode), "BRAF"] <- "T"
muts[as.character(genesToBarcodes(data, genes="KRAS")$KRAS$Tumor_Sample_Barcode), "KRAS"] <- "T"
muts[as.character(genesToBarcodes(data, genes="TP53")$TP53$Tumor_Sample_Barcode), "TP53"] <- "T"
muts[as.character(genesToBarcodes(data, genes="APC")$APC$Tumor_Sample_Barcode), "APC"] <- "T"
write.csv(x=muts, file="../resources/WES/CRC_mutations_pseudobulk.csv")
