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

patient_ids <- c(
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
  "PAT74143",
  "SG00001",
  "SG00002",
  "SG00003",
  "SG00004"
)

# define list of tumor classes
tumor_colors <- c(
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
  "#ffc101", # MSS
  "#ffc101", # MSS
  "#7a4fa3", # MSI
  "#7a4fa3", # MSI
  "#ffc101", # MSS
  "#ffc101" # MSS
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read all patients
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read master WES data object
data <- read.maf(
  paste0("../data/WES/master_depth", depth, "_maftools.maf"),
  clinicalData = clinical
)
# subset to LCM samples only
data <- subsetMaf(
  data,
  clinQuery="Tumor_ID != 'B'"
)
dim(data@data)
data@variants.per.sample

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save oncoplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png(
  filename=paste0(
    "WES_out/LCM_oncoplot_nopathway_depth",
    depth,
    ".png"
  ),
  width=6500,
  height=5000,
  res=400,
  type="cairo"
)
oncoplot(
  maf=data,
  top=50,
  drawRowBar=T,
  gene_mar=7,
  barcode_mar=7,
  showTumorSampleBarcodes=T,
  barcodeSrt=45,
  clinicalFeatures=c("Patient_ID","Subtype"),
  sortByAnnotation=T,
  legend_height=0.7,
); dev.off()

# plot figure with custom altered pathways
png(
  filename=paste0(
    "WES_out/LCM_oncoplot_pathway_depth",
    depth,
    ".png"
  ),
  width=6500,
  height=8000,
  res=400,
  type="cairo"
)
oncoplot(
  maf=data,
  gene_mar=15,
  barcode_mar=6,
  showTumorSampleBarcodes=T,
  barcodeSrt=45,
  clinicalFeatures=c("Patient_ID","Subtype"),
  sortByAnnotation=T,
  legend_height=1,
  pathways=pways,
  colbar_pathway=T,
); dev.off()

# plot figure with top 3 automatically-determined altered pathways
png(
  filename=paste0(
    "WES_out/LCM_oncoplot_autopathway_depth",
    depth,
    ".png"
  ),
  width=6500,
  height=6000,
  res=400,
  type="cairo"
)
oncoplot(
  maf=data,
  gene_mar=12,
  barcode_mar=6,
  showTumorSampleBarcodes=T,
  barcodeSrt=45,
  clinicalFeatures=c("Patient_ID","Subtype"),
  sortByAnnotation=T,
  legend_height=0.5,
  pathways="auto",
  colbar_pathway=T,
); dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MesKit prep
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write MAF file for MesKit analysis
maftools::write.mafSummary(
  data,
  basename=paste0("../data/WES/LCM_depth", depth)
)
# in order to pass to MesKit, edits must be made to the above MAF file:
# 1. rename column 't_ref_count' to 'Ref_allele_depth'
# 2. rename column 't_alt_count' to 'Alt_allele_depth'
t <- read.table(
  paste0("../data/WES/LCM_depth", depth, "_maftools.maf"),
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
  file=paste0("../data/WES/LCM_depth", depth, "_maftools.maf"),
  sep="\t",
  row.names=F,
  quote=F
)
# read in MAF using MesKit
maf2 <- readMaf(
  mafFile=paste0("../data/WES/LCM_depth", depth, "_maftools.maf"),
  clinicalFile="../resources/WES/clinical_labels.txt",
  refBuild = "hg19"
)

# classify shared, private, public mutations
for (pat in patient_ids) {
  print(pat)
  assign(
    paste0(pat, ".mut.class"),
    classifyMut(maf2, class="SP", patient.id=pat)
  )
}

# calculate MATH score of each sample
LCM.math <- mathScore(
  maf2,
  use.tumorSampleLabel = TRUE
)

# save to individual files per patient
for (pat in patient_ids) {
  print(pat)
  write.csv(
    eval(parse(text = paste0("LCM.math$`", pat, "`"))),
    file=paste0("WES_out/", pat, "_MATHscore_depth", depth, ".csv")
  )
}

# save combined MATH scores to master file
write.csv(
  rbind(
    LCM.math$`PAT71397`,LCM.math$`PAT71662`,LCM.math$`PAT73458`,
    LCM.math$`PAT73899`,LCM.math$`PAT74143`,LCM.math$`PAT05785`,
    LCM.math$`PAT01587`,LCM.math$`PAT00222`,LCM.math$`PAT40364`,
    LCM.math$`PAT06439`,LCM.math$`PAT15211`,LCM.math$`PAT01586`,
    LCM.math$`PAT54273`,LCM.math$`PAT33430`,LCM.math$`PAT59460`,
    LCM.math$`PAT59600`,LCM.math$`PAT59667`,LCM.math$`PAT30884`,
    LCM.math$`SG00001`,LCM.math$`SG00002`,LCM.math$`SG00003`,
    LCM.math$`SG00004`
  ),
  file=paste0("WES_out/LCM_MATHscore_depth", depth, ".csv")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cluster mutations on VAF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (pat in patient_ids) {
  print(pat)
  mut.cluster <- mutCluster(
    maf2, patient.id = pat, use.ccf = F, use.tumorSampleLabel = TRUE
  )
  png(
    filename=paste0(
      "WES_out/",
      pat,
      "_MATHclusters_depth",
      depth,
      ".png"
    ),
    width=6000,
    height=4000,
    res=400,
    type="cairo"
  ); print(cowplot::plot_grid(plotlist = mut.cluster$cluster.plot)); dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MesKit oncoplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
driverGene_oncoplot <- grid.grabExpr(plotMutProfile(
  maf2,
  class = "SP",
  geneList = driverGene,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 60
))
hotspot_oncoplot <- grid.grabExpr(plotMutProfile(
  maf2,
  class = "SP",
  geneList = hotspots,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 60
))
COSMIC_oncoplot <- grid.grabExpr(plotMutProfile(
  maf2,
  class = "SP",
  geneList = COSMIC.CRC.genes,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 60
))
png(
  filename=paste0("WES_out/LCM_oncoplots_depth", depth, ".png"),
  width=16000,
  height=5000,
  res=500,
  type="cairo"
)
ggarrange(
  plotlist = list(driverGene_oncoplot, hotspot_oncoplot, COSMIC_oncoplot),
  labels = c("Drivers", "Hotspots (MSKCC)", "COSMIC (CRC)"),
  ncol = 3,
  hjust = 0.5,
  vjust = 0.05,
  widths = c(1, 1, 1)
)+ theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 45, unit = "pt"))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build phylogenetic tree(s)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get master phylo tree
NJ.all <- getPhyloTree(
  maf2,
  patient.id = NULL,
  method = "NJ",
  min.vaf = freq
)

# get phylo tree per patient
for (pat in patient_ids) {
  print(pat)
  assign(
    paste0("NJ.", pat),
    getPhyloTree(
      maf2,
      patient.id = pat,
      method = "NJ",
      min.vaf = freq
    )
  )
}

# combine plots
tree_PAT71397 <- plotPhyloTree(
  NJ.PAT71397,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 300,
)
tree_PAT71397$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT71662 <- plotPhyloTree(
  NJ.PAT71662,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT71662$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT73458 <- plotPhyloTree(
  NJ.PAT73458,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 700,
)
tree_PAT73458$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT73899 <- plotPhyloTree(
  NJ.PAT73899,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT73899$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT74143 <- plotPhyloTree(
  NJ.PAT74143,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT74143$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT05785<- plotPhyloTree(
  NJ.PAT05785,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT05785$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT01587 <- plotPhyloTree(
  NJ.PAT01587,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 50,
  scale.bar.x = -10,
)
tree_PAT01587$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT00222 <- plotPhyloTree(
  NJ.PAT00222,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = -500,
)
tree_PAT00222$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT40364 <- plotPhyloTree(
  NJ.PAT40364,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 300,
)
tree_PAT40364$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT06439 <- plotPhyloTree(
  NJ.PAT06439,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT06439$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT15211 <- plotPhyloTree(
  NJ.PAT15211,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 1000,
  scale.bar.x = -600,
)
tree_PAT15211$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT01586 <- plotPhyloTree(
  NJ.PAT01586,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT01586$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT54273 <- plotPhyloTree(
  NJ.PAT54273,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = -400,
)
tree_PAT54273$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT33430 <- plotPhyloTree(
  NJ.PAT33430,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT33430$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT59460 <- plotPhyloTree(
  NJ.PAT59460,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = -40,
  scale.bar.x = -40,
)
tree_PAT59460$theme$plot.title$colour <- "#7a4fa3" # MSI, SSL="#c4a4e1"

tree_PAT59600 <- plotPhyloTree(
  NJ.PAT59600,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = -40,
  scale.bar.x = -8,
)
tree_PAT59600$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT59667 <- plotPhyloTree(
  NJ.PAT59667,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = -40,
)
tree_PAT59667$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_PAT30884 <- plotPhyloTree(
  NJ.PAT30884,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_PAT30884$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_SG00001 <- plotPhyloTree(
  NJ.SG00001,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 1200,
)
tree_SG00001$theme$plot.title$colour <- "#7a4fa3"  # MSS, TA/TVA="#fee799"

tree_SG00002 <- plotPhyloTree(
  NJ.SG00002,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 800,
  scale.bar.x = -200,
)
tree_SG00002$theme$plot.title$colour <- "#7a4fa3"  # MSS, TA/TVA="#fee799"

tree_SG00003 <- plotPhyloTree(
  NJ.SG00003,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 0,
)
tree_SG00003$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

tree_SG00004 <- plotPhyloTree(
  NJ.SG00004,
  use.tumorSampleLabel = TRUE,
  show.scale.bar = TRUE,
  show.bootstrap = FALSE,
  branchCol =NULL,
  scale.bar.y = 120,
)
tree_SG00004$theme$plot.title$colour <- "#ffc101"  # MSS, TA/TVA="#fee799"

trees <- ggarrange(
  plotlist = list(
    tree_PAT05785, tree_PAT06439, tree_PAT15211, tree_PAT30884, tree_PAT59600,
    tree_PAT59667, tree_PAT71397, tree_PAT71662, tree_PAT73899, tree_PAT74143,
    tree_SG00003, tree_SG00004,
    tree_PAT00222, tree_PAT01586, tree_PAT01587, tree_PAT33430, tree_PAT40364,
    tree_PAT54273, tree_PAT59460, tree_PAT73458,  tree_SG00001, tree_SG00002
  ),
  labels = c("A", "B", "C", "D", "E", "F", "G", "H",
             "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R","S","T","U","V"),
  ncol = 6,
  widths = c(1, 1, 1, 1, 1, 1),
  nrow = 4,
  heights = c(1, 1, 1)
)

onco <- grid.grabExpr(plotMutProfile(
  maf2,
  patient.id = c(
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
  class = "SP",
  geneList = driverGene,
  use.tumorSampleLabel = TRUE,
  topGenesCount = 60
))

# combine into master plot
png(
  filename=paste0(
    "WES_out/LCM_trees_depth",
    depth,
    ".png"
  ),
  width=15000,
  height=5000,
  res=500,
  type="cairo"
)
ggarrange(
  plotlist = list(trees, onco),
  labels = c("", "W"),
  ncol = 2,
  widths = c(6, 4)
)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mutational signatures analysis (COSMIC v2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trimatrix.all <- triMatrix(NJ.all, level = 3)
# Visualize the 96 trinucleodide mutational profile
plotMutSigProfile(trimatrix.all)
fit.all <- fitSignatures(trimatrix.all, signaturesRef = "cosmic_v2")
#fit.all.nat <- fitSignatures(trimatrix.all, signaturesRef = "nature2013")

# plot v2 COSMIC signatures fits
for (pat in patient_ids) {
  print(pat)
  # plot COSMIC signatures fits
  png(
    filename=paste0(
      "WES_out/",
      pat,
      "_COSMICsignaturesfit_depth",
      depth,
      ".png"
    ),
    width=8000,
    height=6000,
    res=400,
    type="cairo"
  ); print(cowplot::plot_grid(
    plotlist = eval(parse(
      text = paste0("plotMutSigProfile(fit.all)$`", pat, "`")
    ))
  )); dev.off()
}

# cosine similarities between mutational profile of each group and COSMIC v2
library(ComplexHeatmap)
for (i in 1:length(patient_ids)) {
  pat <- patient_ids[i]
  print(pat)
  
  tmp <- eval(parse(text = paste0(
    "data.frame(fit.all$`", pat, "`$cosine.similarity)"
  )))
  rownames(tmp) <- unlist(lapply(rownames(tmp), FUN = function(x) {
    clinical[clinical$Tumor_Sample_Barcode == x, "Tumor_Sample_Label"]
  }))
  colnames(tmp) <- str_replace_all(colnames(tmp), "Signature.", "")
  assign(
    paste0("sig.heatmap.", pat),
    grid.grabExpr(draw(ComplexHeatmap::Heatmap(
      tmp,
      name = "Cosine\nsimilarity",
      column_title = pat,
      column_title_gp = gpar(col = tumor_colors[i])
    ))))
}

png(
  filename=paste0(
    "WES_out/LCM_COSMICsignatures_depth",
    depth,
    ".png"
  ),
  width=9000,
  height=6000,
  res=500,
  type="cairo"
)
ggarrange(
  plotlist = list(
    sig.heatmap.PAT00222,
    sig.heatmap.PAT01586,
    sig.heatmap.PAT01587,
    sig.heatmap.PAT05785,
    sig.heatmap.PAT06439,
    sig.heatmap.PAT15211,
    sig.heatmap.PAT30884,
    sig.heatmap.PAT33430,
    sig.heatmap.PAT40364,
    sig.heatmap.PAT54273,
    sig.heatmap.PAT59460,
    sig.heatmap.PAT59600,
    sig.heatmap.PAT59667,
    sig.heatmap.PAT71397,
    sig.heatmap.PAT71662,
    sig.heatmap.PAT73458,
    sig.heatmap.PAT73899,
    sig.heatmap.PAT74143,
    sig.heatmap.SG00001,
    sig.heatmap.SG00002,
    sig.heatmap.SG00003,
    sig.heatmap.SG00004
  ),
  labels = c(
    "T",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    ""
  ),
  hjust = 0.5,
  vjust = 0.05,
  ncol = 4,
  nrow = 6,
  widths = c(1, 1, 1),
  heights = c(1, 1, 1, 1, 1, 1)
) + theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 45, unit = "pt"))
dev.off()

# write master cosine similarity matrix to file
write.csv(
  rbind(
    fit.all$`PAT00222`$cosine.similarity,
    fit.all$`PAT01586`$cosine.similarity,
    fit.all$`PAT01587`$cosine.similarity,
    fit.all$`PAT05785`$cosine.similarity,
    fit.all$`PAT06439`$cosine.similarity,
    fit.all$`PAT15211`$cosine.similarity,
    fit.all$`PAT30884`$cosine.similarity,
    fit.all$`PAT33430`$cosine.similarity,
    fit.all$`PAT40364`$cosine.similarity,
    fit.all$`PAT54273`$cosine.similarity,
    fit.all$`PAT59460`$cosine.similarity,
    fit.all$`PAT59600`$cosine.similarity,
    fit.all$`PAT59667`$cosine.similarity,
    fit.all$`PAT71397`$cosine.similarity,
    fit.all$`PAT71662`$cosine.similarity,
    fit.all$`PAT73458`$cosine.similarity,
    fit.all$`PAT73899`$cosine.similarity,
    fit.all$`PAT74143`$cosine.similarity,
    fit.all$`SG00001`$cosine.similarity,
    fit.all$`SG00002`$cosine.similarity,
    fit.all$`SG00003`$cosine.similarity,
    fit.all$`SG00004`$cosine.similarity
  ),
  file=paste0("WES_out/LCM_COSMICsignatures_depth", depth, ".csv")
)

# make plot per patient with tree
for (pat in patient_ids) {
  print(pat)
  png(
    filename=paste0(
      "WES_out/",
      pat,
      "_tree_depth",
      depth,
      ".png"
    ),
    width=1800,
    height=1800,
    res=700,
    type="cairo"
  )
  print(get(paste0("tree_", pat)))
  dev.off()
}

# make joint plots per patient with tree and oncoplot
for (pat in patient_ids) {
  print(pat)
  # split out gene id
  tmp.mut.class <- cbind(
    get(paste0(pat, ".mut.class")),
    data.frame(
      do.call(
        "rbind",
        strsplit(get(paste0(pat, ".mut.class"))$Mut_ID,":",fixed=T)
      )
    )
  )
  tmp.mut.class$gene_id <- tmp.mut.class$X1
  # oncoplot for panel C
  onco <- grid.grabExpr(plotMutProfile(
    maf2,
    patient.id=pat,
    class = "SP",
    geneList = driverGene,
    use.tumorSampleLabel = TRUE,
    topGenesCount = length(intersect(tmp.mut.class$gene_id, driverGene))
  ))
  # plot phylo tree and oncoplot for panels A and B
  tree.plus.sig <- ggarrange(
    plotlist = list(
      onco,
      get(paste0("tree_", pat))
    ),
    labels = c(
      "A",
      "B"
    ),
    hjust = 0.5,
    vjust = 0.05,
    ncol = 2,
    nrow = 1,
    widths = c(1.5,1)
  ) + theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 45, unit = "pt"))
  # combine into master plot
  png(
    filename=paste0(
      "WES_out/",
      pat,
      "_tree_oncoplot_depth",
      depth,
      ".png"
    ),
    width=4200,
    height=2200,
    res=490,
    type="cairo"
  )
  print(tree.plus.sig)
  dev.off()
}

# make joint plots per patient with tree, COSMIC sigs, and oncoplot
for (pat in patient_ids) {
  print(pat)
  # split out gene id
  tmp.mut.class <- cbind(
    get(paste0(pat, ".mut.class")),
    data.frame(
      do.call(
        "rbind",
        strsplit(get(paste0(pat, ".mut.class"))$Mut_ID,":",fixed=T)
      )
    )
  )
  tmp.mut.class$gene_id <- tmp.mut.class$X1
  # plot phylo tree and COSMIC sigs for panels A and B
  tree.plus.sig <- ggarrange(
    plotlist = list(
      get(paste0("tree_", pat)),
      get(paste0("sig.heatmap.", pat))
    ),
    labels = c(
      "A",
      "B"
    ),
    hjust = 0.5,
    vjust = 0.05,
    ncol = 1,
    nrow = 2,
    heights = c(1, 1)
  ) + theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 45, unit = "pt"))
  # oncoplot for panel C
  onco <- grid.grabExpr(plotMutProfile(
    maf2,
    patient.id=pat,
    class = "SP",
    geneList = driverGene,
    use.tumorSampleLabel = TRUE,
    topGenesCount = length(intersect(tmp.mut.class$gene_id, driverGene))
  ))
  # combine into master plot
  png(
    filename=paste0(
      "WES_out/",
      pat,
      "_tree_COSMICsignatures_depth",
      depth,
      ".png"
    ),
    width=6000,
    height=3500,
    res=500,
    type="cairo"
  )
  print(ggarrange(
    plotlist = list(tree.plus.sig, onco),
    labels = c("", "C"),
    ncol = 2,
    widths = c(6,4)
  ))
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mutational signatures analysis (COSMIC v3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.all.v3 <- fitSignatures(trimatrix.all, signaturesRef = "exome_cosmic_v3")

# now plot v3 COSMIC signatures fits
for (pat in patient_ids) {
  print(pat)
  # plot COSMIC signatures fits
  png(
    filename=paste0(
      "WES_out/",
      pat,
      "_COSMICsignaturesv3fit_depth",
      depth,
      ".png"
    ),
    width=8000,
    height=6000,
    res=400,
    type="cairo"
  ); print(cowplot::plot_grid(
    plotlist = eval(parse(
      text = paste0("plotMutSigProfile(fit.all.v3)$`", pat, "`")
    ))
  )); dev.off()
}

# cosine similarities between mutational profile of each group and COSMICv3
library(ComplexHeatmap)
for (i in 1:length(patient_ids)) {
  pat <- patient_ids[i]
  print(pat)
  
  tmp <- eval(parse(text = paste0(
    "data.frame(fit.all.v3$`", pat, "`$cosine.similarity)"
  )))
  rownames(tmp) <- unlist(lapply(rownames(tmp), FUN = function(x) {
    clinical[clinical$Tumor_Sample_Barcode == x, "Tumor_Sample_Label"]
  }))
  colnames(tmp) <- str_replace_all(colnames(tmp), "Signature.", "")
  assign(
    paste0("sig.heatmap.v3.", pat),
    grid.grabExpr(draw(ComplexHeatmap::Heatmap(
      tmp,
      name = "Cosine\nsimilarity",
      column_title = pat,
      column_title_gp = gpar(col = tumor_colors[i])
    ))))
}

png(
  filename=paste0(
    "WES_out/LCM_COSMICsignaturesv3_depth",
    depth,
    ".png"
  ),
  width=17000,
  height=8000,
  res=500,
  type="cairo"
)
ggarrange(
  plotlist = list(
    sig.heatmap.v3.PAT00222,
    sig.heatmap.v3.PAT01586,
    sig.heatmap.v3.PAT01587,
    sig.heatmap.v3.PAT05785,
    sig.heatmap.v3.PAT06439,
    sig.heatmap.v3.PAT15211,
    sig.heatmap.v3.PAT30884,
    sig.heatmap.v3.PAT33430,
    sig.heatmap.v3.PAT40364,
    sig.heatmap.v3.PAT54273,
    sig.heatmap.v3.PAT59460,
    sig.heatmap.v3.PAT59600,
    sig.heatmap.v3.PAT59667,
    sig.heatmap.v3.PAT71397,
    sig.heatmap.v3.PAT71662,
    sig.heatmap.v3.PAT73458,
    sig.heatmap.v3.PAT73899,
    sig.heatmap.v3.PAT74143
  ),
  labels = c(
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
    "R"
  ),
  hjust = 0.5,
  vjust = 0.05,
  ncol = 3,
  nrow = 6,
  widths = c(1, 1, 1),
  heights = c(1, 1, 1, 1, 1, 1)
) + theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 45, unit = "pt"))
dev.off()

# write master cosine similarity matrix to file
write.csv(
  rbind(
    fit.all.v3$`PAT00222`$cosine.similarity,
    fit.all.v3$`PAT01586`$cosine.similarity,
    fit.all.v3$`PAT01587`$cosine.similarity,
    fit.all.v3$`PAT05785`$cosine.similarity,
    fit.all.v3$`PAT06439`$cosine.similarity,
    fit.all.v3$`PAT15211`$cosine.similarity,
    fit.all.v3$`PAT30884`$cosine.similarity,
    fit.all.v3$`PAT33430`$cosine.similarity,
    fit.all.v3$`PAT40364`$cosine.similarity,
    fit.all.v3$`PAT54273`$cosine.similarity,
    fit.all.v3$`PAT59460`$cosine.similarity,
    fit.all.v3$`PAT59600`$cosine.similarity,
    fit.all.v3$`PAT59667`$cosine.similarity,
    fit.all.v3$`PAT71397`$cosine.similarity,
    fit.all.v3$`PAT71662`$cosine.similarity,
    fit.all.v3$`PAT73458`$cosine.similarity,
    fit.all.v3$`PAT73899`$cosine.similarity,
    fit.all.v3$`PAT74143`$cosine.similarity
  ),
  file=paste0("WES_out/LCM_COSMICsignaturesv3_depth", depth, ".csv")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(
  x=data@variants.per.sample,
  file="WES_out/LCM_variants_per_sample.csv"
)

muts <- data.frame(
  matrix(
    ncol = 4, nrow = length(
      levels(genesToBarcodes(data, genes="BRAF")$BRAF$Tumor_Sample_Barcode)
    )
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
write.csv(x=muts, file="WES_out/CRC_mutations_LCM.csv")
