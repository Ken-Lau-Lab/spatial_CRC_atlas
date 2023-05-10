rm(list = ls())

# SET CORRECT WORKING DIRECTORY
setwd("/home/cody/git/spatial_CRC_atlas/step7/")

# make directory for outputs
dir.create("./ST_out/", showWarnings = FALSE)
dir.create("./ST_out/tradeseq/", showWarnings = FALSE)
dir.create("./ST_out/tradeseq/plots/", showWarnings = FALSE)

library(RColorBrewer)
library(SingleCellExperiment)
library(tradeSeq)
library(UpSetR)
library(dplyr)
library(reshape2)
library(stringr)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))

## ----------------------------------------------------------------------------
# read in data
## ----------------------------------------------------------------------------
tradeseq.data <- read.csv(
  "../step6/ST_out/tradeseq/tradeseq.csv",
  row.names = 1
)

gene.counts <- read.csv(
  "../step6/ST_out/tradeseq/tradeseq_counts.csv",
  row.names = 1
)
gene.counts <- t(gene.counts)

## ----------------------------------------------------------------------------
cellWeights <- tradeseq.data[,c("HM","CIN")]
pseudotime <- tradeseq.data[,c("TMB","CNV.score")]
colnames(pseudotime) <- c("HM","CIN")

set.seed(7)
gene.sce <- fitGAM(
  counts = gene.counts,
  pseudotime = pseudotime,
  cellWeights = cellWeights,
  nknots = 4,
  verbose = T,
  parallel = T
)

## ----------------------------------------------------------------------------
table(rowData(gene.sce)$tradeSeq$converged)

## ----------------------------------------------------------------------------
gene.assocRes <- associationTest(gene.sce, lineages=T, l2fc=log2(2))
HM_genes <- rownames(gene.assocRes)[
  which(gene.assocRes$pvalue_1 <= 0.0001)
]
CIN_genes <- rownames(gene.assocRes)[
  which(gene.assocRes$pvalue_2 <= 0.0001)
]
all_genes <- rownames(gene.assocRes)[
  which(gene.assocRes$pvalue <= 0.0001)
]

UpSetR::upset(
  fromList(
    list(
      CIN = CIN_genes,
      HM = HM_genes,
      all = all_genes
    )
  )
)

# HM genes
o <- order(gene.assocRes$waldStat_1, decreasing = TRUE)
write.csv(
  x=gene.assocRes[o,],
  file=paste0("ST_out/tradeseq/associationTest_HM.csv"),
  row.names=T
)
for (i in seq_len(3)) {
  sigGene <- names(gene.sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/associationTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/associationTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(gene.sce, gene.counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# CIN genes
o <- order(gene.assocRes$waldStat_2, decreasing = TRUE)
write.csv(
  x=gene.assocRes[o,],
  file=paste0("ST_out/tradeseq/associationTest_CIN.csv"),
  row.names=T
)
for (i in seq_len(3)) {
  sigGene <- names(gene.sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/associationTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/associationTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(gene.sce, gene.counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# Immune Exclusion genes
for (i in c("DDR1","TGFBI","PAK4","DPEP1")) {
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/associationTest_",
      i,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/associationTest_",
      i,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(gene.sce, gene.counts, gene = i, alpha = 1, border = TRUE)
  ); dev.off()
}

View(gene.assocRes[c("DDR1","TGFBI","PAK4","DPEP1"),])

## -----------------------------------------------------------------------------
gene.startRes <- startVsEndTest(gene.sce, lineages=TRUE)

# HM genes
o <- order(gene.startRes$logFClineage1, decreasing = TRUE)
write.csv(
  x=gene.startRes[o,],
  file=paste0("ST_out/tradeseq/startVsEndTest_HM.csv"),
  row.names=T
)
for (i in seq_len(3)) {
  sigGene <- names(gene.sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/startVsEndTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/startVsEndTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(gene.sce, gene.counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# CIN genes
o <- order(gene.startRes$logFClineage2, decreasing = TRUE)
write.csv(
  x=gene.startRes[o,],
  file=paste0("ST_out/tradeseq/startVsEndTest_CIN.csv"),
  row.names=T
)
for (i in seq_len(3)) {
  sigGene <- names(gene.sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/startVsEndTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/startVsEndTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(gene.sce, gene.counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

View(gene.startRes[c("DDR1","TGFBI","PAK4","DPEP1"),])

HM_genes_start <- rownames(gene.startRes)[
  which(gene.startRes$logFClineage1 > 2)
]
CIN_genes_start <- rownames(gene.startRes)[
  which(gene.startRes$logFClineage2 > 2)
]

## -----------------------------------------------------------------------------
gene.endRes <- diffEndTest(gene.sce)

# HM genes
o <- order(gene.endRes$logFC1_2, decreasing = TRUE)
write.csv(
  x=gene.endRes[o,],
  file=paste0("ST_out/tradeseq/diffEndTest_HM.csv"),
  row.names=T
)
for (i in seq_len(3)) {
  sigGene <- names(gene.sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/diffEndTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/diffEndTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(gene.sce, gene.counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# CIN genes
o <- order(gene.endRes$logFC1_2, decreasing = FALSE)
write.csv(
  x=gene.endRes[o,],
  file=paste0("ST_out/tradeseq/diffEndTest_CIN.csv"),
  row.names=T
)
for (i in seq_len(3)) {
  sigGene <- names(gene.sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/diffEndTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/diffEndTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(gene.sce, gene.counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

View(gene.endRes[c("DDR1","TGFBI","PAK4","DPEP1"),])

HM_genes_end <- rownames(gene.endRes)[
  which(gene.endRes$logFC1_2 > 1.25)
]
CIN_genes_end <- rownames(gene.endRes)[
  which(gene.endRes$logFC1_2 < -3.35)
]

heat.genes <- union(
  CIN_genes_start,
  union(
    CIN_genes_end,
    union(
      HM_genes_start,
      union(
        HM_genes_end,
        c("DDR1","TGFBI","PAK4","DPEP1")
      )
    )
  )
)
heat.genes <- heat.genes[! heat.genes %in% c("DEPRECATED_ENSG00000211890")]

##------------------------------------------------------------------------------
# define custom gene lists from heat.genes
immune.genes <- c(
  "IGHG1","IGKC","IGHJ6","IGLV3.1","IGKV1D.12","IGKV3OR2.268","IGKV4.1",
  "LY6G6D","CXCL14","HLA.A"
)
ecm.genes <- c(
  "DDR1","PAK4","TGFBI","LUM","FLNA"
)
mucosa.genes <- c(
  "REG1A","SPINK4","REG3A","CRACR2A","FCGBP","MUC4","SERPINA1",
  "FABP1","SLC26A3","PRAP1","ADIRF","KRT7","DPEP1","NOX1","CDX2"
)
nuc.genes <- c(
  "HIST1H2BH","HIST1H3F","HIST1H3G","H2AFJ","SAP18","PUF60"
)
cancer.genes <- c(
  "CTSH","FXYD5","TP53","HMGCS2","MT1G","TGM2","KLK10"
)
gene_df <- rbind(
  data.frame(Function=rep("Immune",length(immune.genes)), row.names = immune.genes),
  data.frame(Function=rep("ECM",length(ecm.genes)), row.names = ecm.genes),
  data.frame(Function=rep("Mucosa",length(mucosa.genes)), row.names = mucosa.genes),
  data.frame(Function=rep("Nuclear/TF",length(nuc.genes)), row.names = nuc.genes),
  data.frame(Function=rep("Cancer",length(cancer.genes)), row.names = cancer.genes)
)

yhatSmooth <- predictSmooth(
  gene.sce,
  gene=rownames(gene_df),
  nPoints=50,
  tidy=FALSE
)
yhatsmooth_immexcl <- yhatSmooth[c("DDR1","TGFBI","PAK4","DPEP1"),]
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
colnames(yhatSmoothScaled) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled), "lineage1", "HM"
)
colnames(yhatSmoothScaled) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled), "lineage2", "CIN"
)

cat_df = data.frame("Tumor class" = c(rep("HM",50),rep("CIN",50)), check.names = F)
rownames(cat_df) = colnames(yhatSmoothScaled)

png(
  filename="ST_out/tradeseq/plots/topgenes_grouped_heatmap.png",
  width=4000,
  height=3000,
  res=500,
  type="cairo"
)
print(pheatmap::pheatmap(
  yhatSmoothScaled,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = cat_df,
  annotation_row = gene_df,
  annotation_colors = list(
    "Tumor class" = c(HM = "#7a4fa3", CIN = "#ffc101")
  ),
  gaps_col = c(50),
  legend = TRUE,
  angle_col=0
))
dev.off()


# define custom gene lists from heat.genes
immune.genes <- c(
  "IGHG1","HLA.A","IGHJ6","IGLV3.1","IGKV1D.12","IGKC","IGKV3OR2.268",
  "LY6G6D","IGKV4.1","CXCL14"
)
ecm.genes <- c(
  "DDR1","TGFBI","PAK4","FLNA","LUM"
)
mucosa.genes <- c(
  "KRT7","REG1A","SPINK4","REG3A","CRACR2A","FCGBP","MUC4","SERPINA1",
  "FABP1","SLC26A3","PRAP1","ADIRF","CDX2","NOX1","DPEP1"
)
nuc.genes <- c(
  "HIST1H2BH","HIST1H3F","SAP18","PUF60","HIST1H3G","H2AFJ"
)
cancer.genes <- c(
  "TGM2","KLK10","CTSH","FXYD5","HMGCS2","MT1G","TP53"
)
gene_df <- rbind(
  data.frame(Function=rep("Mucosa",length(mucosa.genes)), row.names = mucosa.genes),
  data.frame(Function=rep("ECM",length(ecm.genes)), row.names = ecm.genes),
  data.frame(Function=rep("Immune",length(immune.genes)), row.names = immune.genes),
  data.frame(Function=rep("Nuclear/TF",length(nuc.genes)), row.names = nuc.genes),
  data.frame(Function=rep("Cancer",length(cancer.genes)), row.names = cancer.genes)
)
yhatSmooth <- predictSmooth(
  gene.sce,
  gene=row.names(gene_df),
  nPoints=50,
  tidy=FALSE
)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
colnames(yhatSmoothScaled) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled), "lineage1", "HM"
)
colnames(yhatSmoothScaled) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled), "lineage2", "CIN"
)

png(
  filename="ST_out/tradeseq/plots/topgenes_grouped_ordered_heatmap.png",
  width=3500,
  height=3000,
  res=500,
  type="cairo"
)
print(pheatmap::pheatmap(
  yhatSmoothScaled,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  breaks = seq(min(yhatSmoothScaled)-0.8, max(yhatSmoothScaled)-0.7, length.out = 100),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = cat_df,
  annotation_row = gene_df,
  annotation_colors = list(
    "Tumor class" = c(HM = "#7a4fa3", CIN = "#ffc101")
  ),
  labels_row = as.expression(lapply(rownames(yhatSmoothScaled), function(a) bquote(italic(.(a))))),
  gaps_col = c(50),
  legend = TRUE,
  angle_col=0,
  #scale="row",
  gaps_row = c(
    length(mucosa.genes),
    length(mucosa.genes) + length(ecm.genes),
    length(mucosa.genes) + length(ecm.genes) + length(immune.genes),
    length(mucosa.genes) + length(ecm.genes) + length(immune.genes) + length(nuc.genes)
  ),
))
dev.off()

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
# Cell states and gene signatures
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
# read in data
## ----------------------------------------------------------------------------
counts <- read.csv(
  "../step6/ST_out/tradeseq/tradeseq_signatures_cellstates.csv",
  row.names = 1
)
counts <- t(counts)

## ----------------------------------------------------------------------------
cellWeights <- tradeseq.data[,c("HM","CIN")]
pseudotime <- tradeseq.data[,c("TMB","CNV.score")]
colnames(pseudotime) <- c("HM","CIN")

set.seed(7)
sce <- fitGAM(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cellWeights,
  nknots = 4,
  verbose = T,
  parallel=T,
  family = "gaussian"
)

## ----------------------------------------------------------------------------
table(rowData(sce)$tradeSeq$converged)

## ----------------------------------------------------------------------------
assocRes <- associationTest(sce, lineages=T, l2fc=log2(2))

HM_genes <- rownames(assocRes)[
  which(assocRes$pvalue_1 <= 0.05)
]
CIN_genes <- rownames(assocRes)[
  which(assocRes$pvalue_2 <= 0.05)
]
all_genes <- rownames(assocRes)[
  which(assocRes$pvalue <= 0.0001)
]

UpSetR::upset(
  fromList(
    list(
      CIN = CIN_genes,
      HM = HM_genes,
      all = all_genes
    )
  )
)

# HM genes
o <- order(assocRes$waldStat_1, decreasing = TRUE)
write.csv(
  x=assocRes[o,],
  file=paste0("ST_out/tradeseq/associationTest_HM_sigs.csv"),
  row.names=T
)
for (i in seq_len(10)) {
  sigGene <- names(sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/associationTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/associationTest_HM_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(sce, counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# CIN genes
o <- order(assocRes$waldStat_2, decreasing = TRUE)
write.csv(
  x=assocRes[o,],
  file=paste0("ST_out/tradeseq/associationTest_CIN_sigs.csv"),
  row.names=T
)
for (i in seq_len(10)) {
  sigGene <- names(sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/associationTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/associationTest_CIN_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(sce, counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# Immune Exclusion genes
for (i in c("DDR1","PAK4","TGFBI","DPEP1")) {
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/associationTest_",
      i,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/associationTest_",
      i,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(sce, counts, gene = i, alpha = 1, border = TRUE)
  ); dev.off()
}

## -----------------------------------------------------------------------------
startRes <- startVsEndTest(sce, lineages=TRUE)

# HM genes
o <- order(startRes$logFClineage1, decreasing = TRUE)
write.csv(
  x=startRes[o,],
  file=paste0("ST_out/tradeseq/startVsEndTest_HM_sigs.csv"),
  row.names=T
)
for (i in seq_len(10)) {
  sigGene <- names(sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/startVsEndTest_HM_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/startVsEndTest_HM_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(sce, counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# CIN genes
o <- order(startRes$logFClineage2, decreasing = TRUE)
write.csv(
  x=startRes[o,],
  file=paste0("ST_out/tradeseq/startVsEndTest_CIN_sigs.csv"),
  row.names=T
)
for (i in seq_len(10)) {
  sigGene <- names(sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/startVsEndTest_CIN_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/startVsEndTest_CIN_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(sce, counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

HM_genes_start <- rownames(startRes)[
  which(startRes$logFClineage1 > 0.01)
]
CIN_genes_start <- rownames(startRes)[
  which(startRes$logFClineage2 > 0.01)
]

## -----------------------------------------------------------------------------
endRes <- diffEndTest(sce)

# HM genes
o <- order(endRes$logFC1_2, decreasing = TRUE)
write.csv(
  x=endRes[o,],
  file=paste0("ST_out/tradeseq/diffEndTest_HM_sigs.csv"),
  row.names=T
)
for (i in seq_len(10)) {
  sigGene <- names(sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/diffEndTest_HM_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/diffEndTest_HM_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(sce, counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

# CIN genes
o <- order(endRes$logFC1_2, decreasing = FALSE)
write.csv(
  x=endRes[o,],
  file=paste0("ST_out/tradeseq/diffEndTest_CIN_sigs.csv"),
  row.names=T
)
for (i in seq_len(10)) {
  sigGene <- names(sce)[o[i]]
  print(
    paste0(
      "Saving plot: ST_out/tradeseq/plots/diffEndTest_CIN_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    )
  )
  png(
    filename=paste0(
      "ST_out/tradeseq/plots/diffEndTest_CIN_sigs_rank",
      i,
      "_",
      sigGene,
      ".png"
    ),
    res=300,
    width=2000,
    height=1500
  )
  print(
    plotSmoothers(sce, counts, gene = sigGene, alpha = 1, border = TRUE)
  ); dev.off()
}

HM_genes_end <- rownames(endRes)[
  which(endRes$logFC1_2 > 0.1)
]
CIN_genes_end <- rownames(startRes)[
  which(endRes$logFC1_2 < -0.1)
]

heat.genes <- union(
  CIN_genes_start,
  union(
    CIN_genes_end,
    union(
      HM_genes_start,
      HM_genes_end
    )
  )
)

# Fig 5-6 plot
yhatSmooth <- predictSmooth(
  sce,
  gene=make.names(c(
    'Fibrosis','Oxphos','Hypoxia','pEMT','CytoTRACE','iCMS2','CRC2','Stem',
    'FIB2','FIB3','MYE2',
    'T reg suppressive','T cell CD4','TL1','MYE5','T cell CD8','TL2','TL3',
    'Fetal','Metaplasia','iCMS3','GOB','SSC','ABS','CT'
    #'MHC','Bacterial Response','Stress Response','Interferon',
    #'STM','END1','BL1','FIB1','CRC1','MYE1','CRC3','EE1','MYE3','PLA',
    #'MAS','CRC4','TUF','FIB4','END2','EE2','BL2'
  )),#heat.genes,
  nPoints=50,
  tidy=FALSE
)
yhatsmooth_immexcl <- predictSmooth(
  gene.sce,
  gene=c("DDR1","TGFBI","PAK4","DPEP1"),
  nPoints=50,
  tidy=FALSE
)
yhatSmooth <- rbind(yhatsmooth_immexcl, yhatSmooth)
yhatSmooth <- rbind(colMeans(yhatSmooth[1:4,]), yhatSmooth)
rownames(yhatSmooth)[1] <- "IES"
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
colnames(yhatSmoothScaled) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled), "lineage1", "HM"
)
colnames(yhatSmoothScaled) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled), "lineage2", "CIN"
)
rownames(yhatSmoothScaled) <- stringr::str_replace_all(
  rownames(yhatSmoothScaled), "\\.", " "
)

cat_df = data.frame("Tumor class" = c(rep("HM",50),rep("CIN",50)), check.names = F)
rownames(cat_df) = colnames(yhatSmoothScaled)

png(
  filename="ST_out/tradeseq/plots/immexcl_heatmap.png",
  width=4000,
  height=2100,
  res=500,
  type="cairo"
)
print(pheatmap::pheatmap(
  yhatSmoothScaled,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_cols = FALSE,
  border_color=NA,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = cat_df,
  annotation_colors = list(
    "Tumor class" = c(HM = "#7a4fa3", CIN = "#ffc101")
  ),
  labels_row = c(
    "IES", 
    expression(italic("DDR1")), expression(italic("TGFBI")),
    expression(italic("PAK4")), expression(italic("DPEP1")),
    'Fibrosis','Oxphos','Hypoxia','pEMT','CytoTRACE','iCMS2','CRC2','Stem',
    'FIB2','FIB3','MYE2',
    'T reg suppressive','T cell CD4','TL1','MYE5','T cell CD8','TL2','TL3',
    'Fetal','Metaplasia','iCMS3','GOB','SSC','ABS','CT'
  ),
  #annotation_row = state_df,
  #cutree_rows = 4,
  gaps_col = c(50),
  gaps_row = c(1,5),
  #main = "HM",
  legend = TRUE
))
dev.off()

##------------------------------------------------------------------------------

state_df = data.frame("Compartment" = c(
  'Epithelial','Stromal','Immune','Stromal','Epithelial','Immune','Immune','Immune','Epithelial','Epithelial','Epithelial',
  'Epithelial','Epithelial','Immune','Immune','Stromal','Immune','Epithelial','Immune','Immune','Epithelial','Epithelial',
  'Epithelial','Stromal','Stromal','Immune','Stromal','Immune','Epithelial','Immune'
))
rownames(state_df) = c(
  'STM','END1','BL1','FIB1','CRC1','MYE1','TL1','MYE2','CRC2','CT','SSC',
  'CRC3','EE1','MYE3','PLA','FIB2','MYE4','GOB','MAS','MYE5','CRC4','ABS',
  'TUF','FIB3','FIB4','TL2','END2','TL3','EE2','BL2'
)

yhatSmooth <- predictSmooth(
  sce,
  gene=c(
    'ABS','CRC1','CRC2','CRC3','CRC4','CT','EE1','EE2',
    'GOB','SSC','STM','TUF'
  ),#heat.genes,
  nPoints=50,
  tidy=FALSE
)
yhatSmoothScaled.epi <- t(scale(t(yhatSmooth)))
colnames(yhatSmoothScaled.epi) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled.epi), "lineage1", "HM"
)
colnames(yhatSmoothScaled.epi) <- stringr::str_replace_all(
  colnames(yhatSmoothScaled.epi), "lineage2", "CIN"
)

cat_df = data.frame("Tumor class" = c(rep("HM",50),rep("CIN",50)), check.names = F)
rownames(cat_df) = colnames(yhatSmoothScaled.epi)

heatSmooth_Epi <- pheatmap::pheatmap(
  yhatSmoothScaled.epi,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  breaks = seq(min(yhatSmoothScaled.epi), max(yhatSmoothScaled.epi), length.out = 100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = cat_df,
  annotation_colors = list(
    "Tumor class" = c(HM = "#7a4fa3", CIN = "#ffc101")
  ),
  #annotation_row = state_df,
  #cutree_rows = 4,
  gaps_col = c(50),
  #main = "HM",
  legend = TRUE,
  silent = TRUE
)

yhatSmooth <- predictSmooth(
  sce,
  gene=c(
    'END1','END2','FIB1','FIB2',
    'FIB3','FIB4'
  ),#heat.genes,
  nPoints=50,
  tidy=FALSE
)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
heatSmooth_Stroma <- pheatmap::pheatmap(
  yhatSmoothScaled,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  breaks = seq(min(yhatSmoothScaled.epi), max(yhatSmoothScaled.epi), length.out = 100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = cat_df,
  annotation_colors = list(
    "Tumor class" = c(HM = "#7a4fa3", CIN = "#ffc101")
  ),
  #annotation_row = state_df,
  #cutree_rows = 4,
  gaps_col = c(50),
  #main = "HM",
  legend = TRUE,
  silent = TRUE
)

yhatSmooth <- predictSmooth(
  sce,
  gene=c(
    'BL1','BL2','MAS','MYE1','MYE2',
    'MYE3','MYE4','MYE5','PLA','TL1','TL2','TL3'
  ),#heat.genes,
  nPoints=50,
  tidy=FALSE
)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
heatSmooth_Immune <- pheatmap::pheatmap(
  yhatSmoothScaled,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  breaks = seq(min(yhatSmoothScaled.epi), max(yhatSmoothScaled.epi), length.out = 100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = cat_df,
  annotation_colors = list(
    "Tumor class" = c(HM = "#7a4fa3", CIN = "#ffc101")
  ),
  #annotation_row = state_df,
  #cutree_rows = 4,
  gaps_col = c(50),
  #main = "HM",
  legend = TRUE,
  silent = TRUE
)

png(
  filename="ST_out/tradeseq/plots/cellstate_heatmap.png",
  width=5000,
  height=2800,
  res=500,
  type="cairo"
)
print(gridExtra::grid.arrange(
  heatSmooth_Epi[[4]],
  heatSmooth_Stroma[[4]],
  heatSmooth_Immune[[4]],
  heights = c(10,6,10),
  ncol=1
))
dev.off()


all.genes <- c(mucosa.genes, immune.genes, ecm.genes, nuc.genes, cancer.genes)
all.sigs <- c(
  'STM','END1','BL1','FIB1','CRC1','MYE1','TL1','MYE2','CRC2','CT','SSC',
  'CRC3','EE1','MYE3','PLA','FIB2','MYE4','GOB','MAS','MYE5','CRC4','ABS',
  'TUF','FIB3','FIB4','TL2','END2','TL3','EE2','BL2','Fibrosis','Fetal',
  'Oxphos','Hypoxia','pEMT','CytoTRACE','iCMS2','Stem','IES',
  'T.reg.suppressive','T.cell.CD4','T.cell.CD8','Metaplasia','iCMS3'
)

write.csv(
  x=rbind(endRes[all.sigs,], gene.endRes[all.genes,]) %>% tidyr::drop_na(),
  file=paste0("ST_out/tradeseq/diffEndTest_comb.csv"),
  row.names=T
)

write.csv(
  x=rbind(startRes[all.sigs,], gene.startRes[all.genes,]) %>% tidyr::drop_na(),
  file=paste0("ST_out/tradeseq/startVsEndTest_comb.csv"),
  row.names=T
)

write.csv(
  x=rbind(assocRes[all.sigs,], gene.assocRes[all.genes,]) %>% tidyr::drop_na(),
  file=paste0("ST_out/tradeseq/associationTest_comb.csv"),
  row.names=T
)
