# SMART-Seq Single-Cell T-cell Heatmaps - old
# Matthew Buckley

#==================================================================================================
rm(list = ls())
library(tidyverse)
library(reshape2)
library(viridis)
library(Seurat)
library(pvclust)
library(gplots)
library(corrplot)
#install.packages("devtools");devtools::install_github("TomKellyGenetics/heatmap.2x", ref="test")
library(heatmap.2x) 
library(Seurat)

setwd("~/Desktop/Dropbox/Aging\ Single\ Cell\ Paper/Final_Code/")
load("./SmartSeqInputFiles/plot.df.rda")
load("./SmartSeqInputFiles/tcl_v1_tcell.rda")
load("./SmartSeqInputFiles/tcl.markers.rda")

#==================================================================================================
data <- as.matrix(tcl@data)
meta <- tcl@data.info
color <- magma(n = 20, direction = -1)
loc_colors <- c("#C70039", "#188ad0")

# Get genes from seurat markers and use order from seurat's heatmap function.
top20 <- tcl.markers %>% group_by(cluster) %>% top_n(20, avg_diff)
h <- DoHeatmap(tcl, genes.use = top20$gene, do.return = T)
seuratGenes <- rownames(h)

# Custom heatmap
plot_dat <- data[rownames(data) %in% seuratGenes, ]
plot_dat <- plot_dat[, grepl("Old", colnames(plot_dat))]
meta_old <- meta[grepl("Old", rownames(meta)), ]
h2 <- heatmap.2(x = plot_dat,
		  Rowv = T,
		  Colv = T,
		  cexRow = 1.5,
		  scale="none",
		  key = FALSE,
		  hclustfun = hclust,
		  dendrogram = "both",
		  labCol = FALSE,
		  ColSideColors = loc_colors[factor(meta_old$Location)],
		  trace = "none",
		  col=color,
		  main = "SVZ/Blood T-Cell Top Cluster Markers")
h2
