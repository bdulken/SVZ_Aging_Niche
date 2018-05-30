# Figure 2A Receptor tSNE

rm(list=ls())
library(Seurat)
library(ggplot2)

setwd("~/Desktop/Dropbox/Aging\ Single\ Cell\ Paper/Final_Code/")
load("./10xInputFiles/svz_tsne_young_old_tworeps_celllabeled_clusters_markers.Robj")

allGenes <- rownames(svz@data)
data_mat <- as.matrix(svz@data)

df <- data.frame(tsne_x = svz@tsne.rot[,1], tsne_y = svz@tsne.rot[,2])
rownames(df) <- rownames(svz@tsne.rot)

# ifn g receptors
ifngr_genes <- sort(allGenes[grep("Ifngr.", allGenes)])
ifngr_sig <- colSums(data_mat[ifngr_genes, ])
df$ifngr_sig <- ifngr_sig

low_col <- "grey"
high_col <- "darkred"
size_range <- c(1,2)

# Saturate color at 4.5 log-normed expression.
ifngr_sig_cap <- df$ifngr_sig
ifngr_sig_cap[ifngr_sig_cap > 4.5] <- 4.5
df$ifngr_sig_cap <- ifngr_sig_cap

# Plot IFN gamma receptor
q <- ggplot(data = df, aes(x = tsne_x, y = tsne_y,
         alpha = ifngr_sig_cap, color = ifngr_sig_cap))
q <- q + geom_point() + scale_size(range = size_range)
q <- q + scale_colour_gradient(low = low_col, high = high_col)
q <- q + ggtitle("IFN Gamma Receptor") + 
         theme(axis.line = element_blank(), axis.text = element_blank(),
         axis.ticks = element_blank(), axis.title = element_blank())
q
