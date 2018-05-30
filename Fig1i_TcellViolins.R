# 10X Single-Cell T-cell Naive/Memory/Effect Characterization Analysis - Violin Plots
# Created: Jan 30, 2018, Updated: May 27, 2018
# Matthew Buckley

#==================================================================================================

setwd("~/Desktop/Dropbox/Aging\ Single\ Cell\ Paper/Final_Code/")
library(tidyverse)
library(reshape2)
library(viridis)
load(file = "./10xInputFiles/comphensive_df.Robj") # plotting dataframe

#==================================================================================================
# Filter out non-t-cells
tcells <- df[df$celltype == "T_Cells",] # 203 t-cells x 19340 (labels + genes).

# Specify Genes of Interest
total_genes <- colnames(tcells)
genes <- c("Cd3e", "Cd8a", "Cd4",  "Cd62L", "Cd44", "Cd69", "Xcl1","Itgal", "Itga4", "Ifng")

# Rename Sell as Cd62L
colnames(tcells)[colnames(tcells)=="Sell"] <- "Cd62L"

# Subset to genes of interest
tcells_sub <- tcells[, colnames(tcells) %in% genes]

# Create data frame
tcells_col <- tcells_sub
tcells_col$cell <- tcells$id
rownames(tcells_col) <- NULL
tcells_melt <- melt(tcells_col, id.vars = c("cell"), value.name = "expression")
tcells_melt$gene <- tcells_melt$variable; tcells_melt$variable <- NULL
# Split cell column up to get separate age column
tcells_split <- separate(tcells_melt, col = cell, into = c("age", "id"), sep = "_")
# Order genes for plotting
tcells_split$gene <- factor(tcells_split$gene, levels = rev(genes))


#==================================================================================================
# Plot vertical violins.

p <- ggplot(data = tcells_split, aes(x = gene, y = expression)) +
			geom_violin(fill = "#e8000d", alpha = .75, scale = "width") +
			geom_point(alpha = 0.2, position = position_jitter(width=.3, height=0)) +
			theme_classic() +
			labs(y = "Expression") +
			ggtitle("T-cell Phenotype") +
  			theme(plot.title = element_text(hjust = 0.5)) +
  			theme(axis.text.x=element_text(size=12)) +
  			theme(axis.text.y=element_text(size=12, margin = margin(r = 22))) +
  			theme(axis.ticks.y = element_blank()) +
  			theme(axis.title.y = element_text(margin = margin(r = 60))) +
  			theme(axis.title.x = element_text(size = 12)) +
  			coord_flip()
p

#ggsave("plots3_violins/Vert_tCellViolin2.pdf", p, height = 5, width = 4, bg = "transparent")
