# Figure 3B Paired Violins

rm(list=ls())
library(Seurat) # v1.3
library(magrittr)
library(dplyr)
library(devtools)
library(reshape2)
# Mouse gene information. # source("https://bioconductor.org/biocLite.R"); biocLite("org.Mm.eg.db"):
library(org.Mm.eg.db)
# devtools::install_github("stephenturner/msigdf"):
library(msigdf) # Loads useful msigdf.mouse dataframe based on MsigDB
library(annotate) 
library(ggplot2)

setwd("~/Desktop/Dropbox/Aging\ Single\ Cell\ Paper/Final_Code/")
load("./10xInputFiles/svz_tsne_young_old_tworeps_celllabeled_clusters_markers.Robj")

#==================================================================================================
# Refine seurat object and extract data
age_df = data.frame("Age" = gsub("(*.)([12])", "\\1", FetchData(svz, "orig.ident")[,"orig.ident"]))
row.names(age_df) = svz@cell.names
svz <- AddMetaData(object = svz, metadata = age_df, col.name = "Age")
allGenes <- rownames(svz@data)
data_mat <- as.matrix(svz@data) 

# Make dataframe with tSNE coorderinates
df_age <- FetchData(svz, c("Age"))
age_vector <- df_age$Age
age_relevel <- factor(age_vector, levels = c("Young", "Old"))
df <- data.frame(id = colnames(svz@data), celltype = svz@ident, svz@tsne.rot, age = age_relevel)
df$celltype_factor <- factor(df$celltype,  levels=c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts",
                             "Oligodendrocytes","OPC","Endothelial","Pericytes","Microglia",
                             "T_Cells"), ordered = T)

#==================================================================================================
# Setup genesets.

sym <- getSYMBOL(x = as.character(msigdf.mouse$entrez), data="org.Mm.eg.db")
msigdf.mouse$sym <- sym

# Filter for Hallmarks, Kegg, and reactome collections
hall.df <- filter(msigdf.mouse, collection == "hallmark")
hall.kegg.df <- rbind(hall.df, filter(msigdf.mouse, collection == "c2" & grepl("^KEGG", geneset)))
hall.kegg.react.df <- rbind(hall.kegg.df, filter(msigdf.mouse, collection == "c2" & grepl("^REACT", geneset)))

# Create a named list where each name is a geneset and each list contains the corresponding genes.
addGeneSet <- function(gSetName, data) {
	gSet.df <- filter(data, geneset == gSetName)
	hall.kegg.react.list[[gSetName]] = gSet.df$sym
}

hall.kegg.react.list <- list() # Initialize empty list that addGeneSet() will modify.
hall.geneSet.names <- unique(hall.kegg.react.df$geneset) # All gene lists (910)
hall.kegg.react.list <- lapply(hall.geneSet.names, addGeneSet, msigdf.mouse) # Add gene lists
names(hall.kegg.react.list) <- hall.geneSet.names # Add gene list names

# Add signatures to df
for (n in names(hall.kegg.react.list)) {
	subsetMatrix <- data_mat[rownames(data_mat) %in% hall.kegg.react.list[[n]], ] # Only relevant genes
	df[n] <- colSums(subsetMatrix) } # Sum over all relevant genes to get geneset expression value

#==================================================================================================
# Plot paired young/old violins for 8 celltypes depicting IFNg expression.

df8 <- filter(df, celltype %in% c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts",
                             "Oligodendrocytes","Endothelial","Pericytes","Microglia",
                             "T_Cells"))

Plot_Violins <- function(signature, data) {
	p <- ggplot(data = df8, aes_string(x = "age", y = as.name(signature), fill = "age"))
	p <- p + geom_point(color = "black", size = 1, alpha = .3, position = position_jitter(w = 0.35, h = 0))
	p <- p + geom_violin(alpha = .7, draw_quantiles = c(.5))
	p <- p + ggtitle(signature) + scale_fill_manual(values = c("skyblue", "firebrick"))
	p <- p + facet_wrap(~celltype_factor, scales = "free", nrow =2)
	p <- p + theme(axis.text.x = element_blank(),
				   axis.ticks.x = element_blank(),
				   axis.title.x = element_blank(),
				   axis.title.y = element_blank(),
				   axis.text.y = element_text(size = 8))
	p <- p + theme(strip.text.x = element_text(size = 8))
	p
	#ggsave(plot=p, filename = paste0("plots/custom/violin_8L", signature, ".pdf"),
	#width = 6, height = 3, useDingbats=F)
}

interest <- grep("INTERFERON_GAMMA_RESPONSE", colnames(df))
lapply(colnames(df)[interest], Plot_Violins, df)


