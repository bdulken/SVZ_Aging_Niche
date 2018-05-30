# SMART-Seq V5 Single-Cell T-cell Naive/Memory/Effect Characterization Analysis - tSNE
# 4/5/2018
# Matthew Buckley
remove(list=ls())
library(ggplot2)
# old_1.3 <- "https://github.com/satijalab/seurat/archive/v1.3.0.tar.gz"
# install.packages(old_1.3, repos = NULL, type = "source")
# packageVersion("Seurat") # should be 1.3
# x <- ls(getNamespace("Seurat"), all.names=TRUE)
# x[grep("Set", x)] # Check that Setup function is available to confirm v1.3. 
library(Seurat)

setwd("~/Desktop/Dropbox/Aging\ Single\ Cell\ Paper/Final_Code/")
load("./SmartSeqInputFiles/tcl_v1_tcell.rda")


#===================================================================================================
# Custom tSNE
meta <- tcl@data.info

df <- data.frame(id = colnames(tcl@data),
                 Location = meta$Location,
				 Age = meta$Age,
                 Clonal = meta$Clonal,
                 Paired.Clonal = meta$Paired.Clonal,
				 tcl@tsne.rot)

p <- ggplot(data = filter(df, Age == "Old"),
            aes(x = tSNE_1, y = tSNE_2, color = Location)) +
			geom_point(size = 5, alpha = .6) +
			ggtitle("tSNE Projection of T-Cell Transcriptomes") +
			scale_color_manual(values = c("#C70039", "#188ad0"))
p





