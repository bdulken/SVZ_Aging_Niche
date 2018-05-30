# Fig2C Violins

#==================================================================================================
rm(list = ls())
library(tidyverse)
library(reshape2)
library(viridis)
library(Seurat)

setwd("~/Desktop/Dropbox/Aging\ Single\ Cell\ Paper/Final_Code/")
load("./SmartSeqInputFiles/plot.df.rda")

#==================================================================================================

genes <-  c("Pdcd1", "Ifng")

p <- ggplot(data = filter(plot.df, gene %in% genes, Age == "Old"),
			aes(x = gene, y = expression, fill = Location)) +
			geom_violin(position=position_dodge(.8), scale = "width", draw_quantiles = c(.5)) +
			geom_point(alpha = 0.4, size = .4, 
				position = position_jitterdodge(
					jitter.width = NULL,
					jitter.height = 0,
					dodge.width = .8)) +
			theme_classic() +
			labs(y = "Expression") +
  			theme(plot.title = element_text(hjust = 0.5)) +
  			theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
  			theme(axis.title.y = element_text(size = 12)) +
  			theme(axis.title.x = element_blank()) +
  			scale_fill_manual(values = c("#C70039", "#188ad0")) +
  			ggtitle("T cell Phenotype")
p

