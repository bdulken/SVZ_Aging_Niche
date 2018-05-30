# Fig1cd
# Matthew Buckley
# Refactored 5.27.18

rm(list = ls())
library(tidyverse)
library(tibble)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(scales)
library(NCmisc) #v1.1.5
library(Seurat) #v1.3

setwd("~/Desktop/Dropbox/Aging\ Single\ Cell\ Paper/Final_Code/")
load("~/Desktop/Dropbox/BenNSCProject/10X/Seurat/svz_v1.rda") # Input for Fig1C, loads as "svz"
load("./10xInputFiles/df.all.DE2.results.rda") # Input for Fig 1D, loads as "df"

#==================================================================================================
# Fig 1C

new.cluster.ids <- c("Endothelial","Microglia","Oligodendrocytes","Astrocytes_qNSCs",
	"Neuroblasts","aNSCs_NPCs","Pericytes","T_Cells","OPC") 

svzident2<-svz@ident
young1_svzident<-svzident2[grepl("Young1",names(svzident2))]
old1_svzident<-svzident2[grepl("Old1",names(svzident2))]

young2_svzident<-svzident2[grepl("Young2",names(svzident2))]
old2_svzident<-svzident2[grepl("Old2",names(svzident2))]

age_enrich<-c()
age_fac<-c()

for(i in 1:length(new.cluster.ids)){
  
  percent_y1<-length(young1_svzident[grepl(new.cluster.ids[i],young1_svzident)])/length(young1_svzident)
  percent_o1<-length(old1_svzident[grepl(new.cluster.ids[i],old1_svzident)])/length(old1_svzident)
  percent_y2<-length(young2_svzident[grepl(new.cluster.ids[i],young2_svzident)])/length(young2_svzident)
  percent_o2<-length(old2_svzident[grepl(new.cluster.ids[i],old2_svzident)])/length(old2_svzident)
  
  age_enrich<-c(age_enrich,mean((percent_o1/percent_y1),(percent_o2/percent_y2)))
  age_fac<-c(age_fac,new.cluster.ids[i])
}

dataframe<-data.frame(vals = log2(age_enrich), fac = age_fac,col = c(rep(c("#009933"), length(new.cluster.ids))))
dataframe$fac<-as.character(dataframe$fac)
dataframe$fac <- factor(dataframe$fac, levels=c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts","Oligodendrocytes","OPC","Endothelial","Pericytes","Microglia","T_Cells"), ordered = T)

dataframe$celltype2 <- lapply(dataframe$fac, function(x) {gsub("s_", "s/", x)})
dataframe$celltype2 <- lapply(dataframe$celltype2, function(x) {gsub("_", " ", x)})
dataframe$celltype_fac <- factor(dataframe$celltype2, levels=c("Astrocytes/qNSCs","aNSCs/NPCs","Neuroblasts","Oligodendrocytes","OPC","Endothelial","Pericytes","Microglia","T Cells"), ordered = T)

ben_colors <- c("#03c03c","#0054b4","#966fd6","#aec6cf","#ffdf00","#ffb347","#e5aa70","#db7093","#e8000d")
# Plot Fig 1C
p<-ggplot(dataframe)
p<-p+geom_bar(aes(x=dataframe$celltype_fac,y=dataframe$vals,fill=dataframe$celltype_fac),stat="identity",alpha=0.7)
p<-p+labs(x="Cell Types",y="log2 Enrichment")
p<-p+scale_fill_manual(values = ben_colors) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ggtitle("log2 Cell Type Enrichment") +
  theme(axis.title.y = element_text(size = 20, face = "plain", margin = margin(r = 15))) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.ticks.x = element_blank()) +
  theme(plot.title = element_text(size = 20, face = "plain")) +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed") +
  theme(legend.position="none")
p

#==================================================================================================
# Fig 1D

# Adjust colors
ben_colors <- c("#03c03c", "#0054b4", "#966fd6", "#aec6cf",
                "#ffdf00", "#ffb347", "#e5aa70", "#db7093", "#e8000d")
names(ben_colors) <- levels(df$celltype_factor)

# Make custom color column to facilitate grey coloring by threshold.
col <- ben_colors[df$celltype_factor]
col[df$padj > 0.05] <- "#D3D3D3"
df$col <- as.factor(col)

# Plot Fig 1D
q <- ggplot(df, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("Differential Expression Z-score") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
q

#==================================================================================================
# Combine Fig 1CD into one plot

pq <- plot_grid(p, NULL, q, ncol = 1, rel_heights = c(1, 0.15, 1.5))
pq
