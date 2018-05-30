#This script generates the cell type enrichment bar plot shown in figure 1C
#Ben Dulken

rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles")
load("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")

new.cluster.ids <- c("Endothelial","Microglia","Oligodendrocytes","Astrocytes_qNSCs","Neuroblasts","aNSCs_NPCs","Pericytes","T_Cells","OPC") 

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

dataframe<-data.frame(vals=log2(age_enrich),fac=age_fac,col=c(rep(c("#009933"),length(new.cluster.ids))))
dataframe$fac<-as.character(dataframe$fac)
dataframe$fac <- factor(dataframe$fac, levels=c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts","Oligodendrocytes","OPC","Endothelial","Pericytes","Microglia","T_Cells"), ordered = T)
p<-ggplot(dataframe)
p<-p+geom_bar(aes(x=dataframe$fac,y=dataframe$vals,fill=dataframe$fac),stat="identity",alpha=0.7)
p<-p+labs(x="Cell Types",y="log2 enrichment in Old")
p<-p+scale_fill_manual(values=c("#03c03c","#0054b4","#966fd6","#aec6cf","#ffdf00","#ffb347","#e5aa70","#db7093","#e8000d"))
p<-p+theme(legend.position="none")
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
pdf("Perfusion2_celltypeenrich_meanbar_celltypecolored.pdf",height=4,width=6)
print(p)
dev.off()