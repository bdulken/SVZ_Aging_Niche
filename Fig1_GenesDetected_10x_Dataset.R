#This script generates a violin plot showing total unique gene counts for each cell type and each replicate.
#Ben Dulken

rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles")
load("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")

new.cluster.ids <- c("Endothelial","Microglia","Oligodendrocytes","Astrocytes_qNSCs","Neuroblasts","aNSCs_NPCs","Pericytes","T_Cells","OPC")

raw.data.filt<-svz@raw.data[,match(colnames(svz@data),colnames(svz@raw.data))]
genes_detected<-colSums(raw.data.filt!=0)

all_exprs_curr<-svz@data
celltype_fac<-svz@ident
comb_fac_col<-c(rep("Young1",length(all_exprs_curr[1,])))
comb_fac_col[grepl("Old1",colnames(all_exprs_curr))]<-"Old1"
comb_fac_col[grepl("Young2",colnames(all_exprs_curr))]<-"Young2"             
comb_fac_col[grepl("Old2",colnames(all_exprs_curr))]<-"Old2"

dataframe<-data.frame(data=genes_detected,col=comb_fac_col,fac=celltype_fac)
dataframe$col<-as.character(dataframe$col)
dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
dataframe$fac<-as.character(dataframe$fac)
dataframe$fac <- factor(dataframe$fac, levels=c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts","Oligodendrocytes","OPC","Endothelial","Pericytes","Microglia","T_Cells"), ordered = T)
p<-ggplot(dataframe)
#p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
#p<-p+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1],fill=dataframe$col),position = position_jitter(width = .3,height=0),alpha=0.5,color=dataframe$col,size=3)+geom_violin(aes(x=dataframe$fac,y=dataframe[,1],fill=dataframe$col),trim=T,scale="width",alpha=0.5,width=0.5)
p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe[,1],fill=dataframe$col),trim=T,scale="width",alpha=0.5,width=0.5)


p<-p+ theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+theme(axis.text.y=element_text(size=15))
p<-p+theme(axis.text.x=element_text(size=15))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.title=element_text(size=20))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+labs(y="Genes Detected")
p<-p+theme(plot.title = element_text(hjust = 0.5,size=24))
p<-p+scale_fill_manual(values=c("deepskyblue","slateblue","darkorange","firebrick"))
p<-p+theme(axis.title.x=element_blank())
p<-p+labs(title="Genes Detected")

pdf("GenesDetected_allcelltypes_tworeps.pdf",width=18,height=4)
print(p)
dev.off()