#This script generates violin plots in which all cell types are plotted next to each other.
#Ben Dulken
rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles")
load("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")

new.cluster.ids <- c("Endothelial","Microglia","Oligodendrocytes","Astrocytes_qNSCs","Neuroblasts","aNSCs_NPCs","Pericytes","T_Cells","OPC")

ident_curr<-svz@ident[grepl("Microglia",svz@ident)]
all_exprs_curr<-svz@raw.data
all_exprs_curr<-all_exprs_curr[,match(names(ident_curr),colnames(all_exprs_curr))]

all_exprs_curr<-svz@data
celltype_fac<-svz@ident
comb_fac_col<-c(rep("deepskyblue",length(all_exprs_curr[1,])))
comb_fac_col[grepl("Old1",colnames(all_exprs_curr))]<-"darkorange"
comb_fac_col[grepl("Young2",colnames(all_exprs_curr))]<-"slateblue"             
comb_fac_col[grepl("Old2",colnames(all_exprs_curr))]<-"firebrick"


int_genes<-c("Cd8a","Ifng","Cd3e","Cd4")
int_fpkm_glm<-all_exprs_curr[match(int_genes,rownames(all_exprs_curr)),]

for(j in 1:length(int_genes)){
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=t(int_fpkm_glm[j,]),col=comb_fac_col,fac=celltype_fac)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts","Oligodendrocytes","OPC","Endothelial","Pericytes","Microglia","T_Cells"), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),position = position_jitter(width = .3,height=0),alpha=0.5,color=dataframe$col,size=3)+geom_violin(aes(x=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.5,width=0.5)
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
  p<-p+labs(y="Log-normalized counts")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=24))
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=int_genes[j])
  
  pdf(paste(rownames(int_fpkm_glm)[j],"_allcelltypes_tworeps.pdf",sep=""),width=12,height=4)
  print(p)
  dev.off()
}






#Plotting with raw counts for lowly expressd genes

rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles")
load("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")

new.cluster.ids <- c("Endothelial","Microglia","Oligodendrocytes","Astrocytes_qNSCs","Neuroblasts","aNSCs_NPCs","Pericytes","T_Cells","OPC")

all_exprs_curr<-svz@raw.data
all_exprs_curr<-all_exprs_curr[,match(names(svz@ident),colnames(all_exprs_curr))]
celltype_fac<-svz@ident
comb_fac_col<-c(rep("deepskyblue",length(all_exprs_curr[1,])))
comb_fac_col[grepl("Old1",colnames(all_exprs_curr))]<-"darkorange"
comb_fac_col[grepl("Young2",colnames(all_exprs_curr))]<-"slateblue"             
comb_fac_col[grepl("Old2",colnames(all_exprs_curr))]<-"firebrick"


int_genes<-c("Ifnb1","Ifna1")
#int_genes<-rownames(all_exprs_curr)[grep("Ifna",rownames(all_exprs_curr))]
int_fpkm_glm<-all_exprs_curr[match(int_genes,rownames(all_exprs_curr)),]

for(j in 1:length(int_genes)){
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=t(int_fpkm_glm[j,]),col=comb_fac_col,fac=celltype_fac)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts","Oligodendrocytes","OPC","Endothelial","Pericytes","Microglia","T_Cells"), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),position = position_jitter(width = .3,height=0),alpha=0.5,color=dataframe$col,size=3)+geom_violin(aes(x=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.5,width=0.5)
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
  p<-p+labs(y="Log-normalized counts")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=24))
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=int_genes[j])
  
  pdf(paste(rownames(int_fpkm_glm)[j],"_allcelltypes_tworeps_rawcounts.pdf",sep=""),width=12,height=4)
  print(p)
  dev.off()
}

