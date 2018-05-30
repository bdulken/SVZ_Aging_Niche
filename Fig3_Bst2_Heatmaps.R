rm(list=ls())
library(Seurat)
library(ggplot2)
library(viridis)
library(pheatmap)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles/")
load("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")

ident_curr<-c(svz@ident[grepl("Astrocytes_qNSCs",svz@ident)],svz@ident[grepl("aNSCs_NPCs",svz@ident)])
all_exprs_curr<-svz@data[,match(names(ident_curr),colnames(svz@data))]

greaterthan0<-all_exprs_curr>0
greaterthan0sum<-rowSums(greaterthan0)
all_exprs_curr_genefilt<-all_exprs_curr[greaterthan0sum>=10,]

young1_counts<-all_exprs_curr_genefilt[,grepl("Young1",colnames(all_exprs_curr_genefilt))]
old1_counts<-all_exprs_curr_genefilt[,grepl("Old1",colnames(all_exprs_curr_genefilt))]

young2_counts<-all_exprs_curr_genefilt[,grepl("Young2",colnames(all_exprs_curr_genefilt))]
old2_counts<-all_exprs_curr_genefilt[,grepl("Old2",colnames(all_exprs_curr_genefilt))]

comb_fpkm<-cbind(young1_counts,old1_counts,young2_counts,old2_counts)
#Custom genelist

genelist<-c("B2m","Ifi27","Ifitm3","Ifit3","Ifit1","Xaf1","Bst2")
#genelist<-as.vector(read.table("/Volumes/LaCie/09022016_Dhananjay_10X/Analysis/geneCluster.148.txt"))[,1]
genelist_name<-"custom"
fpkm_glm_genes<-comb_fpkm[na.omit(match(tolower(genelist),tolower(rownames(comb_fpkm)))),]

annotation<-data.frame(age=factor(c(rep("Young",length(young1_counts)),
                                    rep("Old",length(old1_counts)),
                                    rep("Young",length(young2_counts)),
                                    rep("Old",length(old2_counts))),levels=c("Young","Old"),labels=c("Young","Old")))
rownames(annotation)<-colnames(fpkm_glm_genes)

Var1 <- c("deepskyblue","#EF3B39")
names(Var1) <- c("Young", "Old")
anno_colors <- list(age = Var1)

pdf(paste(genelist_name,"withastrocytes_global_clustering.pdf",sep=""),width=5,height=3,onefile=F)
pheatmap(fpkm_glm_genes, col = magma(n = 30, direction = -1), cluster_rows=F,cluster_cols=T, annotation=annotation, annotation_colors = anno_colors,show_colnames=F)
dev.off()


#================================================================================================================================================================================

#Microglia Heatmap
ident_curr<-c(svz@ident[grepl("Microglia",svz@ident)])
all_exprs_curr<-svz@data[,match(names(ident_curr),colnames(svz@data))]

greaterthan0<-all_exprs_curr>0
greaterthan0sum<-rowSums(greaterthan0)
all_exprs_curr_genefilt<-all_exprs_curr[greaterthan0sum>=10,]

young1_counts<-all_exprs_curr_genefilt[,grepl("Young1",colnames(all_exprs_curr_genefilt))]
old1_counts<-all_exprs_curr_genefilt[,grepl("Old1",colnames(all_exprs_curr_genefilt))]

young2_counts<-all_exprs_curr_genefilt[,grepl("Young2",colnames(all_exprs_curr_genefilt))]
old2_counts<-all_exprs_curr_genefilt[,grepl("Old2",colnames(all_exprs_curr_genefilt))]

comb_fpkm<-cbind(young1_counts,old1_counts,young2_counts,old2_counts)
#Custom genelist

genelist<-c("B2m","Ifi27","Ifitm3","Ifit3","Ifit1","Xaf1","Bst2")
#genelist<-as.vector(read.table("/Volumes/LaCie/09022016_Dhananjay_10X/Analysis/geneCluster.148.txt"))[,1]
genelist_name<-"custom"
fpkm_glm_genes<-comb_fpkm[na.omit(match(tolower(genelist),tolower(rownames(comb_fpkm)))),]

annotation<-data.frame(age=factor(c(rep("Young",length(young1_counts)),
                                    rep("Old",length(old1_counts)),
                                    rep("Young",length(young2_counts)),
                                    rep("Old",length(old2_counts))),levels=c("Young","Old"),labels=c("Young","Old")))
rownames(annotation)<-colnames(fpkm_glm_genes)

Var1 <- c("deepskyblue","#EF3B39")
names(Var1) <- c("Young", "Old")
anno_colors <- list(age = Var1)

pdf(paste(genelist_name,"microglia_ifn_bst2_clustering.pdf",sep=""),width=6,height=3,onefile=F)
pheatmap(fpkm_glm_genes, col = magma(n = 30, direction = -1), cluster_rows=F,cluster_cols=T, annotation=annotation, annotation_colors = anno_colors,show_colnames=F)
dev.off()



#================================================================================================================================================================================


#Endothelial Heatmap
ident_curr<-c(svz@ident[grepl("Endothelial",svz@ident)])
all_exprs_curr<-svz@data[,match(names(ident_curr),colnames(svz@data))]


greaterthan0<-all_exprs_curr>0
greaterthan0sum<-rowSums(greaterthan0)
all_exprs_curr_genefilt<-all_exprs_curr[greaterthan0sum>=10,]

young1_counts<-all_exprs_curr_genefilt[,grepl("Young1",colnames(all_exprs_curr_genefilt))]
old1_counts<-all_exprs_curr_genefilt[,grepl("Old1",colnames(all_exprs_curr_genefilt))]

young2_counts<-all_exprs_curr_genefilt[,grepl("Young2",colnames(all_exprs_curr_genefilt))]
old2_counts<-all_exprs_curr_genefilt[,grepl("Old2",colnames(all_exprs_curr_genefilt))]

comb_fpkm<-cbind(young1_counts,old1_counts,young2_counts,old2_counts)
#Custom genelist

genelist<-c("B2m","Ifi27","Ifitm3","Ifit3","Ifit1","Xaf1","Bst2")
#genelist<-as.vector(read.table("/Volumes/LaCie/09022016_Dhananjay_10X/Analysis/geneCluster.148.txt"))[,1]
genelist_name<-"custom"
fpkm_glm_genes<-comb_fpkm[na.omit(match(tolower(genelist),tolower(rownames(comb_fpkm)))),]

annotation<-data.frame(age=factor(c(rep("Young",length(young1_counts)),
                                    rep("Old",length(old1_counts)),
                                    rep("Young",length(young2_counts)),
                                    rep("Old",length(old2_counts))),levels=c("Young","Old"),labels=c("Young","Old")))
rownames(annotation)<-colnames(fpkm_glm_genes)

Var1 <- c("deepskyblue","#EF3B39")
names(Var1) <- c("Young", "Old")
anno_colors <- list(age = Var1)

pdf(paste(genelist_name,"endothelial_ifn_Bst2_clustering.pdf",sep=""),width=6,height=3,onefile=F)
pheatmap(fpkm_glm_genes, col = magma(n = 30, direction = -1), cluster_rows=F,cluster_cols=T, annotation=annotation, annotation_colors = anno_colors,show_colnames=F)
dev.off()


