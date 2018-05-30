#=============================================================================================================================================
#IFN High (Fig2b,2c with astrocytes/qNSCs and aNSCs/NPCs separately)


#=============================================================================================================================================
#=============================================================================================================================================
#Astrocytes/qNSCs
rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/")
load("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")
load("hall.kegg.react.list.rda")

ident_curr<-c(svz@ident[grepl("Astrocytes_qNSCs",svz@ident)])
all_exprs_curr<-svz@data[,match(names(ident_curr),colnames(svz@data))]

greaterthan0<-all_exprs_curr>0
greaterthan0sum<-rowSums(greaterthan0)
all_exprs_curr_genefilt<-all_exprs_curr[greaterthan0sum>=5,]

young1_counts<-all_exprs_curr_genefilt[,grepl("Young1",colnames(all_exprs_curr_genefilt))]
old1_counts<-all_exprs_curr_genefilt[,grepl("Old1",colnames(all_exprs_curr_genefilt))]

young2_counts<-all_exprs_curr_genefilt[,grepl("Young2",colnames(all_exprs_curr_genefilt))]
old2_counts<-all_exprs_curr_genefilt[,grepl("Old2",colnames(all_exprs_curr_genefilt))]

comb_fpkm<-cbind(young1_counts,old1_counts,young2_counts,old2_counts)

comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
                 rep("#EF4136",length(old1_counts)),
                 rep("#40BBEC",length(young2_counts)),
                 rep("#EF4136",length(old2_counts)))

gene_sets<-read.table("/Volumes/guacamole/Analyses/NicheAgingSorts/Analysis/GSEALists_Full_01272016.txt",stringsAsFactors=F,header=F,fill=T,sep="\t")

gene_sets_2<-gene_sets[grepl("HALLMARK_INTERFERON",gene_sets[,1]),]
#gene_sets_3<-gene_sets[grepl("IFN",gene_sets[,1]),]
gene_sets<-rbind(gene_sets_2)

#gene_sets<-read.table("/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/mattinterferonlists.gmt.txt",stringsAsFactors=F,header=F,fill=T,sep="\t")
#gene_sets<-rbind.fill(data.frame(c("HALLMARK_INTERFERON_GAMMA_RESPONSE",hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE)))

setwd("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/Pathway_PCA_Violins/New/")

genelist<-hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE
genelist_name<-"HALLMARK_INTERFERON_GAMMA_RESPONSE"
fpkm_glm_genes<-comb_fpkm[na.omit(match(tolower(genelist),tolower(rownames(comb_fpkm)))),]
if(is.null(dim(fpkm_glm_genes))){
  next
}else if(dim(fpkm_glm_genes)[1]<5){
  next
}else{
  #PCA
  comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
                   rep("#EF4136",length(old1_counts)),
                   rep("#40BBEC",length(young2_counts)),
                   rep("#EF4136",length(old2_counts)))
  
  
  #int_fpkm_glm_sum<-colSums(fpkm_glm_genes)
  int_fpkm_glm_sum<-apply(fpkm_glm_genes,2,mean)
  ifn_high<-names(int_fpkm_glm_sum)[int_fpkm_glm_sum>(mean(int_fpkm_glm_sum)+2.2*sd(int_fpkm_glm_sum))]
  
  comb_fpkm_fill<-comb_fpkm_col
  
  comb_fpkm_col[match(ifn_high,colnames(fpkm_glm_genes))]<-"#ff9f00"
  comb_fpkm_col<-factor(comb_fpkm_col)
  
  PCA_int<-prcomp(t(fpkm_glm_genes), scale = T, center = T, retx=T)
  PCA_results<-PCA_int$x
  summa<-summary(PCA_int)
  col<-comb_fpkm_col
  fill<-comb_fpkm_fill
  data_1<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(fpkm_glm_genes),col=col,fill=fill)
  data_1$factors<-as.character(data_1$factors)
  data_1$factors <- factor(data_1$factors, levels=unique(data_1$factors), ordered = T)
  p<-ggplot(data_1)
  p<-p+theme_classic()
  p<-p+geom_point(aes(x=data_1$x,y=data_1$y),fill=data_1$fill,color=data_1$col,size=6,alpha=0.4,shape=21)
  p<- p+ labs(y = paste("PC2  (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[3,1],digits = 2)*100,"% of Variance)", sep = ""))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title.x=element_text(size=24))
  p<-p+theme(axis.text.y=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(size=24))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))
  pdf(paste("Astrocytes_qNSCs_PCA_withprimaryclustering_",genelist_name,".pdf",sep=""),height=5,width=5)
  print(p)
  dev.off()
  
  comb_fpkm_col<-c(rep("Young",length(young1_counts)),
                   rep("Old",length(old1_counts)),
                   rep("Young",length(young2_counts)),
                   rep("Old",length(old2_counts)))
  
  comb_fpkm_col2<-c(rep("#000000",length(comb_fpkm_col)))
  comb_fpkm_col2[match(ifn_high,colnames(fpkm_glm_genes))]<-"#ff9f00"
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=int_fpkm_glm_sum,col=comb_fpkm_col2,fac=comb_fpkm_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_violin(aes(x=dataframe$fac,fill=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.8,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),color=dataframe$col,position = position_jitter(width = .2,height=0),alpha=0.5)
  p<-p+ theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p + theme(legend.position="none")
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="IFN Response Gene Expression")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=NULL)
  p<-p+scale_fill_manual(values=c("#40BBEC","#EF4136"))
  
  pdf(paste(genelist_name,"_Astrocytes_qNSCs_withprimaryclustering_InterferonPathway_Violins.pdf",sep=""),width=4,height=5)
  print(p)
  dev.off()
  
  comb_fpkm_col2<-c(rep("#000000",length(comb_fpkm_col)))
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=int_fpkm_glm_sum,col=comb_fpkm_col2,fac=comb_fpkm_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_violin(aes(x=dataframe$fac,fill=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.8,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),color=dataframe$col,position = position_jitter(width = .2,height=0),alpha=0.5)
  p<-p+ theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p + theme(legend.position="none")
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="IFN Response Gene Expression")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=NULL)
  p<-p+scale_fill_manual(values=c("#40BBEC","#EF4136"))
  
  pdf(paste(genelist_name,"_Astrocytes_qNSCs_withprimaryclustering_InterferonPathway_Violins_noifnhighmarking.pdf",sep=""),width=4,height=5)
  print(p)
  dev.off()
  
}















#=============================================================================================================================================
#=============================================================================================================================================
#aNSCs_NPCs
rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/")
load("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")
load("hall.kegg.react.list.rda")


ident_curr<-c(svz@ident[grepl("aNSCs_NPCs",svz@ident)])
all_exprs_curr<-svz@data[,match(names(ident_curr),colnames(svz@data))]

greaterthan0<-all_exprs_curr>0
greaterthan0sum<-rowSums(greaterthan0)
all_exprs_curr_genefilt<-all_exprs_curr[greaterthan0sum>=5,]

young1_counts<-all_exprs_curr_genefilt[,grepl("Young1",colnames(all_exprs_curr_genefilt))]
old1_counts<-all_exprs_curr_genefilt[,grepl("Old1",colnames(all_exprs_curr_genefilt))]

young2_counts<-all_exprs_curr_genefilt[,grepl("Young2",colnames(all_exprs_curr_genefilt))]
old2_counts<-all_exprs_curr_genefilt[,grepl("Old2",colnames(all_exprs_curr_genefilt))]

comb_fpkm<-cbind(young1_counts,old1_counts,young2_counts,old2_counts)

#comb_fpkm_col<-factor(c(rep("#0099ff",length(young_counts)),rep("#575757",length(old_counts))))
# comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
#                  rep("darkorange",length(old1_counts)),
#                  rep("slateblue",length(young2_counts)),
#                  rep("#EF4136",length(old2_counts)))

comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
                 rep("#EF4136",length(old1_counts)),
                 rep("#40BBEC",length(young2_counts)),
                 rep("#EF4136",length(old2_counts)))


setwd("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/Pathway_PCA_Violins/New/")

genelist<-hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE
genelist_name<-"HALLMARK_INTERFERON_GAMMA_RESPONSE"
fpkm_glm_genes<-comb_fpkm[na.omit(match(tolower(genelist),tolower(rownames(comb_fpkm)))),]


if(is.null(dim(fpkm_glm_genes))){
  next
}else if(dim(fpkm_glm_genes)[1]<5){
  next
}else{
  #PCA
  comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
                   rep("#EF4136",length(old1_counts)),
                   rep("#40BBEC",length(young2_counts)),
                   rep("#EF4136",length(old2_counts)))
  
  
  #int_fpkm_glm_sum<-colSums(fpkm_glm_genes)
  int_fpkm_glm_sum<-apply(fpkm_glm_genes,2,mean)
  ifn_high<-names(int_fpkm_glm_sum)[int_fpkm_glm_sum>(mean(int_fpkm_glm_sum)+2.2*sd(int_fpkm_glm_sum))]
  
  comb_fpkm_fill<-comb_fpkm_col
  
  comb_fpkm_col[match(ifn_high,colnames(fpkm_glm_genes))]<-"#ff9f00"
  comb_fpkm_col<-factor(comb_fpkm_col)
  
  PCA_int<-prcomp(t(fpkm_glm_genes), scale = T, center = T, retx=T)
  PCA_results<-PCA_int$x
  summa<-summary(PCA_int)
  col<-comb_fpkm_col
  fill<-comb_fpkm_fill
  data_1<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(fpkm_glm_genes),col=col,fill=fill)
  data_1$factors<-as.character(data_1$factors)
  data_1$factors <- factor(data_1$factors, levels=unique(data_1$factors), ordered = T)
  p<-ggplot(data_1)
  p<-p+theme_classic()
  p<-p+geom_point(aes(x=data_1$x,y=data_1$y),fill=data_1$fill,color=data_1$col,size=6,alpha=0.4,shape=21)
  p<- p+ labs(y = paste("PC2  (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[3,1],digits = 2)*100,"% of Variance)", sep = ""))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title.x=element_text(size=24))
  p<-p+theme(axis.text.y=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(size=24))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))
  pdf(paste("aNSCs_NPCs_PCA_withprimaryclustering_",genelist_name,".pdf",sep=""),height=5,width=5)
  print(p)
  dev.off()
  
  comb_fpkm_col<-c(rep("Young",length(young1_counts)),
                   rep("Old",length(old1_counts)),
                   rep("Young",length(young2_counts)),
                   rep("Old",length(old2_counts)))
  
  comb_fpkm_col2<-c(rep("#000000",length(comb_fpkm_col)))
  comb_fpkm_col2[match(ifn_high,colnames(fpkm_glm_genes))]<-"#ff9f00"
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=int_fpkm_glm_sum,col=comb_fpkm_col2,fac=comb_fpkm_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_violin(aes(x=dataframe$fac,fill=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.8,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),color=dataframe$col,position = position_jitter(width = .2,height=0),alpha=0.5)
  p<-p+ theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p + theme(legend.position="none")
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="IFN Response Gene Expression")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=NULL)
  p<-p+scale_fill_manual(values=c("#40BBEC","#EF4136"))
  
  pdf(paste(genelist_name,"_aNSCs_NPCs_withprimaryclustering_InterferonPathway_Violins.pdf",sep=""),width=4,height=5)
  print(p)
  dev.off()
  
  comb_fpkm_col2<-c(rep("#000000",length(comb_fpkm_col)))
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=int_fpkm_glm_sum,col=comb_fpkm_col2,fac=comb_fpkm_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_violin(aes(x=dataframe$fac,fill=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.8,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),color=dataframe$col,position = position_jitter(width = .2,height=0),alpha=0.5)
  p<-p+ theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p + theme(legend.position="none")
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="IFN Response Gene Expression")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=NULL)
  p<-p+scale_fill_manual(values=c("#40BBEC","#EF4136"))
  
  pdf(paste(genelist_name,"_aNSCs_NPCs_withprimaryclustering_InterferonPathway_Violins_noifnhighmarking.pdf",sep=""),width=4,height=5)
  print(p)
  dev.off()
  
}





#==================================================
#Two cell types together (aNSCs/NPCS / Astrocytes NPCs)
rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/")
load("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/svz_tsne_young_old_tworeps_celllabeled_clusters.Robj")
load("hall.kegg.react.list.rda")

ident_curr<-c(svz@ident[grepl("Astrocytes_qNSCs",svz@ident)],svz@ident[grepl("aNSCs_NPCs",svz@ident)])
all_exprs_curr<-svz@data[,match(names(ident_curr),colnames(svz@data))]

greaterthan0<-all_exprs_curr>0
greaterthan0sum<-rowSums(greaterthan0)
all_exprs_curr_genefilt<-all_exprs_curr[greaterthan0sum>=5,]

young1_counts<-all_exprs_curr_genefilt[,grepl("Young1",colnames(all_exprs_curr_genefilt))]
old1_counts<-all_exprs_curr_genefilt[,grepl("Old1",colnames(all_exprs_curr_genefilt))]

young2_counts<-all_exprs_curr_genefilt[,grepl("Young2",colnames(all_exprs_curr_genefilt))]
old2_counts<-all_exprs_curr_genefilt[,grepl("Old2",colnames(all_exprs_curr_genefilt))]

comb_fpkm<-cbind(young1_counts,old1_counts,young2_counts,old2_counts)

#comb_fpkm_col<-factor(c(rep("#0099ff",length(young_counts)),rep("#575757",length(old_counts))))
# comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
#                  rep("darkorange",length(old1_counts)),
#                  rep("slateblue",length(young2_counts)),
#                  rep("#EF4136",length(old2_counts)))

comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
                 rep("#EF4136",length(old1_counts)),
                 rep("#40BBEC",length(young2_counts)),
                 rep("#EF4136",length(old2_counts)))

setwd("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/Pathway_PCA_Violins/New/")

genelist<-hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE
genelist_name<-"HALLMARK_INTERFERON_GAMMA_RESPONSE"
fpkm_glm_genes<-comb_fpkm[na.omit(match(tolower(genelist),tolower(rownames(comb_fpkm)))),]

if(is.null(dim(fpkm_glm_genes))){
  next
}else if(dim(fpkm_glm_genes)[1]<5){
  next
}else{
  #PCA
  comb_fpkm_col<-c(rep("#40BBEC",length(young1_counts)),
                   rep("#EF4136",length(old1_counts)),
                   rep("#40BBEC",length(young2_counts)),
                   rep("#EF4136",length(old2_counts)))
  
  
  #int_fpkm_glm_sum<-colSums(fpkm_glm_genes)
  int_fpkm_glm_sum<-apply(fpkm_glm_genes,2,mean)
  ifn_high<-names(int_fpkm_glm_sum)[int_fpkm_glm_sum>(mean(int_fpkm_glm_sum)+2.2*sd(int_fpkm_glm_sum))]
  
  comb_fpkm_fill<-comb_fpkm_col
  
  comb_fpkm_col[match(ifn_high,colnames(fpkm_glm_genes))]<-"#ff9f00"
  comb_fpkm_col<-factor(comb_fpkm_col)
  
  PCA_int<-prcomp(t(fpkm_glm_genes), scale = T, center = T, retx=T)
  PCA_results<-PCA_int$x
  summa<-summary(PCA_int)
  col<-comb_fpkm_col
  fill<-comb_fpkm_fill
  data_1<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(fpkm_glm_genes),col=col,fill=fill)
  data_1$factors<-as.character(data_1$factors)
  data_1$factors <- factor(data_1$factors, levels=unique(data_1$factors), ordered = T)
  p<-ggplot(data_1)
  p<-p+theme_classic()
  p<-p+geom_point(aes(x=data_1$x,y=data_1$y),fill=data_1$fill,color=data_1$col,size=6,alpha=0.4,shape=21)
  p<- p+ labs(y = paste("PC2  (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[3,1],digits = 2)*100,"% of Variance)", sep = ""))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title.x=element_text(size=24))
  p<-p+theme(axis.text.y=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(size=24))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))
  pdf(paste("Astrocytes_qNSCs_aNSCs_NPCs_PCA_withprimaryclustering_",genelist_name,".pdf",sep=""),height=5,width=5)
  print(p)
  dev.off()
  
  comb_fpkm_col<-c(rep("Young",length(young1_counts)),
                   rep("Old",length(old1_counts)),
                   rep("Young",length(young2_counts)),
                   rep("Old",length(old2_counts)))
  
  comb_fpkm_col2<-c(rep("#000000",length(comb_fpkm_col)))
  comb_fpkm_col2[match(ifn_high,colnames(fpkm_glm_genes))]<-"#ff9f00"
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=int_fpkm_glm_sum,col=comb_fpkm_col2,fac=comb_fpkm_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_violin(aes(x=dataframe$fac,fill=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.8,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),color=dataframe$col,position = position_jitter(width = .2,height=0),alpha=0.5)
  p<-p+ theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p + theme(legend.position="none")
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="IFN Response Gene Expression")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=NULL)
  p<-p+scale_fill_manual(values=c("#40BBEC","#EF4136"))
  
  pdf(paste(genelist_name,"_withprimaryclustering_InterferonPathway_Violins.pdf",sep=""),width=4,height=5)
  print(p)
  dev.off()
  
  comb_fpkm_col2<-c(rep("#000000",length(comb_fpkm_col)))
  
  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=int_fpkm_glm_sum,col=comb_fpkm_col2,fac=comb_fpkm_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_violin(aes(x=dataframe$fac,fill=dataframe$fac,y=dataframe[,1]),trim=T,scale="width",alpha=0.8,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),color=dataframe$col,position = position_jitter(width = .2,height=0),alpha=0.5)
  p<-p+ theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p + theme(legend.position="none")
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="IFN Response Gene Expression")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=NULL)
  p<-p+scale_fill_manual(values=c("#40BBEC","#EF4136"))
  
  pdf(paste(genelist_name,"_withprimaryclustering_InterferonPathway_Violins_noifnhighmarking.pdf",sep=""),width=4,height=5)
  print(p)
  dev.off()
}

#==================================================
#Fluidigm C1

rm(list=ls())
library(ggplot2)
library(edgeR)
load("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/hall.kegg.react.list.rda")


#PCA of aging cells with all oligodendrocytes removed and day 4 experiments removed

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Analyses/NicheAgingSorts/Analysis/")
aging_allcells<-read.table("All_highqual_>500genes_agingCells.txt")

setwd("/Volumes/guacamole/Analyses/NicheAgingSorts/Analysis/")
oligos<-as.vector(read.table("Aging_oligos_updated_01202016.txt")[,1])
aging_allcells_nooligo<-aging_allcells[,-na.omit(match(oligos,colnames(aging_allcells)))]
aging_allcells_nooligo<-aging_allcells_nooligo[,!grepl("4_",colnames(aging_allcells_nooligo))]

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-aging_allcells_nooligo>10
greaterthan0sum<-rowSums(greaterthan0)
aging_allcells_nooligo_genefilt<-aging_allcells_nooligo[greaterthan0sum>=5,]
write.table(aging_allcells_nooligo_genefilt,"/Volumes/guacamole/Analyses/NicheAgingSorts/Analysis/FluidigmC1_rawCounts_genefilt_forPaper.txt",quote=F)
aging_allcells_nooligo_genefilt<-aging_allcells_nooligo_genefilt[!grepl("ERCC-",rownames(aging_allcells_nooligo_genefilt)),]

glm <- DGEList(counts=aging_allcells_nooligo_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Annotations/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)


setwd("/Volumes/guacamole/Analyses/NicheAgingSorts/Analysis/")
write.table(fpkm_glm_genefilt_nooligo,"FluidigmC1_fpkm_forPaper.txt",quote=F)
#PCA using differentially expressed genes with age
Young_Old_col<-vector(length=length(fpkm_glm_genefilt_nooligo[1,]),mode="character")
Young_Old_col[grepl("Y",colnames(fpkm_glm_genefilt_nooligo))]<-"#40BBEC"
Young_Old_col[grepl("O",colnames(fpkm_glm_genefilt_nooligo))]<-"#EF4136"

setwd("/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Analysis_twoperfrep/Pathway_PCA_Violins/New/")

genelist<-hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE
genelist_name<-"HALLMARK_INTERFERON_GAMMA_RESPONSE"

fpkm_glm_genes<-fpkm_glm_genefilt_nooligo[na.omit(match(tolower(genelist),tolower(rownames(fpkm_glm_genefilt_nooligo)))),]

PCA_int<-prcomp(t(log(fpkm_glm_genes+1)), scale = T, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
col<-Young_Old_col 
data_1<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(fpkm_glm_genes),col=col)
data_1$factors<-as.character(data_1$factors)
data_1$factors <- factor(data_1$factors, levels=unique(data_1$factors), ordered = T)
p<-ggplot(data_1)
p<-p+theme_classic()
#p<-p+geom_text(aes(x=data_1$x,y=data_1$y),label=data_1$factors,color=data_1$col)
p<-p+geom_point(aes(x=data_1$x,y=data_1$y),color=data_1$col,size=6,alpha=0.6)
p<- p+ labs(y = paste("PC2  (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[3,1],digits = 2)*100,"% of Variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
pdf(paste("05212018_PCA_",genelist_name,".pdf",sep=""),height=5,width=5)
print(p)
dev.off()

Young_Old_col<-vector(length=length(fpkm_glm_genefilt_nooligo[1,]),mode="character")
Young_Old_col[grepl("Y",colnames(fpkm_glm_genefilt_nooligo))]<-"#40BBEC"
Young_Old_col[grepl("O",colnames(fpkm_glm_genefilt_nooligo))]<-"#EF4136"

int_fpkm_glm_sum<-apply(log2(fpkm_glm_genes+1),2,mean)
ifn_high<-names(int_fpkm_glm_sum)[int_fpkm_glm_sum>(mean(int_fpkm_glm_sum)+2.2*sd(int_fpkm_glm_sum))]

Young_Old_fill<-Young_Old_col
Young_Old_col[match(ifn_high,colnames(fpkm_glm_genes))]<-"#ff9f00"
Young_Old_col<-factor(Young_Old_col)

PCA_int<-prcomp(t(log(fpkm_glm_genes+1)), scale = T, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
col<-Young_Old_col 
fill<-Young_Old_fill
data_1<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(fpkm_glm_genes),col=col,fill=fill)
data_1$factors<-as.character(data_1$factors)
data_1$factors <- factor(data_1$factors, levels=unique(data_1$factors), ordered = T)
p<-ggplot(data_1)
p<-p+theme_classic()
#p<-p+geom_text(aes(x=data_1$x,y=data_1$y),label=data_1$factors,color=data_1$col)
p<-p+geom_point(aes(x=data_1$x,y=data_1$y),color=data_1$col,fill=data_1$fill,size=6,alpha=0.6,shape=21)
#p<-p+geom_point(aes(x=data_1$x,y=data_1$y),color=data_1$col,shape=1,size=6,alpha=0.6)
p<- p+ labs(y = paste("PC2  (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[3,1],digits = 2)*100,"% of Variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=20))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
setwd(dir)
pdf(paste("05212018_PCA_ifnmarked",genelist_name,".pdf",sep=""),height=5,width=5)
print(p)
dev.off()

