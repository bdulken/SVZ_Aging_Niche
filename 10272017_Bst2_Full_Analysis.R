rm(list=ls())
library(edgeR)
library(ggplot2)
library(pvclust)

source("/Volumes/guacamole/Software/R_Files_Packages_Functions/KEGG_pathway_enrichment.R")
reps=100000
pvalue_converter<-function(table){
  signed_pval<-c()
  signed_lr<-c()
  for(i in 1:length(table[,1])){
    if(table[i,1]>0){
      signed_pval<-c(signed_pval,(10000-table[i,4]*10000))
      signed_lr<-c(signed_lr,table[i,3])
    }else{
      signed_pval<-c(signed_pval,-(10000-table[i,4]*10000))
      signed_lr<-c(signed_lr,-table[i,3])
    }
  }
  return(cbind(table,signed_pval,signed_lr))
}

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/LaCie/10242017_Bst2_full_RNA_seq/Analysis/")
allpops<-read.table("Bst2_Full_allcounts_new_names.txt")
allpops_exp2<-allpops[,grepl("Exp2",colnames(allpops))]


allpops_genecounts<-colSums(allpops!=0)
#Filtering for expressed genes 
greaterthan0<-allpops>5
greaterthan0sum<-rowSums(greaterthan0)
allpops_genefilt<-allpops[greaterthan0sum>=2,]
allpops_genefilt<-allpops_genefilt[!grepl("ERCC-",rownames(allpops_genefilt)),]

allpops_genecounts<-colSums(allpops_exp2!=0)
#Filtering for expressed genes 
greaterthan0<-allpops_exp2>5
greaterthan0sum<-rowSums(greaterthan0)
allpops_genefilt_exp2<-allpops_exp2[greaterthan0sum>=2,]
allpops_genefilt_exp2<-allpops_genefilt_exp2[!grepl("ERCC-",rownames(allpops_genefilt_exp2)),]


glm <- DGEList(counts=allpops_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/mm_genomes/mm10/mm10.genes.current.gene.sizes.R")
temp<-(as.numeric(exonic.gene.sizes))
names(temp)<-tolower(names(exonic.gene.sizes))
exonic.gene.sizes<-temp
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt <- rpkm(glm.norm, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

glm <- DGEList(counts=allpops_genefilt_exp2)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/mm_genomes/mm10/mm10.genes.current.gene.sizes.R")
temp<-(as.numeric(exonic.gene.sizes))
names(temp)<-tolower(names(exonic.gene.sizes))
exonic.gene.sizes<-temp
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_exp2 <- rpkm(glm.norm, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)


comb_fpkm<-fpkm_glm_genefilt
comb_fac<-c()
for(i in 1:length(colnames(comb_fpkm))){
  temp<-strsplit(colnames(comb_fpkm)[i],split='[_]')[[1]][1]
  comb_fac<-c(comb_fac,temp)
}

comb_fpkm_exp2<-fpkm_glm_genefilt_exp2
comb_fac_exp2<-c()
for(i in 1:length(colnames(comb_fpkm_exp2))){
  temp<-strsplit(colnames(comb_fpkm_exp2)[i],split='[_]')[[1]][1]
  comb_fac_exp2<-c(comb_fac_exp2,temp)
}
#PCA with all detected genes, no ERCC
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = T, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm),col=comb_fac)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y,color=data$col),size=5,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p

pdf("Bst2_Full_Allsamples_PCA.pdf",height=8,width=8)
print(p)
dev.off()

#PCA with all detected genes, no ERCC
PCA_int<-prcomp(t(log2(comb_fpkm_exp2+1)), scale = T, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm_exp2),col=comb_fac_exp2)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y,color=data$col),size=5,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p

pdf("Bst2_Full_Exp2_PCA.pdf",height=8,width=8)
print(p)
dev.off()

comb_fpkm_exp2_noqNSC<-comb_fpkm_exp2[,!grepl("qNSC",colnames(comb_fpkm_exp2))]
comb_fac_exp2_noqNSC<-comb_fac_exp2[!grepl("qNSC",comb_fac_exp2)]
#PCA with all detected genes, no ERCC
PCA_int<-prcomp(t(log2(comb_fpkm_exp2_noqNSC+1)), scale=F, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm_exp2_noqNSC),col=comb_fac_exp2_noqNSC)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y,color=data$col),size=5,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p

pdf("Bst2_noqNSC_Exp2_PCA.pdf",height=8,width=8)
print(p)
dev.off()

col<-c(rep("red",4),rep("green",4),rep("blue",4))
library(scatterplot3d)
pdf("All_miseq_PCA_pdf_3d.pdf",height=8,width=8)
scatterplot3d(PCA_results[,1:3],color=alpha(col,0.5),angle=45,pch=16,cex.symbols=2.0,cex.axis=1.5,cex.lab=2.0,xlab=paste("PC1    (",round(summa$importance[2,1],digits = 3)*100,"% of variance)", sep = ""),
              ylab=paste("PC2    (",round(summa$importance[2,2],digits = 3)*100,"% of variance)", sep = ""),
              zlab=paste("PC3    (",round(summa$importance[2,3],digits = 3)*100,"% of variance)", sep = ""))
dev.off()

#spearman dist function, this is the function that pvclust will call as its distance method
spearman <- function(x, ...) {
  x <- as.matrix(x)
  res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
  res <- as.dist(res)
  attr(res, "method") <- "spearman"
  return(res)
}
# #Cluster with pg and >5 expressing genes
# clust<-pvclust(comb_fpkm,method.dist=spearman,nboot=1000)
# pdf("Bst2_full_cluster_spearman.pdf", width = 10)
# plot(clust)
# dev.off()



#============================================================================================================================================================
#Differential Expression

library(edgeR)
library(locfit)
glm <- DGEList(counts=allpops_genefilt_exp2)
glm.norm <- calcNormFactors(glm,method="TMM")
comb_fac_exp2<-c()
for(i in 1:length(colnames(allpops_genefilt_exp2))){
  temp<-strsplit(colnames(allpops_genefilt_exp2)[i],split='[_]')[[1]][1]
  comb_fac_exp2<-c(comb_fac_exp2,temp)
}
design<-model.matrix(~0+comb_fac_exp2)
y <- estimateDisp(glm.norm, design)

fit <- glmFit(y, design)


#=========================================================================================================
#Old Bst2pos vs Bst2neg comparison
lrt.2vs1 <- glmLRT(fit,contrast=c(-1,1,0,0,0,0))

toptags_sort<-topTags(lrt.2vs1, p.value=0.05, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
write.table(toptags_sort,"OldBst2pos_vs_OldBst2neg.txt")

tests<-decideTestsDGE(lrt.2vs1)

tests_neg<-lrt.2vs1$table[tests==-1,]
tests_neg_sort<-tests_neg[order(tests_neg[,4]),]
tests_pos<-lrt.2vs1$table[tests==1,]
tests_pos_sort<-tests_pos[order(tests_pos[,4]),]

write.table(tests_neg_sort,"OldBst2pos_vs_OldBst2neg_OldBst2Neg_tests.txt")
write.table(tests_pos_sort,"OldBst2pos_vs_OldBst2neg_OldBst2Pos_tests.txt")
write.table(rownames(tests_neg_sort),"OldBst2pos_vs_OldBst2neg_OldBst2Neg_names.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(tests_pos_sort),"OldBst2pos_vs_OldBst2neg_OldBst2Pos_names.txt",quote=F,row.names=F,col.names=F)

alltags_sort<-topTags(lrt.2vs1, p.value=1, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
OldBst2pos.v.OldBst2neg.conv<-pvalue_converter(alltags_sort$table)
OldBst2pos.v.OldBst2neg.conv.order<-OldBst2pos.v.OldBst2neg.conv[order(OldBst2pos.v.OldBst2neg.conv[,7],decreasing=T),]
write.table(OldBst2pos.v.OldBst2neg.conv.order,"OldBst2pos.v.OldBst2neg.order.txt")
rnk<-cbind(toupper(rownames(OldBst2pos.v.OldBst2neg.conv.order)),OldBst2pos.v.OldBst2neg.conv.order[,7])
write.table(rnk,"OldBst2pos.v.OldBst2neg.order.exp2.signed.lr.rnk",quote=F,row.names=F,col.names=F,sep='\t')
OldBst2pos.v.OldBst2neg.conv.order<-OldBst2pos.v.OldBst2neg.conv[order(OldBst2pos.v.OldBst2neg.conv[,6],decreasing=T),]
rnk<-cbind(toupper(rownames(OldBst2pos.v.OldBst2neg.conv.order)),OldBst2pos.v.OldBst2neg.conv.order[,6])
write.table(rnk,"OldBst2pos.v.OldBst2neg.order.exp2.signed.pval.rnk",quote=F,row.names=F,col.names=F,sep='\t')


upper<-toupper(rownames(OldBst2pos.v.OldBst2neg.conv.order))
rownames(OldBst2pos.v.OldBst2neg.conv.order)<-upper
OldBst2pos.v.OldBst2neg.GO<-KEGG.pathway.enrichment(test.results=OldBst2pos.v.OldBst2neg.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c5.all.v5.2.symbols.gmt")
OldBst2pos.v.OldBst2neg.GO.sig<-OldBst2pos.v.OldBst2neg.GO$pathway.results[abs(OldBst2pos.v.OldBst2neg.GO$pathway.results[,2])<0.25,]
write.table(OldBst2pos.v.OldBst2neg.GO.sig,"OldBst2pos.v.OldBst2neg.GO.sig.txt")

upper<-toupper(rownames(OldBst2pos.v.OldBst2neg.conv.order))
rownames(OldBst2pos.v.OldBst2neg.conv.order)<-upper
OldBst2pos.v.OldBst2neg.GO<-KEGG.pathway.enrichment(test.results=OldBst2pos.v.OldBst2neg.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c2.cp.kegg.v5.2.symbols.gmt")
OldBst2pos.v.OldBst2neg.GO.sig<-OldBst2pos.v.OldBst2neg.GO$pathway.results[abs(OldBst2pos.v.OldBst2neg.GO$pathway.results[,2])<0.25,]
write.table(OldBst2pos.v.OldBst2neg.GO.sig,"OldBst2pos.v.OldBst2neg.KEGG.sig.txt")

upper<-toupper(rownames(OldBst2pos.v.OldBst2neg.conv.order))
rownames(OldBst2pos.v.OldBst2neg.conv.order)<-upper
OldBst2pos.v.OldBst2neg.GO<-KEGG.pathway.enrichment(test.results=OldBst2pos.v.OldBst2neg.conv.order,fold.change = "logFC",statistic="signed_lr",num.samples=reps,gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/h.all.v5.2.symbols.gmt")
OldBst2pos.v.OldBst2neg.GO.sig<-OldBst2pos.v.OldBst2neg.GO$pathway.results[abs(OldBst2pos.v.OldBst2neg.GO$pathway.results[,2])<0.25,]
write.table(OldBst2pos.v.OldBst2neg.GO.sig,"OldBst2pos.v.OldBst2neg.Hallmark.sig.txt")

roundsize<-c()
for(i in 1:length(OldBst2pos.v.OldBst2neg.GO.sig[,1])){
  roundsize<-c(roundsize,floor(min(log10(reps),-log10(OldBst2pos.v.OldBst2neg.GO.sig$adj.p.val[i]))))
}
data<-data.frame(names=factor(rownames(OldBst2pos.v.OldBst2neg.GO.sig),
                              levels=rownames(OldBst2pos.v.OldBst2neg.GO.sig)[length(rownames(OldBst2pos.v.OldBst2neg.GO.sig)):1],ordered=T),fc=OldBst2pos.v.OldBst2neg.GO.sig$fold.change,size=roundsize)

p<-ggplot(data)
p<-p+geom_point(aes(x=data$names,y=data$fc,size=data$size))
p<-p+theme_classic()
p<-p+coord_flip()
p<-p+theme(legend.position = "none") 
p<-p+labs(y="Enrichment Score",title=NULL,x=NULL)


pdf("OldBst2pos.v.OldBst2neg.Hallmark.sig.enrichment.pdf")
print(p)
dev.off()


enrichment<-OldBst2pos.v.OldBst2neg.GO.sig<-OldBst2pos.v.OldBst2neg.GO$pathway.results[abs(OldBst2pos.v.OldBst2neg.GO$pathway.results[,3])<0.25,]
color<-vector()
enrich<-vector()
names<-vector()
for(i in 1:length(enrichment[,4])){
  if(enrichment[i,4]<0){
    enrich_temp<-enrichment[i,4]
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#ff3333")
  }else{
    enrich_temp<-enrichment[i,4]
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#00CCFF")
  }
}

namevar1<-gsub("HALLMARK_","",rownames(OldBst2pos.v.OldBst2neg.GO.sig))
namevar2<-gsub("_"," ",namevar1)


data1<-data.frame(names=namevar2,enrich=enrich,color=color)
data1<-data1[order(data1$enrich,decreasing=F),]
data1$names<-as.character(data1$names)
data1$names <- factor(data1$names, levels=unique(data1$names), ordered = T)
data1$color<-as.character(data1$color)
data1$color <- factor(data1$color, levels=unique(data1$color), ordered = T)
p<-ggplot(data1)
p<-p+geom_bar(aes(x=data1$names,y=data1$enrich,fill=data1$color),stat="identity")
p<-p+theme_classic()
p<-p+coord_flip()+labs(y="-log10(Test Statistic)",x=NULL)
p<-p+theme(axis.text.x=element_text(size=15))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=15,face="bold"))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
#p<-p+scale_y_continuous(limits=c(-4,4),breaks=c(-4,-2,0,2,4), labels=c("4", "2", "0","2","4"))
p<-p+scale_fill_manual(values=c("orange","firebrick"))
p<-p+theme(legend.position = "none") 
p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))
p
pdf("OldBst2pos.v.OldBst2neg.Hallmark.sig.enrichment_modpresentation_exp2.pdf",height=6,width=7.2)
print(p)
dev.off()

upper<-toupper(rownames(OldBst2pos.v.OldBst2neg.conv.order))
rownames(OldBst2pos.v.OldBst2neg.conv.order)<-upper
OldBst2pos.v.OldBst2neg.GO<-KEGG.pathway.enrichment(test.results=OldBst2pos.v.OldBst2neg.conv.order,fold.change = "logFC",statistic="signed_lr",num.samples=reps,gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/h.all.v5.2.withmattslists.symbols.gmt")
OldBst2pos.v.OldBst2neg.GO.sig<-OldBst2pos.v.OldBst2neg.GO$pathway.results[abs(OldBst2pos.v.OldBst2neg.GO$pathway.results[,3])<0.25,]
write.table(OldBst2pos.v.OldBst2neg.GO.sig,"OldBst2pos.v.OldBst2neg.Hallmark.withmattslists.sig.txt")


roundsize<-c()
for(i in 1:length(OldBst2pos.v.OldBst2neg.GO.sig[,1])){
  roundsize<-c(roundsize,floor(min(log10(reps),-log10(OldBst2pos.v.OldBst2neg.GO.sig$adj.p.val[i]))))
}
data<-data.frame(names=factor(rownames(OldBst2pos.v.OldBst2neg.GO.sig),
                              levels=rownames(OldBst2pos.v.OldBst2neg.GO.sig)[length(rownames(OldBst2pos.v.OldBst2neg.GO.sig)):1],ordered=T),fc=OldBst2pos.v.OldBst2neg.GO.sig$fold.change,size=roundsize)

p<-ggplot(data)
p<-p+geom_point(aes(x=data$names,y=data$fc,size=data$size))
p<-p+theme_classic()
p<-p+coord_flip()
p<-p+theme(legend.position = "none") 
p<-p+labs(y="Enrichment Score",title=NULL,x=NULL)


pdf("OldBst2pos.v.OldBst2neg.Hallmark.withmattslists.sig.enrichment.pdf")
print(p)
dev.off()

#=========================================================================================================
#Young Bst2pos vs Bst2neg comparison
lrt.2vs1 <- glmLRT(fit,contrast=c(0,0,0,-1,1,0))

toptags_sort<-topTags(lrt.2vs1, p.value=0.05, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
write.table(toptags_sort,"YoungBst2pos_vs_YoungBst2neg.txt")

tests<-decideTestsDGE(lrt.2vs1)

tests_neg<-lrt.2vs1$table[tests==-1,]
tests_neg_sort<-tests_neg[order(tests_neg[,4]),]
tests_pos<-lrt.2vs1$table[tests==1,]
tests_pos_sort<-tests_pos[order(tests_pos[,4]),]

write.table(tests_neg_sort,"YoungBst2pos_vs_YoungBst2neg_YoungBst2Neg_tests.txt")
write.table(tests_pos_sort,"YoungBst2pos_vs_YoungBst2neg_YoungBst2Pos_tests.txt")
write.table(rownames(tests_neg_sort),"YoungBst2pos_vs_YoungBst2neg_YoungBst2Neg_names.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(tests_pos_sort),"YoungBst2pos_vs_YoungBst2neg_YoungBst2Pos_names.txt",quote=F,row.names=F,col.names=F)

alltags_sort<-topTags(lrt.2vs1, p.value=1, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
YoungBst2pos.v.YoungBst2neg.conv<-pvalue_converter(alltags_sort$table)
YoungBst2pos.v.YoungBst2neg.conv.order<-YoungBst2pos.v.YoungBst2neg.conv[order(YoungBst2pos.v.YoungBst2neg.conv[,6],decreasing=T),]
write.table(YoungBst2pos.v.YoungBst2neg.conv.order,"YoungBst2pos.v.YoungBst2neg.order.txt")
rnk<-cbind(toupper(rownames(YoungBst2pos.v.YoungBst2neg.conv.order)),YoungBst2pos.v.YoungBst2neg.conv.order[,6])
write.table(rnk,"YoungBst2pos.v.YoungBst2neg.order.rnk",quote=F,row.names=F,col.names=F,sep='\t')

upper<-toupper(rownames(YoungBst2pos.v.YoungBst2neg.conv.order))
rownames(YoungBst2pos.v.YoungBst2neg.conv.order)<-upper
YoungBst2pos.v.YoungBst2neg.GO<-KEGG.pathway.enrichment(test.results=YoungBst2pos.v.YoungBst2neg.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c5.all.v5.2.symbols.gmt")
YoungBst2pos.v.YoungBst2neg.GO.sig<-YoungBst2pos.v.YoungBst2neg.GO$pathway.results[abs(YoungBst2pos.v.YoungBst2neg.GO$pathway.results[,2])<0.01,]
write.table(YoungBst2pos.v.YoungBst2neg.GO.sig,"YoungBst2pos.v.YoungBst2neg.GO.sig.txt")

upper<-toupper(rownames(YoungBst2pos.v.YoungBst2neg.conv.order))
rownames(YoungBst2pos.v.YoungBst2neg.conv.order)<-upper
YoungBst2pos.v.YoungBst2neg.GO<-KEGG.pathway.enrichment(test.results=YoungBst2pos.v.YoungBst2neg.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c2.cp.kegg.v5.2.symbols.gmt")
YoungBst2pos.v.YoungBst2neg.GO.sig<-YoungBst2pos.v.YoungBst2neg.GO$pathway.results[abs(YoungBst2pos.v.YoungBst2neg.GO$pathway.results[,2])<0.01,]
write.table(YoungBst2pos.v.YoungBst2neg.GO.sig,"YoungBst2pos.v.YoungBst2neg.KEGG.sig.txt")

upper<-toupper(rownames(YoungBst2pos.v.YoungBst2neg.conv.order))
rownames(YoungBst2pos.v.YoungBst2neg.conv.order)<-upper
YoungBst2pos.v.YoungBst2neg.GO<-KEGG.pathway.enrichment(test.results=YoungBst2pos.v.YoungBst2neg.conv.order,fold.change = "logFC",statistic="signed_lr",num.samples=reps,gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/h.all.v5.2.symbols.gmt")
YoungBst2pos.v.YoungBst2neg.GO.sig<-YoungBst2pos.v.YoungBst2neg.GO$pathway.results[abs(YoungBst2pos.v.YoungBst2neg.GO$pathway.results[,2])<0.05,]
write.table(YoungBst2pos.v.YoungBst2neg.GO.sig,"YoungBst2pos.v.YoungBst2neg.Hallmark.sig.txt")

roundsize<-c()
for(i in 1:length(YoungBst2pos.v.YoungBst2neg.GO.sig[,1])){
  roundsize<-c(roundsize,floor(min(log10(reps),-log10(YoungBst2pos.v.YoungBst2neg.GO.sig$adj.p.val[i]))))
}
data<-data.frame(names=factor(rownames(YoungBst2pos.v.YoungBst2neg.GO.sig),
                              levels=rownames(YoungBst2pos.v.YoungBst2neg.GO.sig)[length(rownames(YoungBst2pos.v.YoungBst2neg.GO.sig)):1],ordered=T),fc=YoungBst2pos.v.YoungBst2neg.GO.sig$fold.change,size=roundsize)

p<-ggplot(data)
p<-p+geom_point(aes(x=data$names,y=data$fc,size=data$size))
p<-p+theme_classic()
p<-p+coord_flip()
p<-p+theme(legend.position = "none") 
p<-p+labs(y="Enrichment Score",title=NULL,x=NULL)


pdf("YoungBst2pos.v.YoungBst2neg.Hallmark.sig.enrichment.pdf")
print(p)
dev.off()

#=========================================================================================================
#Young Bst2pos vs Old Bst2pos comparison
lrt.2vs1 <- glmLRT(fit,contrast=c(0,1,0,0,-1,0))

toptags_sort<-topTags(lrt.2vs1, p.value=0.05, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
write.table(toptags_sort,"YoungBst2pos_vs_OldBst2pos.txt")

tests<-decideTestsDGE(lrt.2vs1)

tests_neg<-lrt.2vs1$table[tests==-1,]
tests_neg_sort<-tests_neg[order(tests_neg[,4]),]
tests_pos<-lrt.2vs1$table[tests==1,]
tests_pos_sort<-tests_pos[order(tests_pos[,4]),]

write.table(tests_neg_sort,"YoungBst2pos_vs_OldBst2pos_OldBst2pos_tests.txt")
write.table(tests_pos_sort,"YoungBst2pos_vs_OldBst2pos_YoungBst2Pos_tests.txt")
write.table(rownames(tests_neg_sort),"YoungBst2pos_vs_OldBst2pos_OldBst2pos_names.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(tests_pos_sort),"YoungBst2pos_vs_OldBst2pos_YoungBst2Pos_names.txt",quote=F,row.names=F,col.names=F)

alltags_sort<-topTags(lrt.2vs1, p.value=1, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
YoungBst2pos.v.OldBst2pos.conv<-pvalue_converter(alltags_sort$table)
YoungBst2pos.v.OldBst2pos.conv.order<-YoungBst2pos.v.OldBst2pos.conv[order(YoungBst2pos.v.OldBst2pos.conv[,6],decreasing=T),]
write.table(YoungBst2pos.v.OldBst2pos.conv.order,"YoungBst2pos.v.OldBst2pos.order.txt")
rnk<-cbind(toupper(rownames(YoungBst2pos.v.OldBst2pos.conv.order)),YoungBst2pos.v.OldBst2pos.conv.order[,6])
write.table(rnk,"YoungBst2pos.v.OldBst2pos.order.rnk",quote=F,row.names=F,col.names=F,sep='\t')

upper<-toupper(rownames(YoungBst2pos.v.OldBst2pos.conv.order))
rownames(YoungBst2pos.v.OldBst2pos.conv.order)<-upper
YoungBst2pos.v.OldBst2pos.GO<-KEGG.pathway.enrichment(test.results=YoungBst2pos.v.OldBst2pos.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c5.all.v5.2.symbols.gmt")
YoungBst2pos.v.OldBst2pos.GO.sig<-YoungBst2pos.v.OldBst2pos.GO$pathway.results[abs(YoungBst2pos.v.OldBst2pos.GO$pathway.results[,2])<0.01,]
write.table(YoungBst2pos.v.OldBst2pos.GO.sig,"YoungBst2pos.v.OldBst2pos.GO.sig.txt")

upper<-toupper(rownames(YoungBst2pos.v.OldBst2pos.conv.order))
rownames(YoungBst2pos.v.OldBst2pos.conv.order)<-upper
YoungBst2pos.v.OldBst2pos.GO<-KEGG.pathway.enrichment(test.results=YoungBst2pos.v.OldBst2pos.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c2.cp.kegg.v5.2.symbols.gmt")
YoungBst2pos.v.OldBst2pos.GO.sig<-YoungBst2pos.v.OldBst2pos.GO$pathway.results[abs(YoungBst2pos.v.OldBst2pos.GO$pathway.results[,2])<0.01,]
write.table(YoungBst2pos.v.OldBst2pos.GO.sig,"YoungBst2pos.v.OldBst2pos.KEGG.sig.txt")

upper<-toupper(rownames(YoungBst2pos.v.OldBst2pos.conv.order))
rownames(YoungBst2pos.v.OldBst2pos.conv.order)<-upper
YoungBst2pos.v.OldBst2pos.GO<-KEGG.pathway.enrichment(test.results=YoungBst2pos.v.OldBst2pos.conv.order,fold.change = "logFC",statistic="signed_lr",num.samples=reps,gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/h.all.v5.2.symbols.gmt")
YoungBst2pos.v.OldBst2pos.GO.sig<-YoungBst2pos.v.OldBst2pos.GO$pathway.results[abs(YoungBst2pos.v.OldBst2pos.GO$pathway.results[,2])<0.05,]
write.table(YoungBst2pos.v.OldBst2pos.GO.sig,"YoungBst2pos.v.OldBst2pos.Hallmark.sig.txt")

roundsize<-c()
for(i in 1:length(YoungBst2pos.v.OldBst2pos.GO.sig[,1])){
  roundsize<-c(roundsize,floor(min(log10(reps),-log10(YoungBst2pos.v.OldBst2pos.GO.sig$adj.p.val[i]))))
}
data<-data.frame(names=factor(rownames(YoungBst2pos.v.OldBst2pos.GO.sig),
                              levels=rownames(YoungBst2pos.v.OldBst2pos.GO.sig)[length(rownames(YoungBst2pos.v.OldBst2pos.GO.sig)):1],ordered=T),fc=YoungBst2pos.v.OldBst2pos.GO.sig$fold.change,size=roundsize)

p<-ggplot(data)
p<-p+geom_point(aes(x=data$names,y=data$fc,size=data$size))
p<-p+theme_classic()
p<-p+coord_flip()
p<-p+theme(legend.position = "none") 
p<-p+labs(y="Enrichment Score",title=NULL,x=NULL)


pdf("YoungBst2pos.v.OldBst2pos.Hallmark.sig.enrichment.pdf")
print(p)
dev.off()
#=========================================================================================================
#Young vs Old Bst2neg comparison
lrt.2vs1 <- glmLRT(fit,contrast=c(1,0,0,-1,0,0))

toptags_sort<-topTags(lrt.2vs1, p.value=0.05, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
write.table(toptags_sort,"OldBst2neg_vs_Young.txt")

tests<-decideTestsDGE(lrt.2vs1)

tests_neg<-lrt.2vs1$table[tests==-1,]
tests_neg_sort<-tests_neg[order(tests_neg[,4]),]
tests_pos<-lrt.2vs1$table[tests==1,]
tests_pos_sort<-tests_pos[order(tests_pos[,4]),]

write.table(tests_neg_sort,"OldBst2neg_vs_Young_Young_tests.txt")
write.table(tests_pos_sort,"OldBst2neg_vs_Young_OldBst2neg_tests.txt")
write.table(rownames(tests_neg_sort),"OldBst2neg_vs_Young_Young_names.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(tests_pos_sort),"OldBst2neg_vs_Young_OldBst2neg_names.txt",quote=F,row.names=F,col.names=F)

alltags_sort<-topTags(lrt.2vs1, p.value=1, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
OldBst2neg.v.Young.conv<-pvalue_converter(alltags_sort$table)
OldBst2neg.v.Young.conv.order<-OldBst2neg.v.Young.conv[order(OldBst2neg.v.Young.conv[,6],decreasing=T),]
write.table(OldBst2neg.v.Young.conv.order,"OldBst2neg.v.Young.order.txt")
rnk<-cbind(toupper(rownames(OldBst2neg.v.Young.conv.order)),OldBst2neg.v.Young.conv.order[,6])
write.table(rnk,"OldBst2neg.v.Young.order.exp2.rnk",quote=F,row.names=F,col.names=F,sep='\t')

upper<-toupper(rownames(OldBst2neg.v.Young.conv.order))
rownames(OldBst2neg.v.Young.conv.order)<-upper
OldBst2neg.v.Young.GO<-KEGG.pathway.enrichment(test.results=OldBst2neg.v.Young.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c5.all.v5.2.symbols.gmt")
OldBst2neg.v.Young.GO.sig<-OldBst2neg.v.Young.GO$pathway.results[abs(OldBst2neg.v.Young.GO$pathway.results[,2])<0.01,]
write.table(OldBst2neg.v.Young.GO.sig,"OldBst2neg.v.Young.GO.sig.txt")

upper<-toupper(rownames(OldBst2neg.v.Young.conv.order))
rownames(OldBst2neg.v.Young.conv.order)<-upper
OldBst2neg.v.Young.GO<-KEGG.pathway.enrichment(test.results=OldBst2neg.v.Young.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c2.cp.kegg.v5.2.symbols.gmt")
OldBst2neg.v.Young.GO.sig<-OldBst2neg.v.Young.GO$pathway.results[abs(OldBst2neg.v.Young.GO$pathway.results[,2])<0.01,]
write.table(OldBst2neg.v.Young.GO.sig,"OldBst2neg.v.Young.KEGG.sig.txt")

upper<-toupper(rownames(OldBst2neg.v.Young.conv.order))
rownames(OldBst2neg.v.Young.conv.order)<-upper
OldBst2neg.v.Young.GO<-KEGG.pathway.enrichment(test.results=OldBst2neg.v.Young.conv.order,fold.change = "logFC",statistic="signed_lr",num.samples=reps,gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/h.all.v5.2.symbols.gmt")
OldBst2neg.v.Young.GO.sig<-OldBst2neg.v.Young.GO$pathway.results[abs(OldBst2neg.v.Young.GO$pathway.results[,2])<0.05,]
write.table(OldBst2neg.v.Young.GO.sig,"OldBst2neg.v.Young.Hallmark.sig.txt")



roundsize<-c()
for(i in 1:length(OldBst2neg.v.Young.GO.sig[,1])){
  roundsize<-c(roundsize,floor(min(log10(reps),-log10(OldBst2neg.v.Young.GO.sig$adj.p.val[i]))))
}
data<-data.frame(names=factor(rownames(OldBst2neg.v.Young.GO.sig),
                              levels=rownames(OldBst2neg.v.Young.GO.sig)[length(rownames(OldBst2neg.v.Young.GO.sig)):1],ordered=T),fc=OldBst2neg.v.Young.GO.sig$fold.change,size=roundsize)

p<-ggplot(data)
p<-p+geom_point(aes(x=data$names,y=data$fc,size=data$size))
p<-p+theme_classic()
p<-p+coord_flip()
p<-p+theme(legend.position = "none") 
p<-p+labs(y="Enrichment Score",title=NULL,x=NULL)


pdf("OldBst2neg.v.Young.Hallmark.sig.enrichment.pdf")
print(p)
dev.off()

#=========================================================================================================
#qNSC comparison
lrt.2vs1 <- glmLRT(fit,contrast=c(0,0,1,0,0,-1))

toptags_sort<-topTags(lrt.2vs1, p.value=0.05, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
write.table(toptags_sort,"OldqNSC_vs_YoungqNSC.txt")

tests<-decideTestsDGE(lrt.2vs1)

tests_neg<-lrt.2vs1$table[tests==-1,]
tests_neg_sort<-tests_neg[order(tests_neg[,4]),]
tests_pos<-lrt.2vs1$table[tests==1,]
tests_pos_sort<-tests_pos[order(tests_pos[,4]),]

write.table(tests_neg_sort,"OldqNSC_vs_YoungqNSC_YoungqNSC_tests.txt")
write.table(tests_pos_sort,"OldqNSC_vs_YoungqNSC_OldqNSC_tests.txt")
write.table(rownames(tests_neg_sort),"OldqNSC_vs_YoungqNSC_YoungqNSC_names.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(tests_pos_sort),"OldqNSC_vs_YoungqNSC_OldqNSC_names.txt",quote=F,row.names=F,col.names=F)

alltags_sort<-topTags(lrt.2vs1, p.value=1, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
OldqNSC.v.YoungqNSC.conv<-pvalue_converter(alltags_sort$table)
OldqNSC.v.YoungqNSC.conv.order<-OldqNSC.v.YoungqNSC.conv[order(OldqNSC.v.YoungqNSC.conv[,6],decreasing=T),]
write.table(OldqNSC.v.YoungqNSC.conv.order,"OldqNSC.v.YoungqNSC.order.txt")
rnk<-cbind(toupper(rownames(OldqNSC.v.YoungqNSC.conv.order)),OldqNSC.v.YoungqNSC.conv.order[,6])
write.table(rnk,"OldqNSC.v.YoungqNSC.order.rnk",quote=F,row.names=F,col.names=F,sep='\t')

upper<-toupper(rownames(OldqNSC.v.YoungqNSC.conv.order))
rownames(OldqNSC.v.YoungqNSC.conv.order)<-upper
OldqNSC.v.YoungqNSC.GO<-KEGG.pathway.enrichment(test.results=OldqNSC.v.YoungqNSC.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c5.all.v5.2.symbols.gmt")
OldqNSC.v.YoungqNSC.GO.sig<-OldqNSC.v.YoungqNSC.GO$pathway.results[abs(OldqNSC.v.YoungqNSC.GO$pathway.results[,2])<0.01,]
write.table(OldqNSC.v.YoungqNSC.GO.sig,"OldqNSC.v.YoungqNSC.GO.sig.txt")

upper<-toupper(rownames(OldqNSC.v.YoungqNSC.conv.order))
rownames(OldqNSC.v.YoungqNSC.conv.order)<-upper
OldqNSC.v.YoungqNSC.GO<-KEGG.pathway.enrichment(test.results=OldqNSC.v.YoungqNSC.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c2.cp.kegg.v5.2.symbols.gmt")
OldqNSC.v.YoungqNSC.GO.sig<-OldqNSC.v.YoungqNSC.GO$pathway.results[abs(OldqNSC.v.YoungqNSC.GO$pathway.results[,2])<0.01,]
write.table(OldqNSC.v.YoungqNSC.GO.sig,"OldqNSC.v.YoungqNSC.KEGG.sig.txt")

upper<-toupper(rownames(OldqNSC.v.YoungqNSC.conv.order))
rownames(OldqNSC.v.YoungqNSC.conv.order)<-upper
OldqNSC.v.YoungqNSC.GO<-KEGG.pathway.enrichment(test.results=OldqNSC.v.YoungqNSC.conv.order,fold.change = "logFC",statistic="signed_lr",num.samples=reps,gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/h.all.v5.2.symbols.gmt")
OldqNSC.v.YoungqNSC.GO.sig<-OldqNSC.v.YoungqNSC.GO$pathway.results[abs(OldqNSC.v.YoungqNSC.GO$pathway.results[,2])<0.05,]
write.table(OldqNSC.v.YoungqNSC.GO.sig,"OldqNSC.v.YoungqNSC.Hallmark.sig.txt")



roundsize<-c()
for(i in 1:length(OldqNSC.v.YoungqNSC.GO.sig[,1])){
  roundsize<-c(roundsize,floor(min(log10(reps),-log10(OldqNSC.v.YoungqNSC.GO.sig$adj.p.val[i]))))
}
data<-data.frame(names=factor(rownames(OldqNSC.v.YoungqNSC.GO.sig),
                              levels=rownames(OldqNSC.v.YoungqNSC.GO.sig)[length(rownames(OldqNSC.v.YoungqNSC.GO.sig)):1],ordered=T),fc=OldqNSC.v.YoungqNSC.GO.sig$fold.change,size=roundsize)

p<-ggplot(data)
p<-p+geom_point(aes(x=data$names,y=data$fc,size=data$size))
p<-p+theme_classic()
p<-p+coord_flip()
p<-p+theme(legend.position = "none") 
p<-p+labs(y="Enrichment Score",title=NULL,x=NULL)


pdf("OldqNSC.v.YoungqNSC.Hallmark.sig.enrichment.pdf")
print(p)
dev.off()





#=========================================================================================================

#Young aNSC vs Young qNSC comparison
lrt.2vs1 <- glmLRT(fit,contrast=c(0,0,0,1,0,-1))

toptags_sort<-topTags(lrt.2vs1, p.value=0.05, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
write.table(toptags_sort,"YoungaNSC_vs_YoungqNSC.txt")

tests<-decideTestsDGE(lrt.2vs1)

tests_neg<-lrt.2vs1$table[tests==-1,]
tests_neg_sort<-tests_neg[order(tests_neg[,4]),]
tests_pos<-lrt.2vs1$table[tests==1,]
tests_pos_sort<-tests_pos[order(tests_pos[,4]),]

write.table(tests_neg_sort,"YoungaNSC_vs_YoungqNSC_YoungqNSC_tests.txt")
write.table(tests_pos_sort,"YoungaNSC_vs_YoungqNSC_YoungaNSC_tests.txt")
write.table(rownames(tests_neg_sort),"YoungaNSC_vs_YoungqNSC_YoungqNSC_names.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(tests_pos_sort),"YoungaNSC_vs_YoungqNSC_YoungaNSC_names.txt",quote=F,row.names=F,col.names=F)

alltags_sort<-topTags(lrt.2vs1, p.value=1, sort.by="p.value",n=length(lrt.2vs1$table[,1]))
YoungaNSC.v.YoungqNSC.conv<-pvalue_converter(alltags_sort$table)
YoungaNSC.v.YoungqNSC.conv.order<-YoungaNSC.v.YoungqNSC.conv[order(YoungaNSC.v.YoungqNSC.conv[,6],decreasing=T),]
write.table(YoungaNSC.v.YoungqNSC.conv.order,"YoungaNSC.v.YoungqNSC.order.txt")
rnk<-cbind(toupper(rownames(YoungaNSC.v.YoungqNSC.conv.order)),YoungaNSC.v.YoungqNSC.conv.order[,6])
write.table(rnk,"YoungaNSC.v.YoungqNSC.order.rnk",quote=F,row.names=F,col.names=F,sep='\t')

upper<-toupper(rownames(YoungaNSC.v.YoungqNSC.conv.order))
rownames(YoungaNSC.v.YoungqNSC.conv.order)<-upper
YoungaNSC.v.YoungqNSC.GO<-KEGG.pathway.enrichment(test.results=YoungaNSC.v.YoungqNSC.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c5.all.v5.2.symbols.gmt")
YoungaNSC.v.YoungqNSC.GO.sig<-YoungaNSC.v.YoungqNSC.GO$pathway.results[abs(YoungaNSC.v.YoungqNSC.GO$pathway.results[,2])<0.01,]
write.table(YoungaNSC.v.YoungqNSC.GO.sig,"YoungaNSC.v.YoungqNSC.GO.sig.txt")

upper<-toupper(rownames(YoungaNSC.v.YoungqNSC.conv.order))
rownames(YoungaNSC.v.YoungqNSC.conv.order)<-upper
YoungaNSC.v.YoungqNSC.GO<-KEGG.pathway.enrichment(test.results=YoungaNSC.v.YoungqNSC.conv.order,fold.change = "logFC",statistic="signed_lr",gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/c2.cp.kegg.v5.2.symbols.gmt")
YoungaNSC.v.YoungqNSC.GO.sig<-YoungaNSC.v.YoungqNSC.GO$pathway.results[abs(YoungaNSC.v.YoungqNSC.GO$pathway.results[,2])<0.01,]
write.table(YoungaNSC.v.YoungqNSC.GO.sig,"YoungaNSC.v.YoungqNSC.KEGG.sig.txt")

upper<-toupper(rownames(YoungaNSC.v.YoungqNSC.conv.order))
rownames(YoungaNSC.v.YoungqNSC.conv.order)<-upper
YoungaNSC.v.YoungqNSC.GO<-KEGG.pathway.enrichment(test.results=YoungaNSC.v.YoungqNSC.conv.order,fold.change = "logFC",statistic="signed_lr",num.samples=reps,gene.sets.file="/Volumes/guacamole/Software/R_Files_Packages_Functions/MsigDB_pathways/h.all.v5.2.symbols.gmt")
YoungaNSC.v.YoungqNSC.GO.sig<-YoungaNSC.v.YoungqNSC.GO$pathway.results[abs(YoungaNSC.v.YoungqNSC.GO$pathway.results[,2])<0.05,]
write.table(YoungaNSC.v.YoungqNSC.GO.sig,"YoungaNSC.v.YoungqNSC.Hallmark.sig.txt")



roundsize<-c()
for(i in 1:length(YoungaNSC.v.YoungqNSC.GO.sig[,1])){
  roundsize<-c(roundsize,floor(min(log10(reps),-log10(YoungaNSC.v.YoungqNSC.GO.sig$adj.p.val[i]))))
}
data<-data.frame(names=factor(rownames(YoungaNSC.v.YoungqNSC.GO.sig),
                              levels=rownames(YoungaNSC.v.YoungqNSC.GO.sig)[length(rownames(YoungaNSC.v.YoungqNSC.GO.sig)):1],ordered=T),fc=YoungaNSC.v.YoungqNSC.GO.sig$fold.change,size=roundsize)

p<-ggplot(data)
p<-p+geom_point(aes(x=data$names,y=data$fc,size=data$size))
p<-p+theme_classic()
p<-p+coord_flip()
p<-p+theme(legend.position = "none") 
p<-p+labs(y="Enrichment Score",title=NULL,x=NULL)


pdf("YoungaNSC.v.YoungqNSC.Hallmark.sig.enrichment.pdf")
print(p)
dev.off()










#===========================================================================================================
#Plotting


#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/LaCie/10242017_Bst2_full_RNA_seq/Analysis/")
allpops<-read.table("Bst2_Full_allcounts_new_names.txt")

allpops_genecounts<-colSums(allpops!=0)
#Filtering for expressed genes 

allpops_genefilt<-allpops[!grepl("ERCC-",rownames(allpops)),]
allpops_genefilt_sampfilt<-allpops_genefilt[,grepl("Exp2",colnames(allpops_genefilt))]

glm <- DGEList(counts=allpops_genefilt_sampfilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/mm_genomes/mm10/mm10.genes.current.gene.sizes.R")
temp<-(as.numeric(exonic.gene.sizes))
names(temp)<-tolower(names(exonic.gene.sizes))
exonic.gene.sizes<-temp
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt <- rpkm(glm.norm, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)


comb_fpkm<-fpkm_glm_genefilt

comb_fac<-c()
for(i in 1:length(colnames(comb_fpkm))){
  temp<-strsplit(colnames(comb_fpkm)[i],split='[_]')[[1]][1]
  comb_fac<-c(comb_fac,temp)
}
comb_fac<-factor(comb_fac,levels=c(unique(comb_fac)),ordered=T)

int_genes<-c("Sox2","Ascl1","Bst2","Ifitm3","Usp18","Mki67","Pecam1","Ptprc","Dcx","Dlx2","Jag1","Ccl5","Ifi27","Zfp9",
             "Gsx1","Crmp1","Nfix","Elmo1","Egfr","Cxcl14","Cd9","Thbs4","Prom1","Ntsr2","Atp1a2","Cdk1","Top2a","Prc1",
             "Kif11","Cdk4","Rpl4","Ccnd2","Irf7","Cdkn2a","Htra1","Oas1a","Mog","Mbp","Mag","Vwf")


# int_genes<-c("Cd34","Icam1","Tek")
# 
int_genes<-c("Ascl1","Egfr","Nes","Mki67","Cdkn1b","Myc","Aimp2","Adar","Clu")

int_genes<-c("Gfap","B2m","Plcd4","Irf7","Irf9")
int_genes<-c("Mki67")
int_genes<-c("Dcx","Nrxn3","Dlx2","Dlx1","Dlx6as1")
# 
# int_genes<-c("Trex1")

#int_genes<-rownames(comb_fpkm)[grepl("Ifn",rownames(comb_fpkm))]
#int_genes<-c("Ifit1","Ifit3b","Ifit3","Stat1")

int_genes<-c("Ccl5","Il2","Il18","Il12a","Il2ra","Il12rb1","Il12rb2","Il18r1","Egfr")
int_genes<-c("Dcx","Nrxn3","Dlx2","Dlx1","Dlx6as1")
setwd("/Volumes/LaCie/10242017_Bst2_full_RNA_seq/Analysis/Violins")

match(int_genes,rownames(comb_fpkm))

for(i in int_genes){
  curr<-comb_fpkm[match(i,rownames(comb_fpkm)),]
  
  data<-data.frame(values=curr,fac=comb_fac)
  p<-ggplot(data)
  p<-p+geom_point(aes(x=data$fac,y=data$values,fill=data$fac))
  p<-p+geom_boxplot(aes(x=data$fac,y=data$values),alpha=0.5)
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
  p<-p+labs(x=NULL,y="Expression (FPKM)",title=i)
  p<-p+theme(legend.position = "none") 
  p<-p+theme(axis.text.x=element_text(size=18))
  p<-p+theme(axis.text.y=element_text(size=22))
  p<-p+theme(axis.title.y=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=24))
  pdf(paste(i,"_Bst2Pilot.pdf",sep=""),height=4,width=7)
  print(p)
  dev.off()
  
}
