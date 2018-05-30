#Bst2 Pilot Compiling Counts

rm(list=ls())
library(ggplot2)
library(edgeR)


# All specific pops ---------------------------------------------------------------------

#Reading in all HTseq counts files.
setwd("/Volumes/LaCie/10242017_Bst2_full_RNA_seq/171006_NS500615_0580_AHNLFKBGX3/STAR/HTseqcounts_OnlyOld/")
fileList=list.files(pattern=".counts.txt$")

#Combine the counts from each of the individual files to make one large counts matrix and save the counts matrix to
#a file.
allcounts<-read.table(fileList[1],header=F)
colnames(allcounts)[2]<-strsplit(fileList[1],split="_S")[[1]][1]

for (i in 2:length(fileList)){
  temp = read.table(fileList[i],header=F)
  allcounts=cbind(allcounts,temp[,2])
  colnames(allcounts)[i+1]<-strsplit(fileList[i],split="_S")[[1]][1]
}

allcounts<-allcounts[-c((length(allcounts[,1])-4):length(allcounts[,1])),] #removes HTseq QC values from bottom of counts matrix
allcounts_fin<-allcounts[,2:length(allcounts[1,])]
rownames(allcounts_fin)<-allcounts[,1]

write.table(allcounts_fin,"/Volumes/LaCie/10242017_Bst2_full_RNA_seq/Analysis/Bst2_Full_allcounts_OnlyOld.txt")