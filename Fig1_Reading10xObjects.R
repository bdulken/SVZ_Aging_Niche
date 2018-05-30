
rm(list=ls())
library(Seurat)
library(ggplot2)
library(cellrangerRkit) 

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles")

cellranger_pipestance_path <- "/Volumes/LaCie/03142017_Dhananjay10x_perfusion/Y2-HFYF5BGX2/"
young1_gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)

young1_expression<-exprs(young1_gbm)
young1featureData<-fData(young1_gbm)

cellranger_pipestance_path <- "/Volumes/LaCie/03142017_Dhananjay10x_perfusion/O2-HFYF5BGX2/"
old1_gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)

old1_expression<-exprs(old1_gbm)
old1featureData<-fData(old1_gbm)

cellranger_pipestance_path <- "/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/Y3-H5GLYBGX3/"
young2_gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)

young2_expression<-exprs(young2_gbm)
young2featureData<-fData(young2_gbm)

cellranger_pipestance_path <- "/Volumes/LaCie/05222017_10xGenomics_plusperfusion_v2/O3-H5GLYBGX3/"
old2_gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)

old2_expression<-exprs(old2_gbm)
old2featureData<-fData(old2_gbm)

for(i in 1:length(colnames(young1_expression))){
  colnames(young1_expression)[i]<-paste("Young1_",strsplit(colnames(young1_expression)[i],split='[-]')[[1]][1],sep="")
}


for(i in 1:length(colnames(old1_expression))){
  colnames(old1_expression)[i]<-paste("Old1_",strsplit(colnames(old1_expression)[i],split='[-]')[[1]][1],sep="")
}


for(i in 1:length(colnames(young2_expression))){
  colnames(young2_expression)[i]<-paste("Young2_",strsplit(colnames(young2_expression)[i],split='[-]')[[1]][1],sep="")
}


for(i in 1:length(colnames(old2_expression))){
  colnames(old2_expression)[i]<-paste("Old2_",strsplit(colnames(old2_expression)[i],split='[-]')[[1]][1],sep="")
}

all_exprs<-cbind(young1_expression,young2_expression,old1_expression,old2_expression)

#Aggregate counts by genes
all_exprs_matrix<-as.matrix(all_exprs)
all_exprs_matrix_agg<-aggregate(all_exprs_matrix, by=list(young1featureData[,2]),FUN=sum)

all_exprs_matrix_agg_mod<-all_exprs_matrix_agg[,2:length(all_exprs_matrix_agg[1,])]
rownames(all_exprs_matrix_agg_mod)<-all_exprs_matrix_agg[,1]

#saves counts matrix
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/10xInputFiles")
save(all_exprs_matrix_agg_mod,file="all_exprs_matrix_agg_mod_tworeps.R")

