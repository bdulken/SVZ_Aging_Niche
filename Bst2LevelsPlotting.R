#Bst2 levels
rm(list=ls())
library(ggplot2)
setwd("/Users/bendulken/Documents/Experiments/Interferon/11082017_Bst2_scRNA_seq_NeurosphereFormation/")
values<-read.table("Bst2Levels.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(rownames(values))){
  temp<-strsplit(rownames(values)[i],split='[_]')[[1]][1]
  fac<-c(fac,temp)
}

fac<-factor(fac,levels=c("young","old"),labels=c("Young","Old"), ordered=T)

for(i in 1:length(values[1,])){
  data<-data.frame(values=values[,i],fac=fac)
  p<-ggplot(data)
  p<-p+geom_boxplot(aes(x=data$fac,y=data$values))+geom_point(aes(x=data$fac,y=data$values),alpha=0.7)
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p+labs(x=NULL,y=colnames(values)[i],title=colnames(values)[i])
  p<-p+theme(legend.position = "none") 
  p<-p+theme(axis.text.x=element_text(size=18))
  p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
  p<-p+theme(axis.text.y=element_text(size=16))
  p<-p+theme(axis.title.y=element_text(size=15))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  pdf(paste0(colnames(values)[i],".pdf"),width=8,height=5)
  print(p)
  dev.off()
}

