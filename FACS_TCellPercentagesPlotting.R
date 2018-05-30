rm(list=ls())
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("FACS_TCellFreqs.txt",sep="\t")
values<-values[grepl("svz",rownames(values)),]

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")
fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="_00")[[1]][1])
}

fac<-factor(fac,levels=c("young_svz","old_svz"),ordered=T)

for(i in 1:length(values[1,])){
  data<-data.frame(values=values[,i],fac=fac)
  p<-ggplot(data)
  #p<-p+geom_boxplot(aes(x=data$fac,y=data$values,fill=data$fac))+geom_point(aes(x=data$fac,y=data$values),alpha=0.7)
  p<-p+geom_jitter(aes(x=data$fac,y=data$values,color=data$fac), width=0.1, alpha=0.7,size=4)+ 
    stat_summary(aes(x=data$fac,y=data$values),fun.y=mean, geom = "point",size=9,alpha=0.7,shape=95) + 
    stat_summary(aes(x=data$fac,y=data$values),fun.data = mean_se, geom = "errorbar",width=0.2,size=0.4,alpha=0.7)
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
  p<-p+scale_color_manual(values=c("#40BBEC","#EF4136"))
  pdf(paste0(colnames(values)[i],"_svz_only.pdf"),width=3,height=3)
  print(p)
  dev.off()
}

for(i in 6){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$young_svz,split(values[,i],fac)$old_svz))
}

#Experiment 1 - TCR clonality by PCR
exp1<-values[c(1:2,6,7),]
fac<-c()

for(i in 1:length(exp1[,1])){
  fac<-c(fac,strsplit(rownames(exp1)[i],split="_00")[[1]][1])
}

fac<-factor(fac,levels=c("young_svz","old_svz"),ordered=T)

for(i in 6){
  print(colnames(exp1)[i])
  print(by(exp1[,i],fac,mean))
  print(by(exp1[,i],fac,std.error))
  print(wilcox.test(split(exp1[,i],fac)$young_svz,split(exp1[,i],fac)$old_svz))
}


#Experiment 2 - TCR clonality by PCR
exp2<-values[c(3,4,5,8,9),]
fac<-c()

for(i in 1:length(exp2[,1])){
  fac<-c(fac,strsplit(rownames(exp2)[i],split="_00")[[1]][1])
}

fac<-factor(fac,levels=c("young_svz","old_svz"),ordered=T)

for(i in 6){
  print(colnames(exp2)[i])
  print(by(exp2[,i],fac,mean))
  print(by(exp2[,i],fac,std.error))
  print(wilcox.test(split(exp2[,i],fac)$young_svz,split(exp2[,i],fac)$old_svz))
}

