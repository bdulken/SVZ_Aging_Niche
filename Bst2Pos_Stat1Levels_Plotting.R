#Bst2 Stat1 Levels All Samps
rm(list=ls())
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Bst2Stat1_AllSamps.txt",header=T,sep="\t")


setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("OldBst2Neg","OldBst2Pos"),ordered=T)


data<-data.frame(values=values[,1],fac=fac)
p<-ggplot(data)
#p<-p+geom_boxplot(aes(x=data$fac,y=data$values))+geom_point(aes(x=data$fac,y=data$values),alpha=0.7)
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
p<-p+scale_color_manual(values=c("#EF4136","#9B2828"))
pdf("Stat1allbst2exps.pdf",width=3.5,height=5)
print(p)
dev.off()



for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$OldBst2Neg,split(values[,i],fac)$OldBst2Pos))
}
