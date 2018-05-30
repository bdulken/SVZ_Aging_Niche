#Plotting average Stat1 Fluorescence
#1/9/18
rm(list=ls())
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Stat1_IF_AverageFluorescence.txt",header=T,sep="\t")

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

data<-data.frame(values,fac=fac)
p<-ggplot(data)
#p<-p+geom_boxplot(aes(x=data$fac, y=data[,1],fill=data$fac))+geom_point(aes(x=data$fac, y=data[,1]),alpha=0.9)
p<-p+geom_jitter(aes(x=data$fac,y=data[,1],color=data$fac), width=0.1, alpha=0.7,size=4)+ 
  stat_summary(aes(x=data$fac,y=data[,1]),fun.y=mean, geom = "point",size=9,alpha=0.7,shape=95) + 
  stat_summary(aes(x=data$fac,y=data[,1]),fun.data = mean_se, geom = "errorbar",width=0.2,size=0.4,alpha=0.7)
p<-p+labs(x="Age",y="Stat1 Fluorescence Intensity")
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+theme(legend.position="none")
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.title.x=element_text(size=20))
p<-p+theme(axis.text.y=element_text(size=18))
p<-p+theme(axis.title.y=element_text(size=20))
p<-p+theme(plot.title=element_text(size=18))
p<-p+scale_color_manual(values=c("#40BBEC","#EF4136"))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
#pvalue=0.02

pdf("IF_Stat1_Intensity.pdf",height=5,width=3)
print(p)
dev.off()

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}




#Experiment 1 - Dena/Ashley Mice
exp1<-values[c(1:2,5,6),]
fac<-c()

for(i in c(1,2,5,6)){
  fac<-c(fac,strsplit(rownames(values)[i],split="_0")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

for(i in 1){
  print(by(exp1,fac,mean))
  print(by(exp1,fac,std.error))
  print(wilcox.test(split(exp1,fac)$Young,split(exp1,fac)$Old))
}


#Experiment 1 - Dena/Ashley Mice
exp2<-values[c(3,4,7,8),]
fac<-c()

for(i in c(3,4,7,8)){
  fac<-c(fac,strsplit(rownames(values)[i],split="_0")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

for(i in 1){
  print(by(exp2,fac,mean))
  print(by(exp2,fac,std.error))
  print(wilcox.test(split(exp2,fac)$Young,split(exp2,fac)$Old))
}

