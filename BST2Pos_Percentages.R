#Bst2 Percentages
rm(list=ls())
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Bst2Levels_AllExps.txt",header=T,sep="\t")

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)


data<-data.frame(values,fac=fac)

p<-ggplot(data)
#p<-p+geom_boxplot(aes(x=data$fac, y=data[,2],fill=data$fac))+geom_point(aes(x=data$fac,y=data[,2]),alpha=0.9)
p<-p+geom_jitter(aes(x=data$fac,y=data[,1],color=data$fac), width=0.1, alpha=0.7,size=4)+ 
  stat_summary(aes(x=data$fac,y=data[,1]),fun.y=mean, geom = "point",size=9,alpha=0.7,shape=95) + 
  stat_summary(aes(x=data$fac,y=data[,1]),fun.data = mean_se, geom = "errorbar",width=0.2,size=0.4,alpha=0.7)
p<-p+labs(x="Age",y="Percent BST2 Positive")
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
#p<=0.008
pdf("Bst2Percentages_all.pdf",height=5,width=2.5)
print(p)
dev.off()


for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}


#RNA-seq 09252017
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Bst2Levels_09252017_RNAseq.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



#Bst2 SC 11082017
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Bst2Levels_11082017_sc.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



#Bst2 Edu1
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Bst2Percent_12172017_Edu.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



#Bst2 Edu2
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Bst2PosPercent_02072017_Edu.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}


#Bst2 Edu3
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Bst2Levels_03212018_Edu.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("Young","Old"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}


