#Ki67 levels
rm(list=ls())
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("IC_Stat1_AllSamps.txt",header=T,sep="\t")

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("YoungNSC","OldNSC"),ordered=T)


data<-data.frame(values,fac=fac)

p<-ggplot(data)
#p<-p+geom_boxplot(aes(x=data$fac, y=data[,2],fill=data$fac))+geom_point(aes(x=data$fac,y=data[,2]),alpha=0.9)
p<-p+geom_jitter(aes(x=data$fac,y=data[,1],color=data$fac), width=0.1, alpha=0.7,size=4)+ 
  stat_summary(aes(x=data$fac,y=data[,1]),fun.y=mean, geom = "point",size=9,alpha=0.7,shape=95) + 
  stat_summary(aes(x=data$fac,y=data[,1]),fun.data = mean_se, geom = "errorbar",width=0.2,size=0.4,alpha=0.7)
p<-p+labs(x="Age",y="Percent Stat1 High")
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
pdf("Stat1Values.pdf",height=5,width=2.5)
print(p)
dev.off()


for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}


#Ibu 11162017
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("IC_Stat1values_11162017ibu.txt",header=T,sep="\t")

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("YoungNSC","OldNSC"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



#Ibu 03082017
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("IC_Stat1High_Values_forPaper_03082017ibu.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("YoungNSC","OldNSC"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



#06152016
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("IC_Stat1_HighFreqs_06152016.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("YoungNSC","OldNSC"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



#07152016
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("IC_Stat1_HighFreqs_07152016.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("YoungNSC","OldNSC"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}


#09132016
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("IC_Stat1values_09132016_agegrad.txt",header=T,sep="\t")

fac<-c()
for(i in 1:length(values[,1])){
  fac<-c(fac,strsplit(rownames(values)[i],split="[_0]")[[1]][1])
}

fac<-factor(fac,levels=c("YoungNSC","OldNSC"),ordered=T)

for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}


