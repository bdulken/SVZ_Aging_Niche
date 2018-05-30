#Plotting and stats for T cell quantification and Ki67+Sox2+ positive cells in SVZ.
#Sections are from 4 old mice and 5 young mice.

#For young mice used 3-7 month, for old mice used 20-24 month mice.

#Counts from stains performed on 7/12 and 7/17/17


rm(list=ls())
library(ggplot2)
library(plotrix)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("Immunofluorescence_TCell_counts.txt")
values<-values[grepl("3m",rownames(values))|grepl("7m",rownames(values))|grepl("20m",rownames(values))|grepl("24m",rownames(values)),]

fac<-factor(c(rep("Young",5),rep("Old",4)),
            levels=c("Young","Old"),
            ordered=T)

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")
#Ki67
data<-data.frame(values=values[,1],fac=fac)
p<-ggplot(data)
#p<-p+geom_boxplot(aes(x=data$fac,y=data$values,fill=data$fac))+geom_point(aes(x=data$fac,y=data$values),alpha=0.7)
p<-p+geom_jitter(aes(x=data$fac,y=data$values,color=data$fac), width=0.1, alpha=0.7,size=4)+ 
  stat_summary(aes(x=data$fac,y=data$values),fun.y=mean, geom = "point",size=9,alpha=0.7,shape=95) + 
  stat_summary(aes(x=data$fac,y=data$values),fun.data = mean_se, geom = "errorbar",width=0.2,size=0.4,alpha=0.7)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+labs(x=NULL,y="Sox2+Ki67+ Per SVZ",title="Sox2+Ki67+ Per SVZ")
p<-p+theme(legend.position = "none") 
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.text.y=element_text(size=16))
p<-p+theme(axis.title.y=element_text(size=15))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+scale_color_manual(values=c("#40BBEC","#EF4136"))
pdf("Ki67+_per_svz_extremes.pdf",width=3,height=3)
print(p)
dev.off()

#SVZ t cells
data<-data.frame(values=values[,2],fac=fac)
p<-ggplot(data)
#p<-p+geom_boxplot(aes(x=data$fac,y=data$values,fill=data$fac))+geom_point(aes(x=data$fac,y=data$values),alpha=0.7)
p<-p+geom_jitter(aes(x=data$fac,y=data$values,color=data$fac), width=0.1, alpha=0.7,size=4)+ 
  stat_summary(aes(x=data$fac,y=data$values),fun.y=mean, geom = "point",size=9,alpha=0.7,shape=95) + 
  stat_summary(aes(x=data$fac,y=data$values),fun.data = mean_se, geom = "errorbar",width=0.2,size=0.4,alpha=0.7)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+labs(x=NULL,y="T Cells Per SVZ",title="T Cells Per SVZ")
p<-p+theme(legend.position = "none") 
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.text.y=element_text(size=16))
p<-p+theme(axis.title.y=element_text(size=15))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+scale_color_manual(values=c("#40BBEC","#EF4136"))
pdf("TCells_per_svz_extremes.pdf",width=3,height=3)
print(p)
dev.off()


#total T cells
data<-data.frame(values=values[,3],fac=fac)
p<-ggplot(data)
#p<-p+geom_boxplot(aes(x=data$fac,y=data$values,fill=data$fac))+geom_point(aes(x=data$fac,y=data$values),alpha=0.7)
p<-p+geom_jitter(aes(x=data$fac,y=data$values,color=data$fac), width=0.1, alpha=0.7,size=4)+ 
  stat_summary(aes(x=data$fac,y=data$values),fun.y=mean, geom = "point",size=9,alpha=0.7,shape=95) + 
  stat_summary(aes(x=data$fac,y=data$values),fun.data = mean_se, geom = "errorbar",width=0.2,size=0.4,alpha=0.7)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+labs(x=NULL,y="T Cells Per Section",title="T Cells Per Section")
p<-p+theme(legend.position = "none") 
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.text.y=element_text(size=16))
p<-p+theme(axis.title.y=element_text(size=15))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+scale_color_manual(values=c("#40BBEC","#EF4136"))
pdf("TCells_per_section_extremes.pdf",width=3,height=3)
print(p)
dev.off()

for(i in 1:3){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



#Counting T cells from experiment in which I also performed CD31 staining (Supp Fig 1c)

rm(list=ls())
library(ggplot2)
setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/InputFiles/")
values<-read.table("TCell_Quantification_FromCD31stains.txt")

setwd("/Users/bendulken/Dropbox/Aging Single Cell Paper/Final_Code/OutputFiles/")
fac<-c()
for(i in 1:length(values[,1])){
  temp<-strsplit(rownames(values)[i],split='[.]')[[1]][1]
  fac<-c(fac,temp)
}

fac<-factor(fac,
            levels=c("Young","Old"),
            ordered=T)

data<-data.frame(values=values[,1],fac=fac)
p<-ggplot(data)
p<-p+geom_boxplot(aes(x=data$fac,y=data$values),alpha=0.7,fill="red")
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+labs(x=NULL,y="# CD3+ Per Section")
p<-p+theme(legend.position = "none") 
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.text.y=element_text(size=16))
p<-p+theme(axis.title.y=element_text(size=15))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
pdf("CD3+_per_section.pdf",width=4,height=4)
print(p)
dev.off()


for(i in 1){
  print(colnames(values)[i])
  print(by(values[,i],fac,mean))
  print(by(values[,i],fac,std.error))
  print(wilcox.test(split(values[,i],fac)$Young,split(values[,i],fac)$Old))
}



