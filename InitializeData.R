library(data.table);library(pheatmap)
library(ggplot2);library(hrbrthemes)
library(dplyr);library(RColorBrewer)

library(BiocManager)
library('clusterProfiler')
library(stringr)
library('org.Hs.eg.db')

setwd('E:/ManuscriptI/script')
rm(list = ls())

Protein<-read.csv('./ProteinAbundances.csv',sep = ',')
colnames(Protein)
head(Protein[,1])
rownames(Protein)<-Protein[,1];Protein<-Protein[,-1]
colnames(Protein)
colnames(Protein)<-gsub('Abundance..','',colnames(Protein))
colnames(Protein)<-gsub('..Sample','',colnames(Protein))
colnames(Protein)
Protein<-Protein[-grep('?CON',rownames(Protein)),]
Keratin<-Protein[grep('?Keratin',Protein$Description),]
Protein<-Protein[-grep('?Keratin',Protein$Description),]
Protein<-Protein[-grep('?Trypsin',Protein$Description),]
Description<-data.frame(row.names = rownames(Protein),
                        Description=Protein[,'Description'])
Protein<-dplyr::select(Protein,select=-'Description')

Keratin[is.na(Keratin)]<-0
Keratin<-apply(Keratin[,-c(1,2)],2,sum)
plot(Keratin)
Protein<-dplyr::select(Protein,select=-'F8')

Protein[is.na(Protein)]<-0
###Protein<-Protein[-grep('Q9NYL5|P61247',rownames(Protein)),]
hist(apply(Protein,2,function(x){sum(x>0)}))
Protein<-Protein[,apply(Protein,2,function(x){sum(x>0)>1000})]
colnames(Protein)
hist(apply(Protein,1,function(x){sum(x>0)}))
Original<-Protein[apply(Protein,1,function(x){sum(x>0)>35}),]
Original<-Original[-grep(c('Q9NYL5|P00761'),rownames(Original)),]

head(Original)
Sum<-apply(Original,2,sum)
library(MASS)
fit<-fitdistr(Sum,'normal')
fit
hist(Sum)
min<-fit$estimate[1]-3*fit$estimate[2]
max<-fit$estimate[1]+3*fit$estimate[2]
Sum[Sum<min];Sum[Sum>max]
Original<-dplyr::select(Original,select=-c('F21'))
Log<-log2(Original+1)

### Use All Protein to Normalization
Mean<-apply(Original,2,mean)
Norm<-sweep(Original,2,Mean/mean(Mean),'/')
LogNorm<-log2(Norm+1)
save(Description,
     Original,Log,Norm,LogNorm,
     file='ProteinAbundance.Rdata')

load('ProteinAbundance.Rdata')
write.csv(Norm,file='forCorrelation.csv')

CorMap<-read.csv('ReorderCorrelation.csv')
head(CorMap)
rownames(CorMap)<-CorMap[,1];CorMap<-CorMap[,-1]
cluster<-read.table('CorrelationGroup.txt',sep=' ')
group<-rep('0',dim(CorMap)[1])
for(i in 1:dim(cluster)[1]){
  group[cluster[i,1]:cluster[i,2]]<-rownames(cluster)[i]
}
clus<-group
group<-data.frame(row.names = rownames(CorMap),
                  group=group)
pheatmap(CorMap,
         color = colorRampPalette(c("blue","white","red"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 0.3, cellheight = 0.3,
         show_rownames = F, show_colnames = F,
         annotation_col=group, annotation_row=group)


Cluster<-data.frame(row.names = rownames(group)[group$group>0])
Cluster$group<-group[rownames(Cluster),]
Cluster$Description<-Description[rownames(Cluster),]
for (i in 1:max(as.numeric(group$group))){
  print(i)
  print(Description[rownames(group)[group==i],])
}
sum(group[group!=0])

save(CorMap,group,file = 'CorrelationMap.Rdata')

data<-CorMap
data<-as.data.frame(data)
data$Protein2<-rownames(data)
data<-reshape2::melt(data)
# head(data)rm()
colnames(data)[2]<-'Protein1'
data<-data[data$Protein2!=data$Protein1,]
data$index<-apply(data,1,function(x) 
  paste(sort(x[c(1,2)],decreasing=T)[1],
        sort(x[c(1,2)],decreasing=T)[2],sep='_'))
data<-data[duplicated(data$index)==F,]
data<-data[,-c(4,5)]
Pair<-data
rownames(Pair)<-Pair$index
Pair$index<-paste(Pair$Protein2,Pair$Protein1,sep='_')
rownames(Pair)<-Pair$index
save(Pair,file='CorrelationPair.Rdata')