library(data.table);library(pheatmap)
library(ggplot2);library(hrbrthemes)
library(dplyr);library(RColorBrewer)

library(BiocManager)
library('clusterProfiler')
library(stringr)
library('org.Hs.eg.db')

setwd('E:/ManuscriptI/script')
rm(list = ls())

load('ProteinAbundance.Rdata')
load('CorrelationMap.Rdata')

GN<-gsub('^(.*)GN=','',Description$Description)
GN<-gsub(' (.*)$','',GN)
GN<-data.frame(row.names = rownames(Description),GN=GN)

load('corr_matrix_K562_MALBAC-DT_paper.Rdata')
correlation_matrix<-correlation_matrix[
  -grep('\\.',rownames(correlation_matrix)),
  -grep('\\.',colnames(correlation_matrix))]
correlation_matrix<-correlation_matrix[
  -grep('LINC',rownames(correlation_matrix)),
  -grep('LINC',colnames(correlation_matrix))]
dim(correlation_matrix)
RNA<-rownames(correlation_matrix)[rownames(correlation_matrix)%in%GN$GN]

## F4a
plotdata<-reshape2::melt(RNA)$value
plotdata<-data.frame(plotdata)
plotdata$omic<-'RNA'
temp<-reshape2::melt(CorMap)$value
plotdata<-rbind(data.frame(plotdata=temp,omic='MS'),plotdata)
ggplot(plotdata)+geom_density(aes(x=plotdata,group=omic,fill=omic),alpha=0.7)+
  theme(text=element_text(size=10,family='sans'))+
  xlab('correlation value')+ylab('density')+
  scale_fill_manual(values=c(RNA='#54a965',MS='#dc8350'))

#F4bF4c
ribo<-c("P46779", "P62847", "P62913", "P62829", "P83731", "P35268", "P36578",
        "P18124", "P61353", "P46776", "P15880", "P62244", "P23396", "P46781",
        "P62851", "P62269", "P62277", "Q07020", "P39023", "P62424", "P62701", 
        "P27635", "P61247", "P62753", "P62917", "P13489", "Q02878", "P61313",
        "P18077", "P62888", "P26373", "P62899", "P46777", "P40429", "P46778", 
        "Q02543", "P84098")
Annot_GN<-GN[ribo,'GN']
ribo<-ribo[Annot_GN%in%rownames(RNA)]
Annot_GN<-Annot_GN[Annot_GN%in%rownames(RNA)]
Annot_GN
PlotCorMap<-CorMap[ribo,ribo]
rownames(PlotCorMap)<-Annot_GN;colnames(PlotCorMap)<-Annot_GN
pheatmap(PlotCorMap,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 8, cellheight = 8)
plotdata<-reshape2::melt(CorMap)
plotdata$type<-0            
temp<-reshape::melt(PlotCorMap)
temp$type<-1
plotdata<-rbind(plotdata,temp)
plotdata<-plotdata[plotdata$value<1,]
ggplot(plotdata)+geom_density(aes(x=value,group=type,fill=type,alpha=0.5))+
  guides(fill=guide_legend(title=NULL))

PlotCorMap<-RNA[Annot_GN,Annot_GN]
rownames(PlotCorMap)<-Annot_GN;colnames(PlotCorMap)<-Annot_GN
pheatmap(PlotCorMap,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 8, cellheight = 8)
plotdata<-reshape2::melt(RNA)
plotdata$type<-0            
temp<-reshape::melt(PlotCorMap)
temp$type<-1
colnames(temp)<-colnames(plotdata)
plotdata<-rbind(plotdata,temp)
plotdata<-plotdata[plotdata$value<1,]
ggplot(plotdata)+geom_density(aes(x=value,group=type,fill=type,alpha=0.5))
Annot_GN

## F4d
index<-rownames(group)[group==14]
Annot_GN<-GN[index,'GN']
index<-index[Annot_GN%in%rownames(RNA)==T]
Annot_GN<-GN[index,'GN']
PlotCorMap<-CorMap[index,index]
rownames(PlotCorMap)<-Annot_GN;colnames(PlotCorMap)<-Annot_GN
pheatmap(PlotCorMap,
         color = colorRampPalette(c("blue","white","red"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 8, cellheight = 8)
PlotCorMap<-RNA[Annot_GN,Annot_GN]
rownames(PlotCorMap)<-Annot_GN;colnames(PlotCorMap)<-Annot_GN
pheatmap(PlotCorMap,
         color = colorRampPalette(c("blue","white","red"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 8, cellheight = 8)



#F4e
load('BulkData.Rdata')
plot_CorMap<-CorMap
plot_CorMap<-plot_CorMap[which(duplicated(GN$GN)==F),which(duplicated(GN$GN)==F)]
rownames(plot_CorMap)<-GN[rownames(plot_CorMap),'GN']
colnames(plot_CorMap)<-GN[colnames(plot_CorMap),'GN']

G19<-rownames(group)[group==19]
data<-Bulk[Bulk$Protein_1%in%G19&
             Bulk$Protein_2%in%G19,]
length(unique(c(data$Protein_1,data$Protein_2)))
G19<-G19[G19%in%unique(c(data$Protein_1,data$Protein_2))]
data.2<-data.frame(Protein_1<-data$Protein_2,
                   Protein_2<-data$Protein_1,
                   coregulation_score<-data$coregulation_score)
colnames(data.2)<-colnames(data)
data<-rbind(data,data.2)
rm(data.2)
head(data)
data<-data[duplicated(data[,c(1,2)])==F,]
data<-tidyr::spread(data,Protein_2,coregulation_score)
rownames(data)<-data[,1];data<-data[,-1]
data[is.na(data)]<-0
data<-data[G19,G19]
for(i in 1:dim(data)[1]){
  data[i,i]<-0.3
}
rownames(data)<-GN[rownames(data),'GN']
pheatmap(data,
         color = colorRampPalette(c("white","red"))(101),
         breaks = seq(from = 0, to = 0.3, length.out = 101),
         border_color = NA,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 10, cellheight = 10,
         show_rownames = T, show_colnames = F)
data<-CorMap[G19,G19]
pheatmap(data,
         color = colorRampPalette(c('blue',"white","red"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 10, cellheight = 10,
         show_rownames = T, show_colnames = F)

