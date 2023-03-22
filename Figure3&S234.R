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
load('CPMData.Rdata')
load('293TCorrelationMap.Rdata')
load('OocyteCorrelationMap.Rdata')

GN<-gsub('^(.*)GN=','',Description$Description)
GN<-gsub(' (.*)$','',GN)
GN<-data.frame(row.names = rownames(Description),GN=GN)

#S2a
data<-CorMap; data$Protein2<-rownames(data)
data<-reshape2::melt(data)
data<-data[which(data$Protein2!=data$variable),]
data$group<-'outside CPM'
data[group[data$Protein2,'group']==group[data$variable,'group']&
       group[data$Protein2,'group']>0,'group']<-'within CPM'
data$group<-as.factor(data$group)
ggplot(data)+geom_density(aes(x=value,group=group,fill=group),alpha=0.5)+
  xlab('Pearson correlation coefficient')+
  guides(fill=guide_legend(title=NULL))

#F3a S4a S4b
breaks<-read.table('CorrelationGroup.txt',sep='\t')
breaks<-strsplit(breaks[,1],' ')
breaks<-unlist(breaks)
breaks<-as.numeric(breaks)
index<-apply(group,1,function(x) ifelse(x>0,1,0))
mycolors<-apply(group,1,function(x) ifelse(x>0,'#c44e53','white'))
names(mycolors)<-group$group
ann_colors<-list(group=mycolors)
pheatmap(CorMap,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 0.2, cellheight = 0.2,
         show_rownames = F, show_colnames = F,
         annotation_col=group,annotation_row=group,
         annotation_colors = ann_colors,annotation_legend = F)


#F3b S4d S4e S4f
index<-19
Annot_GN<-GN[rownames(group)[group$group==index],'GN']
PlotCorMap<-CorMap[which(group$group==index),
                   which(group$group==index)]
#S2b S2c
index<-11&57
Annot_GN<-GN[rownames(group)[group$group==37|group$group==57],'GN']
PlotCorMap<-CorMap[which(group$group==37|group$group==57),
                   which(group$group==37|group$group==57)]
rownames(PlotCorMap)<-Annot_GN;colnames(PlotCorMap)<-Annot_GN
svg(filename = paste('E:/ManuscriptI/PlotFigure/CPM',index,'.svg',sep=''),
    width=10,height=10)
pheatmap(PlotCorMap,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 10, cellheight = 10,
         fontfamily='sans',fontsize=10)
graphics.off()

#F2c
load('Enrichment.Rdata')
head(go_stat)
plotdata<-go_stat
plotdata$'PPI Ratio'<-'0'
plotdata[plotdata$PPI>0,'PPI Ratio']<-'1'
plotdata[plotdata$PPI>0.2,'PPI Ratio']<-'2'
plotdata[plotdata$PPI>0.5,'PPI Ratio']<-'2'
plotdata[plotdata$PPI>0.7,'PPI Ratio']<-'3'
plotdata$'GOBP Ratio'<-'0'
plotdata[plotdata$GOBP>0,'GOBP Ratio']<-'1'
plotdata[plotdata$GOBP>0.5,'GOBP Ratio']<-'2'
plotdata[plotdata$GOBP>0.7,'GOBP Ratio']<-'3'
plotdata$'GOCC Ratio'<-'0'
plotdata[plotdata$GOCC>0,'GOCC Ratio']<-'1'
plotdata[plotdata$GOCC>0.5,'GOCC Ratio']<-'2'
plotdata[plotdata$GOCC>0.7,'GOCC Ratio']<-'3'
plotdata$'KEGG Ratio'<-'0'
plotdata[plotdata$KEGG>0,'KEGG Ratio']<-'1'
plotdata[plotdata$KEGG>0.5,'KEGG Ratio']<-'2'
plotdata[plotdata$KEGG>0.7,'KEGG Ratio']<-'3'
plotdata$'Reactome Ratio'<-'0'
plotdata[plotdata$Reactome>0,'Reactome Ratio']<-'1'
plotdata[plotdata$Reactome>0.5,'Reactome Ratio']<-'2'
plotdata[plotdata$Reactome>0.7,'Reactome Ratio']<-'3'
plotdata$Module<-rownames(plotdata)
plotdata<-plotdata[,-c(1:5)]
colnames(plotdata)
plotdata<-reshape2::melt(plotdata[,c('GOBP Ratio','GOCC Ratio','KEGG Ratio',
                                     'Reactome Ratio','PPI Ratio','Module')],id.var='Module')
plotdata<-plotdata[order(plotdata$value),]
plotdata$value<-as.factor(plotdata$value)
mycolors=c('#c44e53','#4d72af','#dc8350','#54a965')
ggplot(plotdata)+geom_bar(aes(fill=value,x=variable),stat='count',
                          position='fill',width=0.5)+
  scale_fill_manual('value',values=mycolors,
                    labels=c('No Enrichment',
                             'Enrichment Ratio >0',
                             'Enrichment Ratio >0.5',
                             'Enrichment Ratio >0.7'))+
  theme(axis.title.y=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.title.x=element_blank(),legend.title = element_blank(),
        text=element_text(size=10,family='sans'))
hist(go_stat$TotalProtein)

#S3
plotdata<-go_stat
plotdata[plotdata>0]<-1
plotdata<-apply(plotdata,2,as.numeric)
heatmap(plotdata,cluster_cols=F,cluster_rows=F)


#S4c
data<-c('AARS1','MARS1','PSAT1','WARS1','ASNS','GARS1','YARS1')
for (i in 1:7){
  print(rownames(Description)[grep(data[i],Description$Description)])}
data<-c('P49588','P56192','Q9Y617','P23381','P08243','P41250','P54577')
HEK<-HEK[data,]
K<-Norm[data,]
plotdata<-cbind(HEK,K)
MEAN<-apply(plotdata,2,sum)
plotdata<-sweep(plotdata,2,MEAN/mean(MEAN),'/')
pheatmap(plotdata,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         cluster_row=F,cluster_cols=F,legend=F,border_color = F,
         show_rownames = F,show_colnames = F
)
breaks = seq(from = -1, to = 1, length.out = 101),
border_color = NA,legend=F,show_colnames = F,show_rownames = F,
cluster_rows = F, cluster_cols = F,
cellwidth = 1, cellheight = 1,
fontfamily='sans',fontsize=10)


##F3.4
load('CorrelationPair.Rdata')
load('PPI.Rdata')

data<-PPI[PPI$Uniprot1!=PPI$Uniprot2,]
step<-c(-0.5,-0.25,0,0.25,0.5,0.75)
# step<-seq(-0.5,1,0.1)
frac<-c()
for(i in step){
  frac<-c(frac,
          sum(data[data$cor>i&data$cor<min(step[step>i]),'combined_score']>700)/
            length(which(Pair$value>i&Pair$value<min(step[step>i]))))
}
plotdata<-data.frame(step=step,frac=frac)
plotdata<-rbind(plotdata,data.frame(step=c(1,1.1),frac=c(0.12,0.47)))  
ggplot(plotdata)+geom_bar(aes(x=step,y=frac),stat = 'identity')

              
