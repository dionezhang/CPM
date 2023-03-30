library(data.table);library(pheatmap)
library(ggplot2);library(hrbrthemes)
library(dplyr);library(RColorBrewer)

library(BiocManager)
library('clusterProfiler')
library(stringr)
library('org.Hs.eg.db')

setwd('E:/project/202106_ReCalculateCorrelation')
rm(list = ls())

load('0308_InitializedProteinAbundance.Rdata')
load('HEK_Abundance.Rdata')

head(HEK)
head(Original)
protein<-intersect(rownames(HEK),rownames(Original))
HEK<-HEK[protein,]
K562<-Original[protein,]
type<-c(rep('HEK',dim(HEK)[2]),rep('K56',dim(K562)[2]))
data<-cbind(HEK,K562)
Mean<-apply(data,2,median)
data<-sweep(data,2,Mean/mean(Mean),'/')

library('FactoMineR')
library('factoextra')
pcap<-PCA(scale(t(data)))
fviz_pca_ind(pcap)
fviz_pca_ind(pcap,col.ind = type,geom.ind='point')
