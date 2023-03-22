library(data.table);library(pheatmap)
library(ggplot2);library(hrbrthemes)
library(dplyr);library(RColorBrewer)
library(reshape2)

library(BiocManager)
library('clusterProfiler')
library(stringr)
library('org.Hs.eg.db')

setwd('E:/ManuscriptI/script')
rm(list = ls())

load('ProteinAbundance.Rdata')
load('CorrelationMap.Rdata')

###PPI Enrichment
# PPI<-read.table('E:/data/9606.protein.links.detailed.v11.0.txt',sep=' ')
# 
# genelist<-rownames(Norm)
# gene.df<-bitr(genelist,fromType='UNIPROT',toType='ENSEMBL',OrgDb=org.Hs.eg.db)
# 
# library(STRINGdb)
# string_db<-STRINGdb$new(version='11',species=9606)
# data_mapped<-string_db$map(gene.df,'ENSEMBL',removeUnmappedRows = T)
# hit<-data_mapped$STRING_id
# data_mapped$STRING_id<-gsub('9606.','',data_mapped$STRING_id)
# head(data_mapped)
# 
# head(PPI)
# PPI<-PPI[PPI$protein1%in%data_mapped$STRING_id&
#            PPI$protein2%in%data_mapped$STRING_id,]
# PPI<-na.omit(PPI)
# length(unique(unique(PPI[,c(1)]),unique(PPI[,c(2)])))
# a<-unique(PPI$protein1)
# index<-0
# for (i in 1:length(a)){
#   index<-index+sum(str_count(PPI$protein1,a[i]))
#   temp<-PPI[index+1:dim(PPI)[1],]
#   temp<-temp[-grep(a[i],temp$protein2),]
#   PPI<-rbind(PPI[1:index,],temp)
#   PPI<-na.omit(PPI)
#   print(index)
# }
# PPI$Uniprot1<-data_mapped[match(PPI$protein1,data_mapped$STRING_id),'UNIPROT']
# PPI$Uniprot2<-data_mapped[match(PPI$protein2,data_mapped$STRING_id),'UNIPROT']
# PPI$cor<-0
# for (i in 1:dim(PPI)[1]){
#   PPI[i,'cor']<-CorMap[PPI[i,'Uniprot1'],PPI[i,'Uniprot2']]
#   print(i)
# }
# hist(PPI$cor)
# save(PPI,file='PPI.Rdata')
load('PPI.Rdata')
PPI<-PPI[PPI$combined_score>=250,]
PPI$group_A<-group[PPI$Uniprot1,'group']
PPI$group_B<-group[PPI$Uniprot2,'group']
table(group)
table(PPI$group_A)
table(PPI$group_A,PPI$group_B)
PPIGroup<-table(PPI$group_A,PPI$group_B)
PPIGroup<-as.data.frame(PPIGroup)
head(PPIGroup)
library(reshape2)
PPIGroup<-as.data.frame(acast(PPIGroup,Var1~Var2))
groupcount<-as.data.frame(table(group))
rownames(groupcount)<-groupcount$group
PPIGroup$TotalCount<-groupcount[rownames(PPIGroup),'Freq']
head(PPIGroup)
PPIRatio<-data.frame(row.names = colnames(PPIGroup))
PPIRatio<-PPIRatio[-grep('Total',rownames(PPIRatio)),]
PPIRatio$PPIRatio<-0;PPIRatio$MeanCor<-0
for (i in 1:67){
  print(i)
  PPIRatio[i,'PPIRatio']<-PPIGroup[i,i]/choose(PPIGroup[i,'TotalCount'],2)
  j<-group$group==rownames(PPIRatio)[i]
  if(length(j)==0){next}
  PPIRatio[i,'MeanCor']<-
    mean(reshape2::melt(CorMap[j,j])[reshape2::melt(CorMap[j,j])$value<1,]$value)
}
save(PPIRatio,file='PPIRatio.Rdata')

go_temp<-enrichGO(rownames(group)[group$group==1],OrgDb = 'org.Hs.eg.db',ont = 'CC',
                  pvalueCutoff = 0.05,keyType = 'UNIPROT' )
go_temp@result$module<-1
go_cc<-go_temp@result
for(i in 2:66){
  go_temp<-enrichGO(rownames(group)[group$group==i],OrgDb = 'org.Hs.eg.db',ont = 'CC',
                    pvalueCutoff = 0.05,keyType = 'UNIPROT' )
  print(i)
  go_temp<-simplify(go_temp)
  go_temp@result$module<-as.numeric(i)
  go_cc<-rbind(go_cc,go_temp@result)
}
go_cc$score<-0
for (i in 1:length(rownames(go_cc))){
  a<-unlist(strsplit(go_cc[i,'GeneRatio'],split='/'))[1]
  b<-unlist(strsplit(go_cc[i,'GeneRatio'],split='/'))[2]
  go_cc[i,'score']<-round(as.numeric(a)/as.numeric(b),2)
}
save(go_cc,file='CCEnrichmentResult.Rdata')

go_temp<-enrichGO(rownames(group)[group$group==1],OrgDb = 'org.Hs.eg.db',
                  ont = 'BP',pvalueCutoff = 0.05,keyType = 'UNIPROT' )
go_temp@result$module<-1
go_bp<-go_temp@result
for(i in 2:66){
  go_temp<-enrichGO(rownames(group)[group$group==i],OrgDb = 'org.Hs.eg.db',
                    ont = 'BP',pvalueCutoff = 0.05,keyType = 'UNIPROT' )
  print(i)
  go_temp<-simplify(go_temp)
  go_temp@result$module<-as.numeric(i)
  go_bp<-rbind(go_bp,go_temp@result)
}
go_bp$score<-0
for (i in 1:length(rownames(go_bp))){
  a<-unlist(strsplit(go_bp[i,'GeneRatio'],split='/'))[1]
  b<-unlist(strsplit(go_bp[i,'GeneRatio'],split='/'))[2]
  go_bp[i,'score']<-round(as.numeric(a)/as.numeric(b),2)
}
save(go_bp,file='BPEnrichmentResult.Rdata')

data<-rownames(Norm)
data<-bitr(rownames(Norm),fromType='UNIPROT',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
data<-data$ENTREZID
ego<-enrichKEGG(gene=data,keyType = 'kegg',organism='hsa',pvalue=0.05)
data$group<-group[data$UNIPROT,'group']
ego_temp<-enrichKEGG(gene=data[data$group==1,'ENTREZID'],
                     keyType = 'kegg',organism='hsa',pvalue=0.05)
ego_temp<-ego_temp@result  
ego_temp$module<-1  
kegg<-ego_temp  
for (i in 2:66){
  print(i)
  ego_temp<-enrichKEGG(gene=data[data$group==i,'ENTREZID'],
                       keyType = 'kegg',organism='hsa',pvalue=0.9)
  ego_temp<-ego_temp@result  
  ego_temp$module<-i  
  kegg<-rbind(kegg,ego_temp)
}  
save(kegg,file='KEGGResult.Rdata')

library(ReactomePA)
id<-bitr(rownames(group),fromType='UNIPROT',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
id<-id[duplicated(id$UNIPROT)==F,]
rownames(id)<-id$UNIPROT
re_temp<-enrichPathway(id[rownames(group)[group$group==1],'ENTREZID']
                          ,organism='human',pvalueCutoff = 0.05 )
re_temp@result$module<-1
Reactome<-re_temp@result
for(i in 2:66){
  re_temp<-enrichPathway(id[rownames(group)[group$group==1],'ENTREZID']
                         ,organism='human',pvalueCutoff = 0.05 )
  print(i)
  re_temp<-simplify(re_temp)
  re_temp@result$module<-as.numeric(i)
  Reactome<-rbind(Reactome,re_temp@result)
}
save(Reactome,file='Reactome.Rdata')

load('CCEnrichmentResult.Rdata')
load('BPEnrichmentResult.Rdata')
load('KEGGResult.Rdata')
load('ReactomeResult.Rdata')
load('PPIRatio.Rdata')

go_bp<-go_bp[go_bp$pvalue<0.01&go_bp$Count>1,]
go_cc<-go_cc[go_cc$pvalue<0.01&go_cc$Count>1,]
kegg<-kegg[kegg$pvalue<0.01&kegg$Count>1,]
Reactome<-Reactome[Reactome$pvalue<0.01&Reactome$Count>1,]

colnames(go_bp)
go_bp$score<-
  as.numeric(apply(go_bp,1,function(x) unlist(strsplit(x[3],split='/'))[1]))/
  as.numeric(apply(go_bp,1,function(x) length(rownames(group)[group$group==as.numeric(x[10])])))
go_bp$score<-
  as.numeric(apply(go_bp,1,function(x) unlist(strsplit(x[3],split='/'))[1]))/
  as.numeric(apply(go_bp,1,function(x) unlist(strsplit(x[3],split='/'))[2]))

colnames(kegg)
kegg$score<-
  as.numeric(apply(kegg,1,function(x) unlist(strsplit(x[3],split='/'))[1]))/
  as.numeric(apply(kegg,1,function(x) length(rownames(group)[group$group==as.numeric(x[10])])))
kegg$score<-
  as.numeric(apply(kegg,1,function(x) unlist(strsplit(x[3],split='/'))[1]))/
  as.numeric(apply(kegg,1,function(x) unlist(strsplit(x[3],split='/'))[2]))

go_cc$score<-
  as.numeric(apply(go_cc,1,function(x) unlist(strsplit(x[3],split='/'))[1]))/
  as.numeric(apply(go_cc,1,function(x) length(rownames(group)[group$group==as.numeric(x[10])])))

Reactome$score<-
  as.numeric(apply(Reactome,1,function(x) unlist(strsplit(x[3],split='/'))[1]))/
  as.numeric(apply(Reactome,1,function(x) length(rownames(group)[group$group==as.numeric(x[10])])))


go_stat<-data.frame(row.names=c(1:66))
colnames(go_cc)
go_stat$GOCC<-0
go_stat$GOBP<-0
go_stat$KEGG<-0
go_stat$Reactome<-0
for (i in c(1:66)){
  go_stat[i,'GOCC']<-ifelse(i%in%go_cc$module,
                           max(go_cc[go_cc$module==i,'score']),0)
  go_stat[i,'GOBP']<-ifelse(i%in%go_bp$module,
                           max(go_bp[go_bp$module==i,'score']),0)
  go_stat[i,'KEGG']<-ifelse(i%in%kegg$module,
                           max(kegg[kegg$module==i,'score']),0)
  go_stat[i,'Reactome']<-ifelse(i%in%Reactome$module,
                               max(Reactome[Reactome$module==i,'score']),0)
}
go_stat$PPI<-round(PPIRatio[-1,'PPIRatio'],2)
go_stat<-round(go_stat,2)
View(go_stat)
save(go_stat,file='Enrichment.Rdata')
