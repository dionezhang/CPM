
# S1a
TMT<-read.csv('K562-TMT-corr-matrix.csv')
rownames(TMT)<-TMT[,1]
TMT<-TMT[,-1]
pheatmap(TMT,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,show_rownames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 1, cellheight = 1,
         fontfamily='sans',fontsize=10)

#S1b
plotdata<-reshape::melt(CorMap)
plotdata$variable<-'LF'
temp<-reshape::melt(TMT)
temp$variable<-'TMT'
plotdata<-rbind(plotdata,temp)
plotdata<-plotdata[plotdata$value<1,]
ggplot(plotdata)+geom_violin(aes(x=variable,y=value,fill=variable))+
  scale_fill_manual(values=c('#c44e53','#4d72af'))+
  scale_y_continuous(breaks=seq(-0.75,0.75,0.25))+
  theme(axis.title.y=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#F2e F2f
head(TMT)
head(CorMap)
Ribosome<-rownames(Description)[grep('ribosomal protein',Description[,1])]
Ribosome<-Ribosome[Ribosome%in%rownames(CorMap)==T]
Ribosome_LF<-CorMap[Ribosome,Ribosome]
Ribosome<-Ribosome[Ribosome%in%rownames(TMT)==T]
Ribosome_TMT<-Ribosome_LF
Ribosome_TMT[Ribosome,Ribosome]<-TMT[Ribosome,Ribosome]
Ribosome_TMT[rownames(Ribosome_TMT)%in%Ribosome==F,c(1:77)]<-0
Ribosome_TMT[c(1:77),rownames(Ribosome_TMT)%in%Ribosome==F]<-0
pheatmap(Ribosome_TMT,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,show_rownames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 3, cellheight = 3,
         fontfamily='sans',fontsize=10)
pheatmap(Ribosome_LF,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,show_rownames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 3, cellheight = 3,
         fontfamily='sans',fontsize=10)
identical(rownames(Ribosome_TMT),rownames(Ribosome_LF))
plotdata<-reshape::melt(Ribosome_LF)
plotdata$variable<-'LF'
Ribosome<-reshape::melt(Ribosome_TMT[Ribosome,Ribosome])
Ribosome$variable<-'TMT'
plotdata<-rbind(plotdata,Ribosome)
plotdata<-plotdata[plotdata$value<1,]
ggplot(plotdata)+geom_violin(aes(x=variable,y=value,fill=variable))+
  scale_fill_manual(values=c('#c44e53','#4d72af'))+
  scale_y_continuous(breaks=seq(-0.75,0.75,0.25))+
  theme(axis.title.y=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

temp<-rownames(CorMap)[rownames(CorMap)%in%rownames(TMT)==F]
temp<-sample(temp,(77-53),replace=F)
Random<-c(temp,sample(rownames(TMT)[rownames(TMT)%in%rownames(CorMap)],53,replace=F))
Random<-sample(Random,77,replace=F)
Random_LF<-CorMap[Random,Random]
Random_TMT<-Random_LF
Random_TMT[Random[Random%in%temp==F],Random[Random%in%temp==F]]<-
  TMT[Random[Random%in%temp==F],Random[Random%in%temp==F]]
Random_TMT[temp,c(1:77)]<-0
Random_TMT[c(1:77),temp]<-0
pheatmap(Random_TMT,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,show_rownames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 3, cellheight = 3,
         fontfamily='sans',fontsize=10)
pheatmap(Random_LF,
         color = colorRampPalette(c("#4d72af","white","#c44e53"))(101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend=F,show_colnames = F,show_rownames = F,
         cluster_rows = F, cluster_cols = F,
         cellwidth = 3, cellheight = 3,
         fontfamily='sans',fontsize=10)
plotdata<-reshape::melt(Random_LF)
plotdata$variable<-'LF'
temp<-reshape::melt(Random_TMT[Random[Random%in%temp==F],Random[Random%in%temp==F]])
temp$variable<-'TMT'
plotdata<-rbind(plotdata,temp)
plotdata<-plotdata[plotdata$value<1,]
ggplot(plotdata)+geom_violin(aes(x=variable,y=value,fill=variable))+
  scale_fill_manual(values=c('#c44e53','#4d72af'))+
  scale_y_continuous(breaks=seq(-0.75,0.75,0.25))+
  theme(axis.title.y=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

