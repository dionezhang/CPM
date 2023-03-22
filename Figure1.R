library(data.table);library(pheatmap)
library(ggplot2);library(hrbrthemes)
library(dplyr);library(RColorBrewer)

setwd('E:/project/202106_ReCalculateCorrelation')
rm(list = ls())


## F1b
plotdata<-data.frame(colnames=c('type','copy','p','error'))
Error<-function(x){
  return(sd(x)/mean(x))
}
RNA<-rpois(1000,17)
RNA_0.2<-rbinom(1000,RNA,0.2)
RNA_0.5<-rbinom(1000,RNA,0.5)
Protein<-rnbinom(100,2.5e2,0.005)
Protein_0.01<-rbinom(100,Protein,0.01)
Protein_0.005<-rbinom(100,Protein,0.005)
plot_RNA<-data.frame(type='RNA',copy=c(RNA_0.5,RNA_0.2),
                     p=c(rep('0.5',1000),rep('0.2',1000)))
plot_Protein<-data.frame(type='Protein',copy=c(Protein_0.01,Protein_0.005),
                         p=c(rep('0.01',100),rep('0.005',100)))
plotdata<-rbind(plot_RNA,plot_Protein)
plotdata$copy<-log10(plotdata$copy)
# plotdata$copy<-round(plotdata$copy,1)
plotdata$p<-factor(plotdata$p,levels=c('0.5','0.2','0.01','0.005'))
Error<-function(x){
  return(sd(x)/mean(x))
}
Error(RNA_0.2)
Error(RNA_0.5)
Error(Protein_0.01)
Error(Protein_0.005)
F1b<-ggplot(plotdata,aes(x=p,y=copy))+geom_violin(aes(fill=type),bw=0.1)+ylim(0,3)+
  geom_boxplot(width=0.05,color='black',alpha=0.2)

## F1c
Simulation_RNA<-function(size,p,round){
  set.seed(123)
  a_real<-rpois(size,20)
  delta<-as.data.frame(0)
  b_real<-matrix(0,nrow=size,ncol=50)
  cor_real<-rep(0,50)
  for (i in 1:50){
    set.seed(i+p*10)
    delta<-rnorm(size,mean=20,sd=(i-1)*2)
    b_real[,i]<-delta+a_real
    cor_real[i]<-cor(a_real,b_real[,i])
  }
  a_real<-round(a_real,0)
  b_real<-round(b_real,0)
  b_real[b_real<0]<-0
  
  set.seed(round)
  a<-rbinom(size,a_real,p)
  b<-matrix(nrow=nrow(b_real),ncol=ncol(b_real))
  cor_test<-rep(0,50)
  for(i in 1:50){
    set.seed(round+i)
    b[,i]<-rbinom(size,b_real[,i],p)
    cor_test[i]<-cor(a,b[,i])
  }
  data.frame(real=cor_real,test=cor_test,
             group=as.factor(size),
             detectability=as.factor(p))
}
Repeat<-function(size,p){
  result<-Simulation_RNA(size,p,100)
  for( i in 1:99){
    result[,c(1:2)]<-result[,c(1:2)]+Simulation_RNA(size,p,i)[,c(1:2)]
  }
  result[,c(1,2)]<-result[,c(1,2)]/100
  result
}
Subsample<-function(result){
  result<-result[order(result$real),]
  return(result[duplicated(round(result$real,1))==F,])
}

result<-Subsample(Repeat(1000,0.1))
result<-rbind(result,Subsample(Repeat(1000,0.3)))
result<-rbind(result,Subsample(Repeat(1000,0.5)))
result<-rbind(result,Subsample(Repeat(1000,1)))
F1c<-ggplot(result)+geom_line(aes(x=real,y=test,group=detectability,
                             color=detectability),size=1)+
  xlab('real correlation')+ylab('measured correlation')+
  theme(panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(limits = c(0,1),breaks=c(0,0.25,0.50,0.75,1))+
  scale_y_continuous(limits = c(0,1),breaks=c(0,0.25,0.50,0.75,1))

## F1d
Simulation_Protein<-function(size,p,round){
  set.seed(123)
  a_real<-rnbinom(size,5e3,0.005)
  delta<-as.data.frame(0)
  b_real<-matrix(0,nrow=size,ncol=50)
  cor_real<-rep(0,50)
  for (i in 1:50){
    set.seed(i+p*100)
    delta<-rnorm(size,mean=2e4,sd=i*2e3)
    b_real[,i]<-delta+a_real
    cor_real[i]<-cor(a_real,b_real[,i])
  }
  a_real<-round(a_real,0)
  b_real<-round(b_real,0)
  b_real[b_real<0]<-0
  
  set.seed(round)
  a<-rbinom(size,a_real,p)
  b<-matrix(nrow=nrow(b_real),ncol=ncol(b_real))
  cor_test<-rep(0,50)
  for(i in 1:50){
    set.seed(round+i)
    b[,i]<-rbinom(size,b_real[,i],p)
    cor_test[i]<-cor(a,b[,i])
  }
  data.frame(real=cor_real,test=cor_test,
             group=as.factor(size),
             detectability=as.factor(p))
}
Repeat<-function(size,p){
  result<-Simulation_Protein(size,p,100)
  for( i in 1:99){
    result[,c(1:2)]<-result[,c(1:2)]+Simulation_Protein(size,p,i)[,c(1:2)]
  }
  result[,c(1,2)]<-result[,c(1,2)]/100
  result
}
Subsample<-function(result){
  result<-result[order(result$real),]
  return(result[duplicated(round(result$real,1))==F,])
}
result<-Subsample(Repeat(100,0.005))
result<-rbind(result,Subsample(Repeat(100,0.01)))
result<-rbind(result,Subsample(Repeat(100,0.05)))
result<-rbind(result,Subsample(Repeat(100,1)))
F1d<-ggplot(result)+geom_line(aes(x=real,y=test,group=detectability,
                             color=detectability),size=1)+
  xlim(0.1,1)+xlab('real correlation')+ylab('measured correlation')+
  theme(panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(limits = c(0,1),breaks=c(0,0.25,0.50,0.75,1))+
  scale_y_continuous(limits = c(0,1),breaks=c(0,0.25,0.50,0.75,1))

