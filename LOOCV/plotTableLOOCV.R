rm(list=ls())

library(rms)
require(ROCR)
set.seed(1)
library(ggplot2)
scale<-"sara"
CL<-2000
metricv<-c("qNet")  #
datasetvec<-c("nanion8rand")  #AP_simulation name. for nanion8rand_AP_simulation is "nanion8rand"
measurevec<-c("AUC1","AUC2","Pairwise Correct Rate","LR+ of Threshold 1", "LR- of Threshold 1",
               "LR+ of Threshold 2", "LR- of Threshold 2", "Mean Classification Error")

all28df<-data.frame(metric=c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error"))
colordf<-data.frame(
   drug = c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "disopyramide", "vandetanib",
             "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","domperidone","astemizole", "pimozide","droperidol","clozapine",
             "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine"),
   classidx=c(rep(2,8),rep(1,11),rep(0,9)),
   coloridx=c(rep(1,8),rep(3,11),rep(2,9)),
   isTraining=c(1,1,1,1,0,0,0,0,
                1,1,1,1,0,0,0,0,0,0,0,
                1,1,1,1,0,0,0,0,0)
 )

for(dataset in datasetvec){
 for(metric in metricv){
   altstr<-paste0(dataset,"_",metric)
   altvec<-0
valdf<-read.delim(paste0("results/uncertainty/",metric,"LOOCV_allprobes.csv"),sep=",",as.is=T)

valdf$classidx<-valdf$CiPA
valdf$class<- 2- valdf$classidx 
#
forROC1<-valdf[,c("inter_prob","classidx","sample", "drug")]
  forROC1$inter_prob<- rowSums(valdf[,c("inter_prob","high_prob")])
  colnames(forROC1)[1]<-"positive_prob"
  forROC1$class <- forROC1$classidx > 0  # now high and intermediate are both "TRUE" or "positive"                           
   forROC1predictions<-do.call(cbind,by(forROC1[,"positive_prob"],forROC1$sample,function(x) x))
   forROC1labels<-do.call(cbind,by(forROC1$class,forROC1$sample,function(x) x))
  ROC1obj<-prediction(forROC1predictions, forROC1labels)        
    AUC1<-performance(ROC1obj,"auc")  #the value is in y.values
    AUC1dist<-formatC(quantile(unlist(AUC1@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(AUC1dist[2]," (",AUC1dist[1]," - ",AUC1dist[3],")"))
    
forROC2<-valdf[,c("high_prob","classidx","sample","drug")]
  colnames(forROC2)[1]<-"positive_prob"
  forROC2$class <- forROC2$classidx > 1                         
   forROC2predictions<-do.call(cbind,by(forROC2[,"positive_prob"],forROC2$sample,function(x) x))  
   forROC2labels<-do.call(cbind,by(forROC2$class,forROC2$sample,function(x) x))
  ROC2obj<-prediction(forROC2predictions, forROC2labels)        
    AUC2<-performance(ROC2obj,"auc")  #the value is in y.values 
     AUC2dist<-formatC(quantile(unlist(AUC2@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(AUC2dist[2]," (",AUC2dist[1]," - ",AUC2dist[3],")"))
	
#pairwise ranking
pairwisefun<- function(fulltable){
cmb <- combn(seq_len(nrow(fulltable)), 2)
mergedtable<-cbind(fulltable[cmb[1,],], fulltable[cmb[2,],])
validpairidx<- (mergedtable[,7]!=mergedtable[,16])&(!mergedtable[,9]|!mergedtable[,18])
correctidx1<- ((mergedtable[,7]>mergedtable[,16])&(mergedtable[,6]<mergedtable[,15]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,6]>mergedtable[,15])) #when predicted class are different
correctidx2<- (mergedtable[,6]==1)&(mergedtable[,15]==1)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,3]>mergedtable[,12]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,3]<mergedtable[,12]))) #when predicted class are both high
correctidx3<- (mergedtable[,6]==3)&(mergedtable[,15]==3)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both low
correctidx4<- (mergedtable[,6]==2)&(mergedtable[,15]==2)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both intermediate
correctidx<- correctidx1|correctidx2|correctidx3|correctidx4
sum(validpairidx&correctidx)/sum(validpairidx)
}
tempdf<-valdf[,c("drug","sample","high_prob","inter_prob","low_prob")]; 
tempdf$pred<- apply(valdf[,5:7],1,which.max) #note here 1 is high, 2 is int, 3 is low
tempdf<-merge(tempdf,colordf,by="drug")
distpairwise<-by(tempdf, valdf$sample, pairwisefun)
  pairwisedist<-formatC(quantile(unlist(distpairwise),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(pairwisedist[2]," (",pairwisedist[1]," - ",pairwisedist[3],")"))
  
sens1<- with(valdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
spec1<-  with(valdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
LRplus1<- sens1/(1-spec1)
LRminus1<- (1-sens1)/spec1

sens2<- with(valdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
spec2<-  with(valdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
LRplus2<- sens2/(1-spec2)
LRminus2<- (1-sens2)/spec2

#try using random sampling to calculate LR
sens1UQ<- by(valdf, valdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
spec1UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
LRplus1UQ<- (sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(1-spec1UQ+rnorm(length(sens1UQ),1e-6,1e-12))
 LRplus1dist<-formatC(quantile(unlist(LRplus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(LRplus1dist[2]," (",LRplus1dist[1]," - ",LRplus1dist[3],")"))
LRminus1UQ<- (1-sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(spec1UQ+rnorm(length(sens1UQ),1e-12,1e-12))
 LRminus1dist<-formatC(quantile(unlist(LRminus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(LRminus1dist[2]," (",LRminus1dist[1]," - ",LRminus1dist[3],")"))

sens2UQ<- by(valdf, valdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
spec2UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
LRplus2UQ<- (sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(1-spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
 LRplus2dist<-formatC(quantile(unlist(LRplus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(LRplus2dist[2]," (",LRplus2dist[1]," - ",LRplus2dist[3],")"))
LRminus2UQ<- (1-sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
 LRminus2dist<-formatC(quantile(unlist(LRminus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(LRminus2dist[2]," (",LRminus2dist[1]," - ",LRminus2dist[3],")"))

    meanerror<- mean(valdf$pred_err)
    lowererror<- meanerror - 1.96* sd(valdf$pred_err)/sqrt(dim(valdf)[1]); uppererror<-meanerror+1.96*sd(valdf$pred_err)/sqrt(dim(valdf)[1])
    errordist<-formatC(c(meanerror,lowererror,uppererror),digits=3,format="g",flag="-",drop0trailing=T)
    altvec<-c(altvec,paste0(errordist[1]," (",errordist[2]," - ",errordist[3],")"))

   altvec<-altvec[-1]
   all28df[,altstr]<-altvec
  }
 }
 
write.table(all28df, "Result_PlotTableLOOCV/al28df.txt",col.names=T,row.names=F,sep="\t",quote=F)
table7<-data.frame(metric=rep(measurevec,each=length(datasetvec)), dataset=rep(datasetvec, length(measurevec)),
                   qNet=NA)
table7[,-1:-2]<- matrix(unlist(t(all28df[,-1])),nrow= length(measurevec)*length(datasetvec),byrow=T)
write.table(table7, "Result_PlotTableLOOCV/Table7.txt",col.names=T,row.names=F,sep="\t",quote=F)
