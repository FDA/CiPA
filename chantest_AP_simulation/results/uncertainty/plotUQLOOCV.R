rm(list=ls())

library(rms)
require(ROCR)
set.seed(1)
library(ggplot2)
scale<-"sara"
CL<-2000
metric<-"APD50_carest"
dataset<-"pureManual4"
outdir<-sprintf("results/%s/scale_%s/CL%d/",dataset,scale,CL)
figdir<-sub("results","figs",outdir)
system(paste0("mkdir -p ",outdir)); system(paste0("mkdir -p ", figdir))

valdf<-read.delim(paste0("../",dataset,"_AP_simulation/results/uncertainty/",metric,"LOOCV_allprobes.csv"),sep=",",as.is=T)
valdf$classidx<-valdf$CiPA
valdf$class<- 2- valdf$classidx #just like before, for plotting

thresholdtable<-read.delim(paste0("../",dataset,"_AP_simulation/results/uncertainty/",metric,"LOOCV_errors.csv"),sep=",",as.is=T)
t1<- with(thresholdtable,-intercept1/slope)
t2<- with(thresholdtable, -intercept2/slope)  #note that because in mycomputTdP_error.R the order is low (0), int(1),high(2) while in plotUQclassification.R the order is reversed, the t1/t2 here is of reversed order compared to those in plotUQclassification.R!
#use formal math
ebeta1<-with(thresholdtable,exp(-intercept1)); ebeta2<- with(thresholdtable,exp(-intercept2)); slope<-with(thresholdtable,slope)
t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/slope
t2 <- log(ebeta2-2*ebeta1)/slope
thresholds<- c(t1,t2)

#AUC of ROC?
forROC1<-valdf[,c("inter_prob","classidx","sample", "drug",metric)]
  forROC1$inter_prob<- rowSums(valdf[,c("inter_prob","high_prob")])
  colnames(forROC1)[1]<-"positive_prob"
  forROC1$class <- forROC1$classidx > 0  # now high and intermediate are both "TRUE" or "positive"                           
   forROC1predictions<-do.call(cbind,by(forROC1[,"positive_prob"],forROC1$sample,function(x) x))
   forROC1labels<-do.call(cbind,by(forROC1$class,forROC1$sample,function(x) x))
  ROC1obj<-prediction(forROC1predictions, forROC1labels)        
    AUC1<-performance(ROC1obj,"auc")  #the value is in y.values
    
forROC2<-valdf[,c("high_prob","classidx","sample","drug",metric)]
  colnames(forROC2)[1]<-"positive_prob"
  forROC2$class <- forROC2$classidx > 1                         
   forROC2predictions<-do.call(cbind,by(forROC2[,"positive_prob"],forROC2$sample,function(x) x))  
   forROC2labels<-do.call(cbind,by(forROC2$class,forROC2$sample,function(x) x))
  ROC2obj<-prediction(forROC2predictions, forROC2labels)        
    AUC2<-performance(ROC2obj,"auc")  #the value is in y.values 

if(1>10){
#decided to perform on the metric directly
  forROC1predictions<-do.call(cbind,by(forROC1[,metric],forROC1$sample,function(x) x))
  if(metric == "qNet" || metric == "model5")
    forROC1predictions<- -forROC1predictions
  ROC1obj<-prediction(forROC1predictions, forROC1labels)        
    AUC1<-performance(ROC1obj,"auc")  #the value is in y.values

 forROC2predictions<-do.call(cbind,by(forROC2[,metric],forROC2$sample,function(x) x))
  if(metric == "qNet" || metric == "model5")
    forROC2predictions<- -forROC2predictions
  ROC2obj<-prediction(forROC2predictions, forROC2labels)        
    AUC2<-performance(ROC2obj,"auc")  #the value is in y.values 
 
}
sens1<- with(valdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
spec1<-  with(valdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
LRplus1<- sens1/(1-spec1)
LRminus1<- (1-sens1)/spec1

sens2<- with(valdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
spec2<-  with(valdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
LRplus2<- sens2/(1-spec2)
LRminus2<- (1-sens2)/spec2

#try using random sampling to calculate LR?
sens1UQ<- by(valdf, valdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
spec1UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
LRplus1UQ<- (sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(1-spec1UQ+rnorm(length(sens1UQ),1e-6,1e-12))
LRminus1UQ<- (1-sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(spec1UQ+rnorm(length(sens1UQ),1e-12,1e-12))
quantile(LRplus1UQ,prob=c(0.025,0.5,0.975))

sens2UQ<- by(valdf, valdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
spec2UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
LRplus2UQ<- (sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(1-spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
LRminus2UQ<- (1-sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
quantile(LRplus2UQ,prob=c(0.025,0.5,0.975)) 


desiredorder<-drugnames<-c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "disopyramide", "vandetanib",
             "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","domperidone","astemizole", "pimozide","droperidol","clozapine",
             "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine")

valdf$drug<- factor(valdf$drug, levels=rev(desiredorder))

figfile<-paste0(figdir,metric,"_violinLOOCV.pdf")
pdf(figfile, width=8, height=4)
p<-ggplot(valdf, aes(x=drug, color=as.character(class), fill=as.character(class)))
#p<-p+geom_violin(aes_string(y=metric), alpha=0.5,width=10, position=position_dodge(width=0.3))
p<-p+geom_violin(aes_string(y=metric), alpha=0.5,scale="width", width=0.5)
#p<-p+geom_violin(aes_string(y=metric), alpha=0.5, adjust=10,width=10)
for(thresh in thresholds)
    p<-p+geom_hline(yintercept=thresh, linetype="dotted")
p<-p+ylab(metric)
p<-p+coord_flip()                      #this is how to set xlim when coord_flip: put "ylim=c()" inside the ()
p<-p+scale_color_brewer(NULL, palette="Set1")
p<-p+scale_fill_brewer(NULL, palette="Set1")
p<-p+theme_bw()
print(p)
dev.off()


